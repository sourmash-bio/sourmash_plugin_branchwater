/// Rust code for pyo3_branchwater.

use pyo3::prelude::*;

use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use anyhow::Result;

use log::{info, error};

pub mod utils;
pub mod manysearch;

use utils::{prepare_query, PrefetchResult, write_prefetch,
            load_sketchlist_filenames, load_sketches,
            load_sketches_above_threshold, consume_query_by_gather,
            is_revindex_database, read_signatures_from_zip,
            build_template};

use manysearch::manysearch;


#[macro_use]
extern crate simple_error;

use sourmash::signature::Signature;
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;
use sourmash::index::revindex::RevIndex;
use sourmash::signature::SigsTrait;

/// Run counter-gather with a query against a list of files.

fn countergather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query_filename: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    gather_output: Option<P>,
    prefetch_output: Option<P>,
) -> Result<()> {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    let query_label = query_filename.to_string();
    eprintln!("Loading query from '{}'", query_label);
    let query = {
        let sigs = Signature::from_path(query_filename)?;

        prepare_query(&sigs, &template)
    };

    // did we find anything matching the desired template?
    let query = match query {
        Some(query) => query,
        None => bail!("No sketch found with scaled={}, k={}", scaled, ksize),
    };

    // build the list of paths to match against.
    eprintln!("Loading matchlist from '{}'", matchlist_filename.as_ref().display());

    let matchlist_paths = load_sketchlist_filenames(&matchlist_filename)?;

    eprintln!("Loaded {} sig paths in matchlist", matchlist_paths.len());

    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes : u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }.try_into()?;

    eprintln!("using threshold overlap: {} {}",
              threshold_hashes, threshold_bp);

    // load a set of sketches, filtering for those with overlaps > threshold
    let result = load_sketches_above_threshold(matchlist_paths,
                                               &template,
                                               &query.minhash,
                                               threshold_hashes)?;
    let matchlist = result.0;
    let skipped_paths = result.1;
    let failed_paths = result.2;

    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} paths - no compatible signatures.",
                  skipped_paths);
    }
    if failed_paths > 0 {
        eprintln!("WARNING: {} signature paths failed to load. See error messages above.",
                  failed_paths);
    }

    if matchlist.is_empty() {
        eprintln!("No matchlist signatures loaded, exiting.");
        return Ok(());
    }

    if prefetch_output.is_some() {
        write_prefetch(query_label.clone(), prefetch_output, &matchlist).ok();
    }

    // run the gather!
    consume_query_by_gather(query.minhash, matchlist, threshold_hashes,
                            gather_output, query_label).ok();
    Ok(())
}

/// Run counter-gather for multiple queries against a list of files.

fn multigather<P: AsRef<Path> + std::fmt::Debug + Clone>(
    query_filenames: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
) -> Result<()> {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    // load the list of query paths
    let querylist_paths = load_sketchlist_filenames(&query_filenames)?;
    println!("Loaded {} sig paths in querylist", querylist_paths.len());

    // build the list of paths to match against.
    println!("Loading matchlist");
    let matchlist_paths = load_sketchlist_filenames(&matchlist_filename)?;
    println!("Loaded {} sig paths in matchlist", matchlist_paths.len());

    let threshold_hashes : u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }.try_into().unwrap();

    println!("threshold overlap: {} {}", threshold_hashes, threshold_bp);

    // Load all the against sketches
    let result = load_sketches(matchlist_paths, &template)?;
    let (sketchlist, skipped_paths, failed_paths) = result;

    eprintln!("Loaded {} sketches to search against.", sketchlist.len());
    if failed_paths > 0 {
        eprintln!("WARNING: {} search paths failed to load. See error messages above.",
                  failed_paths);
    }
    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} search paths - no compatible signatures.",
                  skipped_paths);
    }

    if sketchlist.is_empty() {
        bail!("No sketches loaded to search against!?")
    }

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    querylist_paths
        .par_iter()
        .for_each(|q| {
            // increment counter of # of queries
            let _i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);

            // set query_label to the last path element.
            let query_label = q.clone().into_os_string().into_string().unwrap();
            let query_label = query_label.split('/').last().unwrap().to_string();

            if let Some(query) = {
                // load query from q
                let mut mm = None;
                if let Ok(sigs) = Signature::from_path(dbg!(q)) {
                    mm = prepare_query(&sigs, &template);

                    if mm.is_none() {
                        eprintln!("WARNING: no compatible sketches in path '{}'",
                                  q.display());
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                } else {
                    eprintln!("WARNING: could not load sketches from path '{}'",
                              q.display());
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
                mm
            } {

            // filter first set of matches out of sketchlist
            let matchlist: BinaryHeap<PrefetchResult> = sketchlist
                .par_iter()
                .filter_map(|sm| {
                    let mut mm = None;

                    if let Ok(overlap) = sm.minhash.count_common(&query.minhash, false) {
                        if overlap >= threshold_hashes {
                            let result = PrefetchResult {
                                name: sm.name.clone(),
                                md5sum: sm.md5sum.clone(),
                                minhash: sm.minhash.clone(),
                                overlap,
                            };
                            mm = Some(result);
                        }
                    }
                    mm
                })
                .collect();

            if !matchlist.is_empty() {
                let prefetch_output = format!("{query_label}.prefetch.csv");
                let gather_output = format!("{query_label}.gather.csv");

                // save initial list of matches to prefetch output
                write_prefetch(query_label.clone(), Some(prefetch_output),
                               &matchlist).ok();

                // now, do the gather!
                consume_query_by_gather(query.minhash, matchlist, threshold_hashes,
                                        Some(gather_output), query_label).ok();
            } else {
                println!("No matches to '{}'", query_label);
            }
        } else {
            println!("ERROR loading signature from '{}'", query_label);
        }
        });


    println!("Processed {} queries total.", processed_queries.into_inner());

    let skipped_paths = skipped_paths.into_inner();
    let failed_paths = failed_paths.into_inner();

    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} query paths - no compatible signatures.",
                  skipped_paths);
    }
    if failed_paths > 0 {
        eprintln!("WARNING: {} query paths failed to load. See error messages above.",
                  failed_paths);
    }
        
    Ok(())
}


fn index<P: AsRef<Path>>(
    siglist: P,
    template: Sketch,
    output: P,
    save_paths: bool,
    colors: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut temp_dir = None;
    info!("Loading siglist");

    let index_sigs: Vec<PathBuf>;

    if siglist.as_ref().extension().map(|ext| ext == "zip").unwrap_or(false) {
        let (paths, tempdir) = read_signatures_from_zip(&siglist)?;
        temp_dir = Some(tempdir);
        index_sigs = paths;
    } else {
        index_sigs = load_sketchlist_filenames(&siglist)?;
    }

    // if index_sigs pathlist is empty, bail
    if index_sigs.is_empty() {
        bail!("No signatures to index loaded, exiting.");
    }

    info!("Loaded {} sig paths in siglist", index_sigs.len());
    println!("Loaded {} sig paths in siglist", index_sigs.len());

    // Create or open the RevIndex database with the provided output path and colors flag
    let db = RevIndex::create(output.as_ref(), colors);

    // Index the signatures using the loaded template, threshold, and save_paths option
    db.index(index_sigs, &template, 0.0, save_paths);

    if let Some(temp_dir) = temp_dir {
        temp_dir.close()?;
    }

    Ok(())
}


fn check<P: AsRef<Path>>(index: P, quick: bool) -> Result<(), Box<dyn std::error::Error>> {

    if !is_revindex_database(index.as_ref()) {
        bail!("'{}' is not a valid RevIndex database", index.as_ref().display());
    }

    info!("Opening DB");
    let db = RevIndex::open(index.as_ref(), true);

    info!("Starting check");
    db.check(quick);

    info!("Finished check");
    Ok(())
}


fn mastiff_manysearch<P: AsRef<Path>>(
    queries_file: P,
    index: P,
    template: Sketch,
    minimum_containment: f64,
    output: Option<P>,
    ) -> Result<(), Box<dyn std::error::Error>> {

    if !is_revindex_database(index.as_ref()) {
        bail!("'{}' is not a valid RevIndex database", index.as_ref().display());
    }
    // Open database once
    let db = RevIndex::open(index.as_ref(), true);
    info!("Loaded DB");

    // Load query paths
    let query_paths = load_sketchlist_filenames(&queries_file)?;

    // if query_paths is empty, exit with error
    if query_paths.is_empty() {
        bail!("No query signatures loaded, exiting.");
    }

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query_name,query_md5,match_name,containment,intersect_hashes").unwrap();
        for (query, query_md5, m, cont, overlap) in recv.into_iter() { //m_md5 is missing rn
            writeln!(&mut writer, "\"{}\",{},\"{}\",{},{}",
                        query, query_md5, m, cont, overlap).ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send = query_paths
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load query signature from path:
            if let Ok(query_sig) = Signature::from_path(filename) {
                if let Some(query) = prepare_query(&query_sig, &template) {
                    let query_size = query.minhash.size() as f64;
                    // search mastiff db
                    let counter = db.counter_for_query(&query.minhash);
                    let matches = db.matches_from_counter(counter, minimum_containment as usize);

                    // filter the matches for containment

                    for (path, overlap) in matches {
                        let containment = overlap as f64 / query_size;
                        if containment >= minimum_containment {
                            results.push((query.name.clone(),
                                          query.md5sum.clone(),
                                          path.clone(),
                                          containment,
                                          overlap));
                        }
                    }
                } else {
                    eprintln!("WARNING: no compatible sketches in path '{}'",
                              filename.display());
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            } else {
                let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                eprintln!("WARNING: could not load sketches from path '{}'",
                          filename.display());
            }

            if results.is_empty() {
                None
            } else {
                Some(results)
            }
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    // do some cleanup and error handling -
    if let Err(e) = send {
        error!("Unable to send internal data: {:?}", e);
    }

    if let Err(e) = thrd.join() {
        error!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} search sigs", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} paths - no compatible signatures.",
                  skipped_paths);
    }
    if failed_paths > 0 {
        eprintln!("WARNING: {} signature paths failed to load. See error messages above.",
                  failed_paths);
    }

    Ok(())
}


fn mastiff_manygather<P: AsRef<Path>>(
    queries_file: P,
    index: P,
    template: Sketch,
    threshold_bp: usize,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(index.as_ref()) {
        bail!("'{}' is not a valid RevIndex database", index.as_ref().display());
    }
    // Open database once
    let db = RevIndex::open(index.as_ref(), true);
    info!("Loaded DB");

    // Load query paths
    let query_paths = load_sketchlist_filenames(&queries_file)?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query_name,query_md5,match_name,match_md5,f_match_query,intersect_bp").unwrap();
        for (query, query_md5, m, m_md5, f_match_query, intersect_bp) in recv.into_iter() {
            writeln!(&mut writer, "\"{}\",{},\"{}\",{},{},{}",
                        query, query_md5, m, m_md5, f_match_query, intersect_bp).ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send = query_paths
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load query signature from path:
            if let Ok(query_sig) = Signature::from_path(filename) {
                if let Some(query) = prepare_query(&query_sig, &template) {
                    // let query_size = query.minhash.size() as f64;
                    let threshold = threshold_bp / query.minhash.scaled() as usize;
 
                    // mastiff gather code
                    info!("Building counter");
                    let (counter, query_colors, hash_to_color) = db.prepare_gather_counters(&query.minhash);
                    info!("Counter built");

                    let matches = db.gather(
                        counter,
                        query_colors,
                        hash_to_color,
                        threshold,
                        &query.minhash,
                        &template,
                    );

                    // extract matches from Result
                    if let Ok(matches) = matches {
                        info!("matches: {}", matches.len());
                        for match_ in &matches {
                            results.push((query.name.clone(),
                                      query.md5sum.clone(),
                                      match_.name().clone(),
                                      match_.md5().clone(),
                                      match_.f_match(), // f_match_query
                                      match_.intersect_bp())); // intersect_bp
                        }
                    } else {
                        eprintln!("Error gathering matches: {:?}", matches.err());
                    }

                } else {
                    eprintln!("WARNING: no compatible sketches in path '{}'",
                              filename.display());
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            } else {
                let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                eprintln!("WARNING: could not load sketches from path '{}'",
                          filename.display());
            }

            if results.is_empty() {
                None
            } else {
                Some(results)
            }
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    // do some cleanup and error handling -
    if let Err(e) = send {
        error!("Unable to send internal data: {:?}", e);
    }

    if let Err(e) = thrd.join() {
        error!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} search sigs", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} query paths - no compatible signatures.",
                  skipped_paths);
    }
    if failed_paths > 0 {
        eprintln!("WARNING: {} signature paths failed to load. See error messages above.",
                  failed_paths);
    }

    Ok(())
}


//
// The python wrapper functions.
//

#[pyfunction]
fn do_manysearch(querylist_path: String,
                 siglist_path: String,
                 threshold: f64,
                 ksize: u8,
                 scaled: usize,
                 output_path: Option<String>,
) -> anyhow::Result<u8> {
    // if siglist_path is revindex, run mastiff_manysearch; otherwise run manysearch
    if is_revindex_database(siglist_path.as_ref()) {
        let template = build_template(ksize, scaled);
        match mastiff_manysearch(querylist_path, siglist_path, template, threshold, output_path) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match manysearch(querylist_path, siglist_path, threshold, ksize, scaled,
                     output_path) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
fn do_countergather(query_filename: String,
                    siglist_path: String,
                    threshold_bp: usize,
                    ksize: u8,
                    scaled: usize,
                    output_path_prefetch: Option<String>,
                    output_path_gather: Option<String>,
) -> anyhow::Result<u8> {
    match countergather(query_filename, siglist_path, threshold_bp,
                        ksize, scaled,
                        output_path_prefetch,
                        output_path_gather) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_multigather(query_filenames: String,
                     siglist_path: String,
                     threshold_bp: usize,
                     ksize: u8,
                     scaled: usize,
                     output_path: Option<String>,
) -> anyhow::Result<u8> {
    // if a siglist path is a revindex, run mastiff_manygather. If not, run multigather
    if is_revindex_database(siglist_path.as_ref()) {
        let template = build_template(ksize, scaled);
        match mastiff_manygather(query_filenames, siglist_path, template, threshold_bp, output_path) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match multigather(query_filenames, siglist_path, threshold_bp, ksize, scaled) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
fn set_global_thread_pool(num_threads: usize) -> PyResult<usize> {
    if std::panic::catch_unwind(||
        rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()
    ).is_ok() {
        Ok(rayon::current_num_threads())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Could not set the number of threads. Global thread pool might already be initialized."))
    }
}

#[pyfunction]
fn do_index(siglist: String,
            ksize: u8,
            scaled: usize,
            output: String,
            save_paths: bool,
            colors: bool,
) -> anyhow::Result<u8>{
    // build template from ksize, scaled
    let template = build_template(ksize, scaled);
    match index(siglist, template, output,
                save_paths, colors) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
            }
        }
    }

#[pyfunction]
fn do_check(index: String,
            quick: bool,
    ) -> anyhow::Result<u8>{
    match check(index, quick) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
            }
        }
    }


#[pymodule]
fn pyo3_branchwater(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_manysearch, m)?)?;
    m.add_function(wrap_pyfunction!(do_countergather, m)?)?;
    m.add_function(wrap_pyfunction!(do_multigather, m)?)?;
    m.add_function(wrap_pyfunction!(do_index, m)?)?;
    m.add_function(wrap_pyfunction!(do_check, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    Ok(())
}
