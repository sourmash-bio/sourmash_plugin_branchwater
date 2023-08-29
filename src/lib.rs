/// Rust code for pyo3_branchwater.

use pyo3::prelude::*;

use rayon::prelude::*;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use std::cmp::{PartialOrd, Ordering};

use anyhow::{Result, anyhow};

#[macro_use]
extern crate simple_error;

use log::{error, info};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

/// Track a name/minhash.

struct SmallSignature {
    name: String,
    md5sum: String,
    minhash: KmerMinHash,
}

/// Structure to hold overlap information from comparisons.

struct PrefetchResult {
    name: String,
    md5sum: String,
    minhash: KmerMinHash,
    overlap: u64,
}

impl Ord for PrefetchResult {
    fn cmp(&self, other: &PrefetchResult) -> Ordering {
        self.overlap.cmp(&other.overlap)
    }
}

impl PartialOrd for PrefetchResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for PrefetchResult {
    fn eq(&self, other: &Self) -> bool {
        self.overlap == other.overlap
    }
}

impl Eq for PrefetchResult {}

/// check to see if two KmerMinHash are compatible.
///
/// CTB note: despite the name, downsampling is not performed?
/// Although it checks if they are compatible in one direction...

fn check_compatible_downsample(
    me: &KmerMinHash,
    other: &KmerMinHash,
) -> Result<(), sourmash::Error> {
    /* // ignore num minhashes.
    if self.num != other.num {
        return Err(Error::MismatchNum {
            n1: self.num,
            n2: other.num,
        }
        .into());
    }
    */
    use sourmash::Error;

    if me.ksize() != other.ksize() {
        return Err(Error::MismatchKSizes);
    }
    if me.hash_function() != other.hash_function() {
        // TODO: fix this error
        return Err(Error::MismatchDNAProt);
    }
    if me.max_hash() < other.max_hash() {
        return Err(Error::MismatchScaled);
    }
    if me.seed() != other.seed() {
        return Err(Error::MismatchSeed);
    }
    Ok(())
}


/// Given a vec of search Signatures, each containing one or more sketches,
/// and a template Sketch, return a compatible (& now downsampled)
/// Sketch from the search Signatures..
///
/// CTB note: this will return the first acceptable match, I think, ignoring
/// all others.


fn prepare_query(search_sigs: &[Signature], template: &Sketch) -> Option<SmallSignature> {

    for search_sig in search_sigs.iter() {
        // find exact match for template?
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
            return Some(SmallSignature {
                name: search_sig.name(),
                md5sum: mh.md5sum(),
                minhash: mh.clone()
            });
        } else {
            // no - try to find one that can be downsampled
            if let Sketch::MinHash(template_mh) = template {
                for sketch in search_sig.sketches() {
                    if let Sketch::MinHash(ref_mh) = sketch {
                        if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                            let max_hash = max_hash_for_scaled(template_mh.scaled());
                            let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                            return Some(SmallSignature {
                                name: search_sig.name(),
                                md5sum: ref_mh.md5sum(), // original
                                minhash: mh,             // downsampled
                            });
                        }
                    }
                }
            }
        }
    }
    None
}

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

fn manysearch<P: AsRef<Path>>(
    querylist: P,
    siglist: P,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    output: Option<P>,
) -> Result<()> {
    // construct a MinHash template for loading.
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    // Read in list of query paths.
    eprintln!("Reading list of queries from: '{}'", querylist.as_ref().display());

    // Load all queries into memory at once.
    let querylist_paths = load_sketchlist_filenames(&querylist)?;

    let result = load_sketches(querylist_paths, &template)?;
    let (queries, skipped_paths, failed_paths) = result;

    eprintln!("Loaded {} query signatures", queries.len());
    if failed_paths > 0 {
        eprintln!("WARNING: {} signature paths failed to load. See error messages above.",
                  failed_paths);
    }
    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} paths - no compatible signatures.",
                  skipped_paths);
    }

    if queries.is_empty() {
        bail!("No query signatures loaded, exiting.");
    }

    // Load all _paths_, not signatures, into memory.
    eprintln!("Reading search file paths from: '{}'", siglist.as_ref().display());

    let search_sigs_paths = load_sketchlist_filenames(&siglist)?;
    if search_sigs_paths.is_empty() {
        bail!("No signatures to search loaded, exiting.");
    }

    eprintln!("Loaded {} sig paths to search.", search_sigs_paths.len());

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query_name,query_md5,match_name,match_md5,containment,intersect_hashes").unwrap();
        for (query, query_md5, m, m_md5, cont, overlap) in recv.into_iter() {
            writeln!(&mut writer, "\"{}\",{},\"{}\",{},{},{}",
                     query, query_md5, m, m_md5, cont, overlap).ok();
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

    let send = search_sigs_paths
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load search signature from path:
            if let Ok(search_sigs) = Signature::from_path(filename) {
                if let Some(search_sm) = prepare_query(&search_sigs, &template) {
                    // search for matches & save containment.
                    for q in queries.iter() {
                        let overlap = q.minhash.count_common(&search_sm.minhash, false).unwrap() as f64;
                        let query_size = q.minhash.size() as f64;

                        let containment = overlap / query_size;
                        if containment > threshold {
                            results.push((q.name.clone(),
                                          q.md5sum.clone(),
                                          search_sm.name.clone(),
                                          search_sm.md5sum.clone(),
                                          containment,
                                          overlap))
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

/// Find sketches in 'sketchlist' that overlap with 'query' above
/// specified threshold.

fn prefetch(
    query: &KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
) -> BinaryHeap<PrefetchResult> {
    sketchlist
        .into_par_iter()
        .filter_map(|result| {
            let mut mm = None;
            let searchsig = &result.minhash;
            let overlap = searchsig.count_common(query, false);
            if let Ok(overlap) = overlap {
                if overlap >= threshold_hashes {
                    let result = PrefetchResult {
                        overlap,
                        ..result
                    };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

/// Write list of prefetch matches.

fn write_prefetch<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query_label: String,
    prefetch_output: Option<P>,
    matchlist: &BinaryHeap<PrefetchResult>
) -> Result<()> {
    // Set up a writer for prefetch output
    let prefetch_out: Box<dyn Write> = match prefetch_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(prefetch_out);
    writeln!(&mut writer, "query_filename,match_name,match_md5,intersect_bp").ok();

    for m in matchlist.iter() {
        writeln!(&mut writer, "{},\"{}\",{},{}", query_label,
                 m.name, m.md5sum, m.overlap).ok();
    }

    Ok(())
}

/// Load a list of filenames from a file. Exits on bad lines.

fn load_sketchlist_filenames<P: AsRef<Path>>(sketchlist_filename: &P) ->
    Result<Vec<PathBuf>>
{
    let sketchlist_file = BufReader::new(File::open(sketchlist_filename)?);

    let mut sketchlist_filenames : Vec<PathBuf> = Vec::new();
    for line in sketchlist_file.lines() {
        let line = match line {
            Ok(v) => v,
            Err(_) => return {
                let filename = sketchlist_filename.as_ref().display();
                let msg = format!("invalid line in fromfile '{}'", filename);
                Err(anyhow!(msg))
            },
        };

        if !line.is_empty() {
            let mut path = PathBuf::new();
            path.push(line);
            sketchlist_filenames.push(path);
        }
    }
    Ok(sketchlist_filenames)
}

/// Load a collection of sketches from a file in parallel.

fn load_sketches(sketchlist_paths: Vec<PathBuf>, template: &Sketch) ->
    Result<(Vec<SmallSignature>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let sketchlist : Vec<SmallSignature> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let mut sm = None;

            let filename = m.display();

            if let Ok(sigs) = Signature::from_path(m) {
                sm = prepare_query(&sigs, template);
                if sm.is_none() {
                    // track number of paths that have no matching sigs
                    let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            } else {
                // failed to load from this path - print error & track.
                eprintln!("WARNING: could not load sketches from path '{}'",
                          filename);
                let _i = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
            sm
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);
    Ok((sketchlist, skipped_paths, failed_paths))
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

fn load_sketches_above_threshold(
    sketchlist_paths: Vec<PathBuf>,
    template: &Sketch,
    query: &KmerMinHash,
    threshold_hashes: u64
) ->
    Result<(BinaryHeap<PrefetchResult>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let matchlist: BinaryHeap<PrefetchResult> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m);

            match sigs {
                Ok(sigs) => {
                    let mut mm = None;

                    if let Some(sm) = prepare_query(&sigs, template) {
                        let mh = sm.minhash;
                        if let Ok(overlap) = mh.count_common(query, false) {
                            if overlap >= threshold_hashes {
                                let result = PrefetchResult {
                                    name: sm.name,
                                    md5sum: sm.md5sum,
                                    minhash: mh,
                                    overlap,
                                };
                                mm = Some(result);
                            }
                        }
                    } else {
                        eprintln!("WARNING: no compatible sketches in path '{}'",
                                  m.display());
                        let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    mm
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("WARNING: could not load sketches from path '{}'",
                          m.display());
                    None
                }
            }
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((matchlist, skipped_paths, failed_paths))
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

fn consume_query_by_gather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    mut query: KmerMinHash,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<P>,
    query_label: String
) -> Result<()> {
    // Set up a writer for gather output
    let gather_out: Box<dyn Write> = match gather_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(gather_out);
    writeln!(&mut writer, "query_filename,rank,match_name,match_md5,intersect_bp").ok();

    let mut matching_sketches = matchlist;
    let mut rank = 0;

    while !matching_sketches.is_empty() {
        eprintln!("{} remaining: {} {}", query_label,
                  query.size(), matching_sketches.len());
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query.remove_from(&best_element.minhash)?;

        writeln!(&mut writer, "{},{},\"{}\",{},{}", query_label, rank,
                 best_element.name, best_element.md5sum,
                 best_element.overlap).ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query, matching_sketches, threshold_hashes);
        rank += 1;
    }
    Ok(())
}
                           

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

    write_prefetch(query_label.clone(), prefetch_output, &matchlist).ok();

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

//
// The python wrapper functions.
//

#[pyfunction]
fn do_manysearch(querylist_path: String,
                 siglist_path: String,
                 threshold: f64,
                 ksize: u8,
                 scaled: usize,
                 output_path: String
) -> anyhow::Result<u8> {
    match manysearch(querylist_path, siglist_path, threshold, ksize, scaled,
                     Some(output_path)) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
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
                     scaled: usize
) -> anyhow::Result<u8> {
    match multigather(query_filenames, siglist_path, threshold_bp,
                         ksize, scaled) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn set_global_thread_pool(num_threads: usize) -> PyResult<usize> {
    if let Ok(_) = std::panic::catch_unwind(|| rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()) {
        Ok(rayon::current_num_threads())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("Could not set the number of threads. Global thread pool might already be initialized."))
    }
}

#[pymodule]
fn pyo3_branchwater(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_manysearch, m)?)?;
    m.add_function(wrap_pyfunction!(do_countergather, m)?)?;
    m.add_function(wrap_pyfunction!(do_multigather, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    Ok(())
}
