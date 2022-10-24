use pyo3::prelude::*;
use pyo3::create_exception;
use pyo3::exceptions::PyException;

use rayon::prelude::*;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use std::cmp::{PartialOrd, Ordering};

use log::{error, info};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

create_exception!(pymagsearch, SomeError, pyo3::exceptions::PyException);

impl std::convert::From<SomeError> for PyErr {
    fn from(err: SomeError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}

fn check_compatible_downsample(
    me: &KmerMinHash,
    other: &KmerMinHash,
) -> Result<(), sourmash::Error> {
    /*
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

fn prepare_query(search_sig: &Signature, template: &Sketch) -> Option<KmerMinHash> {
    let mut search_mh = None;
    if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
        search_mh = Some(mh.clone());
    } else {
        // try to find one that can be downsampled
        if let Sketch::MinHash(template_mh) = template {
            for sketch in search_sig.sketches() {
                if let Sketch::MinHash(ref_mh) = sketch {
                    if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                        let max_hash = max_hash_for_scaled(template_mh.scaled());
                        let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                        return Some(mh);
                    }
                }
            }
        }
    }
    search_mh
}

fn search<P: AsRef<Path>>(
    querylist: P,
    siglist: P,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Loading queries");

    let querylist_file = BufReader::new(File::open(querylist)?);

    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    let queries: Vec<(String, KmerMinHash)> = querylist_file
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if !line.is_empty() {
                // skip empty lines
                let mut path = PathBuf::new();
                path.push(line);
                Some(path)
            } else {
                None
            }
        })
        .filter_map(|query| {
            let query_sig = Signature::from_path(query).unwrap();

            let mut query = None;
            for sig in &query_sig {
                if let Some(mh) = prepare_query(sig, &template) {
                    query = Some((sig.name(), mh));
                    break;
                }
            }
            query
        })
        .collect();

    if queries.is_empty() {
        info!("No query signatures loaded, exiting.");
        return Ok(());
    }

    info!("Loaded {} query signatures", queries.len());

    info!("Loading siglist");
    let siglist_file = BufReader::new(File::open(siglist)?);
    let search_sigs: Vec<PathBuf> = siglist_file
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if !line.is_empty() {
                let mut path = PathBuf::new();
                path.push(line);
                Some(path)
            } else {
                None
            }
        })
        .collect();
    info!("Loaded {} sig paths in siglist", search_sigs.len());

    let processed_sigs = AtomicUsize::new(0);

    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // Spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query,Run,containment").unwrap();
        for (query, m, overlap) in recv.into_iter() {
            writeln!(&mut writer, "'{}','{}',{}", query, m, overlap).unwrap();
        }
    });

    let send = search_sigs
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} search sigs", i);
            }

            let mut search_mh = None;
            let search_sig = &Signature::from_path(&filename)
                .unwrap_or_else(|_| panic!("Error processing {:?}", filename))[0];

            if let Some(mh) = prepare_query(search_sig, &template) {
                search_mh = Some(mh);
            }
            let search_mh = search_mh.unwrap();

            let match_fn = filename.clone().into_os_string().into_string().unwrap();
            let mut results = vec![];

            for (name, query) in &queries {
                let overlap =
                    query.count_common(&search_mh, false).unwrap() as f64 / query.size() as f64;
                if overlap > threshold {
                    results.push((name.clone(), match_fn.clone(), overlap))
                }
            }
            if results.is_empty() {
                None
            } else {
                Some(results)
            }
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    if let Err(e) = send {
        error!("Unable to send internal data: {:?}", e);
    }

    if let Err(e) = thrd.join() {
        error!("Unable to join internal thread: {:?}", e);
    }

    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    info!("DONE. Processed {} search sigs", i);
    Ok(())
}

struct SmallSignature {
    name: String,
    minhash: KmerMinHash,
}

struct PrefetchResult {
    name: String,
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

// Find overlapping above specified threshold.

fn filter_overlapping(
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

/// Loads a list of filenames containing sketches from sketchlist_file.

fn load_sketchlist_filenames<P: AsRef<Path>>(sketchlist_file: P) ->
    Result<Vec<PathBuf>, Box<dyn std::error::Error>>
{
    let sketchlist_file = BufReader::new(File::open(sketchlist_file)?);
    let sketchlist_filenames = sketchlist_file
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if !line.is_empty() {
                // skip empty lines
                let mut path = PathBuf::new();
                path.push(line);
                Some(path)
            } else {
                None
            }
        })
        .collect();
    Ok(sketchlist_filenames)
}

/// Load a collection of sketches from a file.

fn load_sketches(sketchlist_paths: Vec<PathBuf>, template: &Sketch) ->
    Result<Vec<SmallSignature>, Box<dyn std::error::Error>>
{
    let sketchlist : Vec<SmallSignature> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m).unwrap();

            let mut sm = None;
            for sig in &sigs {
                if let Some(mh) = prepare_query(sig, &template) {
                    sm = Some(SmallSignature {
                        name: sig.name(),
                        minhash: mh,
                    });
                }
            }
            sm
        })
        .collect();
    Ok(sketchlist)
}

/// Load a collection of sketches from a file, filtering w/query & threshold

fn load_sketches_above_threshold(
    sketchlist_paths: Vec<PathBuf>,
    template: &Sketch,
    query: &KmerMinHash,
    threshold_hashes: u64
) ->
    Result<BinaryHeap<PrefetchResult>, Box<dyn std::error::Error>>
{
    let matchlist: BinaryHeap<PrefetchResult> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m).unwrap();

            let mut mm = None;
            for sig in &sigs {
                if let Some(mh) = prepare_query(sig, &template) {
                    if let Ok(overlap) = mh.count_common(&query, false) {
                        if overlap >= threshold_hashes {
                            let result = PrefetchResult {
                                name: sig.name(),
                                minhash: mh,
                                overlap,
                            };
                            mm = Some(result);
                        }
                    }
                }
            }
            mm
        })
        .collect();

    Ok(matchlist)
}

/// Run counter-gather with a query against a list of files.

fn countergather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query_filename: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    prefetch_output: Option<P>,
    gather_output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    let query_label = query_filename.clone();

    println!("Loading query");
    let query = {
        let sigs = Signature::from_path(dbg!(query_filename)).unwrap();

        let mut mm = None;
        for sig in &sigs {
            if let Some(mh) = prepare_query(sig, &template) {
                mm = Some(mh.clone());
                break;
            }
        }
        mm
    }
    .unwrap();

    // build the list of paths to match against.
    println!("Loading matchlist");
    let matchlist_paths = load_sketchlist_filenames(matchlist_filename).unwrap();
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

    let matchlist = load_sketches_above_threshold(matchlist_paths,
                                                  &template,
                                                  &query,
                                                  threshold_hashes).unwrap();

    if matchlist.is_empty() {
        println!("No matchlist signatures loaded, exiting.");
        return Ok(());
    }

    write_prefetch(query_label.to_string(), prefetch_output, &matchlist)?;

    // now, do the gather!
    consume_query_by_gather(query, matchlist, threshold_hashes, gather_output,
                            query_label.to_string());

    Ok(())
}

fn write_prefetch<P: AsRef<Path> + std::fmt::Debug>(
    query_label: String,
    prefetch_output: Option<P>,
    matches: &BinaryHeap<PrefetchResult>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Write to prefetch output
    let prefetch_out: Box<dyn Write> = match prefetch_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(prefetch_out);
    writeln!(&mut writer, "match,overlap").unwrap();
    for m in matches {
        writeln!(&mut writer, "'{}','{}',{}", query_label, m.name, m.overlap);
    }

    Ok(())
}

fn consume_query_by_gather<P: AsRef<Path> + std::fmt::Debug>(
    query: KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<P>,
    query_label: String,
) -> Result<(), Box<dyn std::error::Error>> {
    // Set up gather output
    let gather_out: Box<dyn Write> = match gather_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(gather_out);
    writeln!(&mut writer, "match,overlap").unwrap();

    let mut query = query;
    let mut sketchlist = sketchlist;

    // loop until no more matching sketches -
    while !sketchlist.is_empty() {
        println!("{}: remaining: {} {}", query_label, query.size(), sketchlist.len());
        let best_element = sketchlist.peek().unwrap();

        // remove!
        println!("{}: removing {}", query_label, best_element.name);
        query.remove_from(&best_element.minhash)?;

        writeln!(&mut writer, "'{}','{}',{}", query_label, best_element.name,
                 best_element.overlap);

        sketchlist = filter_overlapping(&query, sketchlist, threshold_hashes);
    }
    Ok(())
}

/// Run counter-gather for multiple queries against a list of files.

fn countergather2<P: AsRef<Path> + std::fmt::Debug + Clone>(
    query_filenames: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    // load the list of query paths
    let querylist_paths = load_sketchlist_filenames(query_filenames)?;
    println!("Loaded {} sig paths in querylist", querylist_paths.len());

    // build the list of paths to match against.
    println!("Loading matchlist");
    let matchlist_paths = load_sketchlist_filenames(matchlist_filename)?;
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
    let sketchlist = load_sketches(matchlist_paths, &template).unwrap();

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);

    querylist_paths
        .par_iter()
        .for_each(|q| {
            let i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);
            let query_label = q.clone().into_os_string().into_string().unwrap();

            if let Some(query) = {
                // load query from q
                let mut mm = None;
                if let Ok(sigs) = Signature::from_path(dbg!(q)) {
                    for sig in &sigs {
                        if let Some(mh) = prepare_query(sig, &template) {
                            mm = Some(mh);
                            break;
                        }
                    }
                }
                mm
            } {
                // filter first set of matches out of sketchlist
                    let matchlist: BinaryHeap<PrefetchResult> = sketchlist
                    .par_iter()
                    .filter_map(|sm| {
                        let mut mm = None;

                        if let Ok(overlap) = sm.minhash.count_common(&query, false) {
                            if overlap >= threshold_hashes {
                                let result = PrefetchResult {
                                    name: sm.name.clone(),
                                    minhash: sm.minhash.clone(),
                                    overlap,
                                };
                                mm = Some(result);
                            }
                        }
                        mm
                    })
                    .collect();

                if matchlist.len() > 0 {
                    let prefetch_output = format!("prefetch-{i}.csv");
                    let gather_output = format!("gather-{i}.csv");

                    // save initial list of matches to prefetch output
                    write_prefetch(query_label.clone(), Some(prefetch_output),
                                   &matchlist);

                    // now, do the gather!
                    consume_query_by_gather(query, matchlist, threshold_hashes,
                                            Some(gather_output), query_label);
                } else {
                    println!("No matches to '{}'", query_label);
                }
            } else {
                println!("ERROR loading signature from '{}'", query_label);
            }
        });

        
    Ok(())
}

#[pyfunction]
fn do_search(querylist_path: String,
             siglist_path: String,
             threshold: f64,
             ksize: u8,
             scaled: usize,
             output_path: String
) -> PyResult<()> {
    search(querylist_path, siglist_path, threshold, ksize, scaled,
           Some(output_path));
    Ok(())
}

#[pyfunction]
fn do_countergather(query_filename: String,
                    siglist_path: String,
                    threshold_bp: usize,
                    ksize: u8,
                    scaled: usize,
                    output_path_prefetch: String,
                    output_path_gather: String,
) -> PyResult<()> {
    let x = countergather(query_filename, siglist_path, threshold_bp, ksize, scaled,
                  Some(output_path_prefetch), Some(output_path_gather));
    match x {
        Ok(_) => Ok(()),
        Err(error) => Err(SomeError::new_err("fiz")),
    }
}

#[pyfunction]
fn do_countergather2(query_filenames: String,
                    siglist_path: String,
                    threshold_bp: usize,
                    ksize: u8,
                    scaled: usize,
) -> PyResult<()> {
    let x = countergather2(query_filenames, siglist_path, threshold_bp,
                           ksize, scaled);
    match x {
        Ok(_) => Ok(()),
        Err(error) => Err(SomeError::new_err("fiz")),
    }
}

#[pymodule]
fn pymagsearch(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_search, m)?)?;
    m.add_function(wrap_pyfunction!(do_countergather, m)?)?;
    m.add_function(wrap_pyfunction!(do_countergather2, m)?)?;
    m.add("SomeError", py.get_type::<SomeError>())?;
    Ok(())
}
