// TODO:
// * md5sum output by search and countergather are of modified/downsampled,
//   not orig. This is different from sourmash...

use pyo3::prelude::*;
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use std::cmp::{PartialOrd, Ordering};

#[macro_use]
extern crate simple_error;

use log::{error, info};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

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

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
///
/// TODO:
///   - support jaccard as well as containment/overlap
///   - support md5 output columns; other?

fn search<P: AsRef<Path>>(
    querylist: P,
    siglist: P,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    // construct a MinHash template for loading.
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    // Read in list of query paths.
    eprintln!("Reading querylist from: {}", querylist.as_ref().display());
    let querylist_file = BufReader::new(File::open(querylist)?);

    // Load all queries into memory at once.
    let queries: Vec<(String, KmerMinHash)> = querylist_file
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if !line.is_empty() {
                // skip empty lines; load non-empty!
                let mut path = PathBuf::new();
                path.push(line);
                Some(path)
            } else {
                None
            }
        })
        // for non-empty paths, load whichever one matches template.
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
        bail!("No query signatures loaded, exiting.");
    }

    eprintln!("Loaded {} query signatures", queries.len());

    // Load all _paths_, not signatures, into memory.
    eprintln!("Reading search file paths from: {}", siglist.as_ref().display());

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
    if search_sigs.is_empty() {
        bail!("No signatures to search loaded, exiting.");
    }

    eprintln!("Loaded {} sig paths to search.", search_sigs.len());

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query,quer5_md5,match,match_md5,containment").unwrap();
        for (query, query_md5, m, m_md5, overlap) in recv.into_iter() {
            writeln!(&mut writer, "'{}',{},'{}',{},{}",
                     query, query_md5, m, m_md5, overlap).ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);

    let send = search_sigs
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} search sigs", i);
            }

            let mut search_mh = None;
            // load search signature from path:
            let search_sig = &Signature::from_path(&filename)
                .unwrap_or_else(|_| panic!("Error processing {:?}", filename))[0];

            // make sure it is compatible etc.
            if let Some(mh) = prepare_query(search_sig, &template) {
                search_mh = Some(mh);
            }
            // (this will raise an exception if nothing compatible.)
            let search_mh = search_mh.unwrap();

            let mut results = vec![];

            // search for matches & save containment.
            for (name, query) in &queries {
                let overlap =
                    query.count_common(&search_mh, false).unwrap() as f64 / query.size() as f64;
                if overlap > threshold {
                    results.push((name.clone(),
                                  query.md5sum(),
                                  search_sig.name(),
                                  search_sig.md5sum(),
                                  overlap))
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

fn countergather<P: AsRef<Path> + std::fmt::Debug>(
    query_filename: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    gather_output: Option<P>,
    prefetch_output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    eprintln!("Loading query from '{}'", query_filename.as_ref().display());
    let query = {
        let mut mm = None;
        let sigs = Signature::from_path(query_filename)?;

        for sig in &sigs {
            if let Some(mh) = prepare_query(sig, &template) {
                mm = Some(mh.clone());
                break;
            }
        };
        mm
    };

    // did we find one?
    let mut query = match query {
        Some(query) => query,
        None => bail!("No sketch found with scaled={}, k={}", scaled, ksize),
    };

    // build the list of paths to match against.
    eprintln!("Loading matchlist from '{}'", matchlist_filename.as_ref().display());
    let matchlist_paths = load_sketchlist_filenames(matchlist_filename)?;

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
    let matchlist = load_sketches_above_threshold(matchlist_paths,
                                                  &template,
                                                  &query,
                                                  threshold_hashes)?;

    if matchlist.is_empty() {
        eprintln!("No matchlist signatures loaded, exiting.");
        return Ok(());
    }

    // Write all the prefetch matches to prefetch output.
    let prefetch_out: Box<dyn Write> = match prefetch_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(Vec::new()),
    };
    let mut writer = BufWriter::new(prefetch_out);
    writeln!(&mut writer, "match,md5sum,overlap").unwrap();
    for m in &matchlist {
        writeln!(&mut writer, "'{}',{},{}", m.name, m.minhash.md5sum(), m.overlap).ok();
    }
    writer.flush().ok();
    drop(writer);

    // Set up a writer for gather output
    let gather_out: Box<dyn Write> = match gather_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(gather_out);
    writeln!(&mut writer, "rank,match,md5sum,overlap").ok();

    //
    // main loop: loop until no more matching sketches.
    //
    let mut matching_sketches = matchlist;
    let mut rank = 0;

    while !matching_sketches.is_empty() {
        eprintln!("remaining: {} {}", query.size(), matching_sketches.len());
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query.remove_from(&best_element.minhash)?;

        writeln!(&mut writer, "{},'{}',{},{}", rank, best_element.name, best_element.minhash.md5sum(), best_element.overlap).ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query, matching_sketches, threshold_hashes);
        rank = rank + 1;
    }

    Ok(())
}

#[pyfunction]
fn do_search(querylist_path: String,
             siglist_path: String,
             threshold: f64,
             ksize: u8,
             scaled: usize,
             output_path: String
) -> PyResult<u8> {
    match search(querylist_path, siglist_path, threshold, ksize, scaled,
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
) -> PyResult<u8> {
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

#[pymodule]
fn pyo3_branchwater(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_search, m)?)?;
    m.add_function(wrap_pyfunction!(do_countergather, m)?)?;
    Ok(())
}
