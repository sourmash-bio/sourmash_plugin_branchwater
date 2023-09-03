/// Rust code for pyo3_branchwater.

use pyo3::prelude::*;

use rayon::prelude::*;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Seek, Stdout};
use std::path::{Path, PathBuf};

use zip::read::ZipArchive;
use tempfile::tempdir;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::io::Read;

use std::collections::BinaryHeap;

use std::cmp::{PartialOrd, Ordering};

use anyhow::{Result, anyhow};

use needletail::parse_fastx_reader;

#[macro_use]
extern crate simple_error;

use log::{error, info};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;
use sourmash::index::revindex::{RevIndex};
use sourmash::prelude::MinHashOps;
use sourmash::prelude::FracMinHashOps;
use sourmash::cmd::ComputeParameters;

// set up logging to print info messages using env logger
use std::sync::Once;
use env_logger;

static LOGGER_INITIALIZED: Once = Once::new();

fn initialize_logger() {
    LOGGER_INITIALIZED.call_once(|| {
        std::env::set_var("RUST_LOG", "info");
        env_logger::init();
    });
}

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
    template: Sketch,
    threshold: f64,
    output: Option<P>,
) -> Result<()> {

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
    let (send, recv) = std::sync::mpsc::sync_channel::<SearchResult>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output.as_ref());

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
                        let target_size = search_sm.minhash.size() as f64;

                        let containment_query_in_target = overlap / query_size;
                        let containment_in_target = overlap / target_size;
                        let max_containment = containment_query_in_target.max(containment_in_target);
                        let jaccard = overlap / (target_size + query_size - overlap);

                        if containment_query_in_target > threshold {
                            results.push(SearchResult {
                                query_name: q.name.clone(),
                                query_md5: q.md5sum.clone(),
                                match_name: search_sm.name.clone(),
                                containment: containment_query_in_target,
                                intersect_hashes: overlap as usize,
                                match_md5: Some(search_sm.md5sum.clone()),
                                jaccard: Some(jaccard),
                                max_containment: Some(max_containment),
                            });
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

    let mut last_hashes = query.size();
    let mut last_matches = matching_sketches.len();

    eprintln!("{} iter {}: start: query hashes={} matches={}", query_label, rank,
            query.size(), matching_sketches.len());

    while !matching_sketches.is_empty() {
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

        let sub_hashes = last_hashes - query.size();
        let sub_matches = last_matches - matching_sketches.len();

        eprintln!("{} iter {}: remaining: query hashes={}(-{}) matches={}(-{})", query_label, rank,
            query.size(), sub_hashes, matching_sketches.len(), sub_matches);

        last_hashes = query.size();
        last_matches = matching_sketches.len();

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

// mastiff rocksdb functions

fn build_template(ksize: u8, scaled: usize) -> Sketch {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    Sketch::MinHash(template_mh)
}

fn read_signatures_from_zip<P: AsRef<Path>>(
    zip_path: P,
) -> Result<(Vec<PathBuf>, tempfile::TempDir), Box<dyn std::error::Error>> {
    let mut signature_paths = Vec::new();
    let temp_dir = tempdir()?;
    let zip_file = File::open(&zip_path)?;
    let mut zip_archive = ZipArchive::new(zip_file)?;

    for i in 0..zip_archive.len() {
        let mut file = zip_archive.by_index(i)?;
        let mut sig = Vec::new();
        file.read_to_end(&mut sig)?;

        let file_name = Path::new(file.name()).file_name().unwrap().to_str().unwrap();
        if file_name.ends_with(".sig") || file_name.ends_with(".sig.gz") {
            info!("Found signature file: {}", file_name);
            let mut new_file = File::create(temp_dir.path().join(file_name))?;
            new_file.write_all(&sig)?;

            // Push the created path directly to the vector
            signature_paths.push(temp_dir.path().join(file_name));
        }
    }
    println!("wrote {} signatures to temp dir", signature_paths.len());
    Ok((signature_paths, temp_dir))
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

fn is_revindex_database(path: &Path) -> bool {
    // quick file check for Revindex database:
    // is path a directory that contains a file named 'CURRENT'?
    if path.is_dir() {
        let current_file = path.join("CURRENT");
        current_file.exists() && current_file.is_file()
    } else {
        false
    }
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

trait ResultType {
    fn header_fields() -> Vec<&'static str>;
    fn format_fields(&self) -> Vec<String>;
}

struct SearchResult {
    query_name: String,
    query_md5: String,
    match_name: String,
    containment: f64,
    intersect_hashes: usize,
    match_md5: Option<String>,
    jaccard: Option<f64>,
    max_containment: Option<f64>,
}

impl ResultType for SearchResult {
    fn header_fields() -> Vec<&'static str> {
        vec!["query_name", "query_md5", "match_name", "containment", "intersect_hashes", "match_md5", "jaccard", "max_containment"]
    }

    fn format_fields(&self) -> Vec<String> {
        vec![
            format!("\"{}\"", self.query_name),  // Wrap query_name with quotes
            self.query_md5.clone(),
            format!("\"{}\"", self.match_name),  // Wrap match_name with quotes
            self.containment.to_string(),
            self.intersect_hashes.to_string(),
            match &self.match_md5 {
                Some(md5) => md5.clone(),
                None => "".to_string(),
            },
            match &self.jaccard {
                Some(jaccard) => jaccard.to_string(),
                None => "".to_string(),
            },
            match &self.max_containment {
                Some(max_containment) => max_containment.to_string(),
                None => "".to_string(),
            }
        ]
    }
}


enum WriterWrapper {
    File(BufWriter<File>),
    Stdout(Stdout),
}

impl Write for WriterWrapper {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            WriterWrapper::File(writer) => writer.write(buf),
            WriterWrapper::Stdout(stdout) => stdout.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            WriterWrapper::File(writer) => writer.flush(),
            WriterWrapper::Stdout(stdout) => stdout.flush(),
        }
    }
}

impl Seek for WriterWrapper {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match self {
            WriterWrapper::File(writer) => writer.seek(pos),
            WriterWrapper::Stdout(_) => Err(std::io::Error::new(std::io::ErrorKind::Other, "Cannot seek on stdout")),
        }
    }
}

fn open_output_file<P: AsRef<Path>>(output: Option<P>) -> WriterWrapper {
    match output {
        Some(path) => {
            let file = File::create(&path).unwrap_or_else(|e| {
                eprintln!("Error creating output file: {:?}", e);
                std::process::exit(1);
            });
            info!("Writing to file: {}", path.as_ref().display());
            WriterWrapper::File(BufWriter::new(file))
        }
        None => WriterWrapper::Stdout(std::io::stdout()),
    }
}



fn csvwriter_thread<T: ResultType + Send + 'static, P>(
    recv: std::sync::mpsc::Receiver<T>,
    output: Option<P>,
) -> std::thread::JoinHandle<()>
where
    T: ResultType,
    P: Clone + std::convert::AsRef<std::path::Path>,
{
    // create output file
    let out = open_output_file(output);
    // spawn a thread that is dedicated to printing to a buffered output
    std::thread::spawn(move || {
        let mut writer = out;

        let header = T::header_fields();
        if let Err(e) = writeln!(&mut writer, "{}", header.join(",")) {
            error!("Error writing header: {:?}", e);
        }
        writer.flush().unwrap();

        for item in recv.iter() {
            let formatted_fields = item.format_fields();
            if let Err(e) = writeln!(&mut writer, "{}", formatted_fields.join(",")) {
                error!("Error writing item: {:?}", e);
            }
            writer.flush().unwrap();
        }
    })
}


fn sigwriter(
    recv: std::sync::mpsc::Receiver<Signature>,
    output: Option<Box<std::path::PathBuf>>,
) -> std::thread::JoinHandle<()> 
{
    std::thread::spawn(move || {
        if let Some(out_path) = &output {
            if out_path.extension() == Some(std::ffi::OsStr::new("zip")) {
                // ZipWriter 
                let out = open_output_file(Some(&**out_path));
                let options = zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Stored);
                info!("Creating ZipWriter");
                let mut zip = zip::ZipWriter::new(out);
                for sig in recv {
                    write_signature(&sig, &mut Some(&mut zip), Some(options));
                }
            } else {
                // only zip for now
                error!("Output path must be a zip file.");
            }
        } else {
            for sig in recv {
                write_signature(&sig, &mut None::<&mut _>, None);
            }
        }
    })
}


fn write_signature(
    sig: &Signature,
    zip: &mut Option<&mut zip::ZipWriter<WriterWrapper>>,
    zip_options: Option<zip::write::FileOptions>,
) {
    let wrapped_sig = vec![sig];
    let json_bytes = serde_json::to_vec(&wrapped_sig).unwrap();
    info!("Serialized signature to JSON.");
    info!("JSON: {}", String::from_utf8(json_bytes.clone()).unwrap());
    info!("Gzipping signature.");
    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            ).unwrap();
            gz_writer.write_all(&json_bytes).unwrap();
        }
        buffer.into_inner()
    };
    let sig_filename = format!("{}.sig.gz", sig.md5sum());
    match zip {
        Some(zip) => {
            zip.start_file(sig_filename, zip_options.unwrap()).unwrap();
            zip.write_all(&gzipped_buffer).unwrap();
        }
        None => {
            // write gzipped signature to sig_filename
            std::fs::write(&sig_filename, &gzipped_buffer).unwrap();
        }
    }
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
    let (send, recv) = std::sync::mpsc::sync_channel::<SearchResult>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output.as_ref());

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send_result = query_paths
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
                            results.push( SearchResult {
                                query_name: query.name.clone(),
                                query_md5: query.md5sum.clone(),
                                match_name: path.clone(),
                                containment: containment,
                                intersect_hashes: overlap,
                                match_md5: None,
                                jaccard: None,
                                max_containment: None,
                            });
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
          .try_for_each_with(send, |s, results| {
            if let Err(e) = s.send(results) {
                Err(format!("Unable to send internal data: {:?}", e))
            } else {
                Ok(())
            }
        });

    // do some cleanup and error handling -
    if let Err(e) = send_result {
        error!("Error during parallel processing: {}", e);
    }

    // join the writer thread
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


fn manysketch<P: AsRef<Path> + Sync>(
    filelist: P,
    ksize: u8,
    scaled: usize,
    //output as pathbuf
    output: Option<Box<std::path::PathBuf>>,
) -> Result<(), Box<dyn std::error::Error>> {

    // load list of file paths (todo: modify for fasta files)
    let filelist_paths = load_sketchlist_filenames(&filelist)?;

    // if filelist_paths is empty, exit with error
    if filelist_paths.is_empty() {
        bail!("No files to load, exiting.");
    }

    // set up a multi-producer, single-consumer channel that receives Signature
    let (send, recv) = std::sync::mpsc::sync_channel::<Signature>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = sigwriter(recv, output);

    // let thrd = sigwriter(recv, output.as_ref(), false); // set gzip false for now

    // iterate over filelist_paths
    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send_result = filelist_paths
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} fasta files", i);
            }

            // build signature from params
            // future: allow multiple ksizes, moltype, etc
            let cp = ComputeParameters::builder()
            .ksizes(vec![ksize as u32])
            .scaled(scaled as u64)
            .protein(false)
            .dna(true)
            .num_hashes(0)
            .build();

            // print filename so we know what we're working with
            println!("filename: {}", filename.display());
            let mut sig = Signature::from_params(&cp);
            info!("sig initialized");
            let mut data: Vec<u8> = vec![];
            let mut f = File::open(filename).unwrap();
            info!("file opened");
            let _ = f.read_to_end(&mut data);
            info!("file read");

            let mut parser = parse_fastx_reader(&data[..]).unwrap();
            info!("fastx parser created");
            // iterate over records in fasta file
            // todo - after debugging, remove record_count
            let mut record_count = 0;
            while let Some(record_result) = parser.next() {
                match record_result {
                    Ok(record) => {
                        // add sequence to signature
                        sig.add_sequence(&record.seq(), true).unwrap();
                        record_count += 1;
                    }
                    Err(error) => {
                        // Handle the error, if needed
                        println!("Error while processing record: {:?}", error);
                    }
                }
            }
            info!("Processed {} records", record_count);
            info!("sig size: {}", sig.size());
            Some(sig)

        })
        .try_for_each_with(send, |s, sig| s.send(sig));

    // do some cleanup and error handling -
    if let Err(e) = send_result {
        error!("Error during parallel processing: {}", e);
    }

    // join the writer thread
    if let Err(e) = thrd.join() {
        error!("Unable to join internal thread: {:?}", e);
    }

    // todo: collect signature info above and then write to csv manifest in the zip?

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} fasta files", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if skipped_paths > 0 {
        eprintln!("WARNING: skipped {} fasta files.",
                  skipped_paths);
    }

    if failed_paths > 0 {
        eprintln!("WARNING: {} fasta files failed to load. See error messages above.",
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
    let template = build_template(ksize, scaled);
    if is_revindex_database(siglist_path.as_ref()) {
        match mastiff_manysearch(querylist_path, siglist_path, template, threshold, output_path) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match manysearch(querylist_path, siglist_path, template, threshold,
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

#[pyfunction]
fn do_manysketch(filelist: String,
                 ksize: u8,
                 scaled: usize,
                 output: Option<String>,
    ) -> anyhow::Result<u8>{
    initialize_logger(); // initialize logger
    // make output a pathbuf if it exists
    let output = match output {
        Some(path) => Some(Box::new(std::path::PathBuf::from(path))),
        None => None,
    };
    match manysketch(filelist, ksize, scaled, output) {
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
    m.add_function(wrap_pyfunction!(do_manysketch, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    Ok(())
}
