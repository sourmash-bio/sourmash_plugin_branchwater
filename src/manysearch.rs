use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use anyhow::Result;

use crate::utils::{prepare_query,
                   load_sketchlist_filenames, load_sketches};

use log::{error, info};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
///
/// TODO:
///   - support jaccard as well as containment/overlap
///   - support md5 output columns; other?

pub fn manysearch<P: AsRef<Path>>(
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
        writeln!(&mut writer, "query_name,query_md5,match_name,match_md5,containment,jaccard,intersect_hashes").unwrap();
        for (query, query_md5, m, m_md5, cont, jaccard, overlap) in recv.into_iter() {
            writeln!(&mut writer, "\"{}\",{},\"{}\",{},{},{},{}",
                     query, query_md5, m, m_md5, cont, jaccard, overlap).ok();
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

                        let mut merged = q.minhash.clone();
                        merged.merge(&search_sm.minhash).ok();
                        let total_size = merged.size() as f64;

                        let containment = overlap / query_size;
                        let jaccard = overlap / total_size;
                        if containment > threshold {
                            results.push((q.name.clone(),
                                          q.md5sum.clone(),
                                          search_sm.name.clone(),
                                          search_sm.md5sum.clone(),
                                          containment,
                                          jaccard,
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

