/// Run counter-gather for multiple queries against a list of files.

use anyhow::Result;
use rayon::prelude::*;

use sourmash::signature::{Signature, SigsTrait};
use std::path::Path;
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;
use sourmash::prelude::MinHashOps;
use sourmash::prelude::FracMinHashOps;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use crate::utils::{prepare_query, write_prefetch, PrefetchResult,
    load_sketchlist_filenames, load_sketches, consume_query_by_gather};

pub fn fastmultigather<P: AsRef<Path> + std::fmt::Debug + Clone>(
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
            let location = q.clone().into_os_string().into_string().unwrap();
            let location = location.split('/').last().unwrap().to_string();

            let query = match Signature::from_path(dbg!(q)) {
                Ok(sigs) => {
                    let mm = prepare_query(&sigs, &template, &location);

                    if mm.is_none() {
                        eprintln!("WARNING: no compatible sketches in path '{}'",
                                q.display());
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    mm
                },
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!("WARNING: could not load sketches from path '{}'",
                                q.display());
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            };

            if let Some(query) = query {
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
                    let prefetch_output = format!("{location}.prefetch.csv");
                    let gather_output = format!("{location}.gather.csv");

                    // save initial list of matches to prefetch output
                    write_prefetch(&query, Some(prefetch_output), &matchlist).ok();

                    // now, do the gather!
                    consume_query_by_gather(query, matchlist, threshold_hashes,
                                            Some(gather_output)).ok();
                } else {
                    println!("No matches to '{}'", location);
                }
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