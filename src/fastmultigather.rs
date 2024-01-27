/// fastmultigather: Run gather for multiple queries against a list of files.
use anyhow::Result;
use rayon::prelude::*;

use serde::Serialize;
use sourmash::selection::Selection;
use sourmash::sketch::Sketch;
use sourmash::storage::SigStore;
use sourmash::{selection, signature::Signature};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use camino::{Utf8Path, Utf8PathBuf};

use crate::utils::{
    consume_query_by_gather, load_collection, load_sigpaths_from_zip_or_pathlist,
    load_sketches_from_zip_or_pathlist, prepare_query, write_prefetch, PrefetchResult, ReportType,
};

pub fn fastmultigather(
    query_filepath: camino::Utf8PathBuf,
    against_filepath: camino::Utf8PathBuf,
    threshold_bp: usize,
    scaled: usize,
    selection: &Selection,
) -> Result<()> {
    // load the list of query paths
    let query_collection = load_collection(&query_filepath, selection, ReportType::Query)?;
    println!("Loaded {} sig paths in querylist", query_collection.len());

    let threshold_hashes: u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }
    .try_into()
    .unwrap();

    println!("threshold overlap: {} {}", threshold_hashes, threshold_bp);

    // Load all the against sketches
    let against_collection = load_collection(&against_filepath, selection, ReportType::Against)?;
    // load actual signatures
    let mut sketchlist: Vec<SigStore> = vec![];

    for (idx, record) in against_collection.iter() {
        if let Ok(sig) = against_collection.sig_for_dataset(idx) {
            sketchlist.push(sig);
        } else {
            eprintln!("Failed to load 'against' record: {}", record.name());
        }
    }

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    query_collection.par_iter().for_each(|(idx, record)| {
        // increment counter of # of queries. q: could we instead use the index from par_iter()?
        let _i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);
        // Load query sig
        match query_collection.sig_for_dataset(idx) {
            Ok(query_sig) => {
                let prefix = query_sig.name();
                let location = Utf8Path::new(&prefix).file_name().unwrap();
                for sketch in query_sig.iter() {
                    // Access query MinHash
                    if let Sketch::MinHash(query) = sketch {
                        let matchlist: BinaryHeap<PrefetchResult> = sketchlist
                            .par_iter()
                            .filter_map(|sm| {
                                let mut mm = None;
                                // Access against MinHash
                                if let Some(sketch) = sm.sketches().get(0) {
                                    if let Sketch::MinHash(against_sketch) = sketch {
                                        if let Ok(overlap) =
                                            against_sketch.count_common(&query, true)
                                        {
                                            if overlap >= threshold_hashes {
                                                let result = PrefetchResult {
                                                    name: sm.name(),
                                                    md5sum: sm.md5sum().clone(),
                                                    minhash: against_sketch.clone(),
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
                        if !matchlist.is_empty() {
                            let prefetch_output = format!("{}.prefetch.csv", location);
                            let gather_output = format!("{}.gather.csv", location);

                            // Save initial list of matches to prefetch output
                            write_prefetch(&query_sig, Some(prefetch_output), &matchlist).ok();

                            // Now, do the gather!
                            consume_query_by_gather(
                                query_sig.clone(),
                                matchlist,
                                threshold_hashes,
                                Some(gather_output),
                            )
                            .ok();
                        } else {
                            println!("No matches to '{}'", location);
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            record.internal_location()
                        );
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                }
            }
            Err(_) => {
                eprintln!(
                    "WARNING: no compatible sketches in path '{}'",
                    record.internal_location()
                );
                let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
        }
    });

    println!(
        "Processed {} queries total.",
        processed_queries.into_inner()
    );

    let skipped_paths = skipped_paths.into_inner();
    let failed_paths = failed_paths.into_inner();

    if skipped_paths > 0 {
        eprintln!(
            "WARNING: skipped {} query paths - no compatible signatures.",
            skipped_paths
        );
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} query paths failed to load. See error messages above.",
            failed_paths
        );
    }

    Ok(())
}
