/// fastmultigather: Run gather for multiple queries against a list of files.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::selection::Selection;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use camino::Utf8Path as PathBuf;

use crate::utils::{
    consume_query_by_gather, load_collection, load_sketches, write_prefetch, PrefetchResult,
    ReportType,
};

pub fn fastmultigather(
    query_filepath: String,
    against_filepath: String,
    threshold_bp: usize,
    scaled: usize,
    selection: &Selection,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    // load query collection
    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

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

    // load against collection
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;
    // load against sketches into memory, downsampling on the way
    let against = load_sketches(against_collection, selection, ReportType::Against).unwrap();

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    query_collection.par_iter().for_each(|(_idx, record)| {
        // increment counter of # of queries. q: could we instead use the _idx from par_iter(), or will it vary based on thread?
        let _i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);
        // Load query sig (downsampling happens here)
        match query_collection.sig_from_record(record) {
            Ok(query_sig) => {
                let name = query_sig.name();
                let prefix = name.split(' ').next().unwrap_or_default().to_string();
                let location = PathBuf::new(&prefix).file_name().unwrap();
                if let Some(query_mh) = query_sig.minhash() {
                    let matchlist: BinaryHeap<PrefetchResult> = against
                        .iter()
                        .filter_map(|against| {
                            let mut mm: Option<PrefetchResult> = None;
                            if let Ok(overlap) = against.minhash.count_common(query_mh, false) {
                                if overlap >= threshold_hashes {
                                    let result = PrefetchResult {
                                        name: against.name.clone(),
                                        md5sum: against.md5sum.clone(),
                                        minhash: against.minhash.clone(),
                                        overlap,
                                    };
                                    mm = Some(result);
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
                    // different warning here? Could not load sig from record??
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        record.internal_location()
                    );
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            }
            Err(_) => {
                // different warning here? Could not load sig from record??
                eprintln!(
                    "WARNING: no compatible sketches in path '{}'",
                    record.internal_location()
                );
                let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
        }
    });

    println!(
        "DONE. Processed {} queries total.",
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
