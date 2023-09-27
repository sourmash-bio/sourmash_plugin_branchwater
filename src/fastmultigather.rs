/// fastmultigather: Run gather for multiple queries against a list of files.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::signature::Signature;
use sourmash::sketch::Sketch;
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use crate::utils::{
    consume_query_by_gather, load_sigpaths_from_zip_or_pathlist,
    load_sketches_from_zip_or_pathlist, prepare_query, write_prefetch, PrefetchResult, ReportType,
};

pub fn fastmultigather<P: AsRef<Path> + std::fmt::Debug + Clone>(
    query_filenames: P,
    matchlist_filename: P,
    threshold_bp: usize,
    scaled: usize,
    template: Sketch,
) -> Result<()> {
    // load the list of query paths
    let queryfile_name = query_filenames.as_ref().to_string_lossy().to_string();
    let (querylist_paths, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(&query_filenames, &template, ReportType::Query)?;
    println!("Loaded {} sig paths in querylist", querylist_paths.len());

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
    let sketchlist =
        load_sketches_from_zip_or_pathlist(&matchlist_filename, &template, ReportType::Against)?;

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    querylist_paths.par_iter().for_each(|q| {
        // increment counter of # of queries
        let _i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);

        // set query_label to the last path element.
        let location = q.clone().into_os_string().into_string().unwrap();
        let location = location.split('/').last().unwrap().to_string();

        let query = match Signature::from_path(dbg!(q)) {
            Ok(sigs) => {
                let mm = prepare_query(&sigs, &template, &location);

                if mm.is_none() {
                    if !queryfile_name.ends_with(".zip") {
                        eprintln!("WARNING: no compatible sketches in path '{}'", q.display());
                    }
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
                mm
            }
            Err(err) => {
                eprintln!("Sketch loading error: {}", err);
                eprintln!(
                    "WARNING: could not load sketches from path '{}'",
                    q.display()
                );
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
                consume_query_by_gather(query, matchlist, threshold_hashes, Some(gather_output))
                    .ok();
            } else {
                println!("No matches to '{}'", location);
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
