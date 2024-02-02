/// mastiff_manysearch: mastiff-indexed version of manysearch.
use anyhow::Result;
use camino::Utf8PathBuf as PathBuf;
use rayon::prelude::*;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::Sketch;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_collection, ReportType, SearchResult,
};

pub fn mastiff_manysearch(
    queries_path: String,
    index: PathBuf,
    selection: &Selection,
    minimum_containment: f64,
    output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(&index) {
        bail!("'{}' is not a valid RevIndex database", index);
    }
    // Open database once
    let db = RevIndex::open(index, true)?;
    println!("Loaded DB");

    // Load query paths
    let query_collection = load_collection(
        &queries_path,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel::<SearchResult>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send_result = query_collection
        .par_iter()
        .filter_map(|(idx, record)| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];
            match query_collection.sig_for_dataset(idx) {
                Ok(query_sig) => {
                    for sketch in query_sig.iter() {
                        if let Sketch::MinHash(query_mh) = sketch {
                            // let location = query_sig.filename();
                            let query_size = query_mh.size();
                            let counter = db.counter_for_query(&query_mh);
                            let matches =
                                db.matches_from_counter(counter, minimum_containment as usize);

                            // filter the matches for containment
                            for (path, overlap) in matches {
                                let containment = overlap as f64 / query_size as f64;
                                if containment >= minimum_containment {
                                    results.push(SearchResult {
                                        query_name: query_sig.name(),
                                        query_md5: query_sig.md5sum(),
                                        match_name: path.clone(),
                                        containment,
                                        intersect_hashes: overlap,
                                        match_md5: None,
                                        jaccard: None,
                                        max_containment: None,
                                    });
                                }
                            }
                        } else {
                            eprintln!(
                                "WARNING: no compatible sketches in path '{}'",
                                query_sig.filename()
                            );
                            let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                        }
                    }
                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                }
                Err(err) => {
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!(
                        "WARNING: could not load sketches from path '{}'",
                        record.internal_location()
                    );
                    None
                }
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
        eprintln!("Error during parallel processing: {}", e);
    }

    // join the writer thread
    if let Err(e) = thrd.join() {
        eprintln!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} search sigs", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

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
