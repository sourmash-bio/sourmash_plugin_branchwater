/// mastiff_manysearch: mastiff-indexed version of manysearch.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::Sketch;
use std::path::Path;

use sourmash::index::revindex::RevIndex;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_sigpaths_from_zip_or_pathlist, prepare_query,
    ReportType, SearchResult,
};

pub fn mastiff_manysearch<P: AsRef<Path>>(
    queries_file: P,
    index: P,
    template: Sketch,
    minimum_containment: f64,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(index.as_ref()) {
        bail!(
            "'{}' is not a valid RevIndex database",
            index.as_ref().display()
        );
    }
    // Open database once
    let db = RevIndex::open(index.as_ref(), true);
    println!("Loaded DB");

    // Load query paths
    let queryfile_name = queries_file.as_ref().to_string_lossy().to_string();
    let (query_paths, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(&queries_file, &template, ReportType::Query)?;

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
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load query signature from path:
            match Signature::from_path(filename) {
                Ok(query_sig) => {
                    let location = filename.display().to_string();
                    if let Some(query) = prepare_query(&query_sig, &template, &location) {
                        let query_size = query.minhash.size() as f64;
                        // search mastiff db
                        let counter = db.counter_for_query(&query.minhash);
                        let matches =
                            db.matches_from_counter(counter, minimum_containment as usize);

                        // filter the matches for containment
                        for (path, overlap) in matches {
                            let containment = overlap as f64 / query_size;
                            if containment >= minimum_containment {
                                results.push(SearchResult {
                                    query_name: query.name.clone(),
                                    query_md5: query.md5sum.clone(),
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
                        // for reading zips, this is likely not a useful warning and
                        // would show up too often (every sig is stored as individual file).
                        if !queryfile_name.ends_with(".zip") {
                            eprintln!(
                                "WARNING: no compatible sketches in path '{}'",
                                filename.display()
                            );
                        }
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
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
                        filename.display()
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

    // _temp_dir goes out of scope => is deleted.

    Ok(())
}
