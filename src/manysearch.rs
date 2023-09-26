/// multisearch: massively search of many queries against many large subjects.
/// the OG MAGsearch, branchwater, etc.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::Sketch;
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, load_sigpaths_from_zip_or_pathlist, load_sketches_from_zip_or_pathlist,
    prepare_query, ReportType, SearchResult,
};

pub fn manysearch<P: AsRef<Path>>(
    querylist: P,
    siglist: P,
    template: Sketch,
    threshold: f64,
    output: Option<P>,
) -> Result<()> {
    // Read in list of query paths.
    eprintln!(
        "Reading list of queries from: '{}'",
        querylist.as_ref().display()
    );

    // Load all queries into memory at once.
    let queries = load_sketches_from_zip_or_pathlist(querylist, &template, ReportType::Query)?;

    // Load all _paths_, not signatures, into memory.
    let siglist_name = siglist.as_ref().to_string_lossy().to_string();
    let (search_sigs_paths, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(siglist, &template, ReportType::Against)?;

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
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load search signature from path:
            match Signature::from_path(filename) {
                Ok(search_sigs) => {
                    let location = filename.display().to_string();
                    if let Some(search_sm) = prepare_query(&search_sigs, &template, &location) {
                        // search for matches & save containment.
                        for q in queries.iter() {
                            let overlap =
                                q.minhash.count_common(&search_sm.minhash, false).unwrap() as f64;
                            let query_size = q.minhash.size() as f64;
                            let target_size = search_sm.minhash.size() as f64;

                            let containment_query_in_target = overlap / query_size;
                            let containment_in_target = overlap / target_size;
                            let max_containment =
                                containment_query_in_target.max(containment_in_target);
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
                        // for reading zips, this is likely not a useful warning and
                        // would show up too often (every sig is stored as individual file).
                        if !siglist_name.ends_with(".zip") {
                            eprintln!(
                                "WARNING: no compatible sketches in path '{}'",
                                filename.display()
                            );
                        }
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    Some(results)
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
        .try_for_each_with(send, |s, m| s.send(m));

    // do some cleanup and error handling -
    if let Err(e) = send {
        eprintln!("Unable to send internal data: {:?}", e);
    }

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
            "WARNING: skipped {} search paths - no compatible signatures.",
            skipped_paths
        );
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} search paths failed to load. See error messages above.",
            failed_paths
        );
    }

    Ok(())
}
