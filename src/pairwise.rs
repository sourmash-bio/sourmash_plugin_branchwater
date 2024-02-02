/// pairwise: massively parallel in-memory pairwise comparisons.
use anyhow::Result;
use rayon::prelude::*;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, load_collection, load_mh_with_name_and_md5, MultiSearchResult, ReportType,
};
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;

/// Perform pairwise comparisons of all signatures in a list.
///
/// Note: this function loads all _signatures_ into memory.

pub fn pairwise(
    siglist: String,
    threshold: f64,
    selection: &Selection,
    output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all sigs into memory at once.
    let collection = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
    )?;

    if collection.len() <= 1 {
        bail!(
            "Pairwise requires two or more sketches. Check input: '{:?}'",
            &siglist
        )
    }
    let sketches = load_mh_with_name_and_md5(collection, &selection, ReportType::General).unwrap();

    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<MultiSearchResult>(rayon::current_num_threads());

    // // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all signature,
    // Results written to the writer thread above.

    let processed_cmp = AtomicUsize::new(0);

    sketches
        .par_iter()
        .enumerate()
        .for_each(|(idx, (q1, q1_name, q1_md5))| {
            for (j, (q2, q2_name, q2_md5)) in sketches.iter().enumerate().skip(idx + 1) {
                let overlap = q1.count_common(q2, false).unwrap() as f64;
                let query1_size = q1.size() as f64;
                let query2_size = q2.size() as f64;

                let containment_q1_in_q2 = overlap / query1_size;
                let containment_q2_in_q1 = overlap / query2_size;
                let max_containment = containment_q1_in_q2.max(containment_q2_in_q1);
                let jaccard = overlap / (query1_size + query2_size - overlap);

                if containment_q1_in_q2 > threshold || containment_q2_in_q1 > threshold {
                    send.send(
                        (MultiSearchResult {
                            query_name: q1_name.clone(),
                            query_md5: q1_md5.clone(),
                            match_name: q2_name.clone(),
                            match_md5: q2_md5.clone(),
                            containment: containment_q1_in_q2,
                            max_containment: max_containment,
                            jaccard: jaccard,
                            intersect_hashes: overlap,
                        }),
                    )
                    .unwrap();
                }

                let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
                if i % 100000 == 0 {
                    eprintln!("Processed {} comparisons", i);
                }
            }
        });

    // do some cleanup and error handling -
    drop(send); // close the channel

    if let Err(e) = thrd.join() {
        eprintln!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_cmp.load(atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} comparisons", i);

    Ok(())
}
