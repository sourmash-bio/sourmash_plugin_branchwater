/// multisearch: massively parallel in-memory sketch search.
use anyhow::Result;
use rayon::prelude::*;
use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::search_significance::{
    compute_inverse_document_frequency, get_hash_frequencies, get_prob_overlap,
    get_term_frequency_inverse_document_frequency, merge_all_minhashes, Normalization,
};
use crate::utils::multicollection::SmallSignature;
use crate::utils::{csvwriter_thread, load_collection, MultiSearchResult, ReportType};
use sourmash::ani_utils::ani_from_containment;

#[derive(Default, Clone, Debug)]
struct ProbOverlapStats {
    prob_overlap: f64,
    prob_overlap_adjusted: f64,
    containment_adjusted: f64,
    containment_adjusted_log10: f64,
    tf_idf_score: f64,
}

/// Computes probability overlap statistics for a single pair of signatures
fn compute_single_prob_overlap(
    query: &SmallSignature,
    against: &SmallSignature,
    n_comparisons: f64,
    query_merged_frequencies: &HashMap<u64, f64>,
    against_merged_frequencies: &HashMap<u64, f64>,
    query_term_frequencies: &HashMap<String, HashMap<u64, f64>>,
    inverse_document_frequency: &HashMap<u64, f64>,
    containment_query_in_target: f64,
) -> ProbOverlapStats {
    let overlapping_hashvals: Vec<u64> = query
        .minhash
        .intersection(&against.minhash)
        .expect("Intersection of query and against minhashes")
        .0;

    let prob_overlap = get_prob_overlap(
        &overlapping_hashvals,
        query_merged_frequencies,
        against_merged_frequencies,
    );

    let prob_overlap_adjusted = prob_overlap * n_comparisons;
    let containment_adjusted = containment_query_in_target / prob_overlap_adjusted;

    ProbOverlapStats {
        prob_overlap,
        prob_overlap_adjusted,
        containment_adjusted,
        containment_adjusted_log10: containment_adjusted.log10(),
        tf_idf_score: get_term_frequency_inverse_document_frequency(
            &overlapping_hashvals,
            &query_term_frequencies[&query.md5sum],
            inverse_document_frequency,
        ),
    }
}

/// Computes probability overlap statistics for queries and against signatures
/// Estimate probability of overlap between query sig and against sig, using
/// underlying distribution of hashvals for all queries and all againsts
fn compute_prob_overlap_stats(
    queries: &Vec<SmallSignature>,
    againsts: &Vec<SmallSignature>,
) -> (
    f64,
    HashMap<u64, f64>,
    HashMap<u64, f64>,
    HashMap<String, HashMap<u64, f64>>,
    HashMap<u64, f64>,
) {
    let n_comparisons = againsts.len() as f64 * queries.len() as f64;

    // Combine all the queries and against into a single signature each
    eprintln!("Merging queries ...");
    let queries_merged_mh: KmerMinHash =
        merge_all_minhashes(queries).expect("Merging query minhashes");
    eprintln!("\tDone.\n");

    eprintln!("Merging against ...");
    let against_merged_mh: KmerMinHash =
        merge_all_minhashes(againsts).expect("Merging against minhashes");
    eprintln!("\tDone.\n");

    // Compute IDF
    eprintln!("Computing Inverse Document Frequency (IDF) of hashes in all againsts ...");
    let inverse_document_frequency =
        compute_inverse_document_frequency(&against_merged_mh, againsts, Some(true));
    eprintln!("\tDone.\n");

    // Compute frequencies
    eprintln!("Computing frequency of hashvals across all againsts (L1 Norm) ...");
    let against_merged_frequencies =
        get_hash_frequencies(&against_merged_mh, Some(Normalization::L1));
    eprintln!("\tDone.\n");

    eprintln!("Computing frequency of hashvals across all queries (L1 Norm) ...");
    let query_merged_frequencies =
        get_hash_frequencies(&queries_merged_mh, Some(Normalization::L1));
    eprintln!("\tDone.\n");

    // Compute term frequencies
    eprintln!("Computing hashval term frequencies within each query (L2 Norm) ...");
    let query_term_frequencies = queries
        .par_iter()
        .map(|query| {
            (
                query.md5sum.clone(),
                get_hash_frequencies(&query.minhash, Some(Normalization::L2)),
            )
        })
        .collect::<HashMap<String, HashMap<u64, f64>>>();
    eprintln!("\tDone.\n");

    (
        n_comparisons,
        query_merged_frequencies,
        against_merged_frequencies,
        query_term_frequencies,
        inverse_document_frequency,
    )
}

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn multisearch(
    query_filepath: String,
    against_filepath: String,
    threshold: f64,
    selection: Selection,
    allow_failed_sigpaths: bool,
    estimate_ani: bool,
    estimate_prob_overlap: bool,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all queries into memory at once.
    let query_collection = load_collection(
        &query_filepath,
        &selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    let expected_scaled = match selection.scaled() {
        Some(s) => s,
        None => {
            let s = *query_collection.max_scaled().expect("no records!?") as u32;
            eprintln!(
                "Setting scaled={} based on max scaled in query collection",
                s
            );
            s
        }
    };

    let ksize = selection.ksize().unwrap() as f64;

    let mut new_selection = selection;
    new_selection.set_scaled(expected_scaled);

    // update selection with new scaled.
    let query_collection = query_collection.select(&new_selection)?;

    let queries: Vec<SmallSignature> = query_collection.load_sketches()?;

    // Load all against sketches into memory at once.
    let against_collection = load_collection(
        &against_filepath,
        &new_selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    let againsts: Vec<SmallSignature> = against_collection.load_sketches()?;

    let (
        n_comparisons,
        query_merged_frequencies,
        against_merged_frequencies,
        query_term_frequencies,
        inverse_document_frequency,
    ) = if estimate_prob_overlap {
        compute_prob_overlap_stats(&queries, &againsts)
    } else {
        (
            0.0,
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
        )
    };

    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<MultiSearchResult>(rayon::current_num_threads());

    // // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_cmp = AtomicUsize::new(0);

    let send = againsts
        .par_iter()
        .filter_map(|against| {
            let mut results = vec![];
            // search for matches & save containment.
            for query in queries.iter() {
                let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
                if i % 100000 == 0 && i > 0 {
                    eprintln!("Processed {} comparisons", i);
                }

                // be paranoid and check scaled.
                if query.minhash.scaled() != expected_scaled {
                    panic!("different scaled for query");
                }

                if against.minhash.scaled() != expected_scaled {
                    panic!("different scaled for against");
                }

                let overlap = query
                    .minhash
                    .count_common(&against.minhash, false)
                    .expect("cannot compare query and against!?")
                    as f64;
                // use downsampled sizes
                let query_size = query.minhash.size() as f64;
                let target_size = against.minhash.size() as f64;

                let containment_query_in_target = overlap / query_size;

                if containment_query_in_target > threshold {
                    let containment_target_in_query = overlap / target_size;
                    let max_containment =
                        containment_query_in_target.max(containment_target_in_query);
                    let jaccard = overlap / (target_size + query_size - overlap);
                    let mut query_containment_ani = None;
                    let mut match_containment_ani = None;
                    let mut average_containment_ani = None;
                    let mut max_containment_ani = None;
                    let mut prob_overlap: Option<f64> = None;
                    let mut prob_overlap_adjusted: Option<f64> = None;
                    let mut containment_adjusted: Option<f64> = None;
                    let mut containment_adjusted_log10: Option<f64> = None;
                    let mut tf_idf_score: Option<f64> = None;

                    // Compute probability overlap stats if requested
                    if estimate_prob_overlap {
                        let prob_stats = compute_single_prob_overlap(
                            query,
                            against,
                            n_comparisons,
                            &query_merged_frequencies,
                            &against_merged_frequencies,
                            &query_term_frequencies,
                            &inverse_document_frequency,
                            containment_query_in_target,
                        );
                        prob_overlap = Some(prob_stats.prob_overlap);
                        prob_overlap_adjusted = Some(prob_stats.prob_overlap_adjusted);
                        containment_adjusted = Some(prob_stats.containment_adjusted);
                        containment_adjusted_log10 = Some(prob_stats.containment_adjusted_log10);
                        tf_idf_score = Some(prob_stats.tf_idf_score);
                    }

                    // estimate ANI values
                    if estimate_ani {
                        let qani = ani_from_containment(containment_query_in_target, ksize);
                        let mani = ani_from_containment(containment_target_in_query, ksize);
                        query_containment_ani = Some(qani);
                        match_containment_ani = Some(mani);
                        average_containment_ani = Some((qani + mani) / 2.);
                        max_containment_ani = Some(f64::max(qani, mani));
                    }

                    results.push(MultiSearchResult {
                        query_name: query.name.clone(),
                        query_md5: query.md5sum.clone(),
                        match_name: against.name.clone(),
                        match_md5: against.md5sum.clone(),
                        ksize: query.minhash.ksize() as u16,
                        scaled: query.minhash.scaled(),
                        moltype: query.minhash.hash_function().to_string(),
                        containment: containment_query_in_target,
                        max_containment,
                        jaccard,
                        intersect_hashes: overlap,
                        query_containment_ani,
                        match_containment_ani,
                        average_containment_ani,
                        max_containment_ani,
                        prob_overlap,
                        prob_overlap_adjusted,
                        containment_adjusted,
                        containment_adjusted_log10,
                        tf_idf_score,
                    })
                }
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
        eprintln!("Unable to send internal data: {:?}", e);
    }

    if let Err(e) = thrd.join() {
        eprintln!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_cmp.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} comparisons", i);

    Ok(())
}
