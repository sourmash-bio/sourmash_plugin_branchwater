/// multisearch: massively parallel in-memory sketch search.
use anyhow::Result;
use rayon::prelude::*;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::collections::HashMap;


use crate::utils::{
    csvwriter_thread, load_collection, load_sketches, MultiSearchResult, ReportType
};
use crate::search_significance::{
    compute_inverse_document_frequency, 
    get_prob_overlap, 
    get_term_frequency_inverse_document_frequency, 
    merge_all_minhashes,
    get_hash_frequencies,
    Normalization,
};
use sourmash::ani_utils::ani_from_containment;


/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn multisearch(
    query_filepath: String,
    against_filepath: String,
    threshold: f64,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    estimate_ani: bool,
    estimate_prob_overlap: bool,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all queries into memory at once.

    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;
    let queries = load_sketches(query_collection, selection, ReportType::Query).unwrap();

    // Load all against sketches into memory at once.
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    let againsts: Vec<crate::utils::SmallSignature> = load_sketches(against_collection, selection, ReportType::Against).unwrap();


    // if estimate_prob_overlap {
    let n_comparisons: f64 = againsts.len() as f64 * queries.len() as f64;
    // Combine all the queries and against into a single signature each, to get their 
    // underlying distribution of hashes across the whole input
    eprintln!("Merging queries ...");
    let queries_merged_mh: KmerMinHash = merge_all_minhashes(&queries).unwrap();
    eprintln!("\tDone.\n");
    eprintln!("Merging against ...");
    let against_merged_mh: KmerMinHash = merge_all_minhashes(&againsts).unwrap();
    eprintln!("\tDone.\n");       

    // Precompute
    eprintln!("Computing Inverse Document Frequency (IDF) of hashes in all againsts ...");
    let inverse_document_frequency =  compute_inverse_document_frequency(
        &against_merged_mh, 
        &againsts,
        Some(true),
    );
    eprintln!("\tDone.\n");

    eprintln!("Computing frequency of hashvals across all againsts ...");
    let against_merged_frequencies: HashMap<u64, f64> = get_hash_frequencies(
        &against_merged_mh, 
        Some(Normalization::L1), 
    );
    eprintln!("\tDone.\n");

    eprintln!("Computing frequency of hashvals across all queries ...");
    let query_merged_frequencies: HashMap<u64, f64> = get_hash_frequencies(
        &queries_merged_mh, 
        Some(Normalization::L1), 
    );
    eprintln!("\tDone.\n");

    eprintln!("Computing hashval term frequencies within each query ...");
    let query_term_frequencies: HashMap<String, HashMap<u64, f64>> = HashMap::from(
        queries
        .par_iter()
        .map(|query|
             (
                query.md5sum.clone(), 
                get_hash_frequencies(
                    &query.minhash, 
                    Some(Normalization::L1)
                )
            )
        ).collect::<HashMap<String, HashMap<u64, f64>>>()
    );
    eprintln!("\tDone.\n");


    // Can't get below to work because then get "cannot find value `variable_name` in this scope" error:
    // --- begin error ---
    // error[E0425]: cannot find value `queries_merged_mh` in this scope
    //    --> src/multisearch.rs:138:30
    //     |
    // 138 | ...                   &queries_merged_mh, 
    //     |                        ^^^^^^^^^^^^^^^^^ not found in this scope

    // error[E0425]: cannot find value `against_merged_mh` in this scope
    //    --> src/multisearch.rs:139:30
    //     |
    // 139 | ...                   &against_merged_mh, 
    //     |                        ^^^^^^^^^^^^^^^^^ not found in this scope

    // error[E0425]: cannot find value `n_comparisons` in this scope
    //    --> src/multisearch.rs:144:78
    //     |
    // 144 |                     let prob_overlap_adjusted = Some(prob_overlap.unwrap() * n_comparisons);
    //     |                                                                              ^^^^^^^^^^^^^ not found in this scope

    // error[E0425]: cannot find value `inverse_document_frequency` in this scope
    //    --> src/multisearch.rs:154:30
    //     |
    // 154 |   ...                   &inverse_document_frequency
    //     |                          ^^^^^^^^^^^^^^^^^^^^^^^^^^ help: a function with a similar name exists: `compute_inverse_document_frequency`
    //     |
    //    ::: src/search_significance.rs:197:1
    // --- end   error ---
    // } else {
    //     let mut n_comparisons = 0.0;
    //     let mut queries_merged_mh: KmerMinHash = Default::default();
    //     let mut against_merged_mh: KmerMinHash = Default::default();
    // };


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
    let ksize = selection.ksize().unwrap() as f64;
    let logged = false;

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

                let overlap = query.minhash.count_common(&against.minhash, false).unwrap() as f64;
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
                    let mut containment_adjusted: Option<f64>= None;
                    let mut containment_adjusted_log10: Option<f64>= None;
                    let mut tf_idf_score: Option<f64> = None;

                    // if estimate_prob_overlap {
                    let overlapping_hashvals: Vec<u64> = query.minhash.intersection(
                        &against.minhash
                    )
                        .unwrap().0;

                    // Related to getting the probability of overlap between query and target
                    let prob_overlap = Some(
                        get_prob_overlap(
                            &overlapping_hashvals, 
                            &query_merged_frequencies, 
                            &against_merged_frequencies, 
                        )
                    );
                    // Do simple, conservative Bonferroni correction
                    let prob_overlap_adjusted = Some(prob_overlap.unwrap() * n_comparisons);

                    // Rescale the containment based on the probability of overlap, 
                    // such that smaller probability -> higher containment
                    let containment_adjusted = Some(containment_query_in_target / prob_overlap_adjusted.unwrap());
                    let containment_adjusted_log10 = Some(containment_adjusted.unwrap().log10());
                    let tf_idf_score = Some(
                        get_term_frequency_inverse_document_frequency(
                            &overlapping_hashvals, 
                            &query_term_frequencies[&query.md5sum], 
                            &inverse_document_frequency
                        )
                    );
                    // } else {

                    // }


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
