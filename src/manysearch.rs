/// multisearch: massively search of many queries against many large subjects.
/// the OG MAGsearch, branchwater, etc.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
use anyhow::Result;
use rayon::prelude::*;
use stats::{median, stddev};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, load_collection, ManySearchResult, MultiCollection, ReportType,
    SmallSignature,
};
use sourmash::ani_utils::ani_from_containment;
use sourmash::errors::SourmashError;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;

type AbundanceStats = (
    Option<u64>,
    Option<u64>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
);

pub fn manysearch(
    query_filepath: String,
    against_filepath: String,
    selection: Selection,
    threshold: f64,
    output: Option<String>,
    allow_failed_sigpaths: bool,
    ignore_abundance: bool,
    output_all_comparisons: bool,
) -> Result<()> {
    // Load query collection
    let query_collection = load_collection(
        &query_filepath,
        &selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    // Figure out what scaled to use - either from selection, or from query.
    let common_scaled: u32 = if let Some(set_scaled) = selection.scaled() {
        set_scaled
    } else {
        let s = *query_collection.max_scaled().expect("no records!?");
        eprintln!(
            "Setting scaled={} based on max scaled in query collection",
            s
        );
        s
    };

    let mut selection = selection;
    selection.set_scaled(common_scaled);

    // load all query sketches into memory, downsampling on the way
    let query_sketchlist = query_collection.load_sketches()?;

    // Against: Load collection, potentially off disk & not into memory.
    let against_collection = load_collection(
        &against_filepath,
        &selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    let (n_processed, skipped_paths, failed_paths) = manysearch_obj(
        &query_sketchlist,
        &against_collection,
        threshold,
        common_scaled,
        output,
        ignore_abundance,
        output_all_comparisons,
    )?;

    eprintln!("DONE. Processed {} search sigs", n_processed);

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

pub(crate) fn manysearch_obj(
    query_sketchlist: &Vec<SmallSignature>,
    against_collection: &MultiCollection,
    threshold: f64,
    common_scaled: u32,
    output: Option<String>,
    ignore_abundance: bool,
    output_all_comparisons: bool,
) -> Result<(usize, usize, usize)> {
    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<ManySearchResult>(rayon::current_num_threads());

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

    let send = against_collection
        .par_iter()
        .filter_map(|(coll, _idx, record)| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 && i > 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            match coll.sig_from_record(record) {
                Ok(against_sig) => {
                    let against_name = against_sig.name();
                    let against_md5 = against_sig.md5sum();

                    if let Ok(against_mh) =
                        <SigStore as TryInto<KmerMinHash>>::try_into(against_sig)
                    {
                        for query in query_sketchlist.iter() {
                            let sr = calculate_manysearch_result(
                                query,
                                &against_mh,
                                &against_name,
                                &against_md5,
                                threshold,
                                common_scaled,
                                ignore_abundance,
                                output_all_comparisons,
                            );
                            if let Some(sr) = sr {
                                results.push(sr);
                            }
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            record.internal_location()
                        );
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        record.internal_location()
                    );
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            }

            Some(results)
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    send.expect("Unable to send internal data");
    thrd.join().expect("Unable to join internal thread.");

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((i, skipped_paths, failed_paths))
}

// inflate_abundances: "borrow" the abundances from 'against' onto the
// intersection with 'query'.

fn inflate_abundances(
    query: &KmerMinHash,
    against: &KmerMinHash,
) -> Result<AbundanceStats, SourmashError> {
    let abunds: Vec<u64>;
    let sum_weighted: u64;
    let sum_all_abunds: u64 = against.sum_abunds();

    (abunds, sum_weighted) = query.inflated_abundances(against)?;

    let average_abund = sum_weighted as f64 / abunds.len() as f64;
    let median_abund = median(abunds.iter().cloned()).expect("error");
    let std_abund = stddev(abunds.iter().cloned());

    Ok((
        Some(sum_all_abunds),
        Some(sum_weighted),
        Some(average_abund),
        Some(median_abund),
        Some(std_abund),
    ))
}

// calculate_manysearch_result: calculate all the things

fn calculate_manysearch_result(
    query: &SmallSignature,
    against_mh: &KmerMinHash,
    against_name: &str,
    against_md5: &str,
    threshold: f64,
    common_scaled: u32,
    ignore_abundance: bool,
    output_all_comparisons: bool,
) -> Option<ManySearchResult> {
    // be paranoid and confirm scaled match.
    if query.minhash.scaled() != common_scaled {
        panic!("different query scaled");
    }
    if against_mh.scaled() != common_scaled {
        panic!("different against scaled");
    }

    let overlap = query
        .minhash
        .count_common(against_mh, false)
        .expect("incompatible sketches") as f64;

    let query_size = query.minhash.size() as f64;
    let containment_query_in_target = overlap / query_size;

    // only calculate results if we have shared hashes
    if containment_query_in_target > threshold || output_all_comparisons {
        let target_size = against_mh.size() as f64;
        let containment_target_in_query = overlap / target_size;

        let max_containment = containment_query_in_target.max(containment_target_in_query);
        let jaccard = overlap / (target_size + query_size - overlap);

        let qani = ani_from_containment(containment_query_in_target, against_mh.ksize() as f64);
        let mani = ani_from_containment(containment_target_in_query, against_mh.ksize() as f64);
        let query_containment_ani = Some(qani);
        let match_containment_ani = Some(mani);
        let average_containment_ani = Some((qani + mani) / 2.);
        let max_containment_ani = Some(f64::max(qani, mani));

        let calc_abund_stats = against_mh.track_abundance() && !ignore_abundance;
        let (total_weighted_hashes, n_weighted_found, average_abund, median_abund, std_abund) =
            if calc_abund_stats {
                inflate_abundances(&query.minhash, against_mh).ok()?
            } else {
                (None, None, None, None, None)
            };

        let sr = ManySearchResult {
            query_name: query.name.clone(),
            query_md5: query.md5sum.clone(),
            match_name: against_name.to_string(),
            containment: containment_query_in_target,
            intersect_hashes: overlap as u64,
            ksize: query.minhash.ksize() as u16,
            scaled: query.minhash.scaled(),
            moltype: query.minhash.hash_function().to_string(),
            match_md5: Some(against_md5.to_string()),
            jaccard: Some(jaccard),
            max_containment: Some(max_containment),
            average_abund,
            median_abund,
            std_abund,
            query_containment_ani,
            match_containment_ani,
            average_containment_ani,
            max_containment_ani,
            n_weighted_found,
            total_weighted_hashes,
        };
        return Some(sr);
    }
    None
}
