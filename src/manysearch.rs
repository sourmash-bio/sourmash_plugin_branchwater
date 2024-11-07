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

use crate::utils::{csvwriter_thread, load_collection, ReportType, SearchResult};
use sourmash::ani_utils::ani_from_containment;
use sourmash::errors::SourmashError;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;

pub fn manysearch(
    query_filepath: String,
    against_filepath: String,
    selection: Selection,
    threshold: f64,
    output: Option<String>,
    allow_failed_sigpaths: bool,
    ignore_abundance: bool,
) -> Result<()> {
    // Load query collection
    let query_collection = load_collection(
        &query_filepath,
        &selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    // load all query sketches into memory, downsampling on the way
    let query_sketchlist = query_collection.load_sketches(&selection)?;

    // Against: Load collection, potentially off disk & not into memory.
    let against_collection = load_collection(
        &against_filepath,
        &selection,
        ReportType::Against,
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

    let send = against_collection
        .par_iter()
        .filter_map(|(coll, _idx, record)| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 && i > 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // against downsampling happens here
            match coll.sig_from_record(record) {
                Ok(against_sig) => {
                    let against_name = against_sig.name();
                    let against_md5 = against_sig.md5sum();

                    if let Ok(against_mh) = against_sig.try_into() {
                        for query in query_sketchlist.iter() {
                            // avoid calculating details unless there is overlap
                            let overlap = query
                                .minhash
                                .count_common(&against_mh, true)
                                .expect("incompatible sketches")
                                as f64;

                            let query_size = query.minhash.size() as f64;
                            let containment_query_in_target = overlap / query_size;
                            // only calculate results if we have shared hashes
                            if containment_query_in_target > threshold {
                                let target_size = against_mh.size() as f64;
                                let containment_target_in_query = overlap / target_size;

                                let max_containment =
                                    containment_query_in_target.max(containment_target_in_query);
                                let jaccard = overlap / (target_size + query_size - overlap);

                                let qani = ani_from_containment(
                                    containment_query_in_target,
                                    against_mh.ksize() as f64,
                                );
                                let mani = ani_from_containment(
                                    containment_target_in_query,
                                    against_mh.ksize() as f64,
                                );
                                let query_containment_ani = Some(qani);
                                let match_containment_ani = Some(mani);
                                let average_containment_ani = Some((qani + mani) / 2.);
                                let max_containment_ani = Some(f64::max(qani, mani));

                                let calc_abund_stats =
                                    against_mh.track_abundance() && !ignore_abundance;
                                let (
                                    total_weighted_hashes,
                                    n_weighted_found,
                                    average_abund,
                                    median_abund,
                                    std_abund,
                                ) = if calc_abund_stats {
                                    downsample_and_inflate_abundances(&query.minhash, &against_mh)
                                        .ok()?
                                } else {
                                    (None, None, None, None, None)
                                };

                                results.push(SearchResult {
                                    query_name: query.name.clone(),
                                    query_md5: query.md5sum.clone(),
                                    match_name: against_name.clone(),
                                    containment: containment_query_in_target,
                                    intersect_hashes: overlap as u64,
                                    match_md5: Some(against_md5.clone()),
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
                                });
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

fn downsample_and_inflate_abundances(
    query: &KmerMinHash,
    against: &KmerMinHash,
) -> Result<
    (
        Option<u64>,
        Option<u64>,
        Option<f64>,
        Option<f64>,
        Option<f64>,
    ),
    SourmashError,
> {
    let query_scaled = query.scaled();
    let against_scaled = against.scaled();

    let abunds: Vec<u64>;
    let sum_weighted: u64;
    let sum_all_abunds: u64;

    // avoid downsampling if we can
    if against_scaled != query_scaled {
        let against_ds = against
            .clone()
            .downsample_scaled(query.scaled())
            .expect("cannot downsample sketch");
        (abunds, sum_weighted) = query.inflated_abundances(&against_ds)?;
        sum_all_abunds = against_ds.sum_abunds();
    } else {
        (abunds, sum_weighted) = query.inflated_abundances(against)?;
        sum_all_abunds = against.sum_abunds();
    }

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
