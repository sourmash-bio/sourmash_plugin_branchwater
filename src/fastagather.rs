// use crate::utils::buildutils::BuildCollection;
use anyhow::{bail, Result};

use needletail::parse_fastx_file;
// use sourmash::selection::Selection;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::utils::{
    build_selection, consume_query_by_gather, csvwriter_thread, load_collection,
    load_sketches_above_threshold, write_prefetch, BranchwaterGatherResult, ReportType,
};

#[allow(clippy::too_many_arguments)]
pub fn fastagather(
    query_filename: String,
    against_filepath: String,
    input_moltype: String,
    threshold_bp: u64,
    ksize: u32,
    scaled: u32,
    moltype: String,
    prefetch_output: Option<String>,
    gather_output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    // to start, implement straightforward record --> sketch --> gather
    // other ideas/to do:
    // - add full-file (lower resolution) prefetch first, to reduce search space
    // - parallelize and/or batch records?
    // - write function to filter fasta entries for those with matches (or those without)
    // - could use that with this structure for charcoal decontam or other functions
    // - add rocksdb search -- only way this will make sense.

    // Build minhash template based on parsed parameters

    // to do -- use input moltype to check that we can build desired moltype
    let _input_moltype = input_moltype.to_ascii_lowercase();

    let hash_function = HashFunctions::try_from(moltype.as_str())
        .map_err(|_| panic!("Unknown molecule type: {}", moltype))
        .unwrap();

    let mh_template = KmerMinHashBTree::new(scaled, ksize, hash_function, 42, true, 0); // track abund by default
    let selection = build_selection(ksize as u8, Some(scaled), &moltype);

    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes = {
        let x = threshold_bp / scaled as u64;
        if x > 0 {
            x
        } else {
            1
        }
    };

    // load collection to match against.
    let against_collection = load_collection(
        &against_filepath,
        &selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    let failed_records = AtomicUsize::new(0);
    // open file and start iterating through sequences
    // Open fasta file reader
    let mut reader = match parse_fastx_file(query_filename.clone()) {
        Ok(r) => r,
        Err(err) => {
            bail!("Error opening file {}: {:?}", query_filename, err);
        }
    };

    // channel for gather results
    let (send, recv) =
        std::sync::mpsc::sync_channel::<BranchwaterGatherResult>(rayon::current_num_threads());
    let _gather_out_thrd = csvwriter_thread(recv, gather_output);

    // later: can we parallelize across records or sigs? Do we want to batch groups of records for improved gather efficiency?
    while let Some(record_result) = reader.next() {
        // clone sig_templates for use
        // let mut sigcoll = sig_template.clone();
        let mut query_mh = mh_template.clone();
        match record_result {
            Ok(record) => {
                let query_name = std::str::from_utf8(record.id())
                    .expect("record.id() contains invalid UTF-8")
                    .to_string();
                if let Err(err) = query_mh.add_sequence(&record.seq(), true) {
                    eprintln!(
                        "Error building minhash from record: {}, {:?}",
                        query_filename, err
                    );
                    failed_records.fetch_add(1, Ordering::SeqCst);
                }
                let query_md5 = query_mh.md5sum();
                eprintln!("query minhash; {:?}", query_mh);

                // now do prefetch/gather
                let prefetch_result = load_sketches_above_threshold(
                    against_collection.clone(), // can we get rid of this clone??
                    &KmerMinHash::from(query_mh.clone()),
                    threshold_hashes,
                )?;
                let matchlist = prefetch_result.0;
                let _skipped_paths = prefetch_result.1;
                let _failed_paths = prefetch_result.2;

                if prefetch_output.is_some() {
                    write_prefetch(
                        query_filename.clone(),
                        query_name.clone(),
                        query_md5,
                        prefetch_output.clone(),
                        &matchlist,
                    )
                    .ok();
                }

                consume_query_by_gather(
                    query_name.clone(),
                    query_filename.clone(),
                    KmerMinHash::from(query_mh),
                    scaled as u32,
                    matchlist,
                    threshold_hashes,
                    Some(send.clone()), // is this clone ok?
                )
                .ok();
            }
            Err(err) => eprintln!("Error while processing record: {:?}", err),
        }
    }
    Ok(())
}

pub(crate) fn fastmultigather_rocksdb_obj(
    query_collection: &MultiCollection,
    db: &RevIndex,
    selection: &Selection,
    threshold_bp: u32,
    output: Option<String>,
) -> Result<(usize, usize, usize)> {
    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<BranchwaterGatherResult>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all records,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);
    let failed_gathers = AtomicUsize::new(0);

    let send = query_collection
        .par_iter()
        .filter_map(|(coll, _idx, record)| {
            let threshold = threshold_bp / selection.scaled().expect("scaled is not set!?");
            let ksize = selection.ksize().expect("ksize not set!?");

            // query downsampling happens here
            match coll.sig_from_record(record) {
                Ok(query_sig) => {
                    let query_filename = query_sig.filename();
                    let query_name = query_sig.name();
                    let query_md5 = query_sig.md5sum();

                    let mut results = vec![];
                    if let Ok(query_mh) = <SigStore as TryInto<KmerMinHash>>::try_into(query_sig) {
                        let _ = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
                        // Gather!
                        let (counter, query_colors, hash_to_color) =
                            db.prepare_gather_counters(&query_mh);

                        let matches = db.gather(
                            counter,
                            query_colors,
                            hash_to_color,
                            threshold as usize,
                            &query_mh,
                            Some(selection.clone()),
                        );
                        if let Ok(matches) = matches {
                            for match_ in &matches {
                                results.push(BranchwaterGatherResult {
                                    intersect_bp: match_.intersect_bp(),
                                    f_orig_query: match_.f_orig_query(),
                                    f_match: match_.f_match(),
                                    f_unique_to_query: match_.f_unique_to_query(),
                                    f_unique_weighted: match_.f_unique_weighted(),
                                    average_abund: match_.average_abund(),
                                    median_abund: match_.median_abund(),
                                    std_abund: match_.std_abund(),
                                    match_filename: match_.filename().clone(),
                                    match_name: match_.name().clone(),
                                    match_md5: match_.md5().clone(),
                                    f_match_orig: match_.f_match_orig(),
                                    unique_intersect_bp: match_.unique_intersect_bp(),
                                    gather_result_rank: match_.gather_result_rank(),
                                    remaining_bp: match_.remaining_bp(),
                                    query_filename: query_filename.clone(),
                                    query_name: query_name.clone(),
                                    query_md5: query_md5.clone(),
                                    query_bp: query_mh.n_unique_kmers(),
                                    ksize: ksize as u16,
                                    moltype: query_mh.hash_function().to_string(),
                                    scaled: query_mh.scaled(),
                                    query_n_hashes: query_mh.size() as u64,
                                    query_abundance: query_mh.track_abundance(),
                                    query_containment_ani: match_.query_containment_ani(),
                                    match_containment_ani: match_.match_containment_ani(),
                                    average_containment_ani: match_.average_containment_ani(),
                                    max_containment_ani: match_.max_containment_ani(),
                                    n_unique_weighted_found: match_.n_unique_weighted_found(),
                                    sum_weighted_found: match_.sum_weighted_found(),
                                    total_weighted_hashes: match_.total_weighted_hashes(),

                                    query_containment_ani_ci_low: match_
                                        .query_containment_ani_ci_low(),
                                    query_containment_ani_ci_high: match_
                                        .query_containment_ani_ci_high(),
                                    match_containment_ani_ci_low: match_
                                        .match_containment_ani_ci_low(),
                                    match_containment_ani_ci_high: match_
                                        .match_containment_ani_ci_high(),
                                });
                            }
                        } else {
                            eprintln!("Error gathering matches: {:?}", matches.err());
                            let _ = failed_gathers.fetch_add(1, atomic::Ordering::SeqCst);
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            query_filename
                        );
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }

                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                }
                Err(err) => {
                    eprintln!("Error loading sketch: {}", err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            }
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    // do some cleanup and error handling -
    send.expect("Unable to send internal data");
    thrd.join().expect("Unable to join CSV writing thread.");

    // done!
    let n_processed: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);
    let failed_gathers = failed_gathers.load(atomic::Ordering::SeqCst);

    if n_processed == 0 {
        return Err(anyhow::anyhow!("no search sigs found!?"));
    }

    if failed_gathers > 0 {
        return Err(anyhow::anyhow!(
            "{} failed gathers. See error messages above.",
            failed_gathers
        ));
    }

    Ok((n_processed, skipped_paths, failed_paths))
}
