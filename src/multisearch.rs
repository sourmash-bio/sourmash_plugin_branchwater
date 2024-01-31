use anyhow::Result;
/// multisearch: massively parallel in-memory sketch search.
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::Sketch;
use sourmash::storage::SigStore;

use crate::utils::{load_collection, ReportType};

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn multisearch(
    query_filepath: &camino::Utf8PathBuf,
    against_filepath: &camino::Utf8PathBuf,
    threshold: f64,
    selection: &Selection,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all queries into memory at once.

    // let queries = load_sketches_from_zip_or_pathlist(&querylist, &template, ReportType::Query)?;
    let query_collection = load_collection(query_filepath, selection, ReportType::Query)?;
    let mut queries: Vec<SigStore> = vec![];
    for (idx, record) in query_collection.iter() {
        if let Ok(sig) = query_collection.sig_from_record(record)
        // .unwrap()
        // .select(&selection) // if we select here, we downsample and the md5sum changes!
        // ...which means we would lose the original md5sum that is used in the standard gather results.
        {
            queries.push(sig);
        } else {
            eprintln!("Failed to load 'against' record: {}", record.name());
        }
    }

    // Load all against sketches into memory at once.
    // let against = load_sketches_from_zip_or_pathlist(&againstlist, &template, ReportType::Against)?;
    let against_collection = load_collection(against_filepath, selection, ReportType::Against)?;
    let mut against: Vec<SigStore> = vec![];

    for (idx, record) in against_collection.iter() {
        if let Ok(sig) = against_collection.sig_from_record(record)
        // .unwrap()
        // .select(&selection) // if we select here, we downsample and the md5sum changes!
        // ...which means we would lose the original md5sum that is used in the standard gather results.
        {
            against.push(sig);
        } else {
            eprintln!("Failed to load 'against' record: {}", record.name());
        }
    }

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query_name,query_md5,match_name,match_md5,containment,max_containment,jaccard,intersect_hashes").unwrap();
        for (query, query_md5, m, m_md5, cont, max_cont, jaccard, overlap) in recv.into_iter() {
            writeln!(
                &mut writer,
                "\"{}\",{},\"{}\",{},{},{},{},{}",
                query, query_md5, m, m_md5, cont, max_cont, jaccard, overlap
            )
            .ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_cmp = AtomicUsize::new(0);

    let send = against
        .par_iter()
        .filter_map(|target| {
            let mut results = vec![];

            let ds_against_sig = target.clone().select(&selection).unwrap();
            if let Some(against_mh) = ds_against_sig.minhash() {
                // search for matches & save containment.
                for query_sig in queries.iter() {
                    let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
                    if i % 100000 == 0 {
                        eprintln!("Processed {} comparisons", i);
                    }
                    let ds_q = query_sig.clone().select(&selection).unwrap();
                    let query_mh = ds_q.minhash()?;
                    let overlap = query_mh.count_common(&against_mh, false).unwrap() as f64;
                    // use downsampled sizes
                    let query_size = query_mh.size() as f64;
                    let target_size = against_mh.size() as f64;

                    let containment_query_in_target = overlap / query_size;
                    let containment_in_target = overlap / target_size;
                    let max_containment = containment_query_in_target.max(containment_in_target);
                    let jaccard = overlap / (target_size + query_size - overlap);

                    if containment_query_in_target > threshold {
                        results.push((
                            query_sig.name(),
                            query_sig.md5sum(),
                            target.name(),
                            target.md5sum(),
                            containment_query_in_target,
                            max_containment,
                            jaccard,
                            overlap,
                        ))
                    }
                }
                if results.is_empty() {
                    None
                } else {
                    Some(results)
                }
            } else {
                None
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
