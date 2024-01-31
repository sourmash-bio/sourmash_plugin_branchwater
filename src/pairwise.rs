use anyhow::Result;
/// pairwise: massively parallel in-memory pairwise comparisons.
use rayon::prelude::*;
use sourmash::sketch::minhash::KmerMinHash;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::signature::SigsTrait;
use sourmash::sketch::Sketch;

use crate::utils::{load_collection, ReportType};
use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::storage::SigStore;

/// Perform pairwise comparisons of all signatures in a list.
///
/// Note: this function loads all _signatures_ into memory.

pub fn pairwise<P: AsRef<Path>>(
    sigpath: &camino::Utf8PathBuf,
    threshold: f64,
    selection: &Selection,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all sigs into memory at once.
    let collection = load_collection(sigpath, selection, ReportType::Query)?;

    if collection.len() <= 1 {
        bail!(
            "Pairwise requires two or more sketches. Check input: '{:?}'",
            &sigpath
        )
    }

    let mut sketches: Vec<(KmerMinHash, String, String)> = Vec::new();
    for (_idx, record) in collection.iter() {
        if let Ok(sig) = collection.sig_from_record(record) {
            if let Some(ds_mh) = sig.clone().select(&selection)?.minhash().cloned() {
                sketches.push((ds_mh, record.name().to_string(), record.md5().to_string()));
            }
        } else {
            eprintln!("Failed to load record: {}", record.name());
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
                    send.send((
                        q1_name.clone(),
                        q1_md5.clone(),
                        q2_name.clone(),
                        q2_md5.clone(),
                        containment_q1_in_q2,
                        max_containment,
                        jaccard,
                        overlap,
                    ))
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
