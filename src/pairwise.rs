use anyhow::Result;
/// pairwise: massively parallel in-memory sketch search.
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::signature::SigsTrait;
use sourmash::sketch::Sketch;

use crate::utils::{load_sketches_from_zip_or_pathlist, ReportType};

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn pairwise<P: AsRef<Path>>(
    siglist: P,
    threshold: f64,
    template: Sketch,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all queries into memory at once.
    let queries = load_sketches_from_zip_or_pathlist(&siglist, &template, ReportType::Query)?;

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
    // Main loop: iterate (in parallel) over all pairs of search signature paths,
    // loading them individually and comparing them. Stuff results into
    // the writer thread above.
    //

    let processed_cmp = AtomicUsize::new(0);

    queries.par_iter().enumerate().for_each(|(i, q1)| {
        for q2 in &queries[(i + 1)..] {
            let overlap = q1.minhash.count_common(&q2.minhash, false).unwrap() as f64;
            let query1_size = q1.minhash.size() as f64;
            let query2_size = q2.minhash.size() as f64;

            let containment_q1_in_q2 = overlap / query1_size;
            let containment_q2_in_q1 = overlap / query2_size;
            let max_containment = containment_q1_in_q2.max(containment_q2_in_q1);
            let jaccard = overlap / (query1_size + query2_size - overlap);

            if containment_q1_in_q2 > threshold || containment_q2_in_q1 > threshold {
                send.send((
                    q1.name.clone(),
                    q1.md5sum.clone(),
                    q2.name.clone(),
                    q2.md5sum.clone(),
                    containment_q1_in_q2,
                    max_containment,
                    jaccard,
                    overlap,
                )).unwrap();
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
