use anyhow::Result;
/// pairwise: massively parallel in-memory pairwise comparisons.
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::signature::SigsTrait;
use sourmash::sketch::Sketch;

use crate::c_ani::containment_to_distance;

use crate::utils::{load_sketches_from_zip_or_pathlist, ReportType};

/// Perform pairwise comparisons of all signatures in a list.
///
/// Note: this function loads all _signatures_ into memory.

pub fn pairwise<P: AsRef<Path>>(
    siglist: P,
    threshold: f64,
    template: Sketch,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all sigs into memory at once.
    let sigs = load_sketches_from_zip_or_pathlist(&siglist, &template, ReportType::Query)?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(&mut writer, "query_name,query_md5,match_name,match_md5,containment,max_containment,jaccard,intersect_hashes,ani").unwrap();
        for (query, query_md5, m, m_md5, cont, max_cont, jaccard, overlap, ani) in recv.into_iter()
        {
            writeln!(
                &mut writer,
                "\"{}\",{},\"{}\",{},{},{},{},{},{}",
                query, query_md5, m, m_md5, cont, max_cont, jaccard, overlap, ani
            )
            .ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all signature,
    // Results written to the writer thread above.

    let processed_cmp = AtomicUsize::new(0);

    sigs.par_iter().enumerate().for_each(|(i, q1)| {
        for q2 in &sigs[(i + 1)..] {
            let overlap = q1.minhash.count_common(&q2.minhash, false).unwrap() as f64;
            let query1_size = q1.minhash.size() as f64;
            let query2_size = q2.minhash.size() as f64;

            let k_size = q1.minhash.ksize() as u32;
            let scaled = q1.minhash.scaled() as u64;

            let containment_q1_in_q2 = overlap / query1_size;
            let containment_q2_in_q1 = overlap / query2_size;
            let max_containment = containment_q1_in_q2.max(containment_q2_in_q1);
            let jaccard = overlap / (query1_size + query2_size - overlap);

            // TODO: this is equivalent to query1_size and query2_size but un u64
            let sig_1_size = q1.minhash.size() as u64;
            let sig_2_size = q2.minhash.size() as u64;

            let ani_1_in_2 = 1.00
                - containment_to_distance(
                    containment_q1_in_q2,
                    k_size,
                    scaled,
                    Some(sig_2_size),
                    Some(sig_2_size * scaled),
                    Some(0.95),
                    Some(true),
                    Some(1e-3),
                )
                .unwrap()
                .point_estimate as f64;

            let ani_2_in_1 = 1.00
                - containment_to_distance(
                    containment_q2_in_q1,
                    k_size,
                    scaled,
                    Some(sig_1_size),
                    Some(sig_1_size * scaled),
                    Some(0.95),
                    Some(true),
                    Some(1e-3),
                )
                .unwrap()
                .point_estimate as f64;

            let ani = (ani_1_in_2 + ani_2_in_1) / 2.0;

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
                    ani,
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
