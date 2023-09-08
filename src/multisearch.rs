/// multisearch: massively parallel in-memory sketch search.

use rayon::prelude::*;
use anyhow::Result;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;


use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

use crate::utils::{load_sketchlist_filenames, load_sketches, load_sketches_from_zip, report_on_sketch_loading};

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn multisearch<P: AsRef<Path>>(
    querylist: P,
    againstlist: P,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
// construct a MinHash template for loading.
let max_hash = max_hash_for_scaled(scaled as u64);
let template_mh = KmerMinHash::builder()
    .num(0u32)
    .ksize(ksize as u32)
    .max_hash(max_hash)
    .build();
let template = Sketch::MinHash(template_mh);

// Read in list of query paths.
eprintln!("Reading list of queries from: '{}'", querylist.as_ref().display());

// Load all queries into memory at once.
let result = if querylist.as_ref().extension().map(|ext| ext == "zip").unwrap_or(false) {
    load_sketches_from_zip(&querylist, &template)?
} else {
    let querylist_paths = load_sketchlist_filenames(&querylist)?;
    load_sketches(querylist_paths, &template)?
};

let (queries, skipped_paths, failed_paths) = result;
report_on_sketch_loading(&queries, skipped_paths, failed_paths, false)?;


// Read in list of against paths.
eprintln!("Reading list of against paths from: '{}'", againstlist.as_ref().display());

// Load all against sketches into memory at once.
let result = if againstlist.as_ref().extension().map(|ext| ext == "zip").unwrap_or(false) {
    load_sketches_from_zip(&querylist, &template)?
} else {
    let querylist_paths = load_sketchlist_filenames(&againstlist)?;
    load_sketches(querylist_paths, &template)?
};

let (against, skipped_paths, failed_paths) = result;
report_on_sketch_loading(&against, skipped_paths, failed_paths, true)?;
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
        writeln!(&mut writer, "\"{}\",{},\"{}\",{},{},{},{},{}",
                    query, query_md5, m, m_md5, cont, max_cont, jaccard, overlap).ok();
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

        // search for matches & save containment.
        for q in queries.iter() {
            let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 100000 == 0 {
                eprintln!("Processed {} comparisons", i);
            }

            let overlap = q.minhash.count_common(&target.minhash, false).unwrap() as f64;
            let query_size = q.minhash.size() as f64;
            let target_size = target.minhash.size() as f64;

            let containment_query_in_target = overlap / query_size;
            let containment_in_target = overlap / target_size;
            let max_containment = containment_query_in_target.max(containment_in_target);
            let jaccard = overlap / (target_size + query_size - overlap);

            if containment_query_in_target > threshold {
                results.push((q.name.clone(),
                                q.md5sum.clone(),
                                target.name.clone(),
                                target.md5sum.clone(),
                                containment_query_in_target,
                                max_containment,
                                jaccard,
                                overlap))
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