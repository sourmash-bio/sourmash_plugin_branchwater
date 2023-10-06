/// mastiff_manygather: mastiff-indexed version of fastmultigather.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::signature::Signature;
use sourmash::sketch::Sketch;
use std::path::Path;

use sourmash::index::revindex::RevIndex;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::utils::{
    is_revindex_database, load_sigpaths_from_zip_or_pathlist, prepare_query, ReportType,
};

pub fn mastiff_manygather<P: AsRef<Path>>(
    queries_file: P,
    index: P,
    template: Sketch,
    threshold_bp: usize,
    output: Option<P>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(index.as_ref()) {
        bail!(
            "'{}' is not a valid RevIndex database",
            index.as_ref().display()
        );
    }
    // Open database once
    let db = RevIndex::open(index.as_ref(), true);
    println!("Loaded DB");

    // Load query paths
    let queryfile_name = queries_file.as_ref().to_string_lossy().to_string();
    let (query_paths, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(&queries_file, &template, ReportType::Query)?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let out: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let thrd = std::thread::spawn(move || {
        let mut writer = BufWriter::new(out);
        writeln!(
            &mut writer,
            "query_name,query_md5,match_name,match_md5,f_match_query,intersect_bp"
        )
        .unwrap();
        for (query, query_md5, m, m_md5, f_match_query, intersect_bp) in recv.into_iter() {
            writeln!(
                &mut writer,
                "\"{}\",{},\"{}\",{},{},{}",
                query, query_md5, m, m_md5, f_match_query, intersect_bp
            )
            .ok();
        }
    });

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send = query_paths
        .par_iter()
        .filter_map(|filename| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // load query signature from path:
            match Signature::from_path(filename) {
                Ok(query_sig) => {
                    let location = filename.display().to_string();
                    if let Some(query) = prepare_query(&query_sig, &template, &location) {
                        // let query_size = query.minhash.size() as f64;
                        let threshold = threshold_bp / query.minhash.scaled() as usize;

                        // mastiff gather code
                        let (counter, query_colors, hash_to_color) =
                            db.prepare_gather_counters(&query.minhash);

                        let matches = db.gather(
                            counter,
                            query_colors,
                            hash_to_color,
                            threshold,
                            &query.minhash,
                            &template,
                        );

                        // extract matches from Result
                        if let Ok(matches) = matches {
                            for match_ in &matches {
                                results.push((
                                    query.name.clone(),
                                    query.md5sum.clone(),
                                    match_.name().clone(),
                                    match_.md5().clone(),
                                    match_.f_match(), // f_match_query
                                    match_.intersect_bp(),
                                )); // intersect_bp
                            }
                        } else {
                            eprintln!("Error gathering matches: {:?}", matches.err());
                        }
                    } else {
                        if !queryfile_name.ends_with(".zip") {
                            eprintln!(
                                "WARNING: no compatible sketches in path '{}'",
                                filename.display()
                            );
                        }
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                }
                Err(err) => {
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!(
                        "WARNING: could not load sketches from path '{}'",
                        filename.display()
                    );
                    None
                }
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
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} search sigs", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if skipped_paths > 0 {
        eprintln!(
            "WARNING: skipped {} query paths - no compatible signatures.",
            skipped_paths
        );
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} signature paths failed to load. See error messages above.",
            failed_paths
        );
    }

    Ok(())
}
