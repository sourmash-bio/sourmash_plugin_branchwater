/// fastgather: Run gather with a query against a list of files.
use anyhow::Result;

use sourmash::sketch::Sketch;
use sourmash::signature::Signature;
use sourmash::selection::Selection;
use camino;
use std::collections::BinaryHeap;
use crate::utils::PrefetchResult;

use crate::utils::{
    consume_query_by_gather, load_sketches_above_threshold, write_prefetch, ReportType, load_collection
};

pub fn fastgather(
    query_filepath: camino::Utf8PathBuf,
    against_filepath: camino::Utf8PathBuf,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    selection: &Selection,
    gather_output: Option<String>,
    prefetch_output: Option<String>,
) -> Result<()> {

    let query_collection = load_collection(&query_filepath, selection, ReportType::Query)?;

    if query_collection.len() > 1 {
        bail!("Found more than one compatible sketch from '{}'. Fastgather requires a single query sketch.", &query_filepath)
    }

    // build the list of paths to match against.
    eprintln!("Loading matchlist from '{}'", against_filepath);
    let against_collection = load_collection(&against_filepath, selection, ReportType::Against)?; 
    eprintln!("Loaded {} sig paths in matchlist", against_collection.len());
    
    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes: u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }
    .try_into()?;

    eprintln!(
        "using threshold overlap: {} {}",
        threshold_hashes, threshold_bp
    );
    query_collection.iter().for_each(|(idx, record)| {
        // Load query sig
        match query_collection.sig_for_dataset(idx) {
            Ok(query_sig) => {
                let location = query_sig.filename();
                let mut matchlist: BinaryHeap<PrefetchResult> = BinaryHeap::new();
                let mut skipped_paths = 0;
                let mut failed_paths = 0;
    
                for sketch in query_sig.iter() {
                    // Access query MinHash
                    if let Sketch::MinHash(query) = sketch {
                        let result = load_sketches_above_threshold(
                            against_collection,
                            &selection,
                            &query,
                            threshold_hashes,
                        );
    
                        match result {
                            Ok((loaded_matchlist, skipped, failed)) => {
                                matchlist.extend(loaded_matchlist);
                                skipped_paths += skipped;
                                failed_paths += failed;
                            }
                            Err(err) => {
                                eprintln!("Error loading sketches: {:?}", err);
                                failed_paths += 1;
                            }
                        }
                    }
                }
    
                if skipped_paths > 0 {
                    eprintln!(
                        "WARNING: Skipped {} search paths - no compatible signatures.",
                        skipped_paths
                    );
                }
                if failed_paths > 0 {
                    eprintln!(
                        "WARNING: {} search paths failed to load. See error messages above.",
                        failed_paths
                    );
                }
    
                if matchlist.is_empty() {
                    eprintln!("No search signatures loaded for '{}', exiting.", location);
                    return; // Return early if no search signatures loaded
                }
    
                if let Some(prefetch_output) = &prefetch_output {
                    write_prefetch(&query_sig, Some(prefetch_output.clone()), &matchlist).ok();
                }
    
                // Run the gather!
                if let Some(gather_output) = &gather_output {
                    if let Err(err) = consume_query_by_gather(query_sig, matchlist, threshold_hashes, Some(gather_output)) {
                        eprintln!("Error during gather: {:?}", err);
                    }
                }
            }
            Err(_) => {
                eprintln!("WARNING: Could not load query sketch '{}'", record.internal_location());
            }
        }
    });
    Ok(())
}

//     query_collection.iter().for_each(|(idx, record)| {
//         // Load query sig
//         match query_collection.sig_for_dataset(idx) {
//             Ok(query_sig) => {
//                 let location = query_sig.filename();
//                 for sketch in query_sig.iter() {
//                     // Access query MinHash
//                     if let Sketch::MinHash(query) = sketch {
//                         let matchlist: BinaryHeap<PrefetchResult> = sketchlist
//                             .par_iter()
//                             .filter_map(|sm| {
//                                 // Call a function to load sketches above threshold
//                                 let result = load_sketches_above_threshold(
//                                     against_collection,
//                                     &selection,
//                                     &query,
//                                     threshold_hashes,
//                                 )?;
//                             let matchlist = result.0;
//                             let skipped_paths = result.1;
//                             let failed_paths = result.2;
                        
//                             if skipped_paths > 0 {
//                                 eprintln!(
//                                     "WARNING: skipped {} search paths - no compatible signatures.",
//                                     skipped_paths
//                                 );
//                             }
//                             if failed_paths > 0 {
//                                 eprintln!(
//                                     "WARNING: {} search paths failed to load. See error messages above.",
//                                     failed_paths
//                                 );
//                             }
                        
//                             if matchlist.is_empty() {
//                                 eprintln!("No search signatures loaded, exiting.");
//                                 return Ok(());
//                             }
                        
//                             if prefetch_output.is_some() {
//                                 write_prefetch(&query_sig, prefetch_output, &matchlist).ok();
//                             }
                        
//                             // run the gather!
//                             consume_query_by_gather(query_sig, matchlist, threshold_hashes, gather_output).ok();
//                             });
//                         }
//                     }
//                 }
//             }
//             Err(_) => {
//                 eprintln!("WARNING: Could not load query sketch '{}'", record.internal_location());
//             }
//         }
//     });
//     Ok(())
// }
