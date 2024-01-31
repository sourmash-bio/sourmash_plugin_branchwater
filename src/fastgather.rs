/// fastgather: Run gather with a query against a list of files.
use anyhow::Result;

use sourmash::selection::Selection;
// use camino;

use sourmash::prelude::Select;
use sourmash::signature::SigsTrait;

use crate::utils::{
    consume_query_by_gather, load_collection, load_sketches_above_threshold, write_prefetch,
    ReportType,
};

pub fn fastgather(
    query_filepath: &camino::Utf8PathBuf,
    against_filepath: &camino::Utf8PathBuf,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    selection: &Selection,
    gather_output: Option<String>,
    prefetch_output: Option<String>,
) -> Result<()> {
    let query_collection = load_collection(query_filepath, selection, ReportType::Query)?;

    if query_collection.len() != 1 {
        bail!(
            "Fastgather requires a single query sketch. Check input: '{:?}'",
            &query_filepath
        )
    }
    // get single query sig and minhash
    let query_sig = query_collection.sig_for_dataset(0)?; // need original md5sum, etc
                                                          // downsample
    let query_sig_ds = query_sig.clone().select(selection)?;
    let query_mh = match query_sig_ds.minhash() {
        Some(query_mh) => query_mh,
        None => {
            bail!("No query sketch matching selection parameters.");
        }
    };
    // build the list of paths to match against.
    let against_collection = load_collection(against_filepath, selection, ReportType::Against)?;

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

    // load a set of sketches, filtering for those with overlaps > threshold
    let result =
        load_sketches_above_threshold(against_collection, &selection, &query_mh, threshold_hashes)?;
    let matchlist = result.0;
    let skipped_paths = result.1;
    let failed_paths = result.2;
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

    if matchlist.is_empty() {
        eprintln!("No search signatures loaded, exiting.");
        return Ok(());
    }

    if prefetch_output.is_some() {
        write_prefetch(&query_sig, prefetch_output, &matchlist).ok();
    }

    // run the gather!
    consume_query_by_gather(
        query_sig.clone(),
        matchlist,
        threshold_hashes,
        gather_output,
    )
    .ok();
    Ok(())
}
