/// fastgather: Run gather with a query against a list of files.
use anyhow::Result;
use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::sketch::minhash::KmerMinHash;

use crate::utils::{
    consume_query_by_gather, load_collection, load_sketches_above_threshold, write_prefetch,
    ReportType,
};

#[allow(clippy::too_many_arguments)]
pub fn fastgather(
    query_filepath: String,
    against_filepath: String,
    threshold_bp: u64,
    selection: Selection,
    gather_output: Option<String>,
    prefetch_output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    let query_collection = load_collection(
        &query_filepath,
        &selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    if query_collection.len() != 1 {
        bail!(
            "Fastgather requires a single query sketch. Check input: '{:?}'",
            &query_filepath
        )
    }
    // get single query sig and minhash
    let query_sig = query_collection.get_first_sig().expect("no queries!?");

    let query_filename = query_sig.filename();
    let query_name = query_sig.name();
    let query_md5 = query_sig.md5sum();

    // clone here is necessary b/c we use full query_sig in consume_query_by_gather
    let query_sig_ds = query_sig.select(&selection)?; // downsample as needed.
    let query_mh: KmerMinHash = match query_sig_ds.try_into() {
        Ok(query_mh) => query_mh,
        Err(_) => {
            bail!("No query sketch matching selection parameters.");
        }
    };

    let mut against_selection = selection;
    let scaled = query_mh.scaled();
    against_selection.set_scaled(scaled);

    // load collection to match against.
    let against_collection = load_collection(
        &against_filepath,
        &against_selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes = {
        let x = threshold_bp / scaled as u64;
        if x > 0 {
            x
        } else {
            1
        }
    };

    eprintln!(
        "using threshold overlap: {} {}",
        threshold_hashes, threshold_bp
    );

    // load a set of sketches, filtering for those with overlaps > threshold
    let result = load_sketches_above_threshold(against_collection, &query_mh, threshold_hashes)?;
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
        write_prefetch(
            query_filename.clone(),
            query_name.clone(),
            query_md5,
            prefetch_output,
            &matchlist,
        )
        .ok();
    }

    // run the gather!
    consume_query_by_gather(
        query_name,
        query_filename,
        query_mh,
        scaled as u32,
        matchlist,
        threshold_hashes,
        gather_output,
    )
    .ok();
    Ok(())
}
