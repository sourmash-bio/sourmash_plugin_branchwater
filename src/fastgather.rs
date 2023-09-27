/// fastgather: Run gather with a query against a list of files.
use anyhow::Result;

use sourmash::signature::Signature;
use sourmash::sketch::Sketch;
use std::path::Path;

use crate::utils::{
    consume_query_by_gather, load_sigpaths_from_zip_or_pathlist, load_sketches_above_threshold,
    prepare_query, write_prefetch, ReportType,
};

pub fn fastgather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query_filename: P,
    matchlist_filename: P,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    template: Sketch,
    gather_output: Option<P>,
    prefetch_output: Option<P>,
) -> Result<()> {
    let location = query_filename.to_string();
    eprintln!("Loading query from '{}'", location);
    let query = {
        let sigs = Signature::from_path(query_filename)?;

        prepare_query(&sigs, &template, &location)
    };
    // did we find anything matching the desired template?
    let query = match query {
        Some(query) => query,
        None => bail!("No sketch found with scaled={}, k={}", scaled, ksize),
    };

    // build the list of paths to match against.
    eprintln!(
        "Loading matchlist from '{}'",
        matchlist_filename.as_ref().display()
    );

    let matchlist_filename = matchlist_filename.as_ref().to_string_lossy().to_string();
    let (matchlist_paths, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(matchlist_filename, &template, ReportType::Against)?;

    eprintln!("Loaded {} sig paths in matchlist", matchlist_paths.len());

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
    let result = load_sketches_above_threshold(
        matchlist_paths,
        &template,
        &query.minhash,
        threshold_hashes,
    )?;
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
        write_prefetch(&query, prefetch_output, &matchlist).ok();
    }

    // run the gather!
    consume_query_by_gather(query, matchlist, threshold_hashes, gather_output).ok();
    Ok(())
}
