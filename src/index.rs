use sourmash::index::revindex::RevIndex;
use sourmash::sketch::Sketch;
use std::path::Path;

use crate::utils::{load_sigpaths_from_zip_or_pathlist, ReportType};

pub fn index<P: AsRef<Path>>(
    siglist: P,
    template: Sketch,
    output: P,
    save_paths: bool,
    colors: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    let (index_sigs, _temp_dir) =
        load_sigpaths_from_zip_or_pathlist(&siglist, &template, ReportType::Index)?;

    // if index_sigs pathlist is empty, bail
    if index_sigs.is_empty() {
        bail!("No signatures to index loaded, exiting.");
    }

    // Create or open the RevIndex database with the provided output path and colors flag
    let db = RevIndex::create(output.as_ref(), colors);

    // Index the signatures using the loaded template, threshold, and save_paths option
    db.index(index_sigs, &template, 0.0, save_paths);

    Ok(())
}
