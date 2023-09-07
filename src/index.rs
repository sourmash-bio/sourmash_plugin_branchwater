use std::path::{Path, PathBuf};
use sourmash::sketch::Sketch;
use sourmash::index::revindex::RevIndex;


use crate::utils::{read_signatures_from_zip, load_sketchlist_filenames};

pub fn index<P: AsRef<Path>>(
    siglist: P,
    template: Sketch,
    output: P,
    save_paths: bool,
    colors: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut temp_dir = None;
    println!("Loading siglist");

    let index_sigs: Vec<PathBuf>;

    if siglist.as_ref().extension().map(|ext| ext == "zip").unwrap_or(false) {
        let (paths, tempdir) = read_signatures_from_zip(&siglist)?;
        temp_dir = Some(tempdir);
        index_sigs = paths;
    } else {
        index_sigs = load_sketchlist_filenames(&siglist)?;
    }

    // if index_sigs pathlist is empty, bail
    if index_sigs.is_empty() {
        bail!("No signatures to index loaded, exiting.");
    }

    eprintln!("Loaded {} sig paths in siglist", index_sigs.len());

    // Create or open the RevIndex database with the provided output path and colors flag
    let db = RevIndex::create(output.as_ref(), colors);

    // Index the signatures using the loaded template, threshold, and save_paths option
    db.index(index_sigs, &template, 0.0, save_paths);

    if let Some(temp_dir) = temp_dir {
        temp_dir.close()?;
    }

    Ok(())
}