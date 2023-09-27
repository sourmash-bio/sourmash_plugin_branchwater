//use sourmash::index::revindex::RevIndex;
use sourmash::collection::Collection;
use sourmash::index::revindex::RevIndex;
// use sourmash::index::revindex::{prepare_query, RevIndex, RevIndexOps};
use sourmash::manifest::Manifest;
use sourmash::prelude::*;
// use sourmash::signature::{Signature, SigsTrait};
use sourmash::storage::{FSStorage, InnerStorage, ZipStorage};
// use sourmash::sketch::Sketch;
use std::path::Path;
use camino::Utf8PathBuf as PathBuf;

use crate::utils::{load_sigpaths_from_zip_or_pathlist, ReportType};

pub fn index<P: AsRef<Path>>(
    siglist: PathBuf,
    // template: Sketch, 
    manifest: Option<P>,
    selection: Selection,
    output: P,
    save_paths: bool,
    colors: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    // let (index_sigs, _temp_dir) = load_sigpaths_from_zip_or_pathlist(&siglist)?;

    // // if index_sigs pathlist is empty, bail
    // if index_sigs.is_empty() {
    //     bail!("No signatures to index loaded, exiting.");
    // }

    // // Create or open the RevIndex database with the provided output path and colors flag
    // let db = RevIndex::create(output.as_ref(), colors);

    // // Index the signatures using the loaded template, threshold, and save_paths option
    // db.index(index_sigs, &template, 0.0, save_paths);

    let manifest = if let Some(m) = manifest {
        let rdr = std::fs::OpenOptions::new().read(true).open(m.as_ref())?;
        Some(Manifest::from_reader(rdr)?)
    } else {
        None
    };

    let collection = if matches!(&siglist.extension(), Some("zip")) {
        if let Some(m) = manifest {
            let storage = ZipStorage::from_file(siglist)?;
            Collection::new(m, InnerStorage::new(storage))
        } else {
            Collection::from_zipfile(siglist)?
        }
    } else {
        let manifest = manifest.ok_or_else(|| "Need a manifest")?;
        let storage = FSStorage::builder()
            .fullpath("".into())
            .subdir("".into())
            .build();
        Collection::new(manifest, InnerStorage::new(storage))
    };

    RevIndex::create(
        output.as_ref(),
        collection.select(&selection)?.try_into()?,
        colors,
    )?;

    Ok(())
}
