use camino::Utf8PathBuf as PathBuf;
use sourmash::collection::Collection;
use sourmash::index::revindex::RevIndex;
use sourmash::manifest::Manifest;
use sourmash::prelude::*;
use sourmash::storage::{FSStorage, InnerStorage, ZipStorage};
use std::path::Path;

use crate::utils::load_sketchlist_filenames;

pub fn index<P: AsRef<Path>>(
    siglist: PathBuf,
    manifest: Option<P>,
    selection: Selection,
    output: P,
    save_paths: bool,
    colors: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

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
        let manifest = manifest.unwrap_or_else(|| {
            let sig_paths: Vec<_> = load_sketchlist_filenames(&siglist)
                .unwrap_or_else(|_| panic!("Error loading siglist"))
                .into_iter()
                .map(|v| PathBuf::from_path_buf(v).unwrap())
                .collect();
            sig_paths.as_slice().into()
        });
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
