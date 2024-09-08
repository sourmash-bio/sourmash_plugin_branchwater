use anyhow::Context;
use camino::Utf8PathBuf as PathBuf;
use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use sourmash::collection::{Collection, CollectionSet};
use crate::utils::{ load_collection, ReportType };
use crate::utils::multicollection::MultiCollection;

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: &Selection,
    output: P,
    colors: bool,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    let multi: MultiCollection = load_collection(&siglist, selection,
                                                 ReportType::General,
                                                 allow_failed_sigpaths).unwrap();
    if multi.len() == 1 || use_internal_storage {
        let mut collection: CollectionSet;
        if multi.len() == 1 {
            let coll: Collection = Collection::try_from(multi).unwrap();
            collection = coll.select(selection)?.try_into().unwrap();
        } else { // use_internal_storage
            // @CTB warn: loading all the things
            let coll = multi.load_all_sigs(selection).unwrap();
            // @CTB multiple selects...
             collection = coll.select(selection)?.try_into()?;
        }
        eprintln!("Indexing {} sketches.", collection.len());
        let mut index = RevIndex::create(output.as_ref(), collection, colors)?;

        if use_internal_storage {
            index.internalize_storage()?;
        }
        Ok(())
    } else {
        Err(anyhow::anyhow!("Signatures failed to load. Exiting.").into())
    }
}
