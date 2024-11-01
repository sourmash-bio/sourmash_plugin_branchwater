use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;

use crate::utils::{load_collection, ReportType};
use sourmash::collection::{Collection, CollectionSet};

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: Selection,
    output: P,
    colors: bool,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!("Loading sketches from {}", siglist);

    let multi = match load_collection(
        &siglist,
        &selection,
        ReportType::General,
        allow_failed_sigpaths,
    ) {
        Ok(multi) => multi,
        Err(err) => return Err(err.into()),
    };
    eprintln!("Found {} sketches total.", multi.len());

    // Try to convert it into a Collection and then CollectionSet.
    let collection = match Collection::try_from(multi.clone()) {
        // conversion worked!
        Ok(c) => {
            let cs: CollectionSet = c.select(&selection)?.try_into()?;
            Ok(cs)
        }
        // conversion failed; can we/should we load it into memory?
        Err(_) => {
            if use_internal_storage {
                eprintln!("WARNING: loading all sketches into memory in order to index.");
                eprintln!("See 'index' documentation for details.");
                let c: Collection = multi.load_all_sigs(&selection)?;
                let cs: CollectionSet = c.try_into()?;
                Ok(cs)
            } else {
                Err(
                    anyhow::anyhow!("cannot index this type of collection with external storage")
                        .into(),
                )
            }
        }
    };

    match collection {
        Ok(collection) => {
            eprintln!("Indexing {} sketches.", collection.len());
            let mut index = RevIndex::create(output.as_ref(), collection, colors)?;

            if use_internal_storage {
                index.internalize_storage()?;
            }
            Ok(())
        }
        Err(e) => Err(e),
    }
}
