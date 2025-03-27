use anyhow::Result;

use sourmash::index::revindex::disk_revindex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;

use crate::utils::MultiCollection;
use crate::utils::{load_collection, ReportType};
use sourmash::collection::{Collection, CollectionSet};

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: Selection,
    output: P,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<()> {
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

    index_obj(multi, output, use_internal_storage)
}

pub(crate) fn index_obj<P: AsRef<Path>>(
    multi: MultiCollection,
    output: P,
    use_internal_storage: bool,
) -> Result<()> {
    // Try to convert it into a Collection and then CollectionSet.
    let collection = match Collection::try_from(multi.clone()) {
        // conversion worked!
        Ok(coll) => {
            let cs: CollectionSet = coll.try_into()?;
            Ok(cs)
        }
        // conversion failed; can we/should we load it into memory?
        Err(_) => {
            if use_internal_storage {
                eprintln!("WARNING: loading all sketches into memory in order to index.");
                eprintln!("See 'index' documentation for details.");
                let c: Collection = multi.load_all_sigs()?;
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
            let mut index = disk_revindex::DiskRevIndex::create(output.as_ref(), collection)?;

            if use_internal_storage {
                eprintln!("Internalizing storage.");
                index.internalize_storage()?;
            } else {
                eprintln!("Using external storage - not copying sketches.");
            }
            Ok(())
        }
        Err(e) => Err(e),
    }
}
