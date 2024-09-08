use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;

use sourmash::Error;
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
    eprintln!("loaded - {}", multi.len());

    let coll: Result<Collection, &str> = multi.clone().try_into();

    let collection = match coll {
        // if we can convert it, we have a single Collection; use that!
        Ok(coll) => {
            Ok(CollectionSet::from(coll.select(selection).unwrap().try_into()?))
        },
        // alt, our only chance is to load everything into memory.
        Err(_) => {
            if use_internal_storage {
                // @CTB warn: loading all the things
                let coll = multi.load_all_sigs(selection).unwrap();
                // @CTB multiple selects...
                Ok(CollectionSet::from(coll.select(selection).unwrap().try_into()?))
            } else {
                Err(anyhow::anyhow!("failed. Exiting.").into())
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
        },
        Err(e) => Err(e),
    }
}
