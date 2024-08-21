use log::debug;
use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;

use crate::utils::{load_collection, ReportType};
use sourmash::collection::Collection;

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: &Selection,
    output: P,
    colors: bool,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    let multi = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
    )?;

    debug!(
        "loaded multicollection from '{}' with len {}",
        siglist,
        multi.len()
    );

    let sigs = multi.load_sigs()?; // @CTB load into memory :sob:
    let coll = Collection::from_sigs(sigs)?;

    let mut index = RevIndex::create(output.as_ref(), coll.select(selection)?.try_into()?, colors)?;

    if use_internal_storage {
        index.internalize_storage()?;
    }

    Ok(())
}
