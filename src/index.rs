use sourmash::index::revindex::RevIndex;
use sourmash::prelude::*;
use sourmash::index::revindex::RevIndexOps;
use std::path::Path;

use crate::utils::{load_collection, ReportType};

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: &Selection,
    output: P,
    colors: bool,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    let collection = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
    )?;

    let mut index = RevIndex::create(
        output.as_ref(),
        collection.select(selection)?.try_into()?,
        colors,
    )?;

    if use_internal_storage {
        index.internalize_storage()?;
    }

    Ok(())
}
