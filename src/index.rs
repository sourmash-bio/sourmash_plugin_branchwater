use sourmash::index::revindex::RevIndex;
use sourmash::prelude::*;
use std::path::Path;

use crate::utils::{load_collection, ReportType};

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: &Selection,
    output: P,
    save_paths: bool,
    colors: bool,
    allow_failed_sigpaths: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");

    let collection = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
    )?;

    RevIndex::create(
        output.as_ref(),
        collection.select(&selection)?.try_into()?,
        colors,
    )?;

    Ok(())
}
