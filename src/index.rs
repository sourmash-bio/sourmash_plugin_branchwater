use log::debug;
use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;
use std::fs::File;
use std::io::{ BufRead, BufReader };
use anyhow::Context;
use camino::Utf8PathBuf as PathBuf;

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

    let collection = match siglist {
        x if x.ends_with(".zip") => {
            Collection::from_zipfile(x)?
        }
        _ => {
            let file = File::open(siglist.clone()).with_context(|| {
                format!(
                    "Failed to open pathlist file: '{}'",
                    siglist
                )
            })?;

            let reader = BufReader::new(file);

            // load list of paths
            let lines: Vec<_> = reader
                .lines()
                .filter_map(|line| match line {
                    Ok(path) => {
                        let mut filename = PathBuf::new();
                        filename.push(path);
                        Some(filename)
                    }
                    Err(_err) => None,
                })
                .collect();

            Collection::from_paths(&lines)?
        }
    };

    let mut index = RevIndex::create(output.as_ref(),
                                     collection.select(selection)?.try_into()?,
                                     colors)?;

    if use_internal_storage {
        index.internalize_storage()?;
    }

    Ok(())
}
