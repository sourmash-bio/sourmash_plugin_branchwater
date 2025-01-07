use crate::utils::is_revindex_database;
use anyhow::Result;

use sourmash::index::revindex::{RevIndex, RevIndexOps};

pub fn check(index: camino::Utf8PathBuf, quick: bool) -> Result<()> {
    if !is_revindex_database(&index) {
        bail!("'{}' is not a valid RevIndex database", index);
    }

    println!("Opening DB");
    let db = RevIndex::open(index, true, None)?;

    println!("Starting check");
    db.check(quick);

    println!("Finished check");
    Ok(())
}
