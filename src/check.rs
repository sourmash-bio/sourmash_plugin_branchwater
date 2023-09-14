use std::path::Path;

use crate::utils::is_revindex_database;

use sourmash::index::revindex::RevIndex;

pub fn check<P: AsRef<Path>>(index: P, quick: bool) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(index.as_ref()) {
        bail!(
            "'{}' is not a valid RevIndex database",
            index.as_ref().display()
        );
    }

    println!("Opening DB");
    let db = RevIndex::open(index.as_ref(), true);

    println!("Starting check");
    db.check(quick);

    println!("Finished check");
    Ok(())
}
