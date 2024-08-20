//! MultiCollection implementation to handle sketches coming from multiple files.

use rayon::prelude::*;

use anyhow::{anyhow, Context, Result};
use camino::Utf8Path as Path;
use log::debug;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::collection::{Collection, CollectionSet};
use sourmash::encodings::Idx;
use sourmash::errors::SourmashError;
use sourmash::manifest::{Manifest, Record};
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage};

/// A collection of sketches, potentially stored in multiple files.
pub struct MultiCollection {
    collections: Vec<Collection>,
}

impl MultiCollection {
    /// Build from a standalone manifest
    pub fn from_manifest(sigpath: &Path) -> Result<Self> {
        debug!("multi from manifest!");
        let file =
            File::open(sigpath).with_context(|| format!("Failed to open file: '{}'", sigpath))?;

        let reader = BufReader::new(file);
        let manifest = Manifest::from_reader(reader)
            .with_context(|| format!("Failed to read manifest from: '{}'", sigpath))?;

        if manifest.is_empty() {
            Err(anyhow!("could not read as manifest: '{}'", sigpath))
        } else {
            let coll = Collection::new(
                manifest,
                InnerStorage::new(
                    FSStorage::builder()
                        .fullpath("".into())
                        .subdir("".into())
                        .build(),
                ),
            );
            Ok(Self {
                collections: vec![coll],
            })
        }
    }

    /// Load a collection from a .zip file.
    pub fn from_zipfile(sigpath: &Path) -> Result<Self> {
        debug!("multi from zipfile!");
        match Collection::from_zipfile(sigpath) {
            Ok(collection) => Ok(Self {
                collections: vec![collection],
            }),
            Err(_) => bail!("failed to load zipfile: '{}'", sigpath),
        }
    }

    /// Load a collection from a RocksDB.
    pub fn from_rocksdb(sigpath: &Path) -> Result<Self> {
        debug!("multi from rocksdb!");
        match Collection::from_rocksdb(sigpath) {
            Ok(collection) => Ok(Self {
                collections: vec![collection],
            }),
            Err(_) => bail!("failed to load rocksdb: '{}'", sigpath),
        }
    }

    /// Load a collection from a list of paths.
    pub fn from_pathlist(sigpath: &Path) -> Result<(Self, usize)> {
        debug!("multi from pathlist!");
        let file = File::open(sigpath)
            .with_context(|| format!("Failed to open pathlist file: '{}'", sigpath))?;
        let reader = BufReader::new(file);

        // load list of paths
        let lines: Vec<_> = reader
            .lines()
            .filter_map(|line| match line {
                Ok(path) => Some(path),
                Err(_err) => None,
            })
            .collect();

        // load sketches from paths in parallel.
        let n_failed = AtomicUsize::new(0);
        let records: Vec<Record> = lines
            .par_iter()
            .filter_map(|path| match Signature::from_path(path) {
                Ok(signatures) => {
                    let recs: Vec<Record> = signatures
                        .into_iter()
                        .flat_map(|v| Record::from_sig(&v, path))
                        .collect();
                    Some(recs)
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!("WARNING: could not load sketches from path '{}'", path);
                    let _ = n_failed.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            })
            .flatten()
            .collect();

        if records.is_empty() {
            eprintln!("No valid signatures found in pathlist '{}'", sigpath);
        }

        let manifest: Manifest = records.into();
        let collection = Collection::new(
            manifest,
            InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        );
        let n_failed = n_failed.load(atomic::Ordering::SeqCst);

        Ok((
            Self {
                collections: vec![collection],
            },
            n_failed,
        ))
    }

    // Load from a sig file
    pub fn from_signature(sigpath: &Path) -> Result<Self> {
        debug!("multi from signature!");
        let signatures = Signature::from_path(sigpath)
            .with_context(|| format!("Failed to load signatures from: '{}'", sigpath))?;

        let coll = Collection::from_sigs(signatures).with_context(|| {
            format!(
                "Loaded signatures but failed to load as collection: '{}'",
                sigpath
            )
        })?;
        Ok(Self {
            collections: vec![coll],
        })
    }

    pub fn len(&self) -> usize {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        val
    }
    pub fn is_empty(&self) -> bool {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        if val > 0 {
            false
        } else {
            true
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Collection> {
        self.collections.iter()
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (&Collection, Idx, &Record)> {
        // @CTB convert away from vec
        let mut sketchinfo: Vec<(&Collection, Idx, &Record)> = vec![];
        for coll in self.collections.iter() {
            let mut si: Vec<_> = coll
                .par_iter()
                .map(|(_idx, record)| (coll, _idx, record))
                .collect();
            sketchinfo.append(&mut si);
        }
        sketchinfo.into_par_iter()
    }
}

impl Select for MultiCollection {
    fn select(mut self, selection: &Selection) -> Result<Self, SourmashError> {
        self.collections = self
            .collections
            .iter()
            .filter_map(|c| c.clone().select(selection).ok())
            .collect();
        Ok(self)
    }
}

impl TryFrom<MultiCollection> for CollectionSet {
    type Error = SourmashError;

    fn try_from(multi: MultiCollection) -> Result<CollectionSet, SourmashError> {
        let coll = multi.iter().next().unwrap().clone();
        let cs: CollectionSet = coll.try_into()?;
        Ok(cs)
    }
}

/// Track a name/minhash.
pub struct SmallSignature {
    pub collection: Collection,
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}
