//! MultiCollection implementation to handle sketches coming from multiple files.

use rayon::prelude::*;

use anyhow::{anyhow, Context, Result};
use camino::Utf8Path as Path;
use log::debug;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::collections::HashSet;

use sourmash::collection::{Collection, CollectionSet};
use sourmash::encodings::Idx;
use sourmash::errors::SourmashError;
use sourmash::manifest::{Manifest, Record};
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};

/// A collection of sketches, potentially stored in multiple files.
pub struct MultiCollection {
    collections: Vec<Collection>,
}

impl MultiCollection {
    fn new(collections: Vec<Collection>) -> Self {
        Self { collections }
    }

    // Turn a set of paths into list of Collections.
    fn load_set_of_paths(paths: HashSet<String>) -> (Vec<Collection>, usize) {
        let n_failed = AtomicUsize::new(0);

        let colls: Vec<_> = paths
            .par_iter()
            .filter_map(|iloc| match iloc {
                // could just use a variant of load_collection here?
                x if x.ends_with(".zip") => {
                    debug!("loading sigs from zipfile {}", x);
                    Some(Collection::from_zipfile(x).unwrap())
                },
                _ => {
                    debug!("loading sigs from sigfile {}", iloc);
                    let signatures = match Signature::from_path(iloc) {
                        Ok(signatures) => Some(signatures),
                        Err(err) => {
                            eprintln!("Sketch loading error: {}", err);
                            None
                        }
                    };

                    match signatures {
                        Some(signatures) => {
                            Some(Collection::from_sigs(signatures).unwrap())
                        },
                        None => {
                            eprintln!("WARNING: could not load sketches from path '{}'", iloc);
                            let _ = n_failed.fetch_add(1, atomic::Ordering::SeqCst);
                            None
                        }
                    }
                }
            })
            .collect();

        let n_failed = n_failed.load(atomic::Ordering::SeqCst);
        (colls, n_failed)
    }

    /// Build from a standalone manifest
    pub fn from_standalone_manifest(sigpath: &Path) -> Result<Self> {
        debug!("multi from standalone manifest!");
        let file =
            File::open(sigpath).with_context(|| format!("Failed to open file: '{}'", sigpath))?;

        let reader = BufReader::new(file);
        let manifest = Manifest::from_reader(reader)
            .with_context(|| format!("Failed to read manifest from: '{}'", sigpath))?;

        if manifest.is_empty() {
            Err(anyhow!("could not read as manifest: '{}'", sigpath))
        } else {
            let ilocs: HashSet<_> = manifest
                .internal_locations()
                .map(|s| String::from(s))
                .collect();

            let (colls, _n_failed) = MultiCollection::load_set_of_paths(ilocs);
            let colls = colls.into_iter().collect();

            Ok(MultiCollection::new(colls))
        }
    }

    /// Load a collection from a .zip file.
    pub fn from_zipfile(sigpath: &Path) -> Result<Self> {
        debug!("multi from zipfile!");
        match Collection::from_zipfile(sigpath) {
            Ok(collection) => Ok(MultiCollection::new(vec![collection])),
            Err(_) => bail!("failed to load zipfile: '{}'", sigpath),
        }
    }

    /// Load a collection from a RocksDB.
    pub fn from_rocksdb(sigpath: &Path) -> Result<Self> {
        debug!("multi from rocksdb!");
        match Collection::from_rocksdb(sigpath) {
            Ok(collection) => Ok(MultiCollection::new(vec![collection])),
            Err(_) => bail!("failed to load rocksdb: '{}'", sigpath),
        }
    }

    /// Load a collection from a list of paths.
    pub fn from_pathlist(sigpath: &Path) -> Result<(Self, usize)> {
        debug!("multi from pathlist!");
        let file = File::open(sigpath)
            .with_context(|| format!("Failed to open pathlist file: '{}'", sigpath))?;
        let reader = BufReader::new(file);

        // load set of paths
        let lines: HashSet<_> = reader
            .lines()
            .filter_map(|line| match line {
                Ok(path) => Some(path),
                Err(_err) => None,
            })
            .collect();

        let num_to_load = lines.len();

        let (colls, n_failed) = MultiCollection::load_set_of_paths(lines);
        let colls: Vec<_> = colls.into_iter().collect();

        let n_missing = num_to_load - colls.len();

        Ok((MultiCollection::new(colls), n_missing))
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
        Ok(MultiCollection::new(vec![coll]))
    }

    pub fn len(&self) -> usize {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        val
    }
    pub fn is_empty(&self) -> bool {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        val == 0
    }

    pub fn iter(&self) -> impl Iterator<Item = &Collection> {
        self.collections.iter()
    }

    // iterate over tuples
    pub fn item_iter(&self) -> impl Iterator<Item = (&Collection, Idx, &Record)> {
        // CTB: request review by Rust expert pls :). Does this make
        // unnecessary copies??
        let s: Vec<_> = self
            .iter()
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c, _idx, record)))
            .collect();
        s.into_iter()
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (&Collection, Idx, &Record)> {
        // CTB: request review by Rust expert - why can't I use item_iter here?
        // i.e. self.item_iter().into_par_iter()?
        let s: Vec<_> = self
            .iter()
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c, _idx, record)))
            .collect();
        s.into_par_iter()
    }

    pub fn get_first_sig(&self) -> Option<SigStore> {
        if !self.is_empty() {
            let query_item = self.item_iter().next().unwrap();
            let (coll, _, _) = query_item;
            Some(coll.sig_for_dataset(0).ok()?)
        } else {
            None
        }
    }
}

impl Select for MultiCollection {
    fn select(mut self, selection: &Selection) -> Result<Self, SourmashError> {
        // CTB: request review by Rust expert! Is the clone necessary?
        self.collections = self
            .iter()
            .filter_map(|c| c.clone().select(selection).ok())
            .collect();
        Ok(self)
    }
}

impl TryFrom<MultiCollection> for CollectionSet {
    type Error = SourmashError;

    fn try_from(multi: MultiCollection) -> Result<CollectionSet, SourmashError> {
        // CTB: request review by Rust expert! Is the clone necessary?
        let coll = multi.iter().next().unwrap().clone();
        let cs: CollectionSet = coll.try_into()?;
        Ok(cs)
    }
}

/// Track a name/minhash.
pub struct SmallSignature {
    // CTB: request help - can we/should we use references & lifetimes here?
    pub collection: Collection,
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}