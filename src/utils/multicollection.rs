//! MultiCollection implementation to handle sketches coming from multiple files.

use rayon::prelude::*;
use sourmash::prelude::*;

use anyhow::{anyhow, Context, Result};
use camino::Utf8Path as Path;
use camino::Utf8PathBuf;
use log::{debug, trace};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::collection::Collection;
use sourmash::encodings::Idx;
use sourmash::errors::SourmashError;
use sourmash::manifest::{Manifest, Record};
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};
use sourmash::ScaledType;

/// A collection of sketches, potentially stored in multiple files.
#[derive(Clone)]
pub struct MultiCollection {
    collections: Vec<Collection>,
    pub contains_revindex: bool, // track whether one or more Collection is a RevIndex
}

impl MultiCollection {
    fn new(collections: Vec<Collection>, contains_revindex: bool) -> Self {
        Self {
            collections,
            contains_revindex,
        }
    }

    // Try loading a set of paths as JSON files only. Fails on any Err.
    //
    // This is a legacy method that supports pathlists for
    // 'index'. See sourmash-bio/sourmash#3321 for background.
    //
    // Use load_set_of_paths for full generality!
    //
    // CTB NOTE: this could potentially have very poor performance if
    // there are a lot of _good_ files, with one _bad_ one. Look into
    // exiting first loop early.
    fn load_set_of_json_files(paths: &HashSet<String>) -> Result<MultiCollection> {
        // load sketches from paths in parallel.
        let n_failed = AtomicUsize::new(0);
        let records: Vec<Record> = paths
            .par_iter()
            .filter_map(|path| match Signature::from_path(path) {
                Ok(signatures) => {
                    let recs: Vec<Record> = signatures
                        .into_iter()
                        .flat_map(|v| Record::from_sig(&v, path))
                        .collect();
                    Some(recs)
                }
                Err(_) => {
                    let _ = n_failed.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            })
            .flatten()
            .collect();

        let n_failed = n_failed.load(atomic::Ordering::SeqCst);

        if records.is_empty() || n_failed > 0 {
            return Err(anyhow!("cannot load everything as JSON files"));
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
        Ok(MultiCollection::from(collection))
    }

    // Turn a set of paths into list of Collections - works recursively
    // if needed, and can handle paths of any supported type.
    fn load_set_of_paths(paths: &HashSet<String>) -> (MultiCollection, usize) {
        let n_failed = AtomicUsize::new(0);

        // could just use a variant of load_collection here?
        let colls: Vec<MultiCollection> = paths
            .par_iter()
            .filter_map(|iloc| match iloc {
                // load from zipfile
                x if x.ends_with(".zip") => {
                    debug!("loading sigs from zipfile {}", x);
                    let coll = Collection::from_zipfile(x).expect("nothing to load!?");
                    Some(MultiCollection::from(coll))
                }
                // load from CSV
                x if x.ends_with(".csv") => {
                    debug!("vec from pathlist of standalone manifests!");

                    let x: String = x.into();
                    let utf_path: &Path = x.as_str().into();
                    MultiCollection::from_standalone_manifest(utf_path).ok()
                }
                // load from (by default) a sigfile
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
                            let records: Vec<_> = signatures
                                .into_iter()
                                .flat_map(|v| Record::from_sig(&v, iloc))
                                .collect();

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
                            Some(MultiCollection::from(collection))
                        }
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
        (MultiCollection::from(colls), n_failed)
    }

    /// Build from a standalone manifest.  Note: the tricky bit here
    /// is that the manifest may select only a subset of the rows,
    /// using (name, md5) tuples.
    pub fn from_standalone_manifest(sigpath: &Path) -> Result<Self> {
        debug!("multi from standalone manifest!");
        let file =
            File::open(sigpath).with_context(|| format!("Failed to open file: '{}'", sigpath))?;

        let reader = BufReader::new(file);
        let manifest = Manifest::from_reader(reader)
            .with_context(|| format!("Failed to read manifest from: '{}'", sigpath))?;
        debug!("got {} records from standalone manifest", manifest.len());

        if manifest.is_empty() {
            Err(anyhow!("could not read as manifest: '{}'", sigpath))
        } else {
            let ilocs: HashSet<_> = manifest.internal_locations().map(String::from).collect();
            let (mut colls, _n_failed) = MultiCollection::load_set_of_paths(&ilocs);

            colls.intersect_manifest(&manifest);

            Ok(colls)
        }
    }

    /// Load a collection from a .zip file.
    pub fn from_zipfile(sigpath: &Path) -> Result<Self> {
        debug!("multi from zipfile!");
        match Collection::from_zipfile(sigpath) {
            Ok(collection) => Ok(MultiCollection::new(vec![collection], false)),
            Err(_) => bail!("failed to load zipfile: '{}'", sigpath),
        }
    }

    /// Load a collection from a RocksDB.
    pub fn from_rocksdb(sigpath: &Path) -> Result<Self> {
        debug!("multi from rocksdb!");
        // duplicate logic from is_revindex_database
        let path: Utf8PathBuf = sigpath.into();

        let mut is_rocksdb = false;

        if path.is_dir() {
            let current_file = path.join("CURRENT");
            if current_file.exists() && current_file.is_file() {
                is_rocksdb = true;
            }
        }

        if is_rocksdb {
            match Collection::from_rocksdb(sigpath) {
                Ok(collection) => {
                    debug!("...rocksdb successful!");
                    Ok(MultiCollection::new(vec![collection], true))
                }
                Err(_) => bail!("failed to load rocksdb: '{}'", sigpath),
            }
        } else {
            bail!("not a rocksdb: '{}'", sigpath)
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

        let val = MultiCollection::load_set_of_json_files(&lines);

        let (multi, n_failed) = match val {
            Ok(collection) => {
                eprintln!("SUCCEEDED in loading as JSON files, woot woot");
                // CTB note: if any path fails to load,
                // load_set_of_json_files returns Err.
                (collection, 0)
            }
            Err(_) => {
                eprintln!("FAILED to load as JSON files; falling back to general recursive");
                MultiCollection::load_set_of_paths(&lines)
            }
        };

        Ok((multi, n_failed))
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
        Ok(MultiCollection::new(vec![coll], false))
    }

    pub fn len(&self) -> usize {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        val
    }

    pub fn is_empty(&self) -> bool {
        let val: usize = self.collections.iter().map(|c| c.len()).sum();
        val == 0
    }

    pub fn max_scaled(&self) -> Option<&ScaledType> {
        self.item_iter().map(|(_, _, record)| record.scaled()).max()
    }

    // iterate over tuples
    pub fn item_iter(&self) -> impl Iterator<Item = (&Collection, Idx, &Record)> {
        let s: Vec<_> = self
            .collections
            .iter()
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c, _idx, record)))
            .collect();
        s.into_iter()
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (&Collection, Idx, &Record)> {
        // first create a Vec of all triples (Collection, Idx, Record)
        let s: Vec<_> = self
            .collections
            .iter() // CTB: are we loading things into memory here? No...
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c, _idx, record)))
            .collect();
        // then return a parallel iterator over the Vec.
        s.into_par_iter()
    }

    pub fn get_first_sig(&self) -> Option<SigStore> {
        if !self.is_empty() {
            let query_item = self.item_iter().next()?;
            let (coll, _, _) = query_item;
            Some(coll.sig_for_dataset(0).ok()?)
        } else {
            None
        }
    }

    // Load all sketches into memory, using SmallSignature to track original
    // signature metadata.
    pub fn load_sketches(self, selection: &Selection) -> Result<Vec<SmallSignature>> {
        if self.contains_revindex {
            eprintln!("WARNING: loading all sketches from a RocksDB into memory!");
        }
        let sketchinfo: Vec<_> = self
            .par_iter()
            .filter_map(|(coll, _idx, record)| match coll.sig_from_record(record) {
                Ok(sig) => {
                    trace!(
                        "MultiCollection load sketch: from:{} idx:{} loc:{}",
                        coll.storage().spec(),
                        _idx,
                        record.internal_location()
                    );

                    let sig_name = sig.name();
                    let sig_md5 = sig.md5sum();
                    let selected_sig = sig.select(selection).ok()?;
                    let minhash = selected_sig.try_into().expect("cannot extract sketch");

                    Some(SmallSignature {
                        location: record.internal_location().to_string(),
                        name: sig_name,
                        md5sum: sig_md5,
                        minhash,
                    })
                }
                Err(_) => {
                    eprintln!(
                        "FAILED to load sketch from '{}'",
                        record.internal_location()
                    );
                    None
                }
            })
            .collect();

        Ok(sketchinfo)
    }

    fn intersect_manifest(&mut self, manifest: &Manifest) {
        for coll in self.collections.iter_mut() {
            coll.intersect_manifest(manifest);
        }
    }

    // Load all sketches into memory, producing an in-memory Collection.
    pub fn load_all_sigs(self, selection: &Selection) -> Result<Collection> {
        let all_sigs: Vec<Signature> = self
            .par_iter()
            .filter_map(|(coll, _idx, record)| match coll.sig_from_record(record) {
                Ok(sig) => {
                    let sig = sig.select(selection).ok()?;
                    Some(Signature::from(sig))
                }
                Err(_) => {
                    eprintln!(
                        "FAILED to load sketch from '{}'",
                        record.internal_location()
                    );
                    None
                }
            })
            .collect();
        Ok(Collection::from_sigs(all_sigs)?)
    }
}

impl Select for MultiCollection {
    fn select(self, selection: &Selection) -> Result<Self, SourmashError> {
        let collections = self
            .collections
            .into_iter()
            .filter_map(|c| c.select(selection).ok())
            .collect();

        Ok(MultiCollection::new(collections, self.contains_revindex))
    }
}

// Convert a single Collection into a MultiCollection
impl From<Collection> for MultiCollection {
    fn from(coll: Collection) -> Self {
        // CTB: how can we check if revindex?
        MultiCollection::new(vec![coll], false)
    }
}

// Merge a bunch of MultiCollection structs into one
impl From<Vec<MultiCollection>> for MultiCollection {
    fn from(multi: Vec<MultiCollection>) -> Self {
        let mut x: Vec<Collection> = vec![];
        let mut contains_revindex = false;
        for mc in multi.into_iter() {
            for coll in mc.collections.into_iter() {
                x.push(coll);
            }
            contains_revindex = contains_revindex || mc.contains_revindex;
        }
        MultiCollection::new(x, contains_revindex)
    }
}

// Extract a single Collection from a MultiCollection, if possible
impl TryFrom<MultiCollection> for Collection {
    type Error = &'static str;

    fn try_from(multi: MultiCollection) -> Result<Self, Self::Error> {
        if multi.collections.len() == 1 {
            // this must succeed b/c len > 0
            Ok(multi.collections.into_iter().next().unwrap())
        } else {
            Err("More than one Collection in this MultiCollection; cannot convert")
        }
    }
}

/// Track a name/minhash.
pub struct SmallSignature {
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}
