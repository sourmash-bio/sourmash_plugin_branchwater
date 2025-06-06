//! MultiCollection implementation to handle sketches coming from multiple files.


// Challenges/design considerations:
// * we need to keep original manifest (or collection?) around to
//   support proper reporting. If Collection, this is _pre_ selection...
// * we need ~CollectionSet around for searching and indexing.
// * So, we need a struct that keeps original stuff around, as well as
//   CollectionSet/RevIndex. This latter stuff needs to be an enum so we
//   can track. May also want to track in-mem vs on-disk for CollectionSet.
// * Do we want to track this as an Option, as in, no select has been called?

// Things that seem OK:
// * prefetch creates a RevIndex + CounterGather + orig manifest object; the
//   RevIndex is needed for CounterGather to retrieve matches.
// * here, we suffer a DiskRevIndex.clone() if needed; CounterGather is
//   new.
// * Manifest could be/should be a reference, though.


// OR...
// * accept that MultiCollection and MultiCollectionSet are different things.
// * have members of MultiCollectionSet take a reference to MultiCollection?

// OR...
// * look into defining my own storage?

use rayon::prelude::*;
use sourmash::prelude::*;

use anyhow::{anyhow, Context, Result};
use camino::Utf8Path as Path;
use camino::Utf8PathBuf;
use log::{debug, trace};
use std::collections::HashSet;
use std::fs::{metadata, File};
use std::io::{BufRead, BufReader};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::collection::{ Collection, CollectionSet };
use sourmash::encodings::Idx;
use sourmash::errors::SourmashError;
use sourmash::index::linear::LinearIndex;
use sourmash::index::revindex::{ CounterGather, RevIndex, RevIndexOps };
use sourmash::index::revindex::mem_revindex::MemRevIndex;
use sourmash::manifest::{Manifest, Record};
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};
use sourmash::ScaledType;

pub struct PrefetchItem<'a> {
    pub revindex: RevIndex,     // collection for search results
    pub cg: CounterGather,      // search results
    pub mf: &'a Manifest,           // original manifest
}

impl PrefetchItem<'_> {
    pub fn found_hashes(&self, template_mh: &KmerMinHash) -> KmerMinHash {
        self.cg.found_hashes(template_mh)
    }
}


pub struct PrefetchContainer<'a> {
    pub matchlists: Vec<PrefetchItem<'a>>,
}

impl PrefetchContainer<'_> {
    // sum the lengths of the contents
    pub fn len(&self) -> usize {
        let mut l = 0;
        for item in self.matchlists.iter() {
            l += item.cg.len();
        }
        l
    }

    // true if all contained matchlists are false
    pub fn is_empty(&self) -> bool {
        for item in self.matchlists.iter() {
            if !item.cg.is_empty() { return false; };
        }
        return true;
    }

    // find best match across all CounterGather objects. Return
    // matching Signature + original Record.
    pub fn peek(&self, threshold_hashes: u64) -> Option<(Signature, &Record)> {
        let mut best_idx = None;
        let mut best_overlap = 0;
        let mut best_item = None;
        
        for item in self.matchlists.iter() {
            if !item.cg.is_empty() { // @CTB test else
                if let Some((idx, overlap)) = item.cg.peek(threshold_hashes as usize) {
                    if overlap > best_overlap {
                        best_idx = Some(idx);
                        best_item = Some(item);
                        best_overlap = overlap;
                    }
                }
            }
        }

        // did we find something?
        if let Some(item) = best_item {
            let idx = best_idx.unwrap();
            let match_sig: Signature = item.revindex.collection().sig_for_dataset(idx).expect("cannot load").into();
            let orig_record: &Record = item.mf.get_record(idx).expect("cannot retrieve original record!?");
            Some((match_sig, orig_record))
        } else {
            None
        }
    }

    // consume the best match across all CounterGathers.
    pub fn consume(self, intersect_mh: &KmerMinHash) -> Self {
        let mut updated: Vec<PrefetchItem> = vec![];
        
        for mut item in self.matchlists.into_iter() {
            item.cg.consume(intersect_mh);
            updated.push(item);
        }
        PrefetchContainer { matchlists : updated }
    }

    pub fn found_hashes(&self, template_mh: &KmerMinHash) -> Result<KmerMinHash> {
        let mut new_mh = template_mh.clone();
        new_mh.clear();

        for item in self.matchlists.iter() {
            let found = item.found_hashes(&template_mh);
            new_mh.merge(&found)?;
        }
        Ok(new_mh)
    }
}


trait Searchable {
    fn prefetch(&self, query: &KmerMinHash, threshold_hashes: u64) ->
        Result<(RevIndex, CounterGather, &Manifest, usize, usize)>;
    fn prefetch_iter(&self, query: &KmerMinHash, threshold_hashes: u64) ->
        impl Iterator<Item = Idx>;
    fn len(&self) -> usize;
    fn iter(&self) -> impl Iterator<Item = (Idx, &Record)>;
    fn collection(&self) -> &CollectionSet;
    fn intersect_manifest(&mut self, manifest: &Manifest);
    fn get_orig_manifest(&self) -> &Manifest;
    fn get_signature(&self, idx: Idx) -> Result<Signature>;
    fn get_record(&self, idx: Idx) -> Option<&Record>;
    fn load(self) -> Result<Self> where Self: Sized;
}


// @CTB do we want to split disk/mem collection?
#[derive(Clone)]
enum SearchContainer<'a> {
    InvertedIndex(RevIndex, &'a Manifest),
    LinearCollection(CollectionSet, &'a Manifest),
}

impl Searchable for SearchContainer<'_> {
    /// find all overlapping sketches; return RevIndex and CounterGather
    /// structs containing the results.
    ///
    /// For pre-existing RevIndexes, this will return a clone of the revindex;
    /// for other collections, it will return a MemRevIndex.

    fn prefetch(&self, query: &KmerMinHash, threshold_hashes: u64) ->
        Result<(RevIndex, CounterGather, &Manifest, usize, usize)> {
            match self {
                SearchContainer::InvertedIndex(revindex, mf) => {
                    let cg = revindex.prepare_gather_counters(query, None);
                    // @CTB clone
                    Ok((revindex.clone(), cg, mf, 0, 0))
                },
                SearchContainer::LinearCollection(coll, mf) => {
                    let (revindex, cg, skip, fail) = 
                        load_sketches_above_threshold_sigs_XXX(coll,
                                                               query,
                                                               threshold_hashes)?;
                    // @CTB clone
                    Ok((revindex, cg, mf, skip, fail))
                                                           
                },
            }
        }

    fn prefetch_iter(&self, query: &KmerMinHash, threshold_hashes: u64) -> impl Iterator<Item = Idx> {
        match self {
            SearchContainer::InvertedIndex(revindex, _mf) => {
                let counter = revindex.counter_for_query(&query, None);
                counter
                    .most_common()
                    .into_iter()
                    .filter_map(move |(dataset_id, size)| {
                        if size as u64 >= threshold_hashes {
                            Some(dataset_id)
                        } else {
                            None
                        }
                    })
            }
            SearchContainer::LinearCollection(_coll, _mf) => {
                panic!("foo");
            }
        }
    }

    fn len(&self) -> usize {
        match self {
            SearchContainer::InvertedIndex(revindex, _mf) => revindex.len(),
            SearchContainer::LinearCollection(coll, _mf) => coll.len(),
        }
    }

    /// iterate over (Idx, &Record)
    fn iter(&self) -> impl Iterator<Item = (Idx, &Record)> {
        match self {
            SearchContainer::InvertedIndex(revindex, _mf) => {
                revindex.collection().iter()
            }
            SearchContainer::LinearCollection(coll, _mf) => {
                coll.iter()
            }
        }
    }

    fn collection(&self) -> &CollectionSet {
        match self {
            SearchContainer::InvertedIndex(revindex, _mf) => {
                revindex.collection()
            }
            SearchContainer::LinearCollection(coll, _mf) => {
                coll
            }
        }
    }

    /// retrieve the original manifest for this SearchContainer, prior
    /// to any downsampling/indexing.
    fn get_orig_manifest(&self) -> &Manifest {
        match self {
            SearchContainer::InvertedIndex(_, mf) => {
                &mf
            }
            SearchContainer::LinearCollection(_, mf) => {
                &mf
            }
        }
    }

    // @CTB this will need to be updated in tricky ways.
    fn intersect_manifest(&mut self, manifest: &Manifest) {
        match self {
            SearchContainer::InvertedIndex(_revindex, _mf) => {
                panic!("foo 3");
            }
            SearchContainer::LinearCollection(coll, _mf) => {
                coll.intersect_manifest(manifest);
            }
        }
    }

    fn get_signature(&self, idx: Idx) -> Result<Signature> {
        let sig = self.collection().sig_for_dataset(idx)?;
        Ok(sig.into())
    }

    fn get_record(&self, idx: Idx) -> Option<&Record> {
        self.collection().manifest().get_record(idx)
    }

    fn load(self) -> Result<Self> {
        match self {
            // in an index? leave alone.
            SearchContainer::InvertedIndex(revindex, mf) =>
                Ok(SearchContainer::InvertedIndex(revindex, mf)),
            // load collection into memory.
            // @CTB what if it is already loaded??
            SearchContainer::LinearCollection(coll, mf) => {
                let coll = coll.into_inner();
                let coll = coll.load_into_memory()?;
                Ok(SearchContainer::LinearCollection(coll.try_into()?, mf))
            }
        }
    }
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap. SIGNATURES VERSION 2 @CTB.
/// @CTB can we refactor to use linear?

pub fn load_sketches_above_threshold_sigs_XXX(
    collection: &Collection,
    query: &KmerMinHash,
    threshold_hashes: u64,
) -> Result<(RevIndex, CounterGather, usize, usize)> {
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let collection: CollectionSet = collection.clone().try_into()?;
    let selection = collection.selection();

    let matchlist: Vec<Signature> = collection
        .par_iter()
        .filter_map(|(_idx, against_record)| {
            let mut results = Vec::new();
            // Load against into memory
            if let Ok(against_sig) = collection.sig_from_record(against_record) {
                let against_filename = against_sig.filename();
                let orig_sig = against_sig.clone();
                let against_mh: KmerMinHash = against_sig.try_into().expect("cannot get sketch");

                let against_mh_ds = against_mh
                    .downsample_scaled(query.scaled())
                    .expect("cannot downsample sketch");

                // good? ok, store as candidate from prefetch.
                if let Ok(overlap) = against_mh_ds.count_common(query, false) {
                    if overlap > 0 && overlap >= threshold_hashes {
                        results.push(orig_sig.into());
                    }
                } else {
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        against_filename
                    );
                    let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            } else {
                // this shouldn't happen here anymore -- likely would happen at load_collection
                eprintln!(
                    "WARNING: could not load sketches for record '{}'",
                    against_record.internal_location()
                );
                let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
            if results.is_empty() {
                None
            } else {
                Some(results)
            }
        })
        .flatten()
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    let revindex = MemRevIndex::new_with_sigs(matchlist,
                                              &selection,
                                              threshold_hashes as usize,
                                              None)?;

    let cg = revindex.prepare_gather_counters(query, None);
    

    Ok((revindex, cg, skipped_paths, failed_paths))
}

// @CTB enum_dispatch
#[derive(Clone)]
enum LoadedDatabase {
    InvertedIndex(RevIndex),
    LinearCollection(Collection),
}

impl LoadedDatabase {
    fn len(&self) -> usize {
        match self {
            LoadedDatabase::InvertedIndex(revindex) => revindex.len(),
            LoadedDatabase::LinearCollection(coll) => coll.len(),
        }
    }

    fn collection(&self) -> &Collection {
        match self {
            LoadedDatabase::InvertedIndex(revindex) => {
                revindex.collection()
            }
            LoadedDatabase::LinearCollection(coll) => {
                coll
            }
        }
    }
    fn manifest(&self) -> &Manifest {
        self.collection().manifest()
    }
}

#[derive(Clone)]
pub struct MultiCollection {
    dbs: Vec<LoadedDatabase>
}

// A collection of databases, including indexes, on-disk collections, and
// in-memory collections.
impl MultiCollection {
    fn new(dbs: Vec<LoadedDatabase>) -> Self {
        Self { dbs }
    }

    /// top level load function; tries to load anything and everything passed
    /// in.
    pub fn load(sigpath: &Path) -> Result<(Self, usize)> {
        let mut last_error = None;

        let collection = if sigpath.extension().map_or(false, |ext| ext == "zip") {
            match MultiCollection::from_zipfile(&sigpath) {
                Ok(coll) => Some((coll, 0)),
                Err(e) => {
                    last_error = Some(e);
                    None
                }
            }
        } else {
            None
        };

        let collection = collection.or_else(|| match MultiCollection::from_rocksdb(&sigpath) {
            Ok(coll) => Some((coll, 0)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

        // we support RocksDB directory paths, but nothing else, unlike sourmash.
        if collection.is_none() {
            let path_metadata = metadata(sigpath).expect("getting path metadata failed");
            if path_metadata.is_dir() {
                bail!("arbitrary directories are not supported as input");
            }
        }

        let collection =
            collection.or_else(
                || match MultiCollection::from_standalone_manifest(&sigpath) {
                    Ok(coll) => Some((coll, 0)),
                    Err(e) => {
                        last_error = Some(e);
                        None
                    }
                },
            );

        let collection = collection.or_else(|| match MultiCollection::from_signature(&sigpath) {
            Ok(coll) => Some((coll, 0)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

        let collection = collection.or_else(|| match MultiCollection::from_pathlist(&sigpath) {
            Ok((coll, n_failed)) => Some((coll, n_failed)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

        if let Some((collection, n_failed)) = collection {
            Ok((collection, n_failed))
        } else {
            let err = last_error.expect("no collection loaded but no error set!?");
            Err(err)
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

            // @CTB colls.intersect_manifest(&manifest);

            Ok(colls)
        }
    }

    /// Load a collection from a .zip file.
    pub fn from_zipfile(sigpath: &Path) -> Result<Self> {
        debug!("multi from zipfile!");
        match Collection::from_zipfile(sigpath) {
            Ok(collection) => {
                Ok(MultiCollection::new(vec![LoadedDatabase::LinearCollection(collection)]))
            },
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
            let db = match RevIndex::open(sigpath, true, None) {
                Ok(db) => db,
                Err(e) => {
                    return Err(anyhow::anyhow!(
                        "cannot open RocksDB database. Error is: {}",
                        e
                    ))
                }
            };
            Ok(MultiCollection::new(vec![LoadedDatabase::InvertedIndex(db)]))
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
        Ok(MultiCollection::new(vec![LoadedDatabase::LinearCollection(coll)]))
    }

    pub fn len(&self) -> usize {
        let val: usize = self.dbs.iter().map(|c| c.len()).sum();
        val
    }

/*
    pub fn is_empty(&self) -> bool {
        let val: usize = self.dbs.iter().map(|c| c.len()).sum();
        val == 0
    }
*/
    pub fn max_scaled(&self) -> Option<&ScaledType> {
        //self.item_iter().map(|(_, _, record)| record.scaled()).max()
        self.dbs.iter().filter_map(|db| {
            if let Some((_, max_scaled)) = db.collection().min_max_scaled() {
                Some(max_scaled)
            } else {
                None
            }
        }).max()
    }

    pub fn select<'a>(&'a self, selection: &Selection) -> Result<MultiCollectionSet<'a>, SourmashError> {
        let collections = self
            .dbs
            .iter()
            .map(|c| {
                match c {
                    LoadedDatabase::LinearCollection(coll) => {
                        let coll = coll.clone();
                        let coll = coll.select(selection).expect("failed select");
                        let cs: CollectionSet = coll.try_into().expect("incomplete selection!?");
                        let mf = c.manifest();
                        SearchContainer::LinearCollection(cs, mf)
                    },
                    LoadedDatabase::InvertedIndex(revindex) => {
                        let new_ri = revindex.clone();
                        let mf = revindex.collection().manifest();
                        SearchContainer::InvertedIndex(new_ri, mf)
                    }
                }
            })
            .collect();

        Ok(MultiCollectionSet { collections })
    }
}
    
/// A collection of sketches, potentially stored in multiple files.
#[derive(Clone)]
pub struct MultiCollectionSet<'a> {
    collections: Vec<SearchContainer<'a>>,
}

impl<'a> MultiCollectionSet<'a> {
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

    pub fn selection(&self) -> Selection {
        // @CTB
        // turn into a collection set, if possible, and then extract first
        // selection. Better if it makes an error I think?
        let coll = self.collections.iter().next().expect("empty?!");
        let cs: CollectionSet = coll.collection().clone().try_into().expect("err");
        cs.selection()
    }

    // iterate over tuples
    pub fn item_iter(&self) -> impl Iterator<Item = (&CollectionSet, Idx, &Record)> {
        let s: Vec<_> = self
            .collections
            .iter()
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c.collection(), _idx, record)))
            .collect();
        s.into_iter()
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (&CollectionSet, Idx, &Record)> {
        // first create a Vec of all triples (Collection, Idx, Record)
        let s: Vec<_> = self
            .collections
            .iter() // CTB: are we loading things into memory here? No...
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c.collection(), _idx, record)))
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
    // @CTB refactor / use Self
    pub fn load_sketches(self) -> Result<Vec<SmallSignature>> {
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
                    let sig_md5 = record.md5().clone();
                    let minhash: KmerMinHash = sig.try_into().expect("cannot extract sketch");

                    Some(SmallSignature {
                        location: record.internal_location().to_string(),
                        name: sig_name,
                        md5sum: sig_md5,
                        minhash,
                    })
                }
                Err(_) => {
                    eprintln!(
                        "FAILED to load sketch from '{}' (1)",
                        record.internal_location()
                    );
                    None
                }
            })
            .collect();

        Ok(sketchinfo)
    }

    // Load all sketches into memory, return new MultiCollection.
    // @CTB rename
    pub fn load_sketches2(self) -> Result<Self> {
        let n_failed = AtomicUsize::new(0);

        let new_coll: Vec<_> = self
            .collections
            .into_iter()
            .map(|s| s.load().expect("couldn't load collection"))
            .collect();

        Ok(Self { collections: new_coll })
    }

    fn intersect_manifest(&mut self, manifest: &Manifest) {
        for coll in self.collections.iter_mut() {
            coll.intersect_manifest(manifest);
        }
    }

    pub fn prefetch(&self, query: &KmerMinHash, threshold_hashes: u64) -> Result<(PrefetchContainer, usize, usize)> {
        let mut skipped_paths = 0;
        let mut failed_paths = 0;

        let mut res = PrefetchContainer { matchlists: vec![] };

        for searchable in self.collections.iter() {
            let (revindex, cg, mf, skip, fail) = searchable.prefetch(query, threshold_hashes)?;
            skipped_paths += skip;
            failed_paths += fail;

            res.matchlists.push(PrefetchItem { revindex, cg, mf });
        }
        Ok((res, skipped_paths, failed_paths))
    }
/*
    pub fn prefetch_iter(&self, query: &KmerMinHash, threshold_hashes: u64) -> impl Iterator<Item = (&'a SearchContainer, Idx)> + '_ {
        self.collections
            .iter()
            .flat_map(move |s| s
                      .prefetch_iter(query, threshold_hashes)
                      .map(move |idx| (s, idx))
                      )
    }
*/      
    pub fn prefetch_consume(self,
                            _query: &KmerMinHash,
                            _threshold_hashes: u64,
    ) -> Result<(Vec<(RevIndex, CounterGather)>, usize, usize)> {
/*        let pairs: Vec<(RevIndex, CounterGather)> = self
            .collections
            .iter()
            .map(|c| {
                ;
            }
        for collection in collections.iter() {
            
            }
                */
        panic!("foo");
    }

    // Load all sketches into memory, producing an in-memory Collection.
    // @CTB refactor?
    pub fn load_all_sigs(self) -> Result<Collection> {
        let all_sigs: Vec<Signature> = self
            .par_iter()
            .filter_map(|(coll, _idx, record)| match coll.sig_from_record(record) {
                Ok(sig) => Some(Signature::from(sig)),
                Err(_) => {
                    eprintln!(
                        "FAILED to load sketch from '{}' (3)",
                        record.internal_location()
                    );
                    None
                }
            })
            .collect();
        Ok(Collection::from_sigs(all_sigs)?)
    }
}

// Convert a single Collection into a MultiCollection
impl From<Collection> for MultiCollection {
    fn from(coll: Collection) -> Self {
        MultiCollection::new(vec![LoadedDatabase::LinearCollection(coll)])
    }
}

// Merge a bunch of MultiCollection structs into one
impl From<Vec<MultiCollection>> for MultiCollection {
    fn from(multi: Vec<MultiCollection>) -> Self {
        let mut x: Vec<LoadedDatabase> = vec![];
        for mc in multi.into_iter() {
            for coll in mc.dbs.into_iter() {
                x.push(coll);
            }
        }
        MultiCollection::new(x)
    }
}

impl TryFrom<MultiCollectionSet<'_>> for CollectionSet {
    type Error = &'static str;

    fn try_from(multi: MultiCollectionSet) -> Result<Self, Self::Error> {
        if multi.collections.len() == 1 {
            // this must succeed b/c len > 0
            let cont = multi.collections.into_iter().next().unwrap();

            let cs = cont.collection().clone();
            Ok(cs)
        } else {
            Err("More than one CollectionSet in this MultiCollectionSet; cannot convert")
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
