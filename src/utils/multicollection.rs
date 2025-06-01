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

use sourmash::collection::{ Collection, CollectionSet };
use sourmash::encodings::Idx;
use sourmash::errors::SourmashError;
use sourmash::index::revindex::{ CounterGather, RevIndex, RevIndexOps };
use sourmash::index::revindex::mem_revindex::MemRevIndex;
use sourmash::manifest::{Manifest, Record};
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};
use sourmash::ScaledType;


pub struct PrefetchContainer {
    pub matchlists: Vec<(RevIndex, CounterGather, Manifest)>,
}

impl PrefetchContainer {
    // sum the lengths of the contents
    pub fn len(&self) -> usize {
        let mut l = 0;
        for (_, cg, _) in self.matchlists.iter() {
            l += cg.len();
        }
        l
    }

    // true if all contained matchlists are false
    pub fn is_empty(&self) -> bool {
        for (_, cg, _) in self.matchlists.iter() {
            if !cg.is_empty() { return false; };
        }
        return true;
    }

    // find best match across all CounterGather objects
    pub fn peek(&self, threshold_hashes: u64) -> Option<(Signature, &Record)> {
        let mut best_idx = None;
        let mut best_overlap = 0;
        let mut best_revindex = None;
        let mut best_mf = None;
        
        for (revindex, cg, mf) in self.matchlists.iter() {
            if !cg.is_empty() {
                let (idx, overlap) = cg.peek(threshold_hashes as usize).expect("empty?!");
                if overlap > best_overlap {
                    best_idx = Some(idx);
                    best_overlap = overlap;
                    best_revindex = Some(revindex);
                    best_mf = Some(mf);
                }
            }
        }

        // did we find something?
        if let Some(revindex) = best_revindex {
            let idx = best_idx.unwrap();
            let mf = best_mf.unwrap();
            let match_sig: Signature = revindex.collection().sig_for_dataset(idx).expect("cannot load").into();
            let orig_record: &Record = mf.get_record(idx).expect("cannot retrieve original record!?");
            Some((match_sig, orig_record))
        } else {
            None
        }
    }

    // consume the best match across all CounterGathers.
    pub fn consume(self, intersect_mh: &KmerMinHash) -> Self {
        let mut updated: Vec<(RevIndex, CounterGather, Manifest)> = vec![];
        
        for (revindex, mut cg, mf) in self.matchlists.into_iter() {
            cg.consume(intersect_mh);
            updated.push((revindex, cg, mf));
        }
        PrefetchContainer { matchlists : updated }
    }
}


trait Searchable {
    fn prefetch(&self, query: &KmerMinHash, threshold_hashes: u64) ->
        Result<(RevIndex, CounterGather, Manifest, usize, usize)>;
    fn len(&self) -> usize;
    fn select(&self, selection: &Selection) -> Result<Self, SourmashError> where Self: Sized;
    fn iter(&self) -> impl Iterator<Item = (Idx, &Record)>;
    fn collection(&self) -> &Collection;
    fn intersect_manifest(&mut self, manifest: &Manifest);
    fn get_orig_manifest(&self) -> &Manifest;
}


#[derive(Clone)]
enum SearchContainer {
    InvertedIndex(RevIndex, Manifest),
    LinearCollection(Collection, Manifest),
}

impl Searchable for SearchContainer {
    /// find all overlapping sketches; return RevIndex and CounterGather
    /// structs containing the results.
    ///
    /// For pre-existing RevIndexes, this will return a clone of the revindex;
    /// for other collections, it will return a MemRevIndex.

    fn prefetch(&self, query: &KmerMinHash, threshold_hashes: u64) ->
        Result<(RevIndex, CounterGather, Manifest, usize, usize)> {
            match self {
                SearchContainer::InvertedIndex(revindex, mf) => {
                    let cg = revindex.prepare_gather_counters(query, None);
                    Ok((revindex.clone(), cg, mf.clone(), 0, 0))
                },
                SearchContainer::LinearCollection(coll, mf) => {
                    let (revindex, cg, skip, fail) = 
                        load_sketches_above_threshold_sigs_XXX(coll,
                                                               query,
                                                               threshold_hashes)?;
                    Ok((revindex, cg, mf.clone(), skip, fail))
                                                           
                },
            }
        }

    fn len(&self) -> usize {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => revindex.len(),
            SearchContainer::LinearCollection(coll, mf) => coll.len(),
        }
    }

    // @CTB this will need to be updated for orig_manifest.
    fn select(&self, selection: &Selection) -> Result<Self, SourmashError> {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => Ok(SearchContainer::InvertedIndex(revindex.clone(), mf.clone())), // @CTB
            SearchContainer::LinearCollection(coll, mf) => {
                let c2 = coll.clone().select(selection)?;
                Ok(SearchContainer::LinearCollection(c2, mf.clone()))
            }
        }
    }

    /// iterate over (Idx, &Record)
    fn iter(&self) -> impl Iterator<Item = (Idx, &Record)> {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => {
                revindex.collection().iter()
            }
            SearchContainer::LinearCollection(coll, mf) => {
                coll.iter()
            }
        }
    }

    fn collection(&self) -> &Collection {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => {
                revindex.collection()
            }
            SearchContainer::LinearCollection(coll, mf) => {
                coll
            }
        }
    }

    /// retrieve the original manifest for this SearchContainer, prior
    /// to any downsampling/indexing.
    fn get_orig_manifest(&self) -> &Manifest {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => {
                &mf
            }
            SearchContainer::LinearCollection(coll, mf) => {
                &mf
            }
        }
    }

    // @CTB this will need to be updated in tricky ways.
    fn intersect_manifest(&mut self, manifest: &Manifest) {
        match self {
            SearchContainer::InvertedIndex(revindex, mf) => {
                panic!("foo 3");
            }
            SearchContainer::LinearCollection(coll, mf) => {
                coll.intersect_manifest(manifest);
            }
        }
    }
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap. SIGNATURES VERSION 2 @CTB.

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
                let against_md5 = against_record.md5().clone(); // keep original md5sum

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

/// A collection of sketches, potentially stored in multiple files.
#[derive(Clone)]
pub struct MultiCollection {
    collections: Vec<SearchContainer>,
}

impl MultiCollection {
    fn new(collections: Vec<SearchContainer>) -> Self {
        Self {
            collections,
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
            Ok(collection) => {
                let mf = collection.manifest().clone();
                Ok(MultiCollection::new(vec![SearchContainer::LinearCollection(collection, mf)]))
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
            let mf = db.collection().manifest().clone();
            Ok(MultiCollection::new(vec![SearchContainer::InvertedIndex(db, mf)]))
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
        let mf = coll.manifest().clone();
        Ok(MultiCollection::new(vec![SearchContainer::LinearCollection(coll, mf)]))
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

    pub fn selection(&self) -> Selection {
        // @CTB
        // turn into a collection set, if possible, and then extract first
        // selection. Better if it makes an error I think?
        let coll = self.collections.iter().next().expect("empty?!");
        let cs: CollectionSet = coll.collection().clone().try_into().expect("err");
        cs.selection()
    }

    // iterate over tuples
    pub fn item_iter(&self) -> impl Iterator<Item = (&Collection, Idx, &Record)> {
        let s: Vec<_> = self
            .collections
            .iter()
            .flat_map(|c| c.iter().map(move |(_idx, record)| (c.collection(), _idx, record)))
            .collect();
        s.into_iter()
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (&Collection, Idx, &Record)> {
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

    // Load all sketches into a MemRevIndex, return new MultiCollection.
    pub fn load_sketches_revindex(self) -> Result<Self> {
        let mut n_failed = AtomicUsize::new(0);
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
                    Some(sig.into())
                }
                Err(_) => {
                    eprintln!(
                        "FAILED to load sketch from '{}' (2)",
                        record.internal_location()
                    );
                    let _ = n_failed.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            })
            .collect();

        let n_failed = n_failed.load(atomic::Ordering::SeqCst);

        if n_failed > 0 {
            return Err(anyhow::anyhow!("{n_failed} sketches failed to load. See error messages above."));
        }

        let combined_mf: Vec<Record> = self
            .collections
            .iter()
            .flat_map(|c| c.get_orig_manifest()
                      .iter()
                      .map(|record| record.clone()))
            .collect();
        let combined_mf: Manifest = combined_mf.into();

        let selection = self.selection(); // @CTB should we need this?
        let revindex = MemRevIndex::new_with_sigs(sketchinfo,
                                                  &selection,
                                                  0,
                                                  None)?;

        Ok(MultiCollection::new(vec![ SearchContainer::InvertedIndex(revindex, combined_mf) ]))
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

            res.matchlists.push((revindex, cg, mf));
        }
        Ok((res, skipped_paths, failed_paths))
    }
        
    pub fn prefetch_consume(self,
                            query: &KmerMinHash,
                            threshold_hashes: u64,
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

impl Select for MultiCollection {
    fn select(self, selection: &Selection) -> Result<Self, SourmashError> {
        let collections = self
            .collections
            .into_iter()
            .filter_map(|c| c.clone().select(selection).ok())
            .collect();

        Ok(MultiCollection::new(collections))
    }
}

// Convert a single Collection into a MultiCollection
impl From<Collection> for MultiCollection {
    fn from(coll: Collection) -> Self {
        let mf = coll.manifest().clone();
        MultiCollection::new(vec![SearchContainer::LinearCollection(coll, mf)])
    }
}

// Merge a bunch of MultiCollection structs into one
impl From<Vec<MultiCollection>> for MultiCollection {
    fn from(multi: Vec<MultiCollection>) -> Self {
        let mut x: Vec<SearchContainer> = vec![];
        for mc in multi.into_iter() {
            for coll in mc.collections.into_iter() {
                x.push(coll);
            }
        }
        MultiCollection::new(x)
    }
}

// Extract a single Collection from a MultiCollection, if possible
impl TryFrom<MultiCollection> for Collection {
    type Error = &'static str;

    fn try_from(multi: MultiCollection) -> Result<Self, Self::Error> {
        if multi.collections.len() == 1 {
            // this must succeed b/c len > 0
            Ok(multi.collections.into_iter().next().unwrap().collection().clone())
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
