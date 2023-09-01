/// Utility functions for pyo3_branchwater.

use rayon::prelude::*;

use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use zip::read::ZipArchive;
use tempfile::tempdir;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use anyhow::{Result, anyhow};

use std::cmp::{PartialOrd, Ordering};

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;
use sourmash::prelude::MinHashOps;
use sourmash::prelude::FracMinHashOps;

/// Track a name/minhash.

pub struct SmallSignature {
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}

/// Structure to hold overlap information from comparisons.

pub struct PrefetchResult {
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
    pub overlap: u64,
}

impl Ord for PrefetchResult {
    fn cmp(&self, other: &PrefetchResult) -> Ordering {
        self.overlap.cmp(&other.overlap)
    }
}

impl PartialOrd for PrefetchResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for PrefetchResult {
    fn eq(&self, other: &Self) -> bool {
        self.overlap == other.overlap
    }
}

impl Eq for PrefetchResult {}

/// check to see if two KmerMinHash are compatible.
///
/// CTB note: despite the name, downsampling is not performed?
/// Although it checks if they are compatible in one direction...

fn check_compatible_downsample(
    me: &KmerMinHash,
    other: &KmerMinHash,
) -> Result<(), sourmash::Error> {
    /* // ignore num minhashes.
    if self.num != other.num {
        return Err(Error::MismatchNum {
            n1: self.num,
            n2: other.num,
        }
        .into());
    }
    */
    use sourmash::Error;

    if me.ksize() != other.ksize() {
        return Err(Error::MismatchKSizes);
    }
    if me.hash_function() != other.hash_function() {
        // TODO: fix this error
        return Err(Error::MismatchDNAProt);
    }
    if me.max_hash() < other.max_hash() {
        return Err(Error::MismatchScaled);
    }
    if me.seed() != other.seed() {
        return Err(Error::MismatchSeed);
    }
    Ok(())
}


/// Given a vec of search Signatures, each containing one or more sketches,
/// and a template Sketch, return a compatible (& now downsampled)
/// Sketch from the search Signatures..
///
/// CTB note: this will return the first acceptable match, I think, ignoring
/// all others.


pub fn prepare_query(search_sigs: &[Signature], template: &Sketch) -> Option<SmallSignature> {

    for search_sig in search_sigs.iter() {
        // find exact match for template?
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
            return Some(SmallSignature {
                name: search_sig.name(),
                md5sum: mh.md5sum(),
                minhash: mh.clone()
            });
        } else {
            // no - try to find one that can be downsampled
            if let Sketch::MinHash(template_mh) = template {
                for sketch in search_sig.sketches() {
                    if let Sketch::MinHash(ref_mh) = sketch {
                        if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                            let max_hash = max_hash_for_scaled(template_mh.scaled());
                            let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                            return Some(SmallSignature {
                                name: search_sig.name(),
                                md5sum: ref_mh.md5sum(), // original
                                minhash: mh,             // downsampled
                            });
                        }
                    }
                }
            }
        }
    }
    None
}

/// Write list of prefetch matches.

pub fn write_prefetch<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query_label: String,
    prefetch_output: Option<P>,
    matchlist: &BinaryHeap<PrefetchResult>
) -> Result<()> {
    // Set up a writer for prefetch output
    let prefetch_out: Box<dyn Write> = match prefetch_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(prefetch_out);
    writeln!(&mut writer, "query_filename,match_name,match_md5,intersect_bp").ok();

    for m in matchlist.iter() {
        writeln!(&mut writer, "{},\"{}\",{},{}", query_label,
                 m.name, m.md5sum, m.overlap).ok();
    }

    Ok(())
}

/// Load a list of filenames from a file. Exits on bad lines.

pub fn load_sketchlist_filenames<P: AsRef<Path>>(sketchlist_filename: &P) ->
    Result<Vec<PathBuf>>
{
    let sketchlist_file = BufReader::new(File::open(sketchlist_filename)?);

    let mut sketchlist_filenames : Vec<PathBuf> = Vec::new();
    for line in sketchlist_file.lines() {
        let line = match line {
            Ok(v) => v,
            Err(_) => return {
                let filename = sketchlist_filename.as_ref().display();
                let msg = format!("invalid line in fromfile '{}'", filename);
                Err(anyhow!(msg))
            },
        };

        if !line.is_empty() {
            let mut path = PathBuf::new();
            path.push(line);
            sketchlist_filenames.push(path);
        }
    }
    Ok(sketchlist_filenames)
}

/// Load a collection of sketches from a file in parallel.

pub fn load_sketches(sketchlist_paths: Vec<PathBuf>, template: &Sketch) ->
    Result<(Vec<SmallSignature>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let sketchlist : Vec<SmallSignature> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let mut sm = None;

            let filename = m.display();

            if let Ok(sigs) = Signature::from_path(m) {
                sm = prepare_query(&sigs, template);
                if sm.is_none() {
                    // track number of paths that have no matching sigs
                    let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            } else {
                // failed to load from this path - print error & track.
                eprintln!("WARNING: could not load sketches from path '{}'",
                          filename);
                let _i = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
            sm
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);
    Ok((sketchlist, skipped_paths, failed_paths))
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

pub fn load_sketches_above_threshold(
    sketchlist_paths: Vec<PathBuf>,
    template: &Sketch,
    query: &KmerMinHash,
    threshold_hashes: u64
) ->
    Result<(BinaryHeap<PrefetchResult>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let matchlist: BinaryHeap<PrefetchResult> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m);

            match sigs {
                Ok(sigs) => {
                    let mut mm = None;

                    if let Some(sm) = prepare_query(&sigs, template) {
                        let mh = sm.minhash;
                        if let Ok(overlap) = mh.count_common(query, false) {
                            if overlap >= threshold_hashes {
                                let result = PrefetchResult {
                                    name: sm.name,
                                    md5sum: sm.md5sum,
                                    minhash: mh,
                                    overlap,
                                };
                                mm = Some(result);
                            }
                        }
                    } else {
                        eprintln!("WARNING: no compatible sketches in path '{}'",
                                  m.display());
                        let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    mm
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("WARNING: could not load sketches from path '{}'",
                          m.display());
                    None
                }
            }
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((matchlist, skipped_paths, failed_paths))
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

pub fn consume_query_by_gather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    mut query: KmerMinHash,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<P>,
    query_label: String
) -> Result<()> {
    // Set up a writer for gather output
    let gather_out: Box<dyn Write> = match gather_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(gather_out);
    writeln!(&mut writer, "query_filename,rank,match_name,match_md5,intersect_bp").ok();

    let mut matching_sketches = matchlist;
    let mut rank = 0;

    while !matching_sketches.is_empty() {
        eprintln!("{} remaining: {} {}", query_label,
                  query.size(), matching_sketches.len());
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query.remove_from(&best_element.minhash)?;

        writeln!(&mut writer, "{},{},\"{}\",{},{}", query_label, rank,
                 best_element.name, best_element.md5sum,
                 best_element.overlap).ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query, matching_sketches, threshold_hashes);
        rank += 1;
    }
    Ok(())
}
                           

/// Find sketches in 'sketchlist' that overlap with 'query' above
/// specified threshold.

fn prefetch(
    query: &KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
) -> BinaryHeap<PrefetchResult> {
    sketchlist
        .into_par_iter()
        .filter_map(|result| {
            let mut mm = None;
            let searchsig = &result.minhash;
            let overlap = searchsig.count_common(query, false);
            if let Ok(overlap) = overlap {
                if overlap >= threshold_hashes {
                    let result = PrefetchResult {
                        overlap,
                        ..result
                    };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

// mastiff rocksdb functions

pub fn build_template(ksize: u8, scaled: usize) -> Sketch {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    Sketch::MinHash(template_mh)
}

pub fn read_signatures_from_zip<P: AsRef<Path>>(
    zip_path: P,
) -> Result<(Vec<PathBuf>, tempfile::TempDir), Box<dyn std::error::Error>> {
    let mut signature_paths = Vec::new();
    let temp_dir = tempdir()?;
    let zip_file = File::open(&zip_path)?;
    let mut zip_archive = ZipArchive::new(zip_file)?;

    for i in 0..zip_archive.len() {
        let mut file = zip_archive.by_index(i)?;
        let mut sig = Vec::new();
        file.read_to_end(&mut sig)?;

        let file_name = Path::new(file.name()).file_name().unwrap().to_str().unwrap();
        if file_name.ends_with(".sig") || file_name.ends_with(".sig.gz") {
            let new_path = temp_dir.path().join(file_name);
            let mut new_file = File::create(&new_path)?;
            new_file.write_all(&sig)?;
            signature_paths.push(new_path);
        }
    }
    println!("wrote {} signatures to temp dir", signature_paths.len());
    Ok((signature_paths, temp_dir))
}

pub fn is_revindex_database(path: &Path) -> bool {
    // quick file check for Revindex database:
    // is path a directory that contains a file named 'CURRENT'?
    if path.is_dir() {
        let current_file = path.join("CURRENT");
        current_file.exists() && current_file.is_file()
    } else {
        false
    }
}