/// Utility functions for sourmash_plugin_branchwater.
use rayon::prelude::*;
use sourmash::encodings::HashFunctions;
use sourmash::manifest::Manifest;
use sourmash::selection::Select;

use anyhow::{anyhow, Context, Result};
use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use csv::Writer;
use serde::ser::Serializer;
use serde::Serialize;
use std::cmp::{Ordering, PartialOrd};
use std::collections::BinaryHeap;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::panic;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::collection::Collection;
use sourmash::manifest::Record;
use sourmash::selection::Selection;
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};

/// Track a name/minhash.

pub struct SmallSignature {
    pub location: String,
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

/// Find sketches in 'sketchlist' that overlap with 'query' above
/// specified threshold.

pub fn prefetch(
    query_mh: &KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
) -> BinaryHeap<PrefetchResult> {
    sketchlist
        .into_par_iter()
        .filter_map(|result| {
            let mut mm = None;
            let searchsig = &result.minhash;
            // downsample within count_common
            let overlap = searchsig.count_common(query_mh, true);
            if let Ok(overlap) = overlap {
                if overlap >= threshold_hashes {
                    let result = PrefetchResult { overlap, ..result };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

/// Write list of prefetch matches.
pub fn write_prefetch(
    query: &SigStore,
    prefetch_output: Option<String>,
    matchlist: &BinaryHeap<PrefetchResult>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Define the writer to stdout by default
    let mut writer: Box<dyn Write> = Box::new(std::io::stdout());

    if let Some(output_path) = &prefetch_output {
        // Account for potential missing dir in output path
        let directory_path = Path::new(output_path).parent();

        // If a directory path exists in the filename, create it if it doesn't already exist
        if let Some(dir) = directory_path {
            create_dir_all(dir)?;
        }

        let file = File::create(output_path)?;
        writer = Box::new(BufWriter::new(file));
    }

    writeln!(
        &mut writer,
        "query_filename,query_name,query_md5,match_name,match_md5,intersect_bp"
    )
    .ok();

    for m in matchlist.iter() {
        writeln!(
            &mut writer,
            "{},\"{}\",{},\"{}\",{},{}",
            query.filename(),
            query.name(),
            query.md5sum(),
            m.name,
            m.md5sum,
            m.overlap
        )
        .ok();
    }

    Ok(())
}

pub fn load_fasta_fromfile(sketchlist_filename: String) -> Result<Vec<(String, PathBuf, String)>> {
    let mut rdr = csv::Reader::from_path(sketchlist_filename)?;

    // Check for right header
    let headers = rdr.headers()?;
    if headers.len() != 3
        || headers.get(0).unwrap() != "name"
        || headers.get(1).unwrap() != "genome_filename"
        || headers.get(2).unwrap() != "protein_filename"
    {
        return Err(anyhow!(
            "Invalid header. Expected 'name,genome_filename,protein_filename', but got '{}'",
            headers.iter().collect::<Vec<_>>().join(",")
        ));
    }

    let mut results = Vec::new();

    let mut row_count = 0;
    let mut genome_count = 0;
    let mut protein_count = 0;
    // Create a HashSet to keep track of processed rows.
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.records() {
        let record = result?;

        // Skip duplicated rows
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());
        row_count += 1;
        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();

        let genome_filename = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'genome_filename' field"))?;
        if !genome_filename.is_empty() {
            results.push((
                name.clone(),
                PathBuf::from(genome_filename),
                "dna".to_string(),
            ));
            genome_count += 1;
        }

        let protein_filename = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'protein_filename' field"))?;
        if !protein_filename.is_empty() {
            results.push((name, PathBuf::from(protein_filename), "protein".to_string()));
            protein_count += 1;
        }
    }
    // Print warning if there were duplicated rows.
    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!(
        "Loaded {} rows in total ({} genome and {} protein files)",
        row_count, genome_count, protein_count
    );
    Ok(results)
}

// Load all compatible minhashes from a collection into memory
// also store sig name and md5 alongside, as we usually need those
pub fn load_sketches(
    collection: Collection,
    selection: &Selection,
    report_type: ReportType,
) -> Result<Vec<SmallSignature>> {
    let mut sketchinfo: Vec<SmallSignature> = Vec::new();
    for (_idx, record) in collection.iter() {
        if let Ok(sig) = collection.sig_from_record(record) {
            if let Some(minhash) = sig.clone().select(selection)?.minhash().cloned() {
                sketchinfo.push(SmallSignature {
                    location: record.internal_location().to_string(),
                    name: sig.name(),
                    md5sum: sig.md5sum(),
                    minhash,
                })
            }
        } else {
            bail!(
                "Error: Failed to load {} record: {}",
                report_type,
                record.name()
            );
        }
    }
    Ok(sketchinfo)
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

pub fn load_sketches_above_threshold(
    against_collection: Collection,
    query: &KmerMinHash,
    threshold_hashes: u64,
) -> Result<(BinaryHeap<PrefetchResult>, usize, usize)> {
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let matchlist: BinaryHeap<PrefetchResult> = against_collection
        .par_iter()
        .filter_map(|(_idx, against_record)| {
            let mut results = Vec::new();
            // Load against into memory
            if let Ok(against_sig) = against_collection.sig_from_record(against_record) {
                if let Some(against_mh) = against_sig.minhash() {
                    // if let Some(against_mh) = against_sig.select(&selection).unwrap().minhash() { // downsample via select
                    // currently downsampling here to avoid changing md5sum
                    if let Ok(overlap) = against_mh.count_common(query, true) {
                        //downsample via count_common
                        if overlap >= threshold_hashes {
                            let result = PrefetchResult {
                                name: against_record.name().to_string(),
                                md5sum: against_mh.md5sum(),
                                minhash: against_mh.clone(),
                                overlap,
                            };
                            results.push(result);
                        }
                    }
                } else {
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        against_sig.filename()
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

    Ok((matchlist, skipped_paths, failed_paths))
}

pub enum ReportType {
    Query,
    Against,
    General,
}

impl std::fmt::Display for ReportType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let description = match self {
            ReportType::Query => "query",
            ReportType::Against => "search",
            ReportType::General => "signature",
        };
        write!(f, "{}", description)
    }
}

pub fn collection_from_zipfile(sigpath: &Path, report_type: &ReportType) -> Result<Collection> {
    match Collection::from_zipfile(sigpath) {
        Ok(collection) => Ok(collection),
        Err(_) => bail!("failed to load {} zipfile: '{}'", report_type, sigpath),
    }
}

fn collection_from_manifest(
    sigpath: &Path,
    report_type: &ReportType,
) -> Result<Collection, anyhow::Error> {
    let file = File::open(sigpath)
        .with_context(|| format!("Failed to open {} file: '{}'", report_type, sigpath))?;

    let reader = BufReader::new(file);
    let manifest = Manifest::from_reader(reader).with_context(|| {
        format!(
            "Failed to read {} manifest from: '{}'",
            report_type, sigpath
        )
    })?;

    if manifest.is_empty() {
        // If the manifest is empty, return an error constructed with the anyhow! macro
        Err(anyhow!("could not read as manifest: '{}'", sigpath))
    } else {
        // If the manifest is not empty, proceed to create and return the Collection
        Ok(Collection::new(
            manifest,
            InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        ))
    }
}

fn collection_from_pathlist(
    sigpath: &Path,
    report_type: &ReportType,
) -> Result<(Collection, usize), anyhow::Error> {
    let file = File::open(sigpath).with_context(|| {
        format!(
            "Failed to open {} pathlist file: '{}'",
            report_type, sigpath
        )
    })?;
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
        .filter_map(|path| match Signature::from_path(&path) {
            Ok(signatures) => {
                let recs: Vec<Record> = signatures
                    .into_iter()
                    .flat_map(|v| Record::from_sig(&v, &path))
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
        eprintln!(
            "No valid signatures found in {} pathlist '{}'",
            report_type, sigpath
        );
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

    Ok((collection, n_failed))
}

fn collection_from_signature(sigpath: &Path, report_type: &ReportType) -> Result<Collection> {
    let signatures = Signature::from_path(sigpath).with_context(|| {
        format!(
            "Failed to load {} signatures from: '{}'",
            report_type, sigpath
        )
    })?;

    Collection::from_sigs(signatures).with_context(|| {
        format!(
            "Loaded {} signatures but failed to load as collection: '{}'",
            report_type, sigpath
        )
    })
}

pub fn load_collection(
    siglist: &String,
    selection: &Selection,
    report_type: ReportType,
    allow_failed: bool,
) -> Result<Collection> {
    let sigpath = PathBuf::from(siglist);

    if !sigpath.exists() {
        bail!("No such file or directory: '{}'", &sigpath);
    }

    // disallow rocksdb input here
    if is_revindex_database(&sigpath) {
        bail!("Cannot load {} signatures from a 'rocksdb' database. Please use sig, zip, or pathlist.", report_type);
    }

    eprintln!("Reading {}(s) from: '{}'", report_type, &siglist);
    let mut last_error = None;

    let collection = if sigpath.extension().map_or(false, |ext| ext == "zip") {
        match collection_from_zipfile(&sigpath, &report_type) {
            Ok(coll) => Some((coll, 0)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        }
    } else {
        None
    };

    let collection =
        collection.or_else(|| match collection_from_manifest(&sigpath, &report_type) {
            Ok(coll) => Some((coll, 0)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

    let collection =
        collection.or_else(|| match collection_from_signature(&sigpath, &report_type) {
            Ok(coll) => Some((coll, 0)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

    let collection =
        collection.or_else(|| match collection_from_pathlist(&sigpath, &report_type) {
            Ok((coll, n_failed)) => Some((coll, n_failed)),
            Err(e) => {
                last_error = Some(e);
                None
            }
        });

    match collection {
        Some((coll, n_failed)) => {
            let n_total = coll.len();
            let selected = coll.select(selection)?;
            let n_skipped = n_total - selected.len();
            report_on_collection_loading(
                &selected,
                n_skipped,
                n_failed,
                report_type,
                allow_failed,
            )?;
            Ok(selected)
        }
        None => {
            if let Some(e) = last_error {
                Err(e)
            } else {
                // Should never get here
                Err(anyhow!(
                    "Unable to load the collection for an unknown reason."
                ))
            }
        }
    }
}

/// Uses the output of collection loading function to report the
/// total number of sketches loaded, as well as the number of files,
/// if any, that failed to load or contained no compatible sketches.
/// If no sketches were loaded, bail.
///
/// # Arguments
///
/// * `sketchlist` - A slice of loaded `SmallSignature` sketches.
/// * `skipped_paths` - # paths that contained no compatible sketches.
/// * `failed_paths` - # paths that failed to load.
/// * `report_type` - ReportType Enum (Query or Against). Used to specify
///                   which sketch input this information pertains to.
///
/// # Returns
///
/// Returns `Ok(())` if at least one signature was successfully loaded.
/// Returns an error if no signatures were loaded.
///
/// # Errors
///
/// Returns an error if:
/// * No signatures were successfully loaded.
pub fn report_on_collection_loading(
    collection: &Collection,
    skipped_paths: usize,
    failed_paths: usize,
    report_type: ReportType,
    allow_failed: bool,
) -> Result<()> {
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} {} paths failed to load. See error messages above.",
            failed_paths, report_type
        );
        if !allow_failed {
            bail! {"Signatures failed to load. Exiting."}
        }
    }
    if skipped_paths > 0 {
        eprintln!(
            "WARNING: skipped {} {} paths - no compatible signatures.",
            skipped_paths, report_type
        );
    }

    // Validate sketches
    if collection.is_empty() {
        eprintln!("No {} signatures loaded, exiting.", report_type);
        return Ok(());
    }
    eprintln!("Loaded {} {} signature(s)", collection.len(), report_type);
    Ok(())
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

pub fn consume_query_by_gather(
    query: SigStore,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<String>,
) -> Result<()> {
    // Define the writer to stdout by default
    let mut writer: Box<dyn Write> = Box::new(std::io::stdout());

    if let Some(output_path) = &gather_output {
        // Account for potential missing dir in output path
        let directory_path = Path::new(output_path).parent();

        // If a directory path exists in the filename, create it if it doesn't already exist
        if let Some(dir) = directory_path {
            create_dir_all(dir)?;
        }

        let file = File::create(output_path)?;
        writer = Box::new(BufWriter::new(file));
    }
    writeln!(
        &mut writer,
        "query_filename,rank,query_name,query_md5,match_name,match_md5,intersect_bp"
    )
    .ok();

    let mut matching_sketches = matchlist;
    let mut rank = 0;

    let mut last_matches = matching_sketches.len();

    // let location = query.location;
    let location = query.filename(); // this is different (original fasta filename) than query.location was (sig name)!!

    let orig_query_mh = query.minhash().unwrap();
    let mut query_mh = orig_query_mh.clone();
    let mut last_hashes = orig_query_mh.size();

    eprintln!(
        "{} iter {}: start: query hashes={} matches={}",
        location,
        rank,
        orig_query_mh.size(),
        matching_sketches.len()
    );

    while !matching_sketches.is_empty() {
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query_mh.remove_from(&best_element.minhash)?;

        writeln!(
            &mut writer,
            "{},{},\"{}\",{},\"{}\",{},{}",
            location,
            rank,
            query.name(),
            query.md5sum(),
            best_element.name,
            best_element.md5sum,
            best_element.overlap
        )
        .ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query_mh, matching_sketches, threshold_hashes);
        rank += 1;

        let sub_hashes = last_hashes - query_mh.size();
        let sub_matches = last_matches - matching_sketches.len();

        eprintln!(
            "{} iter {}: remaining: query hashes={}(-{}) matches={}(-{})",
            location,
            rank,
            query_mh.size(),
            sub_hashes,
            matching_sketches.len(),
            sub_matches
        );

        last_hashes = query_mh.size();
        last_matches = matching_sketches.len();
    }
    Ok(())
}

pub fn build_selection(ksize: u8, scaled: usize, moltype: &str) -> Selection {
    let hash_function = match moltype {
        "dna" => HashFunctions::Murmur64Dna,
        "protein" => HashFunctions::Murmur64Protein,
        "dayhoff" => HashFunctions::Murmur64Dayhoff,
        "hp" => HashFunctions::Murmur64Hp,
        _ => panic!("Unknown molecule type: {}", moltype),
    };
    // let hash_function = HashFunctions::try_from(moltype)
    //     .map_err(|_| panic!("Unknown molecule type: {}", moltype))
    //     .unwrap();

    Selection::builder()
        .ksize(ksize.into())
        .scaled(scaled as u32)
        .moltype(hash_function)
        .build()
}

pub fn is_revindex_database(path: &camino::Utf8PathBuf) -> bool {
    // quick file check for Revindex database:
    // is path a directory that contains a file named 'CURRENT'?
    if path.is_dir() {
        let current_file = path.join("CURRENT");
        current_file.exists() && current_file.is_file()
    } else {
        false
    }
}

#[derive(Serialize)]
pub struct SearchResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub containment: f64,
    pub intersect_hashes: usize,
    pub match_md5: Option<String>,
    pub jaccard: Option<f64>,
    pub max_containment: Option<f64>,
}

#[derive(Serialize)]
pub struct BranchwaterGatherResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub match_md5: String,
    pub f_match_query: f64,
    pub intersect_bp: usize,
}

#[derive(Serialize)]
pub struct MultiSearchResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub match_md5: String,
    pub containment: f64,
    pub max_containment: f64,
    pub jaccard: f64,
    pub intersect_hashes: f64,
}

#[derive(Serialize)]
pub struct ManifestRow {
    pub md5: String,
    pub md5short: String,
    pub ksize: u32,
    pub moltype: String,
    pub num: u32,
    pub scaled: u64,
    pub n_hashes: usize,
    pub with_abundance: BoolPython,
    pub name: String,
    pub filename: String,
    pub internal_location: String,
}

// A wrapper type for booleans to customize serialization
pub struct BoolPython(bool);

impl Serialize for BoolPython {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self.0 {
            true => serializer.serialize_str("True"),
            false => serializer.serialize_str("False"),
        }
    }
}

pub fn make_manifest_row(
    sig: &Signature,
    filename: &Path,
    internal_location: &str,
    scaled: u64,
    num: u32,
    abund: bool,
    is_dna: bool,
    is_protein: bool,
) -> ManifestRow {
    if is_dna && is_protein {
        panic!("Both is_dna and is_protein cannot be true at the same time.");
    } else if !is_dna && !is_protein {
        panic!("Either is_dna or is_protein must be true.");
    }
    let moltype = if is_dna {
        "DNA".to_string()
    } else {
        "protein".to_string()
    };
    let sketch = &sig.sketches()[0];
    let ksize: u32 = if is_dna {
        sketch.ksize() as u32
    } else {
        sketch.ksize() as u32 / 3
    };
    ManifestRow {
        internal_location: internal_location.to_string(),
        md5: sig.md5sum(),
        md5short: sig.md5sum()[0..8].to_string(),
        ksize: ksize,
        moltype,
        num,
        scaled,
        n_hashes: sketch.size(),
        with_abundance: BoolPython(abund),
        name: sig.name().to_string(),
        filename: filename.to_string(),
    }
}

pub fn open_stdout_or_file(output: Option<String>) -> Box<dyn Write + Send + 'static> {
    // if output is a file, use open_output_file
    if let Some(path) = output {
        let outpath: PathBuf = path.into();
        Box::new(open_output_file(&outpath))
    } else {
        Box::new(std::io::stdout())
    }
}

pub fn open_output_file(output: &PathBuf) -> BufWriter<File> {
    let file = File::create(output).unwrap_or_else(|e| {
        eprintln!("Error creating output file: {:?}", e);
        std::process::exit(1);
    });
    BufWriter::new(file)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    pub ksize: u32,
    pub track_abundance: bool,
    pub num: u32,
    pub scaled: u64,
    pub seed: u32,
    pub is_protein: bool,
    pub is_dna: bool,
}
use std::hash::Hash;
use std::hash::Hasher;

impl Hash for Params {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.track_abundance.hash(state);
        self.num.hash(state);
        self.scaled.hash(state);
        self.seed.hash(state);
        self.is_protein.hash(state);
        self.is_dna.hash(state);
    }
}

pub enum ZipMessage {
    SignatureData(Vec<Signature>, Vec<Params>, PathBuf),
    WriteManifest,
}

pub fn sigwriter(
    recv: std::sync::mpsc::Receiver<ZipMessage>,
    output: String,
) -> std::thread::JoinHandle<Result<()>> {
    std::thread::spawn(move || -> Result<()> {
        // cast output as pathbuf
        let outpath: PathBuf = output.into();

        let file_writer = open_output_file(&outpath);

        let options = zip::write::FileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .large_file(true);
        let mut zip = zip::ZipWriter::new(file_writer);
        let mut manifest_rows: Vec<ManifestRow> = Vec::new();
        // keep track of md5sum occurrences to prevent overwriting duplicates
        let mut md5sum_occurrences: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();

        while let Ok(message) = recv.recv() {
            match message {
                ZipMessage::SignatureData(sigs, params, filename) => {
                    if sigs.len() != params.len() {
                        bail!("Mismatched lengths of signatures and parameters");
                    }
                    for (sig, param) in sigs.iter().zip(params.iter()) {
                        let md5sum_str = sig.md5sum();
                        let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
                        *count += 1;
                        let sig_filename = if *count > 1 {
                            format!("signatures/{}_{}.sig.gz", md5sum_str, count)
                        } else {
                            format!("signatures/{}.sig.gz", md5sum_str)
                        };
                        write_signature(sig, &mut zip, options, &sig_filename);
                        manifest_rows.push(make_manifest_row(
                            sig,
                            &filename,
                            &sig_filename,
                            param.scaled,
                            param.num,
                            param.track_abundance,
                            param.is_dna,
                            param.is_protein,
                        ));
                    }
                }
                ZipMessage::WriteManifest => {
                    println!("Writing manifest");
                    // Start the CSV file inside the zip
                    zip.start_file("SOURMASH-MANIFEST.csv", options).unwrap();
                    // write manifest version line
                    writeln!(&mut zip, "# SOURMASH-MANIFEST-VERSION: 1.0").unwrap();
                    // scoped block for csv writing
                    {
                        let mut csv_writer = Writer::from_writer(&mut zip);

                        for row in &manifest_rows {
                            if let Err(e) = csv_writer.serialize(row) {
                                eprintln!("Error writing item: {:?}", e);
                            }
                        }
                        //  CSV writer must be manually flushed to ensure all data is written
                        if let Err(e) = csv_writer.flush() {
                            eprintln!("Error flushing CSV writer: {:?}", e);
                        }
                    } // drop csv writer here

                    // Properly finish writing to the ZIP file
                    if let Err(e) = zip.finish() {
                        eprintln!("Error finalizing ZIP file: {:?}", e);
                    }
                }
            }
        }
        Ok(())
    })
}

pub fn csvwriter_thread<T: Serialize + Send + 'static>(
    recv: std::sync::mpsc::Receiver<T>,
    output: Option<String>,
) -> std::thread::JoinHandle<()> {
    // create output file
    let out = open_stdout_or_file(output);
    // spawn a thread that is dedicated to printing to a buffered output
    std::thread::spawn(move || {
        let mut writer = Writer::from_writer(out);

        for res in recv.iter() {
            if let Err(e) = writer.serialize(res) {
                eprintln!("Error writing item: {:?}", e);
            }
        }
        writer.flush().expect("Failed to flush writer.");
    })
}

pub fn write_signature(
    sig: &Signature,
    zip: &mut zip::ZipWriter<BufWriter<File>>,
    zip_options: zip::write::FileOptions,
    sig_filename: &str,
) {
    let wrapped_sig = vec![sig];
    let json_bytes = serde_json::to_vec(&wrapped_sig).unwrap();

    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )
            .unwrap();
            gz_writer.write_all(&json_bytes).unwrap();
        }
        buffer.into_inner()
    };

    zip.start_file(sig_filename, zip_options).unwrap();
    zip.write_all(&gzipped_buffer).unwrap();
}
