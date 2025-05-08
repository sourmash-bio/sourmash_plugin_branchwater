//! Utility functions for `sourmash_plugin_branchwater`.
use rayon::prelude::*;

use sourmash::encodings::HashFunctions;
use sourmash::selection::Select;
use sourmash::ScaledType;

use anyhow::{anyhow, Result};
use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use csv::Writer;
use glob::glob;
use serde::{Deserialize, Serialize};
// use rust_decimal::{MathematicalOps, Decimal};
use std::cmp::{Ordering, PartialOrd};
use std::collections::BinaryHeap;
use std::fs::{create_dir_all, metadata, File};
use std::io::{BufWriter, Write};
use std::panic;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::sync::mpsc::{Receiver, SyncSender};
use std::thread::JoinHandle;
use zip::write::{FileOptions, ZipWriter};
use zip::CompressionMethod;

use sourmash::ani_utils::{ani_ci_from_containment, ani_from_containment};
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use stats::{median, stddev};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

pub mod multicollection;
pub use multicollection::{MultiCollection, SmallSignature};

pub mod buildutils;
use buildutils::{BuildCollection, BuildManifest};

/// Structure to hold overlap information from comparisons.
pub struct PrefetchResult {
    pub name: String,
    pub md5sum: String,
    pub location: String,
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
            let overlap = searchsig.count_common(query_mh, false);
            if let Ok(overlap) = overlap {
                if overlap > 0 && overlap >= threshold_hashes {
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
    query_filename: String,
    query_name: String,
    query_md5: String,
    prefetch_output: Option<String>,
    matchlist: &BinaryHeap<PrefetchResult>,
) -> Result<()> {
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
            query_filename, query_name, query_md5, m.name, m.md5sum, m.overlap
        )
        .ok();
    }

    Ok(())
}

pub struct FastaData {
    pub name: String,
    pub paths: Vec<PathBuf>,
    pub input_type: String, // to do - this could use moltype enum?
}

enum CSVType {
    Assembly,
    Reads,
    Prefix,
    Unknown,
}

fn detect_csv_type(headers: &csv::StringRecord) -> CSVType {
    if headers.len() == 3
        && headers.get(0).unwrap() == "name"
        && headers.get(1).unwrap() == "genome_filename"
        && headers.get(2).unwrap() == "protein_filename"
    {
        CSVType::Assembly
    } else if headers.len() == 3
        && headers.get(0).unwrap() == "name"
        && headers.get(1).unwrap() == "read1"
        && headers.get(2).unwrap() == "read2"
    {
        CSVType::Reads
    } else if headers.len() == 4
        && headers.get(0).unwrap() == "name"
        && headers.get(1).unwrap() == "input_moltype"
        && headers.get(2).unwrap() == "prefix"
        && headers.get(3).unwrap() == "exclude"
    {
        CSVType::Prefix
    } else {
        CSVType::Unknown
    }
}

pub fn load_fasta_fromfile(
    sketchlist_filename: String,
    force: bool,
) -> Result<(Vec<FastaData>, usize)> {
    let mut rdr = csv::Reader::from_path(sketchlist_filename)?;

    // Check for right header
    let headers = rdr.headers()?;

    match detect_csv_type(headers) {
        CSVType::Assembly => process_assembly_csv(rdr),
        CSVType::Reads => process_reads_csv(rdr),
        CSVType::Prefix => process_prefix_csv(rdr, force),
        CSVType::Unknown => Err(anyhow!(
            "Invalid header. Expected 'name,genome_filename,protein_filename', 'name,read1,read2', or 'name,input_moltype,prefix,exclude', but got '{}'",
            headers.iter().collect::<Vec<_>>().join(",")
        )),
    }
}

fn process_assembly_csv(mut rdr: csv::Reader<std::fs::File>) -> Result<(Vec<FastaData>, usize)> {
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

        // Handle optional genome_filename
        if let Some(genome_filename) = record.get(1) {
            if !genome_filename.is_empty() {
                results.push(FastaData {
                    name: name.clone(),
                    paths: vec![PathBuf::from(genome_filename)],
                    input_type: "dna".to_string(),
                });
                genome_count += 1;
            }
        }
        // Handle optional protein_filename
        if let Some(protein_filename) = record.get(2) {
            if !protein_filename.is_empty() {
                results.push(FastaData {
                    name: name.clone(),
                    paths: vec![PathBuf::from(protein_filename)],
                    input_type: "protein".to_string(),
                });
                protein_count += 1;
            }
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
    let n_fastas = genome_count + protein_count;
    Ok((results, n_fastas))
}

fn process_reads_csv(mut rdr: csv::Reader<std::fs::File>) -> Result<(Vec<FastaData>, usize)> {
    let mut results = Vec::new();
    let mut processed_rows = std::collections::HashSet::new();
    let mut read1_count = 0;
    let mut read2_count = 0;
    let mut duplicate_count = 0;

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());

        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();
        let read1 = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'read1' field"))?;
        read1_count += 1;
        let mut paths = vec![PathBuf::from(read1)];
        // allow missing read2
        let read2 = record
            .get(2)
            .and_then(|r2| if r2.is_empty() { None } else { Some(r2) });
        if let Some(r2) = read2 {
            paths.push(PathBuf::from(r2));
            read2_count += 1;
        }
        results.push(FastaData {
            name: name.clone(),
            paths,
            input_type: "dna".to_string(),
        });
    }

    println!("Found 'reads' CSV, assuming all files are DNA.");
    println!(
        "Loaded {} rows in total ({} with read1 and {} with read2), {} duplicates skipped.",
        processed_rows.len(),
        read1_count,
        read2_count,
        duplicate_count
    );

    let n_fastas = read1_count + read2_count;

    Ok((results, n_fastas))
}

fn process_prefix_csv(
    mut rdr: csv::Reader<std::fs::File>,
    force: bool,
) -> Result<(Vec<FastaData>, usize)> {
    let mut results = Vec::new();
    let mut dna_count = 0;
    let mut protein_count = 0;
    let mut processed_rows = HashSet::new();
    let mut duplicate_count = 0;
    let mut all_paths = HashSet::new(); // track FASTA in use
    let mut duplicate_paths_count = HashMap::new();

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());

        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();

        let moltype = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'input_moltype' field"))?
            .to_string();

        // Validate moltype
        match moltype.as_str() {
            "protein" | "dna" | "DNA" => (),
            _ => return Err(anyhow!("Invalid 'input_moltype' field value: {}", moltype)),
        }

        // For both prefix and exclude, automatically append wildcard for expected "prefix" matching
        let prefix = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'prefix' field"))?
            .to_string()
            + "*";

        // optional exclude pattern
        let exclude = record.get(3).map(|s| s.to_string() + "*");

        // Use glob to find and collect all paths that match the prefix
        let included_paths = glob(&prefix)
            .expect("Failed to read glob pattern for included paths")
            .filter_map(Result::ok)
            .map(|path| PathBuf::from(path.to_str().expect("Path is not valid UTF-8")))
            .collect::<HashSet<PathBuf>>();

        // Use glob to find and collect all paths that match the exclude_prefix, if any
        let excluded_paths = if let Some(ref exclude_pattern) = exclude {
            glob(exclude_pattern)
                .expect("Failed to read glob pattern for excluded paths")
                .filter_map(Result::ok)
                .map(|path| PathBuf::from(path.to_str().expect("Path is not valid UTF-8")))
                .collect::<HashSet<PathBuf>>()
        } else {
            HashSet::new()
        };

        // Exclude the excluded_paths from included_paths
        let filtered_paths: Vec<PathBuf> = included_paths
            .difference(&excluded_paths)
            .cloned()
            .collect();

        // Track duplicates among filtered paths
        for path in &filtered_paths {
            if !all_paths.insert(path.clone()) {
                *duplicate_paths_count.entry(path.clone()).or_insert(0) += 1;
            }
        }

        if !filtered_paths.is_empty() {
            match moltype.as_str() {
                "dna" | "DNA" => dna_count += filtered_paths.len(),
                "protein" => protein_count += filtered_paths.len(),
                _ => {} // should not get here b/c validated earlier
            }
            results.push(FastaData {
                name: name.clone(),
                paths: filtered_paths.to_vec(),
                input_type: moltype.clone(),
            });
        }
    }

    let total_duplicate_paths: usize = duplicate_paths_count.values().sum();

    println!("Found 'prefix' CSV. Using 'glob' to find files based on 'prefix' column.");
    if total_duplicate_paths > 0 {
        eprintln!("Found identical FASTA paths in more than one row!");
        eprintln!("Duplicated paths:");
        for path in duplicate_paths_count.keys() {
            eprintln!("{:?}", path);
        }
        if !force {
            return Err(anyhow!(
                "Duplicated FASTA files found. Please use --force to bypass this check."
            ));
        } else {
            eprintln!("--force is set. Continuing...")
        }
    }
    println!(
        "Loaded {} rows in total ({} DNA FASTA and {} protein FASTA), {} duplicate rows skipped.",
        processed_rows.len(),
        dna_count,
        protein_count,
        duplicate_count,
    );

    let n_fastas = dna_count + protein_count;

    Ok((results, n_fastas))
}

/////////

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

pub fn load_sketches_above_threshold(
    against_collection: MultiCollection,
    query: &KmerMinHash,
    threshold_hashes: u64,
) -> Result<(BinaryHeap<PrefetchResult>, usize, usize)> {
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    if against_collection.contains_revindex {
        eprintln!("WARNING: loading all sketches from a RocksDB into memory!");
    }
    let matchlist: BinaryHeap<PrefetchResult> = against_collection
        .par_iter()
        .filter_map(|(coll, _idx, against_record)| {
            let mut results = Vec::new();
            // Load against into memory
            if let Ok(against_sig) = coll.sig_from_record(against_record) {
                let against_filename = against_sig.filename();
                let against_mh: KmerMinHash = against_sig.try_into().expect("cannot get sketch");
                let against_md5 = against_record.md5().clone(); // keep original md5sum

                let against_mh_ds = against_mh
                    .downsample_scaled(query.scaled())
                    .expect("cannot downsample sketch");

                // good? ok, store as candidate from prefetch.
                if let Ok(overlap) = against_mh_ds.count_common(query, false) {
                    if overlap > 0 && overlap >= threshold_hashes {
                        let result = PrefetchResult {
                            name: against_record.name().to_string(),
                            md5sum: against_md5,
                            minhash: against_mh_ds,
                            location: against_record.internal_location().to_string(),
                            overlap,
                        };
                        results.push(result);
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
            ReportType::General => "analysis",
        };
        write!(f, "{}", description)
    }
}

/// Load a multi collection from a path - this is the new top-level load function.

pub fn load_collection(
    siglist: &String,
    selection: &Selection,
    report_type: ReportType,
    allow_failed: bool,
) -> Result<MultiCollection> {
    let sigpath = PathBuf::from(siglist);

    if !sigpath.exists() {
        bail!("No such file or directory: '{}'", &sigpath);
    }

    eprintln!("Reading {}(s) from: '{}'", report_type, &siglist);
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
        let path_metadata = metadata(sigpath.clone()).expect("getting path metadata failed");
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
/// * `collection` - A MultiCollection.
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
    collection: &MultiCollection,
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
        bail!("No {} signatures loaded, exiting.", report_type);
    }
    eprintln!("Loaded {} {} signature(s)", collection.len(), report_type);
    Ok(())
}

// branchwater version that allows using PrefetchResult
#[allow(clippy::too_many_arguments)]
pub fn branchwater_calculate_gather_stats(
    orig_query: &KmerMinHash,
    query: &KmerMinHash,
    // these are separate in PrefetchResult, so just pass them separately in here
    match_mh: &KmerMinHash,
    match_name: String,
    match_md5: String,
    match_size: u64,
    match_filename: String,
    gather_result_rank: u64,
    sum_weighted_found: u64,
    total_weighted_hashes: u64,
    calc_abund_stats: bool,
    calc_ani_ci: bool,
    confidence: Option<f64>,
) -> Result<InterimGatherResult> {
    // bp remaining in subtracted query
    let remaining_bp = (query.size() as u64 - match_size) * query.scaled() as u64;

    // stats for this match vs original query
    let (intersect_orig, _) = match_mh.intersection_size(orig_query).unwrap();
    let intersect_bp = match_mh.scaled() as u64 * intersect_orig;
    let f_orig_query = intersect_orig as f64 / orig_query.size() as f64;
    let f_match_orig = intersect_orig as f64 / match_mh.size() as f64;

    // stats for this match vs current (subtracted) query
    let f_match = match_size as f64 / match_mh.size() as f64;
    let unique_intersect_bp = match_mh.scaled() as u64 * match_size;
    let f_unique_to_query = match_size as f64 / orig_query.size() as f64;

    // // get ANI values
    let ksize = match_mh.ksize() as f64;
    let query_containment_ani = ani_from_containment(f_orig_query, ksize);
    let match_containment_ani = ani_from_containment(f_match_orig, ksize);
    let mut query_containment_ani_ci_low = None;
    let mut query_containment_ani_ci_high = None;
    let mut match_containment_ani_ci_low = None;
    let mut match_containment_ani_ci_high = None;

    if calc_ani_ci {
        let n_unique_kmers = match_mh.n_unique_kmers();
        let (qani_low, qani_high) = ani_ci_from_containment(
            f_unique_to_query,
            ksize,
            match_mh.scaled(),
            n_unique_kmers,
            confidence,
        )?;
        query_containment_ani_ci_low = Some(qani_low);
        query_containment_ani_ci_high = Some(qani_high);

        let (mani_low, mani_high) = ani_ci_from_containment(
            f_match,
            ksize,
            match_mh.scaled(),
            n_unique_kmers,
            confidence,
        )?;
        match_containment_ani_ci_low = Some(mani_low);
        match_containment_ani_ci_high = Some(mani_high);
    }

    let average_containment_ani = (query_containment_ani + match_containment_ani) / 2.0;
    let max_containment_ani = f64::max(query_containment_ani, match_containment_ani);

    // set up non-abundance weighted values
    let mut f_unique_weighted = f_unique_to_query;
    let mut average_abund = 1.0;
    let mut median_abund = 1.0;
    let mut std_abund = 0.0;
    // should these default to the unweighted numbers?
    let mut n_unique_weighted_found = 0;
    let mut sum_total_weighted_found = 0;

    // If abundance, calculate abund-related metrics (vs current query)
    if calc_abund_stats {
        // take abunds from subtracted query
        let (abunds, unique_weighted_found) = match match_mh.inflated_abundances(query) {
            Ok((abunds, unique_weighted_found)) => (abunds, unique_weighted_found),
            Err(e) => return Err(e.into()),
        };

        n_unique_weighted_found = unique_weighted_found;
        sum_total_weighted_found = sum_weighted_found + n_unique_weighted_found;
        f_unique_weighted = n_unique_weighted_found as f64 / total_weighted_hashes as f64;

        average_abund = n_unique_weighted_found as f64 / abunds.len() as f64;

        // todo: try to avoid clone for these?
        median_abund = median(abunds.iter().cloned()).expect("cannot calculate median");
        std_abund = stddev(abunds.iter().cloned());
    }

    let result = InterimGatherResult {
        intersect_bp,
        f_orig_query,
        f_match,
        f_unique_to_query,
        f_unique_weighted,
        average_abund,
        median_abund,
        std_abund,
        match_filename,
        match_name,
        match_md5,
        f_match_orig,
        unique_intersect_bp,
        gather_result_rank,
        remaining_bp,
        n_unique_weighted_found,
        query_containment_ani,
        query_containment_ani_ci_low,
        query_containment_ani_ci_high,
        match_containment_ani_ci_low,
        match_containment_ani_ci_high,
        match_containment_ani,
        average_containment_ani,
        max_containment_ani,
        sum_weighted_found: sum_total_weighted_found,
        total_weighted_hashes,
    };
    Ok(result)
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

pub fn consume_query_by_gather(
    query_name: String,
    query_filename: String,
    orig_query_mh: KmerMinHash,
    scaled: u32,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<SyncSender<BranchwaterGatherResult>>,
) -> Result<()> {
    let mut matching_sketches = matchlist;
    let mut rank = 0;

    let mut last_matches = matching_sketches.len();

    let query_bp = orig_query_mh.n_unique_kmers();
    let query_n_hashes = orig_query_mh.size() as u64;
    let mut query_moltype = orig_query_mh.hash_function().to_string();
    if query_moltype.to_lowercase() == "dna" {
        query_moltype = query_moltype.to_uppercase();
    }
    let query_md5sum: String = orig_query_mh.md5sum().clone();
    let query_scaled = orig_query_mh.scaled();

    let total_weighted_hashes = orig_query_mh.sum_abunds();
    let ksize = orig_query_mh.ksize() as u16;
    let calc_abund_stats = orig_query_mh.track_abundance();
    let orig_query_size = orig_query_mh.size();
    let mut last_hashes = orig_query_size;

    // this clone is necessary because we iteratively change things!
    // to do == use this to subtract hashes instead
    // let mut query_mh = KmerMinHashBTree::from(orig_query_mh.clone());
    let mut query_mh = orig_query_mh.clone();

    let mut orig_query_ds = orig_query_mh.downsample_scaled(scaled)?;

    // track for full gather results
    let mut sum_weighted_found = 0;

    // set some bools
    let calc_ani_ci = false;
    let ani_confidence_interval_fraction = None;

    eprintln!(
        "{} iter {}: start: query hashes={} matches={}",
        query_filename,
        rank,
        orig_query_size,
        matching_sketches.len()
    );

    while !matching_sketches.is_empty() {
        let best_element = matching_sketches.peek().unwrap();

        query_mh = query_mh.downsample_scaled(best_element.minhash.scaled())?;

        // CTB: won't need this if we do not allow multiple scaleds;
        // see sourmash-bio/sourmash#2951
        orig_query_ds = orig_query_ds
            .downsample_scaled(best_element.minhash.scaled())
            .expect("cannot downsample");

        //calculate full gather stats
        let match_ = branchwater_calculate_gather_stats(
            &orig_query_ds,
            &query_mh,
            &best_element.minhash,
            best_element.name.clone(),
            best_element.md5sum.clone(),
            best_element.overlap,
            best_element.location.clone(),
            rank,
            sum_weighted_found,
            total_weighted_hashes,
            calc_abund_stats,
            calc_ani_ci,
            ani_confidence_interval_fraction,
        )?;

        // build full gather result, then write
        let gather_result = BranchwaterGatherResult {
            intersect_bp: match_.intersect_bp,
            f_orig_query: match_.f_orig_query,
            f_match: match_.f_match,
            f_unique_to_query: match_.f_unique_to_query,
            f_unique_weighted: match_.f_unique_weighted,
            average_abund: match_.average_abund,
            median_abund: match_.median_abund,
            std_abund: match_.std_abund,
            match_filename: match_.match_filename.clone(), // to do: get match filename
            match_name: match_.match_name.clone(),
            match_md5: match_.match_md5.clone(),
            f_match_orig: match_.f_match_orig,
            unique_intersect_bp: match_.unique_intersect_bp,
            gather_result_rank: match_.gather_result_rank as u32,
            remaining_bp: match_.remaining_bp,
            query_filename: query_filename.clone(),
            query_name: query_name.clone(),
            query_md5: query_md5sum.clone(),
            query_bp,
            ksize,
            moltype: query_moltype.clone(),
            scaled: query_scaled,
            query_n_hashes,
            query_abundance: query_mh.track_abundance(),
            query_containment_ani: match_.query_containment_ani,
            match_containment_ani: match_.match_containment_ani,
            average_containment_ani: match_.average_containment_ani,
            max_containment_ani: match_.max_containment_ani,
            n_unique_weighted_found: match_.n_unique_weighted_found,
            sum_weighted_found: match_.sum_weighted_found,
            total_weighted_hashes: match_.total_weighted_hashes,

            query_containment_ani_ci_low: match_.query_containment_ani_ci_low,
            query_containment_ani_ci_high: match_.query_containment_ani_ci_high,
            match_containment_ani_ci_low: match_.match_containment_ani_ci_low,
            match_containment_ani_ci_high: match_.match_containment_ani_ci_high,
        };
        sum_weighted_found = gather_result.sum_weighted_found;
        // send result to channel => CSV file.
        if let Some(ref s) = gather_output {
            s.send(gather_result)?;
        }

        // remove!
        query_mh.remove_from(&best_element.minhash)?;
        // to do -- switch to KmerMinHashTree, for faster removal.
        //query.remove_many(best_element.iter_mins().copied())?; // from sourmash core

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query_mh, matching_sketches, threshold_hashes);
        rank += 1;

        let sub_hashes = last_hashes - query_mh.size();
        let sub_matches = last_matches - matching_sketches.len();

        eprintln!(
            "{} iter {}: remaining: query hashes={}(-{}) matches={}(-{})",
            query_filename,
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

pub fn build_selection(ksize: u8, scaled: Option<u32>, moltype: &str) -> Selection {
    let hash_function = match moltype {
        "DNA" => HashFunctions::Murmur64Dna,
        "protein" => HashFunctions::Murmur64Protein,
        "dayhoff" => HashFunctions::Murmur64Dayhoff,
        "hp" => HashFunctions::Murmur64Hp,
        "skipm1n3" => HashFunctions::Murmur64Skipm1n3,
        "skipm2n3" => HashFunctions::Murmur64Skipm2n3,
        _ => panic!("Unknown molecule type: {}", moltype),
    };
    // let hash_function = HashFunctions::try_from(moltype)
    //     .map_err(|_| panic!("Unknown molecule type: {}", moltype))
    //     .unwrap();

    if let Some(scaled) = scaled {
        Selection::builder()
            .ksize(ksize.into())
            .scaled(scaled)
            .moltype(hash_function)
            .build()
    } else {
        Selection::builder()
            .ksize(ksize.into())
            .moltype(hash_function)
            .build()
    }
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
pub struct ManySearchResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub containment: f64,
    pub intersect_hashes: u64,
    pub ksize: u16,
    pub scaled: u32,
    pub moltype: String,
    pub match_md5: Option<String>,
    pub jaccard: Option<f64>,
    pub max_containment: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_abund: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub median_abund: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub std_abund: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub match_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub n_weighted_found: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_weighted_hashes: Option<u64>,
}

pub struct InterimGatherResult {
    intersect_bp: u64,
    f_orig_query: f64,
    f_match: f64,
    f_unique_to_query: f64,
    f_unique_weighted: f64,
    average_abund: f64,
    median_abund: f64,
    std_abund: f64,
    match_filename: String,
    match_name: String,
    match_md5: String,
    f_match_orig: f64,
    unique_intersect_bp: u64,
    gather_result_rank: u64,
    remaining_bp: u64,
    n_unique_weighted_found: u64,
    total_weighted_hashes: u64,
    sum_weighted_found: u64,
    query_containment_ani: f64,
    query_containment_ani_ci_low: Option<f64>,
    query_containment_ani_ci_high: Option<f64>,
    match_containment_ani: f64,
    match_containment_ani_ci_low: Option<f64>,
    match_containment_ani_ci_high: Option<f64>,
    average_containment_ani: f64,
    max_containment_ani: f64,
}

#[derive(Serialize)]
pub struct BranchwaterGatherResult {
    pub intersect_bp: u64,
    pub f_orig_query: f64,
    pub f_match: f64,
    pub f_unique_to_query: f64,
    pub f_unique_weighted: f64,
    pub average_abund: f64,
    pub median_abund: f64,
    pub std_abund: f64,
    pub match_filename: String,
    pub match_name: String,
    pub match_md5: String,
    pub f_match_orig: f64,
    pub unique_intersect_bp: u64,
    pub gather_result_rank: u32,
    pub remaining_bp: u64,
    pub query_filename: String,
    pub query_name: String,
    pub query_md5: String,
    pub query_bp: u64,
    pub ksize: u16,
    pub moltype: String,
    pub scaled: u32,
    pub query_n_hashes: u64,
    pub query_abundance: bool,
    pub query_containment_ani: f64,
    pub match_containment_ani: f64,
    pub average_containment_ani: f64,
    pub max_containment_ani: f64,
    pub n_unique_weighted_found: u64,
    pub sum_weighted_found: u64,
    pub total_weighted_hashes: u64,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_containment_ani_ci_low: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_containment_ani_ci_high: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub match_containment_ani_ci_low: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub match_containment_ani_ci_high: Option<f64>,
}

#[derive(Serialize, Deserialize)]
pub struct MultiSearchResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub match_md5: String,
    pub containment: f64,
    pub max_containment: f64,
    pub jaccard: f64,
    pub intersect_hashes: f64,
    pub ksize: u16,
    pub scaled: u32,
    pub moltype: String,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub match_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_containment_ani: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_containment_ani: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub prob_overlap: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub prob_overlap_adjusted: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    // max_containment / prob_overlap -> Bigger means less likely to be random
    pub containment_adjusted: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    // logged version is easier to plot/prioritize
    pub containment_adjusted_log10: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub tf_idf_score: Option<f64>,
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
    pub scaled: ScaledType,
    pub seed: u32,
    pub is_protein: bool,
    pub is_dayhoff: bool,
    pub is_hp: bool,
    pub is_dna: bool,
}

impl Hash for Params {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.track_abundance.hash(state);
        self.num.hash(state);
        self.scaled.hash(state);
        self.seed.hash(state);
        self.is_protein.hash(state);
        self.is_dayhoff.hash(state);
        self.is_hp.hash(state);
        self.is_dna.hash(state);
    }
}

pub fn zipwriter_handle(
    recv: Receiver<Option<BuildCollection>>,
    output: String,
) -> JoinHandle<Result<()>> {
    std::thread::spawn(move || -> Result<()> {
        // Convert output to PathBuf
        let outpath: PathBuf = output.into();
        let file_writer = open_output_file(&outpath);

        let options = FileOptions::default()
            .compression_method(CompressionMethod::Stored)
            .unix_permissions(0o644)
            .large_file(true);

        let mut zip = ZipWriter::new(file_writer);
        let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();
        let mut zip_manifest = BuildManifest::new();

        // Process each incoming Option<BuildCollection>
        while let Ok(message) = recv.recv() {
            match message {
                Some(mut build_collection) => {
                    // Use BuildCollection's method to write signatures to the zip file
                    match build_collection.write_sigs_to_zip(
                        &mut zip,
                        &mut md5sum_occurrences,
                        &options,
                    ) {
                        Ok(_) => {
                            zip_manifest.extend_from_manifest(&build_collection.manifest);
                        }
                        Err(e) => {
                            let error = e.context("Error processing signature in BuildCollection");
                            eprintln!("Error: {}", error);
                            return Err(error);
                        }
                    }
                }
                None => {
                    // Finalize and write the manifest when None is received
                    println!("Writing manifest");
                    zip_manifest.write_manifest_to_zip(&mut zip, &options)?;
                    zip.finish()?;
                    break;
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
