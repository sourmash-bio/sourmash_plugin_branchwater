/// sigcat: concatenate signatures into a single sourmash zip file
use anyhow::Result;
use pyo3::Python;
use sourmash::{collection::Collection, selection::Selection};

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8PathBuf as PathBuf;
use sourmash::manifest::{Manifest, Record};
use sourmash::prelude::Select;
use zip::write::{FileOptions, ZipWriter};
use zip::CompressionMethod;
use piz::ZipArchive;
use rayon::ThreadPoolBuilder;


use crate::utils::{load_collection, open_output_file, write_signature, ReportType, zipwriter_handle};


// Use piz to parallel-read signatures from a zip file and send BuildCollection objects to writer
pub fn zipreader_spawn(zip_path: &Utf8PathBuf, tx: Sender<Option<BuildCollection>>) -> Result<()> {
    let file = File::open(zip_path)?;
    let archive = ZipArchive::new(file)?;
    let pool = ThreadPoolBuilder::new().build()?;

    archive.entries_parallel(&pool, |name, mut entry| {
        if name.ends_with(".sig.gz") && name.starts_with("signatures/") {
            let mut buf = vec![];
            entry.read_to_end(&mut buf)?;
            let (_fmt, decompressed) = niffler::decompress(&buf)?;
            let sigs: Vec<Signature> = serde_json::from_slice(&decompressed)?;
            for sig in sigs {
                let record = BuildRecord::from(&sig);
                let mut build_coll = BuildCollection::new();
                build_coll.sigs.push(sig);
                build_coll.manifest.add_record(record);
                tx.send(Some(build_coll)).expect("send failed");
            }
        }
        Ok(())
    })?;

    tx.send(None).expect("failed to send final None");
    Ok(())
}


pub fn sig_cat(
    py: Python,
    input_sigs: String,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    output_path: String,
    verbose: bool,
) -> Result<()> {
    // Check if output_path ends with ".zip"
    if !output_path.ends_with(".zip") {
        return Err(anyhow::anyhow!("Output file must end with '.zip'"));
    }

    // Convert input string to path set
    let sig_paths: std::collections::HashSet<String> =
        input_sigs.split_whitespace().map(|s| s.to_string()).collect();

    // Load all input collections
    let (mut multi, _nloaded) = MultiCollection::load_set_of_paths(&sig_paths);

    if !selection.is_empty() {
        multi = multi.select(selection)?;
    }

    if multi.is_empty() {
        bail!("No signatures to concatenate, exiting.");
    }

    // Setup for writing zip
    let outpath: PathBuf = output_path.into();
    let writer = open_output_file(&outpath);
    let options = FileOptions::default()
        .compression_method(CompressionMethod::Stored)
        .unix_permissions(0o644)
        .large_file(true);
    let mut zip = ZipWriter::new(writer);

    // Track md5sums to handle duplicates
    let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();
    let mut build_manifest = BuildManifest::new();

    // info for progress reporting
    let processed = AtomicUsize::new(0);
    let total_sigs = multi.len();

    // Iterate over all items in MultiCollection
    for (coll, _idx, record) in multi.item_iter() {
        py.check_signals();

        // progress reporting
        let p = processed.load(Ordering::SeqCst);
        if verbose || total_sigs <= 100 {
            eprintln!("processed {} of {}", p, total_sigs);
        } else if processed % (total_sigs / 100).max(1) == 0 {
            eprintln!("processed {} of {} ({}%)", p, total_sigs, (p * 100) / total_sigs);
        }

        let sig = coll.sig_from_record(record).context("failed to load signature")?;
        // do we need to do a select here?? I think maybe yes? .select(selection)?
        // i think this is required for downsampling, right?

        let md5 = sig.md5sum();
        let count = md5sum_occurrences.entry(md5.clone()).or_insert(0);
        *count += 1;

        let sig_filename = if *count > 1 {
            format!("signatures/{}_{}.sig.gz", md5, count)
        } else {
            format!("signatures/{}.sig.gz", md5)
        };

        // Write signature to zip
        write_signature(&sig, &mut zip, options.clone(), &sig_filename)?;
         
        // Build new BuildRecord and update fields
        let mut brec = BuildRecord::from_record(record);
        brec.set_internal_location(Some(sig_filename.clone().into()));
        brec.set_md5(Some(md5.clone()));
        brec.set_md5short(Some(md5[..8].to_string()));
        brec.set_n_hashes(Some(sig.get_sketch()?.size()));
        brec.set_name(sig.name().map(|s| s.to_string()));
        brec.set_filename(sig.filename().map(|s| s.to_string()));
        brec.sequence_added = true;

        build_manifest.add_record(brec);
        processed.fetch_add(1, Ordering::SeqCst);
    }

    // Write manifest
    println!("Writing manifest...");
    build_manifest
        .write_manifest_to_zip(&mut zip, &options)
        .context("failed to write manifest")?;

    zip.finish()?;

    let total = collected_sigs.load(Ordering::SeqCst);
    eprintln!("Concatenated {} signatures into '{}'.", total, outpath);
    
    Ok(())
}
