/// sigcat: concatenate signatures into a single sourmash zip file
use anyhow::Result;
use camino::Utf8PathBuf;
use pyo3::Python;
use sourmash::{collection::Collection, selection::Selection};

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::utils::buildutils::{BuildCollection, BuildManifest, BuildRecord};
use crate::utils::multicollection::MultiCollection;
use rayon::iter::ParallelIterator;
use sourmash::prelude::Select;
use sourmash::signature::Signature;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;
use zip::write::FileOptions;
use zip::ZipWriter;

use std::sync::{Once, OnceLock};

static TERMINATE_FLAG: OnceLock<Arc<AtomicBool>> = OnceLock::new();
static CTRL_HANDLER_INIT: Once = Once::new();

pub fn setup_ctrlc_handler() -> Arc<AtomicBool> {
    let flag = TERMINATE_FLAG.get_or_init(|| Arc::new(AtomicBool::new(false)));

    CTRL_HANDLER_INIT.call_once(|| {
        let flag_clone = Arc::clone(flag);
        let _ = ctrlc::set_handler(move || {
            eprintln!("Received Ctrl-C, requesting shutdown...");
            flag_clone.store(true, Ordering::SeqCst);
        });
    });

    Arc::clone(flag)
}

pub fn precompressed_zipwriter_handle(
    recv: Receiver<Option<Vec<CompressedSig>>>,
    output: Utf8PathBuf,
    cancel: Arc<AtomicBool>,
) -> thread::JoinHandle<Result<()>> {
    thread::spawn(move || -> Result<()> {
        let outpath = output.clone();
        let incomplete_path = outpath.with_extension("zip.incomplete");
        let file_writer = std::fs::File::create(&incomplete_path)?;
        let mut zip = ZipWriter::new(file_writer);

        let options = FileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .unix_permissions(0o644)
            .large_file(true);

        let mut zip_manifest = BuildManifest::new();
        let mut wrote_any_sigs = false;

        while let Ok(message) = recv.recv() {
            if cancel.load(Ordering::SeqCst) {
                eprintln!("Termination requested, exiting early...");
                return Ok(()); // or early return / cleanup
            }
            match message {
                Some(batch) => {
                    if !wrote_any_sigs {
                        wrote_any_sigs = true;
                    }
                    for compressed in batch {
                        zip.start_file(compressed.filename, options)?;
                        zip.write_all(&compressed.data)?;
                        zip_manifest.add_record(compressed.record);
                    }
                }
                None => {
                    if wrote_any_sigs {
                        zip_manifest.write_manifest_to_zip(&mut zip, &options)?;
                        zip.finish()?;
                        std::fs::rename(&incomplete_path, &outpath)?;
                    } else {
                        drop(zip);
                        std::fs::remove_file(&incomplete_path).ok();
                    }
                    break;
                }
            }
        }

        Ok(())
    })
}

pub struct CompressedSig {
    pub filename: String,
    pub data: Vec<u8>,
    pub record: BuildRecord,
}

pub fn compress_batch(
    mut build_collection: BuildCollection,
    md5sum_occurrences: &Arc<Mutex<HashMap<String, usize>>>,
) -> Result<Vec<CompressedSig>> {
    let mut output = Vec::new();
    for (record, sig) in &build_collection {
        if !record.sequence_added {
            continue;
        }

        let compressed = compress_sig(record.clone(), sig, md5sum_occurrences)?;
        output.push(compressed);
    }

    Ok(output)
}

pub fn compress_sig(
    mut record: BuildRecord,
    sig: &Signature,
    md5sum_occurrences: &Arc<Mutex<HashMap<String, usize>>>,
) -> Result<CompressedSig> {
    let md5sum_str = sig.md5sum();
    let sig_filename = {
        let mut md5sums = md5sum_occurrences.lock().unwrap();
        let count = md5sums.entry(md5sum_str.clone()).or_insert(0);
        *count += 1;

        if *count > 1 {
            format!("signatures/{}_{}.sig.gz", md5sum_str, count)
        } else {
            format!("signatures/{}.sig.gz", md5sum_str)
        }
    };

    record.set_internal_location(Some(sig_filename.clone().into()));

    let wrapped_sig = vec![sig.clone()];
    let json_bytes = serde_json::to_vec(&wrapped_sig)?;

    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )?;
            gz_writer.write_all(&json_bytes)?;
            gz_writer.flush()?;
        }
        buffer.into_inner()
    };

    Ok(CompressedSig {
        filename: sig_filename,
        data: gzipped_buffer,
        record,
    })
}

// batched reading and sending to the writer thread
pub fn zipreader_spawn(
    zip_path: &Utf8PathBuf,
    tx: SyncSender<Option<Vec<CompressedSig>>>,
    selection: &Selection,
    batch_size: usize,
    verbose: bool,
    md5sum_occurrences: Arc<Mutex<HashMap<String, usize>>>,
    cancel: Arc<AtomicBool>,
) -> Result<usize> {
    let collection = Collection::from_zipfile(zip_path.clone())?;
    let manifest = collection.manifest().clone();
    let selected = manifest.select(selection)?;

    let total = selected.iter().count();
    let processed = AtomicUsize::new(0);
    let mut batch = BuildCollection::new();

    let mut final_count = 0;
    let mut scope_result: Result<()> = Ok(());

    rayon::scope(|s| {
        for record in selected.iter() {
            if cancel.load(Ordering::SeqCst) {
                eprintln!("Termination requested, exiting early...");
                scope_result = Ok(()); // or early return / cleanup
                return;
            }

            let sig = match collection.sig_from_record(record) {
                Ok(s) => s,
                Err(e) => {
                    scope_result = Err(e.into());
                    return;
                }
            };

            let build_rec = BuildRecord::from_record(record);

            batch.sigs.push(sig.into());
            batch.manifest.add_record(build_rec);

            let count = processed.fetch_add(1, Ordering::SeqCst) + 1;
            final_count = count;

            if verbose || total <= 100 || count % (total / 100).max(1) == 0 {
                println!(
                    "{zip_path}: processed {} of {} ({}%)",
                    count,
                    total,
                    (count * 100) / total
                );
            }

            if batch.sigs.len() >= batch_size {
                let to_compress = std::mem::take(&mut batch);
                let md5sums = Arc::clone(&md5sum_occurrences);
                let tx_clone = tx.clone();
                // spawn a new thread to compress and send to writer thread
                rayon::spawn_fifo(move || match compress_batch(to_compress, &md5sums) {
                    Ok(compressed) => {
                        let _ = tx_clone.send(Some(compressed));
                    }
                    Err(e) => eprintln!("Compression failed: {e}"),
                });
            }
        }

        if !batch.sigs.is_empty() {
            if let Ok(compressed) = compress_batch(batch, &md5sum_occurrences) {
                let _ = tx.send(Some(compressed));
            }
        }
    });

    scope_result?;
    eprintln!(
        "finished reading {}: found {} matching signatures",
        zip_path, total
    );
    Ok(final_count)
}

// Handle non-zip inputs using MultiCollection and rayon
pub fn multicollection_reader(
    input_paths: &[Utf8PathBuf],
    tx: SyncSender<Option<Vec<CompressedSig>>>,
    selection: &Selection,
    batch_size: usize,
    verbose: bool,
    md5sum_occurrences: Arc<Mutex<HashMap<String, usize>>>,
    cancel: Arc<AtomicBool>,
) -> Result<usize> {
    let pathset: HashSet<String> = input_paths.iter().map(|p| p.to_string()).collect();
    let (mut multi, _nfailed) = MultiCollection::load_set_of_paths(&pathset);
    multi = multi.select(selection)?;

    let total = multi.len();
    let processed = AtomicUsize::new(0);
    let batch_accumulator = Arc::new(Mutex::new(Vec::with_capacity(batch_size)));

    multi.par_iter().for_each(|(coll, _idx, record)| {
        if cancel.load(Ordering::SeqCst) {
            return;
        }

        let sig = match coll.sig_from_record(record) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Signature read error: {e}");
                return;
            }
        };

        let build_rec = BuildRecord::from_record(record);

        let compressed = match compress_sig(build_rec, &sig, &md5sum_occurrences) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("Compression error: {e}");
                return;
            }
        };

        let mut batch = batch_accumulator.lock().unwrap();
        batch.push(compressed);

        let count = processed.fetch_add(batch_size, Ordering::SeqCst);
        if verbose || total <= 100 || count % (total / 100).max(1) == 0 {
            println!(
                "non-zips: processed {} of {} ({}%)",
                count,
                total,
                (count * 100) / total
            );
        }

        if batch.len() >= batch_size {
            let to_send = std::mem::take(&mut *batch);
            let _ = tx.send(Some(to_send));
        }
    });

    // After all: send any leftovers
    if let Ok(mut leftover) = batch_accumulator.lock() {
        if !leftover.is_empty() {
            let _ = tx.send(Some(std::mem::take(&mut *leftover)));
        }
    }

    Ok(total)
}

// Expand pathlists (.txt files) into a flat list of paths
fn expand_input_paths(paths: Vec<String>) -> Result<Vec<Utf8PathBuf>> {
    let mut expanded = Vec::new();
    for path in paths {
        let path = Utf8PathBuf::from(path);
        if path.extension().map_or(false, |e| e == "txt") {
            let content = std::fs::read_to_string(&path)?;
            for line in content.lines().filter(|l| !l.is_empty()) {
                expanded.push(Utf8PathBuf::from(line));
            }
        } else {
            expanded.push(path);
        }
    }
    Ok(expanded)
}

pub fn expand_and_partition_inputs(
    inputs: Vec<String>,
) -> Result<(Vec<Utf8PathBuf>, Vec<Utf8PathBuf>)> {
    let paths = expand_input_paths(inputs)?;
    let (zip_inputs, other_inputs): (Vec<_>, Vec<_>) = paths
        .into_iter()
        .partition(|p| p.extension().map_or(false, |e| e == "zip"));
    Ok((zip_inputs, other_inputs))
}

pub fn sig_cat(
    py: Python,
    inputs: Vec<String>,
    output: String,
    selection: &Selection,
    batch_size: usize,
    verbose: bool,
) -> Result<()> {
    // Check if output_path ends with ".zip"
    if !output.ends_with(".zip") {
        return Err(anyhow::anyhow!("Output file must end with '.zip'"));
    }

    // init channels and writer thread
    let (tx, rx): (
        SyncSender<Option<Vec<CompressedSig>>>,
        Receiver<Option<Vec<CompressedSig>>>,
    ) = sync_channel(rayon::current_num_threads());
    // let writer_handle = zipwriter_handle(rx, output.clone());
    let cancel_flag = setup_ctrlc_handler();
    let writer_handle =
        precompressed_zipwriter_handle(rx, output.clone().into(), cancel_flag.clone());

    let total_written = std::sync::Arc::new(AtomicUsize::new(0));
    // flatten input paths and split into zip / non-zip
    let (zip_inputs, other_inputs) = expand_and_partition_inputs(inputs)?;

    eprintln!(
        "Found {} zip files and {} other files to process.",
        zip_inputs.len(),
        other_inputs.len()
    );

    py.check_signals()?;

    // set up writer stuff
    let md5sum_occurrences = Arc::new(Mutex::new(HashMap::new()));

    // spawn processing
    rayon::scope(|s| {
        for zip_path in &zip_inputs {
            let tx = tx.clone();
            let selection = selection.clone();
            let total_written = total_written.clone();
            let md5sums = Arc::clone(&md5sum_occurrences);
            let cancel = cancel_flag.clone();
            s.spawn(move |_| {
                if let Ok(n) = zipreader_spawn(
                    &zip_path, tx, &selection, batch_size, verbose, md5sums, cancel,
                ) {
                    total_written.fetch_add(n, Ordering::SeqCst);
                }
            });
        }

        if !other_inputs.is_empty() {
            let tx = tx.clone();
            let selection = selection.clone();
            let total_written = total_written.clone();
            let md5sums = Arc::clone(&md5sum_occurrences);
            let cancel = cancel_flag.clone();
            s.spawn(move |_| {
                if let Ok(n) = multicollection_reader(
                    &other_inputs,
                    tx,
                    &selection,
                    batch_size,
                    verbose,
                    md5sums,
                    cancel,
                ) {
                    total_written.fetch_add(n, Ordering::SeqCst);
                }
            });
        }
    });

    // After all reading threads finish, send None to signal completion (and write the manifest)
    tx.send(None).expect("failed to send final None");
    // Now wait for the writer thread to finish
    writer_handle.join().expect("writer thread panicked")?;

    if total_written.load(Ordering::SeqCst) == 0 {
        bail!("No signatures could be written to the output file.");
    }

    eprintln!(
        "Concatenated {} signatures into '{}'.",
        total_written.load(Ordering::SeqCst),
        output
    );

    Ok(())
}
