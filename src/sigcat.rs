/// sigcat: concatenate signatures into a single sourmash zip file
use anyhow::Result;
use camino::Utf8PathBuf;
use pyo3::Python;
use sourmash::{collection::Collection, selection::Selection};

use std::sync::atomic::{AtomicUsize, Ordering};

// use sourmash::manifest::{Manifest, Record};
use crate::utils::buildutils::{BuildCollection, BuildRecord};
use crate::utils::multicollection::MultiCollection;
use rayon::iter::ParallelIterator;
use sourmash::prelude::Select;
use std::sync::mpsc::{sync_channel, SyncSender};

use crate::utils::zipwriter_handle;

// batched reading and sending to the writer thread
pub fn zipreader_spawn(
    zip_path: &Utf8PathBuf,
    tx: SyncSender<Option<BuildCollection>>,
    selection: &Selection,
    batch_size: usize,
    verbose: bool,
) -> Result<()> {
    let collection = Collection::from_zipfile(zip_path)?;
    let manifest = collection.manifest().clone();
    let selected = manifest.select(selection)?;

    let total = selected.iter().count();
    let processed = AtomicUsize::new(0);
    let mut batch = BuildCollection::new();

    for record in selected.iter() {
        let sig = collection.sig_from_record(record)?;
        let build_rec = BuildRecord::from_record(record);

        // do we need to select on the sig to downsample?? if so, need to update build rec?

        batch.sigs.push(sig.into());
        batch.manifest.add_record(build_rec);

        let count = processed.fetch_add(1, Ordering::SeqCst) + 1;
        if verbose || total <= 100 {
            eprintln!("processed {} of {}", count, total);
        } else if count % (total / 100).max(1) == 0 {
            let percent = (count * 100) / total;
            eprintln!("... {}% done", percent);
        }

        if batch.sigs.len() >= batch_size {
            tx.send(Some(std::mem::take(&mut batch)))?;
        }
    }

    if !batch.sigs.is_empty() {
        tx.send(Some(batch))?;
    }

    Ok(())
}

// Handle non-zip inputs using MultiCollection and rayon
pub fn multicollection_reader(
    input_paths: &[Utf8PathBuf],
    tx: SyncSender<Option<BuildCollection>>,
    selection: &Selection,
    batch_size: usize,
    verbose: bool,
) -> Result<()> {
    let pathset: std::collections::HashSet<String> =
        input_paths.iter().map(|p| p.to_string()).collect();

    let (mut multi, _nloaded) = MultiCollection::load_set_of_paths(&pathset);
    multi = multi.select(selection)?;

    let total = multi.len();
    let processed = AtomicUsize::new(0);
    let batch = std::sync::Mutex::new(BuildCollection::new());

    multi
        .par_iter()
        .try_for_each(|(coll, _idx, record)| -> Result<()> {
            let sig = coll.sig_from_record(record)?;
            let record = BuildRecord::from_record(record);
            let count = processed.fetch_add(1, Ordering::SeqCst) + 1;

            {
                let mut batch = batch.lock().unwrap();
                batch.sigs.push(sig.into());
                batch.manifest.add_record(record);

                if batch.sigs.len() >= batch_size {
                    tx.send(Some(std::mem::take(&mut *batch)))
                        .expect("send failed");
                }
            }

            if verbose || total <= 100 {
                eprintln!("processed {} of {}", count, total);
            } else if count % (total / 100).max(1) == 0 {
                let percent = (count * 100) / total;
                eprintln!("... {}% done", percent);
            }

            Ok::<(), anyhow::Error>(())
        })?;
    let mut batch = batch.lock().unwrap();
    if !batch.sigs.is_empty() {
        tx.send(Some(std::mem::take(&mut *batch)))
            .expect("send failed");
    }

    Ok(())
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
    let (tx, rx) = sync_channel::<Option<BuildCollection>>(rayon::current_num_threads());
    let writer_handle = zipwriter_handle(rx, output.clone());

    eprintln!("spawned writer thread...");

    // flatten input paths and split into zip / non-zip
    let (zip_inputs, other_inputs) = expand_and_partition_inputs(inputs)?;

    py.check_signals()?;

    // did we spawn a thread?
    let did_spawn = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));

    // spawn processing
    rayon::scope(|s| {
        for zip_path in &zip_inputs {
            let tx = tx.clone();
            let selection = selection.clone();
            s.spawn(move |_| {
                if let Err(e) = zipreader_spawn(zip_path, tx, &selection, batch_size, verbose) {
                    eprintln!("Error in zipreader_spawn: {}", e);
                }
            });
        }

        if !other_inputs.is_empty() {
            let tx = tx.clone();
            let selection = selection.clone();
            s.spawn(move |_| {
                if let Err(e) =
                    multicollection_reader(&other_inputs, tx, &selection, batch_size, verbose)
                {
                    eprintln!("Error in multicollection_reader: {}", e);
                }
            });
        }
    });

    // After all reading threads finish, send None to signal completion (and write the manifest)
    tx.send(None).expect("failed to send final None");
    // drop(tx);
    // Now wait for the writer thread to finish
    writer_handle.join().expect("writer thread panicked")?;

    Ok(())
}
