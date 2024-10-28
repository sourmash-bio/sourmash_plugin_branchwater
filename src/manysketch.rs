/// manysketch: massively parallel sketching of sequence files.
use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::utils::{load_fasta_fromfile, parse_params_str, sigwriter, Params};
use camino::Utf8Path as Path;
use needletail::parse_fastx_file;
use sourmash::cmd::ComputeParameters;
use sourmash::signature::Signature;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

pub fn build_siginfo(params: &[Params], input_moltype: &str) -> Vec<Signature> {
    let mut sigs = Vec::new();

    for param in params.iter().cloned() {
        match input_moltype {
            // if dna, only build dna sigs. if protein, only build protein sigs, etc
            "dna" | "DNA" if !param.is_dna => continue,
            "protein" if !param.is_protein && !param.is_dayhoff && !param.is_hp => continue,
            _ => (),
        }

        // Adjust ksize value based on the is_protein flag
        let adjusted_ksize = if param.is_protein || param.is_dayhoff || param.is_hp {
            param.ksize * 3
        } else {
            param.ksize
        };

        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.is_protein)
            .dna(param.is_dna)
            .dayhoff(param.is_dayhoff)
            .hp(param.is_hp)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        let sig = Signature::from_params(&cp);
        sigs.push(sig);
    }

    sigs
}

pub fn manysketch(
    filelist: String,
    param_str: String,
    output: String,
    singleton: bool,
    force: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let (fileinfo, n_fastas) = match load_fasta_fromfile(filelist, force) {
        Ok((file_info, n_fastas)) => (file_info, n_fastas),
        Err(e) => bail!("Could not load fromfile csv. Underlying error: {}", e),
    };

    // if no files to process, exit with error
    if n_fastas == 0 {
        bail!("No files to load, exiting.");
    }

    // if output doesn't end in zip, bail
    if Path::new(&output)
        .extension()
        .map_or(true, |ext| ext != "zip")
    {
        bail!("Output must be a zip file.");
    }

    // set up a multi-producer, single-consumer channel that receives Signature
    let (send, recv) =
        std::sync::mpsc::sync_channel::<Option<Vec<Signature>>>(rayon::current_num_threads());
    // need to use Arc so we can write the manifest after all sigs have written
    // let send = std::sync::Arc::new(send);

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = sigwriter(recv, output);

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };

    // iterate over filelist_paths
    let processed_fastas = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);
    let skipped_paths: AtomicUsize = AtomicUsize::new(0);

    // set reporting threshold at every 5% or every 1 fasta, whichever is larger)
    let reporting_threshold = std::cmp::max(n_fastas / 20, 1);

    let send_result = fileinfo
        .par_iter()
        .filter_map(|fastadata| {
            let name = &fastadata.name;
            let filenames = &fastadata.paths;
            let moltype = &fastadata.input_type;
            // build sig templates for these sketches from params, check if there are sigs to build
            let sig_templates = build_siginfo(&params_vec, moltype);
            // if no sigs to build, skip this iteration
            if sig_templates.is_empty() {
                skipped_paths.fetch_add(filenames.len(), atomic::Ordering::SeqCst);
                processed_fastas.fetch_add(1, atomic::Ordering::SeqCst);
                return None;
            }

            let mut sigs = sig_templates.clone();
            // have name / filename been set for each sig yet?
            let mut set_name = false;
            // if merging multiple files, sourmash sets filename as last filename
            let last_filename = filenames.last().unwrap();

            for filename in filenames {
                // increment processed_fastas counter; make 1-based for % reporting
                let i = processed_fastas.fetch_add(1, atomic::Ordering::SeqCst);
                // progress report at threshold
                if (i + 1) % reporting_threshold == 0 {
                    let percent_processed = (((i + 1) as f64 / n_fastas as f64) * 100.0).round();
                    eprintln!(
                        "Starting file {}/{} ({}%)",
                        (i + 1),
                        n_fastas,
                        percent_processed
                    );
                }

                // Open fasta file reader
                let mut reader = match parse_fastx_file(filename) {
                    Ok(r) => r,
                    Err(err) => {
                        eprintln!("Error opening file {}: {:?}", filename, err);
                        failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                        return None;
                    }
                };

                // parse fasta and add to signature
                while let Some(record_result) = reader.next() {
                    match record_result {
                        Ok(record) => {
                            // do we need to normalize to make sure all the bases are consistently capitalized?
                            // let norm_seq = record.normalize(false);
                            sigs.iter_mut().for_each(|sig| {
                                if singleton {
                                    let record_name = std::str::from_utf8(record.id())
                                        .expect("could not get record id");
                                    sig.set_name(record_name);
                                    sig.set_filename(filename.as_str());
                                } else if !set_name {
                                    sig.set_name(name);
                                    // sourmash sets filename to last filename if merging fastas
                                    sig.set_filename(last_filename.as_str());
                                };
                                if moltype == "protein" {
                                    sig.add_protein(&record.seq())
                                        .expect("Failed to add protein");
                                } else {
                                    sig.add_sequence(&record.seq(), true)
                                        .expect("Failed to add sequence");
                                    // if not force, panics with 'N' in dna sequence
                                }
                            });
                            if !set_name {
                                set_name = true;
                            }
                        }
                        Err(err) => eprintln!("Error while processing record: {:?}", err),
                    }
                    if singleton {
                        // write sigs immediately to avoid memory issues
                        if let Err(e) = send.send(Some(sigs.clone())) {
                            eprintln!("Unable to send internal data: {:?}", e);
                            return None;
                        }
                        sigs = sig_templates.clone();
                    }
                }
            }
            // if singleton sketches, they have already been written; only write aggregate sketches
            if singleton {
                None
            } else {
                Some(sigs)
            }
        })
        .try_for_each_with(
            send.clone(),
            |s: &mut std::sync::mpsc::SyncSender<Option<Vec<Signature>>>, sigs| {
                if let Err(e) = s.send(Some(sigs)) {
                    Err(format!("Unable to send internal data: {:?}", e))
                } else {
                    Ok(())
                }
            },
        );

    // Send None to sigwriter to signal completion + write manifest
    if let Err(e) = send.send(None) {
        eprintln!("Unable to send completion signal: {:?}", e);
    }
    // do some cleanup and error handling -
    if let Err(e) = send_result {
        eprintln!("Error during parallel processing: {}", e);
    }

    // join the writer thread
    if let Err(e) = thrd
        .join()
        .unwrap_or_else(|e| Err(anyhow!("Thread panicked: {:?}", e)))
    {
        eprintln!("Error in sigwriter thread: {:?}", e);
    }

    // done!
    let i: usize = processed_fastas.load(atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} fasta files", i);

    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if failed_paths == i {
        bail!("Could not load fasta files: no signatures created.");
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} fasta files failed to load. See error messages above.",
            failed_paths
        );
    }

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    if skipped_paths == i {
        bail!("No fasta files compatible with provided sketch parameters: no signatures created.");
    }
    if skipped_paths > 0 {
        eprintln!(
            "WARNING: {} fasta files skipped - no compatible signatures.",
            skipped_paths
        );
    }

    Ok(())
}
