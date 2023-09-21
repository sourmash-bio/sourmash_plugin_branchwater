/// manysketch: massively parallel sketching of sequence files.
use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::utils::{load_fasta_fromfile, sigwriter, Params, ZipMessage};
use needletail::parse_fastx_file;
use sourmash::cmd::ComputeParameters;
use sourmash::signature::Signature;
use std::path::Path;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

fn parse_params_str(params_strs: String) -> Result<Vec<Params>, String> {
    let mut unique_params: std::collections::HashSet<Params> = std::collections::HashSet::new();

    // split params_strs by _ and iterate over each param
    for p_str in params_strs.split('_').collect::<Vec<&str>>().iter() {
        let items: Vec<&str> = p_str.split(',').collect();

        let mut ksizes = Vec::new();
        let mut track_abundance = false;
        let mut num = 0;
        let mut scaled = 1000;
        let mut seed = 42;
        let mut is_protein = false;
        let mut is_dna = true;

        for item in items.iter() {
            match *item {
                _ if item.starts_with("k=") => {
                    let k_value = item[2..]
                        .parse()
                        .map_err(|_| format!("cannot parse k='{}' as a number", &item[2..]))?;
                    ksizes.push(k_value);
                }
                "abund" => track_abundance = true,
                "noabund" => track_abundance = false,
                _ if item.starts_with("num=") => {
                    num = item[4..]
                        .parse()
                        .map_err(|_| format!("cannot parse num='{}' as a number", &item[4..]))?;
                }
                _ if item.starts_with("scaled=") => {
                    scaled = item[7..]
                        .parse()
                        .map_err(|_| format!("cannot parse scaled='{}' as a number", &item[7..]))?;
                }
                _ if item.starts_with("seed=") => {
                    seed = item[5..]
                        .parse()
                        .map_err(|_| format!("cannot parse seed='{}' as a number", &item[5..]))?;
                }
                "protein" => {
                    is_protein = true;
                    is_dna = false;
                }
                "dna" => {
                    is_protein = false;
                    is_dna = true;
                }
                _ => return Err(format!("unknown component '{}' in params string", item)),
            }
        }

        for &k in &ksizes {
            let param = Params {
                ksize: k,
                track_abundance,
                num,
                scaled,
                seed,
                is_protein,
                is_dna,
            };
            unique_params.insert(param);
        }
    }

    Ok(unique_params.into_iter().collect())
}

fn build_siginfo(
    params: &[Params],
    moltype: &str,
    name: &str,
    filename: &Path,
) -> (Vec<Signature>, Vec<Params>) {
    let mut sigs = Vec::new();
    let mut params_vec = Vec::new();

    for param in params.iter().cloned() {
        match moltype {
            // if dna, only build dna sigs. if protein, only build protein sigs
            "dna" if !param.is_dna => continue,
            "protein" if !param.is_protein => continue,
            _ => (),
        }

        // Adjust ksize value based on the is_protein flag
        let adjusted_ksize = if param.is_protein {
            param.ksize * 3
        } else {
            param.ksize
        };

        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.is_protein)
            .dna(param.is_dna)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        // let sig = Signature::from_params(&cp); // cant set name with this
        let template = sourmash::cmd::build_template(&cp);
        let sig = Signature::builder()
            .hash_function("0.murmur64")
            .name(Some(name.to_string()))
            .filename(Some(filename.to_string_lossy().into_owned()))
            .signatures(template)
            .build();
        sigs.push(sig);

        params_vec.push(param);
    }

    (sigs, params_vec)
}

pub fn manysketch<P: AsRef<Path> + Sync>(
    filelist: P,
    param_str: String,
    output: String,
) -> Result<(), Box<dyn std::error::Error>> {
    let fileinfo = match load_fasta_fromfile(&filelist) {
        Ok(result) => result,
        Err(e) => bail!("Could not load fromfile csv. Underlying error: {}", e),
    };

    // if no files to process, exit with error
    let n_fastas = fileinfo.len();
    if n_fastas == 0 {
        bail!("No files to load, exiting.");
    }

    // if output doesnt end in zip, bail
    if Path::new(&output)
        .extension()
        .map_or(true, |ext| ext != "zip")
    {
        bail!("Output must be a zip file.");
    }

    // set up a multi-producer, single-consumer channel that receives Signature
    let (send, recv) = std::sync::mpsc::sync_channel::<ZipMessage>(rayon::current_num_threads());
    // need to use Arc so we can write the manifest after all sigs have written
    let send = std::sync::Arc::new(send);

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = sigwriter::<&str>(recv, output);

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
        .filter_map(|(name, filename, moltype)| {
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

            // build sig templates from params
            let (mut sigs, sig_params) = build_siginfo(&params_vec, moltype, name, filename);
            // if no sigs to build, skip
            if sigs.is_empty() {
                let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                return None;
            }

            // Open fasta file reader
            let mut reader = match parse_fastx_file(filename) {
                Ok(r) => r,
                Err(err) => {
                    eprintln!("Error opening file {}: {:?}", filename.display(), err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    return None;
                }
            };
            // parse fasta and add to signature
            while let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        // do we need to normalize to make sure all the bases are consistently capitalized?
                        // let norm_seq = record.normalize(false);
                        for sig in &mut sigs {
                            if moltype == "protein" {
                                sig.add_protein(&record.seq()).unwrap();
                            } else {
                                sig.add_sequence(&record.seq(), true).unwrap();
                                // if not force, panics with 'N' in dna sequence
                            }
                        }
                    }
                    Err(err) => {
                        eprintln!("Error while processing record: {:?}", err);
                    }
                }
            }
            Some((sigs, sig_params, filename))
        })
        .try_for_each_with(
            send.clone(),
            |s: &mut std::sync::Arc<std::sync::mpsc::SyncSender<ZipMessage>>,
             (sigs, sig_params, filename)| {
                if let Err(e) = s.send(ZipMessage::SignatureData(
                    sigs,
                    sig_params,
                    filename.clone(),
                )) {
                    Err(format!("Unable to send internal data: {:?}", e))
                } else {
                    Ok(())
                }
            },
        );

    // After the parallel work, send the WriteManifest message
    std::sync::Arc::try_unwrap(send)
        .unwrap()
        .send(ZipMessage::WriteManifest)
        .unwrap();

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
