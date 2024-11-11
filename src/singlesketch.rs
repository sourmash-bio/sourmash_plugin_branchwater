use crate::utils::parse_params_str;
use anyhow::{bail, Result};
use camino::Utf8Path as Path;
use needletail::{parse_fastx_file, parse_fastx_reader};
use std::fs::File;
use std::io::{self, BufWriter, Write};

pub fn singlesketch(
    input_filename: String,
    param_str: String,
    output: String,
    name: String,
) -> Result<()> {
    // Parse parameter string into params_vec
    let param_result = parse_params_str(param_str.clone());
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };

    let moltype = if param_str.contains("dna") {
        "dna"
    } else if param_str.contains("protein") {
        "protein"
    } else if param_str.contains("dayhoff") {
        "dayhoff"
    } else if param_str.contains("hp") {
        "hp"
    } else {
        bail!("Unrecognized molecule type in params string");
    };

    // Build signature templates based on parsed parameters and detected moltype
    let mut sigs = crate::manysketch::build_siginfo(&params_vec, moltype);

    if sigs.is_empty() {
        bail!("No signatures to build for the given parameters.");
    }

    // Open FASTA file reader
    let mut reader = if input_filename == "-" {
        let stdin = std::io::stdin();
        parse_fastx_reader(stdin)?
    } else {
        parse_fastx_file(&input_filename)?
    };

    // Counter for the number of sequences processed (u64)
    let mut sequence_count: u64 = 0;

    // Parse FASTA and add to signature
    while let Some(record_result) = reader.next() {
        match record_result {
            Ok(record) => {
                sigs.iter_mut().for_each(|sig| {
                    if moltype == "protein" {
                        sig.add_protein(&record.seq())
                            .expect("Failed to add protein");
                    } else {
                        sig.add_sequence(&record.seq(), true)
                            .expect("Failed to add sequence");
                    }
                });
                sequence_count += 1;
            }
            Err(err) => eprintln!("Error while processing record: {:?}", err),
        }
    }

    // Set name and filename for signatures
    sigs.iter_mut().for_each(|sig| {
        sig.set_name(&name); // Use the provided name
        sig.set_filename(&input_filename);
    });

    // Check if the output is stdout or a file
    if output == "-" {
        // Write signatures to stdout
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        serde_json::to_writer(&mut handle, &sigs)?;
        handle.flush()?;
    } else {
        // Write signatures to output file
        let outpath = Path::new(&output);
        let file = File::create(outpath)?;
        let mut writer = BufWriter::new(file);

        // Write in JSON format
        serde_json::to_writer(&mut writer, &sigs)?;
    }

    eprintln!(
        "calculated {} signatures for {} sequences in {}",
        sigs.len(),
        sequence_count,
        input_filename
    );

    Ok(())
}
