use crate::utils::parse_params_str;
/// sketch: sketch a single sequence file.
use anyhow::{bail, Result};
use camino::Utf8Path as Path;
use needletail::{parse_fastx_file, parse_fastx_reader};
use std::fs::File;
use std::io::BufWriter;

pub fn singlesketch(
    input_filename: String,
    moltype: String,
    param_str: String,
    output: String,
) -> Result<()> {
    // Parse parameter string into params_vec
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };

    // Build signature templates based on parameters and molecule type
    let sig_templates = crate::manysketch::build_siginfo(&params_vec, &moltype);

    if sig_templates.is_empty() {
        bail!("No signatures to build for the given parameters.");
    }

    let mut sigs = sig_templates.clone();

    // Open FASTA file reader
    let mut reader = if input_filename == "-" {
        let stdin = std::io::stdin();
        parse_fastx_reader(stdin)?
    } else {
        parse_fastx_file(&input_filename)?
    };

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
            }
            Err(err) => eprintln!("Error while processing record: {:?}", err),
        }
    }

    // Set name and filename for signatures
    sigs.iter_mut().for_each(|sig| {
        sig.set_name(&input_filename);
        sig.set_filename(&input_filename);
    });

    // Write signatures to output file
    let outpath = Path::new(&output);
    let file = File::create(outpath)?;
    let mut writer = BufWriter::new(file);

    // Write in JSON format
    serde_json::to_writer(&mut writer, &sigs)?;

    Ok(())
}
