use anyhow::{bail, Context, Result};
use camino::Utf8Path as Path;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::utils::buildutils::BuildCollection;

pub fn singlesketch(
    input_filename: String,
    input_moltype: String,
    param_str: String,
    output: String,
    name: String,
) -> Result<()> {
    // parse params --> signature templates
    let sig_template_result = BuildCollection::from_param_str(param_str.as_str());
    let mut sigs = match sig_template_result {
        Ok(sigs) => sigs,
        Err(e) => {
            bail!("Failed to parse params string: {}", e);
        }
    };

    let input_moltype = input_moltype.to_ascii_lowercase();

    // Build signature templates based on parsed parameters and detected moltype
    if sigs.is_empty() {
        bail!("No signatures to build for the given parameters.");
    }

    // to do --> actually handle the Result/Err from these
    // also to do --> add sequence counting so can print?
    let sequence_count =
        sigs.build_sigs_from_file_or_stdin(&input_moltype, name, input_filename.clone())?;

    eprintln!(
        "calculated {} signatures for {} sequences in {}",
        sigs.size(),
        sequence_count,
        input_filename
    );
    // Counter for the number of sequences processed (u64)
    // let mut sequence_count: u64 = 0;

    // Check if the output is stdout or a file
    if output == "-" {
        // Write signatures to stdout
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        sigs.write_sigs_as_json(&mut handle)
            .context("Failed to write signatures to stdout")?;
        handle.flush().context("Failed to flush stdout")?;
    } else {
        // Write signatures to an output file
        let outpath = Path::new(&output);
        let file = File::create(outpath).context(format!("Failed to create file: {}", output))?;
        let mut writer = BufWriter::new(file);

        sigs.write_sigs_as_json(&mut writer)
            .context(format!("Failed to write signatures to file: {}", output))?;
    }

    Ok(())
}
