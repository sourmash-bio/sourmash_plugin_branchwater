use crate::utils::buildutils::BuildCollection;
use anyhow::{bail, Result};

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

    let sequence_count =
        sigs.build_sigs_from_file_or_stdin(&input_moltype, name, input_filename.clone())?;

    eprintln!(
        "calculated {} signatures for {} sequences in {}",
        sigs.size(),
        sequence_count,
        input_filename
    );

    // Write signatures to stdout or output file
    sigs.write_sigs(&output)?;

    Ok(())
}
