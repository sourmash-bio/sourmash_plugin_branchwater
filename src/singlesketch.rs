use crate::utils::buildutils::BuildCollection;
use anyhow::{bail, Result};

pub fn singlesketch(
    input_filenames: Vec<String>,
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
    // print sig templates to build
    let _params = sigs.summarize_params();

    let input_moltype = input_moltype.to_ascii_lowercase();

    // Build signature templates based on parsed parameters and detected moltype
    if sigs.is_empty() {
        bail!("No signatures to build for the given parameters.");
    }

    let mut sequence_count = 0;
    for input_filename in input_filenames.iter() {
        sequence_count += sigs.build_sigs_from_file_or_stdin(
            &input_moltype,
            name.clone(),
            input_filename.clone(),
        )?;
    }

    eprintln!(
        "calculated {} signatures for {} sequences in {} files",
        sigs.size(),
        sequence_count,
        input_filenames.len(),
    );

    // Write signatures to stdout or output file
    sigs.write_sigs(&output)?;

    Ok(())
}
