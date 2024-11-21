use crate::utils::buildutils::BuildCollection;
use anyhow::{bail, Result};

#[allow(clippy::too_many_arguments)]
pub fn fastagather(
    input_filename: String,
    input_moltype: String,
    threshold_bp: u64,
    selection: &Selection,
    gather_output: Option<String>,
    prefetch_output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    
    // to start, implement straightforward record --> sketch --> gather
    // other ideas:
    // - add full-file (lower resolution) prefetch first, to reduce search space
    // - parallelize and/or batch records

    // Build signature templates based on parsed parameters
    let sig_template_result = BuildCollection::from_selection(selection);
    let mut sig_templates = match sig_template_result {
        Ok(sig_templates) => sig_templates,
        Err(e) => {
            bail!("Failed to build template signatures: {}", e);
        }
    };

    if sigs.is_empty() {
        bail!("No signatures to build for the given parameters.");
    }

    let input_moltype = input_moltype.to_ascii_lowercase();

    let mut against_selection = selection;
    let scaled = query_mh.scaled();
    against_selection.set_scaled(scaled);

    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes = {
        let x = threshold_bp / scaled as u64;
        if x > 0 {
            x
        } else {
            1
        }
    };

    // load collection to match against.
    let against_collection = load_collection(
        &against_filepath,
        &against_selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;
    
    let failed_records = AtomicUsize::new(0);
    // open file and start iterating through sequences
    // Open fasta file reader
    let mut reader = match parse_fastx_file(filename) {
        Ok(r) => r,
        Err(err) => {
            bail!("Error opening file {}: {:?}", filename, err);
        }
    };
    
    // later: can we parallelize across records or sigs? Do we want to batch groups of records for improved gather efficiency?
    while let Some(record_result) = reader.next() {
        // clone sig_templates for use 
        sigs = sig_templates.clone();
        match record_result {
            Ok(record) => {
                if let Err(err) = sigs.build_singleton_sigs(
                    record,
                    input_moltype,
                    filename.to_string(),
                ) {
                    eprintln!(
                        "Error building signatures from file: {}, {:?}",
                        filename, err
                    );
                    failed_records.fetch_add(1, atomic::Ordering::SeqCst);
                }
                for (rec, query_sig) in sigs.iter(){
                    let query_md5 = query_sig.md5sum();
                    let query_mh = sig.minhash().expect("could not get minhash from sig");
                    let query_name = query_sig.name(); // this is actually just record.id --> so maybe don't get it from sig here?

                    // now do prefetch/gather
                    let prefetch_result = load_sketches_above_threshold(against_collection, &query_mh, threshold_hashes)?;
                    let matchlist = prefetch_result.0;
                    let skipped_paths = prefetch_result.1;
                    let failed_paths = prefetch_result.2;

                    if prefetch_output.is_some() {
                        write_prefetch(
                            query_filename.clone(),
                            query_name.clone(),
                            query_md5,
                            prefetch_output,
                            &matchlist,
                        )
                        .ok();
                    }

                    consume_query_by_gather(
                        query_name,
                        query_filename,
                        query_mh,
                        scaled as u32,
                        matchlist,
                        threshold_hashes,
                        gather_output,
                    )
                    .ok();
                }
            }
            Err(err) => eprintln!("Error while processing record: {:?}", err),
        }
    }
    Ok(())
}
