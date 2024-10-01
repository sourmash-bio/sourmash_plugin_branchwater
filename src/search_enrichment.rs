use sourmash::Signature as SourmashSignature;

pub fn get_inverse_document_frequency(hashval: u64, signatures: Vec<SourmashSignature>) {
    // Implementation of tf-idf for hashvals and signatures
    // https://en.wikipedia.org/wiki/Tf%E2%80%93idf

    // Total number of documents in the corpus
    let n_signatures = signatures.len();

    // Number of documents where term t appears
    let n_sigs_with_hashval = signatures.par_iter().map(|&sig| {
        match {
            sig.mins().contains(hashval) => 1.0,
            _ => 0.0
        }
    }).collect().sum();

    let inverse_document_frequency = n_signatures / n_sigs_with_hashval;
}