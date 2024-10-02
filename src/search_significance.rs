// Functions to compute statisical signifiance of search results

use rayon::prelude::*;

use sourmash::ffi::signature::SourmashSignature;
use sourmash::sketch::minhash::KmerMinHash;
use logsumexp::LogSumExp;
use std::collections::{HashMap, HashSet};
use sourmash::Error;

use crate::utils::{
    SmallSignature
};

// #[cfg(feature = "maths")]
pub fn get_hash_frequencies<'a>(
    hashvals: &'a Vec<u64>,
    minhash: &KmerMinHash,
    logged: bool,
    // Output hashmap borrows the hashvalues from intersection input
) -> HashMap<&'a u64, f64> {

// pub fn get_hash_frequency(
// hashvals: &Vec<u64>,
// minhash: &KmerMinHash,
// logged: bool,
// // Output hashmap borrows the hashvalues from hashvals input
// ) -> HashMap<u64, f64> {
    
    let sum_abunds: f64 = minhash.sum_abunds() as f64;
    let minhash_abunds: HashMap<u64, f64> = minhash.to_vec_abunds().into_iter().map(|(hashval, abund)| (hashval, abund as f64)).collect();

    eprintln!("--- hashvals, length: {} ---", hashvals.len());
    // for value in intersection {
    //     eprintln!("{}", value);
    // }

    eprintln!("--- sum_abunds: {} ---", sum_abunds);

    eprintln!("--- minhash_abunds ---");
    for (key, value) in &minhash_abunds {
        if *value != 1.0 {
            // Print only abundances that are greater than 1
            eprintln!("{}:\t{}", key, value);
        }
    }

    let mut frequencies: HashMap<&u64, f64> = HashMap::from(hashvals
        .par_iter()
        .map(|hashval| 
            // TODO: add a match statement here to error out properly if the hashval was not found 
            // in the minhash_abunds for some reason (shouldn't happen but ... computers be crazy)
            (hashval, 
                minhash_abunds[hashval] / sum_abunds
            )
        ).collect::<HashMap<&u64, f64>>()
    );

    if logged {
        let _ = frequencies
            .values_mut()
            .map(| freq| 
                freq.ln()
            );
    } 

    return frequencies;
}

// #[cfg(feature = "maths")]
pub fn get_prob_overlap(intersection: &Vec<u64>, queries_merged_mh: &KmerMinHash, database_merged_mh: &KmerMinHash, logged: bool) -> f64 {
    let query_frequencies: HashMap<&u64, f64> = get_hash_frequencies(intersection, queries_merged_mh, logged);
    let database_frequencies: HashMap<&u64, f64> = get_hash_frequencies(intersection, database_merged_mh, logged);

    // It's not guaranteed to me that the MinHashes from the query and database are in the same order, so iterate over one of them
    // and use a hashmap to retrieve the frequency value of the other
    let mut prob_overlap: f64;
    if logged {
        prob_overlap = query_frequencies
            .iter()
            .map(|(hashval, freq)| freq + database_frequencies[hashval])
            .ln_sum_exp();
        // ln_sum_exp uses natural log
        // -> change of base to log 10 for interpretability
        prob_overlap = prob_overlap.ln() / 10_f64.ln();
    } else {
        prob_overlap = query_frequencies
            .par_iter()
            .map(|(hashval, freq)| freq * database_frequencies[hashval]).sum();
    }

    return prob_overlap;
}

// TODO: How to accept SourmashSignature objects? Signature.minhash is Option<&KmerMinHash>, 
// so it's not guaranteed for a SourmashSignature to have a minhash object. Is there a way to 
// only accept SourmashSignature objects that have `.minhash` present?
pub fn merge_all_minhashes(sigs: &Vec<SmallSignature>) -> Result<KmerMinHash, Error> {
    if sigs.is_empty() {
        eprintln!("Signature list is empty");
        std::process::exit(1);
    }

    eprintln!("\n--- in merge_all_minhashes ---");

    let first_sig = &sigs[0];

    // Use the first signature to instantiate the merging of all minhashes
    let mut combined_mh = KmerMinHash::new(
        first_sig.minhash.scaled().try_into().unwrap(),
        first_sig.minhash.ksize().try_into().unwrap(),
        first_sig.minhash.hash_function().clone(),
        // accessing first_sig.minhash.seed is private -> hardcode instead
        42u64,
        first_sig.minhash.abunds().is_some(),
        first_sig.minhash.num(),
    );

    let hashes_with_abund: Vec<(u64, u64)> = sigs
        .par_iter()
        .map(|sig| sig.minhash.to_vec_abunds())
        .flatten()
        .collect();

    // Maybe should be using add_many_with_abund?
    // `merge` is a pretty heavy operation
    for sig in sigs.iter() {
        // Merging of signatures happens in place with KmerMinHash objects
        // Rust Question: Does par_iter() make sense here or does that 
        // mess with the combined_mh state, since the merging happens in-place?
        eprintln!("Adding signature, sig.name: {}, sig.md5sum: {}, sig.location: {}", sig.name, sig.md5sum, sig.location);
        eprintln!("combined_mh.n_unique_kmers(): {}", combined_mh.n_unique_kmers());
        let _ = combined_mh.merge(&sig.minhash);
    }
    eprintln!("Final combined_mh.n_unique_kmers(): {}", combined_mh.n_unique_kmers());
    Ok(combined_mh)
}

pub fn get_inverse_document_frequency(hashval: u64, signatures: Vec<SmallSignature>) {
    // Implementation of tf-idf for hashvals and signatures
    // https://en.wikipedia.org/wiki/Tf%E2%80%93idf

    // Total number of documents in the corpus
    let n_signatures = signatures.len();

    // Number of documents where term t appears
    let n_sigs_with_hashval: f64 = signatures.par_iter().map(|&sig| -> f64 {
        match sig.minhash.mins().contains(&hashval) {
            true => 1.0,
            false => 0.0,
        }
    }).collect()<Vec<f64>>.sum();

    let inverse_document_frequency = n_signatures / n_sigs_with_hashval;
}