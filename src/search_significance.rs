// Functions to compute statisical signifiance of search results

use rayon::prelude::*;

use sourmash::sketch::minhash::KmerMinHash;
use logsumexp::LogSumExp;
use std::collections::HashMap;
use sourmash::Error;
use sourmash::signature::SigsTrait;


pub enum Normalization {
    // L1 norm is the equivalent of frequencies/probabilities, as the counts
    // are divided by the length of the Vec<> object, or mathematically, the
    // number of items in the Vec object, assuming unit length for each
    L1,

    // L2 norm divides the counts by the sum of squares of all counts
    // L2 norm is the Euclidean distance of the counts in N-dimensional vector space
    // When track_abundance=False, L1 and L2 norms are equivalent, since all "counts" are 1
    L2
}

use crate::utils::SmallSignature;

pub fn value_count(hashes_with_abund: Vec<(u64, u64)>) -> HashMap<u64, u64> {
    // Coun the number of times the abundances appear in a hashes with abundance output
    let mut value_counts: HashMap<u64, u64>  = HashMap::new();

    for (_key, value) in hashes_with_abund.iter() {
        // Count number of times we see each value
        *value_counts.entry(*value).or_insert(0) += 1;
    }

    for (value, count) in value_counts.iter() {
        eprintln!("value: {}, count: {}", value, count);
    }

    return value_counts;
}

// #[cfg(feature = "maths")]
pub fn get_hash_frequencies<'a>(
    hashvals: &'a Vec<u64>,
    minhash: &KmerMinHash,
    normalization: Option<Normalization>,
    logged: Option<bool>,
    // Output hashmap borrows the hashvalues from intersection input
) -> HashMap<&'a u64, f64> {

// pub fn get_hash_frequency(
// hashvals: &Vec<u64>,
// minhash: &KmerMinHash,
// logged: bool,
// // Output hashmap borrows the hashvalues from hashvals input
// ) -> HashMap<u64, f64> {
    
    let minhash_abunds: HashMap<u64, f64> = minhash
        .to_vec_abunds()
        .into_iter()
        .map(|(hashval, abund)| (hashval, abund as f64))
        .collect();

    let abund_normalization = match normalization {
        Some(Normalization::L1) => minhash.sum_abunds(),
        Some(Normalization::L2) => { 
            minhash_abunds
                .values()
                .map(|abund| abund * abund )
                .collect<Vec<f64>>()
                .sum() as f64
            }
        _ => eprintln!("Invalid Normalization")
    };

    eprintln!("--- hashvals, length: {} ---", hashvals.len());
    // for value in intersection {
    //     eprintln!("{}", value);
    // }

    eprintln!("--- abund_normalization: {} ---", abund_normalization);

    eprintln!("--- minhash_abunds ---");
    // for (key, value) in &minhash_abunds {
    //     if *value != 1.0 {
    //         // Print only abundances that are greater than 1
    //         eprintln!("{}:\t{}", key, value);
    //     } 
    // }

    let mut frequencies: HashMap<&u64, f64> = HashMap::from(hashvals
        .par_iter()
        .map(|hashval| 
            // TODO: add a match statement here to error out properly if the hashval was not found 
            // in the minhash_abunds for some reason (shouldn't happen but ... computers be crazy)
            (hashval, 
                minhash_abunds[hashval] / abund_normalization
            )
        ).collect::<HashMap<&u64, f64>>()
    );

    match logged {
        Some(true) => {
            let _ = frequencies
            .values_mut()
            .map(| freq| 
                freq.ln()
            );
        },
        _ => ()

    };

    return frequencies;
}

// #[cfg(feature = "maths")]
pub fn get_prob_overlap(
    hashvals: &Vec<u64>, 
    queries_merged_mh: &KmerMinHash, 
    database_merged_mh: &KmerMinHash, 
    logged: Option<bool>
) -> f64 {
    let query_frequencies: HashMap<&u64, f64> = get_hash_frequencies(
        hashvals, 
        queries_merged_mh,
        Some(Normalization::L1), 
        logged
    );
    let database_frequencies: HashMap<&u64, f64> = get_hash_frequencies(hashvals, database_merged_mh, Some(Normalization::L1), logged);

    // It's not guaranteed to me that the MinHashes from the query and database are in the same order, so iterate over one of them
    // and use a hashmap to retrieve the frequency value of the other
    let prob_overlap = match logged {
        Some(true) => {
            let mut p_overlap = query_frequencies
                .iter()
                .map(|(hashval, freq)| freq + database_frequencies[hashval])
                .ln_sum_exp();
            // ln_sum_exp uses natural log
            // -> change of base to log 10 for interpretability
            p_overlap = p_overlap.ln() / 10_f64.ln();
            p_overlap
        },
        Some(false) => {
            let p_overlap = query_frequencies
                .par_iter()
                .map(|(hashval, freq)| freq * database_frequencies[hashval]).sum();
            p_overlap
        },
        _ => 0.0,
    };

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
        first_sig.minhash.hash_function(),
        // accessing first_sig.minhash.seed is private -> hardcode instead
        first_sig.minhash.seed(),
        first_sig.minhash.track_abundance(),
        first_sig.minhash.num(),
    );

    let hashes_with_abund: Vec<(u64, u64)> = sigs
        .par_iter()
        .map(|sig| sig.minhash.to_vec_abunds())
        .flatten()
        .collect();

    eprintln!("combined_mh.track_abundance(): {}", combined_mh.track_abundance());
    eprintln!("Before adding any hashes: combined_mh.n_unique_kmers(): {}", combined_mh.n_unique_kmers());
    _ = combined_mh.add_many_with_abund(&hashes_with_abund);

    // // Maybe should be using add_many_with_abund?
    // // `merge` is a pretty heavy operation
    // for sig in sigs.iter() {
    //     // Merging of signatures happens in place with KmerMinHash objects
    //     // Rust Question: Does par_iter() make sense here or does that 
    //     // mess with the combined_mh state, since the merging happens in-place?
    //     eprintln!("Adding signature, sig.name: {}, sig.md5sum: {}, sig.location: {}", sig.name, sig.md5sum, sig.location);
    //     eprintln!("combined_mh.n_unique_kmers(): {}", combined_mh.n_unique_kmers());
    //     let _ = combined_mh.merge(&sig.minhash);
    // }
    // combined_mh.reset_md5sum();

    // eprintln!("Final combined_mh.n_unique_kmers(): {}", combined_mh.n_unique_kmers());
    eprintln!("Final combined_mh.mins().len(): {}", combined_mh.mins().len());

    let _  = value_count(combined_mh.to_vec_abunds());

    Ok(combined_mh)
}

pub fn get_inverse_document_frequency(hashval: u64, signatures: &Vec<SmallSignature>, smooth_idf: Option<bool>) -> f64 {
    // Inverse document frequency tells us how unique this hashval is to the query database
    // When the value is near 0, then this hashval appears in all signatures
    // When the value is very large, equal to the number of signatures, then the hashval is 
    // unique to a single signature

    // Total number of documents in the corpus
    let n_signatures = signatures.len() as f64;

    // Number of documents where term t appears
    let n_sigs_with_hashval: f64 = signatures.par_iter().map(|sig| -> f64 {
        match sig.minhash.mins().contains(&hashval) {
            true => 1.0,
            false => 0.0,
        }
    }).sum();

    let inverse_document_frequency = match smooth_idf {
        // Add 1 to not totally ignore terms that appear in all documents
        // scikit-learn documentation (assumed to implement best practices for document classification): 
        // > "The effect of adding “1” to the idf in the equation above is that terms with zero idf,
        // > i.e., terms that occur in all documents in a training set, will not be entirely ignored."
        // Source: https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
        Some(true) => ( (1.0 + n_signatures) / (n_sigs_with_hashval) ).ln() + 1.0,
        Some(false) => (n_signatures / (n_sigs_with_hashval) ).ln() + 1.0,
        _ => 1.0
    };

    return inverse_document_frequency;
}

pub fn get_term_frequency_inverse_document_frequency(
    hashvals: &Vec<u64>, 
    query: &SmallSignature, 
    againsts: &Vec<SmallSignature>,
    smooth_idf: Option<bool>,
) -> f64 {
    // Implementation of tf-idf for hashvals and signatures
    // https://en.wikipedia.org/wiki/Tf%E2%80%93idf
    // Square the abundances to use an L2 norm -> why?
    // Because this is the default setting in scikit-learn's battle-tested tf-idf methods:
    // https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
    // https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfVectorizer.html
    let term_frequencies: HashMap<&u64, f64> = get_hash_frequencies(
        hashvals, 
        &query.minhash, 
        Some(Normalization::L2), 
        Some(false),
    );

    let inverse_document_frequencies: HashMap<&u64, f64> = HashMap::from(hashvals
        .par_iter()
        .map(|hashval| 
            (hashval, get_inverse_document_frequency(*hashval, againsts, smooth_idf))
        ).collect::<HashMap<&u64, f64>>()
    );
        
    // Multiply each hashval's term frequency and inverse document frequency, and sum the products
    let tf_idf_score: f64 = term_frequencies.par_iter().map(
        |(hashval, term_frequency)| 
        term_frequency * inverse_document_frequencies[hashval]
    ).sum();

    return tf_idf_score;
}