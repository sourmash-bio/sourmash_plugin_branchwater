// Functions to compute statisical signifiance of search results

use rayon::prelude::*;

use sourmash::sketch::minhash::KmerMinHash;
use logsumexp::LogSumExp;
use std::collections::{HashMap, HashSet};
use sourmash::Error;
use sourmash::signature::SigsTrait;
use std::fmt::{self, Display, Formatter};
use crate::utils::SmallSignature;

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

impl Display for Normalization {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::L1 => write!(f, "L1 Normalization"),
            Self::L2 => write!(f, "L2 Normalization"),
        }
    }
}



pub fn print_value_count(hashes_with_abund: Vec<(u64, u64)>) -> HashMap<u64, u64> {
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

pub fn print_freq_count(frequencies: HashMap<&u64, f64>, freq_type: &str, query_name: &str, n_hashes: usize) -> HashMap<String, u64> {
    // Coun the number of times the abundances appear in a hashes with abundance output
    let mut freq_counts: HashMap<String, u64>  = HashMap::new();

    for (_key, value) in frequencies.iter() {
        // Count number of times we see each value
        *freq_counts.entry(value.to_string()).or_insert(0) += 1;
    }

    for (freq, count) in freq_counts.iter() {
        eprintln!("{}\t{}\tn_hashes: {}\tfreq: {}, count: {}", freq_type, query_name, n_hashes, freq, count);
    }

    return freq_counts;
}

// #[cfg(feature = "maths")]
pub fn get_hash_frequencies<'a>(
    minhash: &KmerMinHash,
    normalization: Option<Normalization>,
) -> HashMap<u64, f64> {

    let minhash_abunds: HashMap<u64, f64> = minhash
        .to_vec_abunds()
        .into_par_iter()
        .map(|(hashval, abund)| (hashval, abund as f64))
        .collect();

    let abund_normalization: f64 = match normalization {
        Some(Normalization::L1) => minhash.sum_abunds() as f64,
        Some(Normalization::L2) => { 
            minhash_abunds
                .par_iter()
                .map(|(hashval, abund)| abund * abund )
                .sum::<f64>() as f64
            }
        // TODO: this should probably be an error
        _ => 0.0,
    };

    // eprintln!("abund_normalization ({}): {}", normalization.unwrap(), abund_normalization);

    let mut frequencies: HashMap<u64, f64> = HashMap::from(
        minhash_abunds
        .par_iter()
        .map(|(hashval, abund)| 
            // TODO: add a match statement here to error out properly if the hashval was not found 
            // in the minhash_abunds for some reason (shouldn't happen but ... computers be crazy)
            (
                *hashval, 
                abund / abund_normalization
            )
        ).collect::<HashMap<u64, f64>>()
    );

    return frequencies;
}

// #[cfg(feature = "maths")]
pub fn get_prob_overlap(
    hashvals: &Vec<u64>, 
    query_frequencies: &HashMap<u64, f64>,
    against_frequencies: &HashMap<u64, f64>
) -> f64 {

    // It's not guaranteed to me that the MinHashes from the query and database are in the same order, so iterate over one of them
    // and use a hashmap to retrieve the frequency value of the other
    let prob_overlap =  hashvals
                .par_iter()
                .map(|hashval| 
                    query_frequencies[hashval] * against_frequencies[hashval]
                ).sum();

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

    _ = combined_mh.add_many_with_abund(&hashes_with_abund);

    Ok(combined_mh)
}

pub fn compute_inverse_document_frequency(
    against_merged_mh: &KmerMinHash,
    againsts: &Vec<SmallSignature>, 
    smooth_idf: Option<bool>
) -> HashMap<u64, f64> {
    // Compute inverse document frequency (IDF) of all 
    // Inverse document frequency tells us how unique this hashval is to the query database
    // When the value is near 0, then this hashval appears in all signatures
    // When the value is very large, equal to the number of signatures, then the hashval is 
    // unique to a single signature

    // Total number of documents in the corpus
    let n_signatures = againsts.len() as f64;

    let againsts_hashes: Vec<HashSet<&u64>> = againsts
        .par_iter()
        .map(|sig| 
            HashSet::from_iter(sig.minhash.iter_mins())
        )
        .collect::<Vec<HashSet<&u64>>>();

    // Number of documents where hashvals appear
    // hashmap of: { hashval: n_sigs_with_hashval }
    let document_frequency: HashMap<&u64, f64> = HashMap::from(
        against_merged_mh
            .iter_mins()
            .map(|hashval| 
                (hashval, againsts_hashes
                    .par_iter()
                    .map(|hashset| 
                        f64::from(u32::from(hashset.contains(&hashval)))
                    ).sum()
                )
            ).collect::<HashMap<&u64, f64>>()
        );

    let inverse_document_frequency: HashMap<u64, f64> = HashMap::from(
        document_frequency
            .par_iter()
            .map(|(hashval, n_sigs_with_hashval)|
                (**hashval, match smooth_idf {
                    // Add 1 to not totally ignore terms that appear in all documents
                    // scikit-learn documentation (assumed to implement best practices for document classification): 
                    // > "The effect of adding “1” to the idf in the equation above is that terms with zero idf,
                    // > i.e., terms that occur in all documents in a training set, will not be entirely ignored."
                    // Source: https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
                    Some(true) => ( (1.0 + n_signatures) / (1.0 + n_sigs_with_hashval) ).ln() + 1.0,
                    Some(false) => (n_signatures / (n_sigs_with_hashval) ).ln() + 1.0,
                    _ => 1.0
                }
            )
        ).collect::<HashMap<u64, f64>>()
    );

    //     match smooth_idf {
    //     // Add 1 to not totally ignore terms that appear in all documents
    //     // scikit-learn documentation (assumed to implement best practices for document classification): 
    //     // > "The effect of adding “1” to the idf in the equation above is that terms with zero idf,
    //     // > i.e., terms that occur in all documents in a training set, will not be entirely ignored."
    //     // Source: https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
    //     Some(true) => ( (1.0 + n_signatures) / (1.0 + n_sigs_with_hashval) ).ln() + 1.0,
    //     Some(false) => (n_signatures / (n_sigs_with_hashval) ).ln() + 1.0,
    //     _ => 1.0
    // };

    return inverse_document_frequency;
}


pub fn get_hashval_inverse_document_frequency(hashval: u64, signatures: &Vec<SmallSignature>, smooth_idf: Option<bool>) -> f64 {
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
        Some(true) => ( (1.0 + n_signatures) / (1.0 + n_sigs_with_hashval) ).ln() + 1.0,
        Some(false) => (n_signatures / (n_sigs_with_hashval) ).ln() + 1.0,
        _ => 1.0
    };

    return inverse_document_frequency;
}

pub fn get_term_frequency_inverse_document_frequency(
    hashvals: &Vec<u64>, 
    query_term_frequencies: &HashMap<u64, f64>,
    inverse_document_frequency: &HashMap<u64, f64>,
) -> f64 {
    // Implementation of tf-idf for hashvals and signatures
    // https://en.wikipedia.org/wiki/Tf%E2%80%93idf
    // Square the abundances to use an L2 norm -> why?
    // Because this is the default setting in scikit-learn's battle-tested tf-idf methods:
    // https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
    // https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfVectorizer.html
    
    // let _ = print_freq_count(term_frequencies.clone(), &"tf", query_name, hashvals.len());

    // eprintln!("- inverse_document_frequencies: {} -", query.name);
    // let _ = print_freq_count(inverse_document_frequencies.clone(), &"idf", query_name, hashvals.len());
        
    // Multiply each hashval's term frequency and inverse document frequency, and sum the products
    let tf_idf: HashMap<&u64, f64> = HashMap::from(
        hashvals
        .par_iter()
        .map(
        |hashval| 
            (
                hashval, 
                query_term_frequencies[hashval] * inverse_document_frequency[hashval]
            )
        ).collect::<HashMap<&u64, f64>>()
    );

    // let _ = print_freq_count(tf_idf.clone(), &"tf-idf", query_name, tf_idf.len());

    let tf_idf_score: f64 = tf_idf.values().sum();

    // eprintln!("\tquery: {}\tn_hashes: {}\ttf_idf_score: {}", query_name, hashvals.len(), tf_idf_score);

    return tf_idf_score;
}