// Functions to compute statisical signifiance of search results

use rayon::prelude::*;

use crate::utils::multicollection::SmallSignature;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::Error;
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Display, Formatter};

pub enum Normalization {
    // L1 norm is the equivalent of frequencies/probabilities, as the counts
    // are divided by the length of the Vec<> object, or mathematically, the
    // number of items in the Vec object, assuming unit length for each
    L1,

    // L2 norm divides the counts by the sum of squares of all counts
    // L2 norm is the Euclidean distance of the counts in N-dimensional vector space
    // When track_abundance=False, L1 and L2 norms are equivalent, since all "counts"
    // are 1, and even if you square it, 1^2 = 1
    L2,
}

impl Display for Normalization {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::L1 => write!(f, "L1 Normalization"),
            Self::L2 => write!(f, "L2 Normalization"),
        }
    }
}

pub fn get_hash_frequencies(
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
        Some(Normalization::L2) => minhash_abunds
            .par_iter()
            .map(|(_hashval, abund)| abund * abund)
            .sum::<f64>(),
        // TODO: this should probably be an error
        _ => 0.0,
    };

    let frequencies: HashMap<u64, f64> = minhash_abunds
        .par_iter()
        .map(|(hashval, abund)|
            // TODO: add a match statement here to error out properly if the hashval was not found 
            // in the minhash_abunds for some reason (shouldn't happen but ... computers be crazy)
            (
                *hashval,
                abund / abund_normalization
            ))
        .collect::<HashMap<u64, f64>>();

    frequencies
}

// #[cfg(feature = "maths")]
pub fn get_prob_overlap(
    hashvals: &Vec<u64>,
    query_frequencies: &HashMap<u64, f64>,
    against_frequencies: &HashMap<u64, f64>,
) -> f64 {
    // It's not guaranteed to me that the MinHashes from the query and database are in the same order, so iterate over one of them
    // and use a hashmap to retrieve the frequency value of the other
    let prob_overlap = hashvals
        .par_iter()
        .map(|hashval| query_frequencies[hashval] * against_frequencies[hashval])
        .sum();

    prob_overlap
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
        first_sig.minhash.scaled(),
        first_sig.minhash.ksize() as u32,
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
    smooth_idf: Option<bool>,
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
        .map(|sig| HashSet::from_iter(sig.minhash.iter_mins()))
        .collect::<Vec<HashSet<&u64>>>();

    // Number of documents where hashvals appear
    // hashmap of: { hashval: n_sigs_with_hashval }
    let document_frequency: HashMap<&u64, f64> = against_merged_mh
        .iter_mins()
        .par_bridge()
        .map(|hashval| {
            (
                hashval,
                againsts_hashes
                    .par_iter()
                    .map(|hashset| f64::from(u32::from(hashset.contains(&hashval))))
                    .sum(),
            )
        })
        .collect::<HashMap<&u64, f64>>();

    let inverse_document_frequency: HashMap<u64, f64> = document_frequency
        .par_iter()
        .map(|(hashval, n_sigs_with_hashval)| {
            (
                **hashval,
                match smooth_idf {
                    // Add 1 to not totally ignore terms that appear in all documents
                    // scikit-learn documentation (assumed to implement best practices for document classification):
                    // > "The effect of adding “1” to the idf in the equation above is that terms with zero idf,
                    // > i.e., terms that occur in all documents in a training set, will not be entirely ignored."
                    // Source: https://scikit-learn.org/1.5/modules/generated/sklearn.feature_extraction.text.TfidfTransformer.html
                    Some(true) => ((1.0 + n_signatures) / (1.0 + n_sigs_with_hashval)).ln() + 1.0,
                    Some(false) => (n_signatures / (n_sigs_with_hashval)).ln() + 1.0,
                    _ => 1.0,
                },
            )
        })
        .collect::<HashMap<u64, f64>>();

    inverse_document_frequency
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

    // Multiply each hashval's term frequency and inverse document frequency, and sum the products
    let tf_idf: HashMap<&u64, f64> = hashvals
        .par_iter()
        .map(|hashval| {
            (
                hashval,
                query_term_frequencies[hashval] * inverse_document_frequency[hashval],
            )
        })
        .collect::<HashMap<&u64, f64>>();

    let tf_idf_score: f64 = tf_idf.values().sum();

    tf_idf_score
}
