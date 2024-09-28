use sourmash::sketch::minhash::KmerMinHash;
use crate::Error;

pub fn get_prob_overlap(
    query_mh: &KmerMinHash,
    database_mh: &KmerMinHash,
) -> Result<f64, Error> {
    let query_intersection = query_mh.intersection(database_mh);

    return 0.0
}