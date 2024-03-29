use crate::utils;
use crate::models;

use log::debug;
use rand::distributions::WeightedIndex;
use utils::neat_rng::NeatRng;
use models::nucleotides::Nuc;

#[derive(Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum SequencingErrorType {
    SNPError,
    IndelError,
}

pub struct SequencingError {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it.
    //
    // chromosome: uniquely identifies the contig from the fasta file.
    // position: the position within the reference (0-indexed) where the first ref base is
    //      located.
    // reference: the vector of nucelotides on chromosome starting at position which are altered
    //      for this variant
    // alternate: the vector of nucleotides that take the place of reference.
    chromosome: String,
    position: usize,
    error_type: SequencingErrorType,
    reference: Vec<Nuc>,
    alternate: Vec<Nuc>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use lib::nucleotides::Nuc::*;

    #[test]
    fn test_get_len() {
        let my_variant = SequencingError {
            chromosome: "chr1".to_string(),
            position: 10,
            error_type: SequencingErrorType::IndelError,
            reference: vec![A, A, A, G, G, T],
            alternate: vec![A],
        };

        assert_eq!(my_variant.get_length(), 5)
    }
}