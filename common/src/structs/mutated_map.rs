//! Note sure if these are structs or models or something in between
//! The point of this is to store a map of valid regions, including
//! the ref versus alt version at the positions where mutations occur.
//! This allows us to choose when we write the fastq, whether to use
//! the ref or alt version.
use crate::structs::{
    fasta_map::{FastaMapError, decode_block_filename}, 
    nucleotides::Nucleotide, 
    variants::{Genotypes, Variant, VariantError},
};
use std::collections::HashMap;
use simple_rng::NeatRng;
use thiserror::Error;
use std::path::PathBuf;

#[derive(Debug, Error)]
pub enum MutatedMapError {
    #[error("Error while mutating position {0}")]
    MutatePositionError(usize),
    #[error("Error retrieving location from variant")]
    VariantLocationError(#[from] VariantError),
    #[error("Error from FastaMap: {0}")]
    FastaMapError(#[from] FastaMapError),
}

#[derive(Debug, Clone)]
pub struct MutatedMap {
    // This is the index of positions for the variants
    pub flagged_positions: Vec<usize>,
    // This is the filename of the sequence block we are mutating
    pub sequence_block: PathBuf,
    // start and end points of block
    pub block_interval: (usize, usize),
    // The variant map, maps positions to the variant
    pub variant_map: HashMap<usize, Variant>,
}

impl MutatedMap {
    pub fn new(
        sequence_block: PathBuf,
        variant_vec: Vec<Variant>
    ) -> Result<Self, MutatedMapError> {
        // Reconstruct the range from the filename
        let block_interval = decode_block_filename(&sequence_block)?;
        // Build the variant_map
        let mut variant_map: HashMap<usize, Variant> = HashMap::new();
        let mut flagged_positions: Vec<usize> = Vec::new();
        for variant in variant_vec {
            let location = variant.get_loc()?;
            variant_map.insert(location, variant);
            flagged_positions.push(location);
        }

        Ok(MutatedMap {
            flagged_positions,
            sequence_block,
            block_interval,
            variant_map
        })
    }

    pub fn is_flagged(&self, position: &usize) -> bool {
        // Checks if this position is flagged
        self.flagged_positions.contains(position)
    }

    pub fn contains(&self, (x, y): (usize, usize)) -> Result<bool, MutatedMapError> {
        let first_clause = x >= self.block_interval.0;
        let second_clause = y <= self.block_interval.1;
        Ok(first_clause && second_clause)
    }

    pub fn mutate_position(
        &self,
        position: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<Nucleotide>, MutatedMapError> {
        // This method searches for a variant at the given position and decides whether to return
        // the ref, the alt, or the original seq. Note that this function assumes that 
        // the calling function has already checked thath position is within variant_map, or else
        // this will throw an index out of range error.
        match self.variant_map[&position].genotype {
            Genotypes::Homozygous => return Ok(self.variant_map[&position].alternate.clone()),
            Genotypes::Heterozygous => {
                // 50/50 chance we mutate
                if rng.gen_bool(0.5).unwrap() {
                    return Ok(self.variant_map[&position].reference.clone())
                } else {
                    return Ok(self.variant_map[&position].alternate.clone())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, G};
    use crate::structs::variants::VariantType;

    #[test]
    fn test_is_flagged() {
        let sequence_block = PathBuf::from("chr1_0001000_0002000.fa");

        let variant = Variant::new(
            VariantType::SNP,
            1003,
            &vec![A],
            &vec![G],
            &mut vec![1, 0],
        ).unwrap();

        let map = MutatedMap::new(
            sequence_block,
            vec![variant],
        ).unwrap();

        assert_eq!(map.is_flagged(&1003), true);
        assert_eq!(map.is_flagged(&1007), false);
    }

    #[test]
    fn test_contains() {
        let sequence_block = PathBuf::from("chr1_0001000_0002000.fa");

        let variant = Variant::new(
            VariantType::SNP,
            1003,
            &vec![A],
            &vec![G],
            &mut vec![1, 0],
        ).unwrap();

        let map = MutatedMap::new(
            sequence_block,
            vec![variant],
        ).unwrap();

        assert_eq!(map.contains((1003, 1300)).unwrap(), true);
        assert_eq!(map.contains((1900, 2100)).unwrap(), false);
    }

    #[test]
    fn test_mutate_position_homozygous_always_returns_alt() {
        let sequence_block = PathBuf::from("chr1_0001000_0002000.fa");
        let variant = Variant::new(
            VariantType::SNP,
            1003,
            &vec![A],
            &vec![G],
            &mut vec![1, 1], // homozygous
        ).unwrap();
        let map = MutatedMap::new(sequence_block, vec![variant]).unwrap();
        let mut rng = simple_rng::NeatRng::new_from_seed(&vec![
            "test".to_string()
        ]).unwrap();
        // Homozygous must always return the alternate
        for _ in 0..10 {
            let result = map.mutate_position(1003, &mut rng).unwrap();
            assert_eq!(result, vec![G]);
        }
    }

    #[test]
    fn test_mutate_position_heterozygous_returns_ref_or_alt() {
        let sequence_block = PathBuf::from("chr1_0001000_0002000.fa");
        let variant = Variant::new(
            VariantType::SNP,
            1003,
            &vec![A],
            &vec![G],
            &mut vec![1, 0], // heterozygous
        ).unwrap();
        let map = MutatedMap::new(sequence_block, vec![variant]).unwrap();
        let mut rng = simple_rng::NeatRng::new_from_seed(&vec![
            "test".to_string()
        ]).unwrap();
        // Heterozygous must return either ref or alt, never anything else
        let mut saw_ref = false;
        let mut saw_alt = false;
        for _ in 0..20 {
            let result = map.mutate_position(1003, &mut rng).unwrap();
            assert!(result == vec![A] || result == vec![G]);
            if result == vec![A] { saw_ref = true; }
            if result == vec![G] { saw_alt = true; }
        }
        assert!(saw_ref, "expected at least one ref result in 20 trials");
        assert!(saw_alt, "expected at least one alt result in 20 trials");
    }
}