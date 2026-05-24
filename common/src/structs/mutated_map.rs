use crate::rng::NeatRng;
use crate::structs::{
    nucleotides::Nucleotide,
    variants::{Genotype, Variant, VariantError},
};
use log::debug;
use std::collections::HashMap;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum MutatedMapError {
    #[error("Error while mutating position {0}")]
    MutatePositionError(usize),
    #[error("Error retrieving location from variant")]
    VariantLocationError(#[from] VariantError),
}

#[derive(Debug, Clone)]
pub struct MutatedMap {
    pub flagged_positions: Vec<usize>,
    pub block_interval: (usize, usize),
    pub variant_map: HashMap<usize, Variant>,
}

impl MutatedMap {
    pub fn from_interval(
        start: usize,
        end: usize,
        variant_vec: Vec<Variant>,
    ) -> Result<Self, MutatedMapError> {
        let mut variant_map: HashMap<usize, Variant> = HashMap::new();
        let mut flagged_positions: Vec<usize> = Vec::new();
        for variant in variant_vec {
            let location = variant.get_loc()?;
            if variant_map.contains_key(&location) {
                debug!(
                    "Two mutations sampled at position {}; keeping first, discarding second",
                    location
                );
                continue;
            }
            variant_map.insert(location, variant);
            flagged_positions.push(location);
        }
        Ok(MutatedMap {
            flagged_positions,
            block_interval: (start, end),
            variant_map,
        })
    }

    pub fn is_flagged(&self, position: &usize) -> bool {
        self.flagged_positions.contains(position)
    }

    pub fn contains(&self, (x, y): (usize, usize)) -> bool {
        x >= self.block_interval.0 && y <= self.block_interval.1
    }

    pub fn mutate_position(
        &self,
        position: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<Nucleotide>, MutatedMapError> {
        match self.variant_map[&position].genotype {
            Genotype::Homozygous => Ok(self.variant_map[&position].alternate.get_vec().unwrap().clone()),
            Genotype::Heterozygous => {
                if rng.gen_bool(0.5).unwrap() {
                    Ok(self.variant_map[&position].reference.clone())
                } else {
                    Ok(self.variant_map[&position].alternate.get_vec().unwrap().clone())
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

    fn make_map(start: usize, end: usize, variant: Variant) -> MutatedMap {
        MutatedMap::from_interval(start, end, vec![variant]).unwrap()
    }

    #[test]
    fn test_is_flagged() {
        let variant =
            Variant::new(VariantType::SNP, 1003, &vec![A], &vec![G], &mut vec![1, 0]).unwrap();
        let map = make_map(1000, 2000, variant);
        assert!(map.is_flagged(&1003));
        assert!(!map.is_flagged(&1007));
    }

    #[test]
    fn test_contains() {
        let variant =
            Variant::new(VariantType::SNP, 1003, &vec![A], &vec![G], &mut vec![1, 0]).unwrap();
        let map = make_map(1000, 2000, variant);
        assert!(map.contains((1003, 1300)));
        assert!(!map.contains((1900, 2100)));
    }

    #[test]
    fn test_mutate_position_homozygous_always_returns_alt() {
        let variant =
            Variant::new(VariantType::SNP, 1003, &vec![A], &vec![G], &mut vec![1, 1]).unwrap();
        let map = make_map(1000, 2000, variant);
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        for _ in 0..10 {
            assert_eq!(map.mutate_position(1003, &mut rng).unwrap(), vec![G]);
        }
    }

    #[test]
    fn test_mutate_position_heterozygous_returns_ref_or_alt() {
        let variant =
            Variant::new(VariantType::SNP, 1003, &vec![A], &vec![G], &mut vec![1, 0]).unwrap();
        let map = make_map(1000, 2000, variant);
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        let mut saw_ref = false;
        let mut saw_alt = false;
        for _ in 0..20 {
            let result = map.mutate_position(1003, &mut rng).unwrap();
            assert!(result == vec![A] || result == vec![G]);
            if result == vec![A] {
                saw_ref = true;
            }
            if result == vec![G] {
                saw_alt = true;
            }
        }
        assert!(saw_ref, "expected at least one ref result in 20 trials");
        assert!(saw_alt, "expected at least one alt result in 20 trials");
    }
}
