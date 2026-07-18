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
    // Symbolic / structural variants (CNV, DEL, DUP, INS, INV, BND, ...) live here
    // separately from variant_map: they shouldn't flag positions for per-base
    // mutation during read generation (their ALT is a tag like `<DEL>`, not a
    // literal base sequence), but they still need to round-trip to the output VCF.
    pub sv_records: Vec<Variant>,
}

/// Per-variant allelic-depth accumulator populated by the gen-reads fragment
/// loop. The key is the variant's absolute per-contig reference position
/// (matching [`Variant::location`]); the value is `(ref_count, alt_count)` —
/// how many simulated reads carried each allele at that locus across the
/// current contig. Drives `FORMAT/AD`, `FORMAT/DP`, and `FORMAT/AF` in the
/// golden VCF (see issue #176).
///
/// Only literal SNP / insertion / deletion variants are tracked here.
/// Symbolic / structural variants use span-based, not point-based, depth
/// semantics — they emit `.` placeholders on the AD/DP/AF fields.
pub type AdCounter = HashMap<usize, (u32, u32)>;

impl MutatedMap {
    pub fn from_interval(
        start: usize,
        end: usize,
        variant_vec: Vec<Variant>,
    ) -> Result<Self, MutatedMapError> {
        let mut variant_map: HashMap<usize, Variant> = HashMap::new();
        let mut flagged_positions: Vec<usize> = Vec::new();
        let mut sv_records: Vec<Variant> = Vec::new();
        for variant in variant_vec {
            if variant.alternate.is_symbolic() {
                sv_records.push(variant);
                continue;
            }
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
        // Keep sorted so `is_flagged` and the read-window overlap query in
        // `write_block_fastq` can binary-search instead of scanning every
        // variant. With dense models (e.g. the corpus-aggregated COSMIC tumor
        // rate at ~5.5e-3 → ~220k SNVs on chr22) the old linear scans made
        // read writing O(fragments × variants); sorting once here makes the
        // per-read lookup O(log V + hits).
        flagged_positions.sort_unstable();
        Ok(MutatedMap {
            flagged_positions,
            block_interval: (start, end),
            variant_map,
            sv_records,
        })
    }

    pub fn is_flagged(&self, position: &usize) -> bool {
        // flagged_positions is kept sorted by from_interval, so binary-search
        // rather than a linear scan (called per-base in the chimeric subseq path).
        self.flagged_positions.binary_search(position).is_ok()
    }

    pub fn contains(&self, (x, y): (usize, usize)) -> bool {
        x >= self.block_interval.0 && y <= self.block_interval.1
    }

    pub fn mutate_position(
        &self,
        position: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<Nucleotide>, MutatedMapError> {
        let variant = &self.variant_map[&position];
        // Defensive: symbolic / structural ALTs shouldn't be in variant_map
        // (from_interval routes them to sv_records), so mutate_position should
        // never be called for one. If a future change leaks one through, fall
        // back to emitting the reference instead of panicking on as_literal().
        let alt_bases = match variant.alternate.as_literal() {
            Some(bases) => bases.to_vec(),
            None => return Ok(variant.reference.clone()),
        };
        // An explicit allele_fraction (input VCF, #398) emits alt with probability f.
        // Otherwise the Genotype default: homozygous always alt (no rng draw, so the
        // default path stays byte-identical), het a 0.5 coin flip.
        match variant.allele_fraction {
            Some(f) => {
                if rng.gen_bool(f).unwrap() {
                    Ok(alt_bases)
                } else {
                    Ok(variant.reference.clone())
                }
            }
            None => match variant.genotype {
                Genotype::Homozygous => Ok(alt_bases),
                Genotype::Heterozygous => {
                    if rng.gen_bool(0.5).unwrap() {
                        Ok(variant.reference.clone())
                    } else {
                        Ok(alt_bases)
                    }
                }
            },
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

    #[test]
    fn test_from_interval_routes_symbolic_to_sv_records() {
        use crate::structs::variants::{AlternateType, Provenance, SvData, SvType};
        let snp =
            Variant::new(VariantType::SNP, 50, &vec![A], &vec![G], &mut vec![1, 1]).unwrap();
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 100,
            reference: vec![A],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        let map = MutatedMap::from_interval(0, 200, vec![snp, sv]).unwrap();
        // SNP flags position 50; symbolic SV does NOT flag position 100.
        assert!(map.is_flagged(&50));
        assert!(!map.is_flagged(&100));
        assert_eq!(map.variant_map.len(), 1);
        assert_eq!(map.sv_records.len(), 1);
        assert_eq!(map.sv_records[0].location, 100);
    }

    #[test]
    fn test_mutate_position_returns_reference_for_symbolic_leak() {
        // Defensive: from_interval routes symbolic ALTs to sv_records, but if
        // a future change leaks one into variant_map, mutate_position must not
        // panic on as_literal().unwrap() — it should fall back to reference.
        use crate::structs::variants::{AlternateType, Provenance, SvData, SvType};
        use std::collections::HashMap;
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 75,
            reference: vec![A],
            alternate: AlternateType::Symbolic(SvData::new("<DUP>", SvType::Dup)),
            genotype_str: "1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        let mut variant_map = HashMap::new();
        variant_map.insert(75usize, sv);
        let map = MutatedMap {
            flagged_positions: vec![75],
            block_interval: (0, 200),
            variant_map,
            sv_records: Vec::new(),
        };
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        assert_eq!(map.mutate_position(75, &mut rng).unwrap(), vec![A]);
    }
}
