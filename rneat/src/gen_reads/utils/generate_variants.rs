use crate::{
    common::{
        models::mutation_model::MutationModel,
        structs::{sequence_block::SequenceBlock, variants::Variant},
    },
    gen_reads::errors::GenerateReadsError,
};
use common::rng::NeatRng;
use log::debug;

/// Generates random variants for a contig.
///
/// `rate_segments` is a sorted, non-overlapping list of `(start, end, rate)` triples
/// covering the mutable positions of the contig. N-regions and excluded positions
/// (e.g. positions already occupied by input variants) must be absent from the list —
/// the caller is responsible for building a compact segment representation.
pub fn generate_variants(
    sequence_block: &SequenceBlock,
    rate_segments: &[(usize, usize, f64)],
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsError> {
    let mut block_variants: Vec<Variant> = Vec::new();

    if rate_segments.is_empty() {
        return Ok(Some(block_variants));
    }

    // Cumulative weights for weighted segment selection via bisect.
    let mut cumulative: Vec<f64> = Vec::with_capacity(rate_segments.len());
    let mut acc = 0.0_f64;
    for &(start, end, rate) in rate_segments {
        acc += (end - start) as f64 * rate;
        cumulative.push(acc);
    }
    let total_weight = acc;

    for _ in 0..num_mutations {
        let target = rng.random()? * total_weight;
        let seg_idx = cumulative
            .partition_point(|&c| c <= target)
            .min(rate_segments.len() - 1);
        let (seg_start, seg_end, seg_rate) = rate_segments[seg_idx];
        let weight_before = if seg_idx > 0 {
            cumulative[seg_idx - 1]
        } else {
            0.0
        };
        let weight_in_seg = (target - weight_before).max(0.0);
        let pos_in_seg = (weight_in_seg / seg_rate) as usize;
        let location = (seg_start + pos_in_seg).min(seg_end - 1);

        debug!("location: {}", location);
        block_variants.push(mutation_model.generate_mutation(
            &sequence_block.sequence,
            location,
            ploidy,
            rng,
        )?);
    }

    Ok(Some(block_variants))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{
        models::mutation_model::MutationModel,
        structs::sequence_block::{RegionType, SequenceBlock, SequenceMap},
    };
    use common::rng::NeatRng;
    use common::structs::nucleotides::Nucleotide;

    fn make_block(len: usize) -> SequenceBlock {
        let sequence: Vec<Nucleotide> = (0..len)
            .map(|i| match i % 4 {
                0 => Nucleotide::A,
                1 => Nucleotide::C,
                2 => Nucleotide::G,
                _ => Nucleotide::T,
            })
            .collect();
        SequenceBlock {
            contig: "test".to_string(),
            ref_start: 0,
            ref_end: len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, len)],
        }
    }

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap()
    }

    #[test]
    fn test_generate_variants_count() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &segments, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn test_generate_variants_positions_in_range() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let variants = generate_variants(&block, &segments, &model, 10, 2, &mut rng)
            .unwrap()
            .unwrap();
        for v in variants {
            assert!(
                v.location < 200,
                "Variant location {} out of range",
                v.location
            );
        }
    }

    #[test]
    fn test_generate_variants_empty_segments_returns_empty() {
        let block = make_block(200);
        let segments: Vec<(usize, usize, f64)> = vec![];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &segments, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert!(
            result.is_empty(),
            "Empty segments should produce no variants"
        );
    }

    #[test]
    fn test_generate_variants_deterministic() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let result1 = generate_variants(&block, &segments, &model, 5, 2, &mut make_rng())
            .unwrap()
            .unwrap();
        let result2 = generate_variants(&block, &segments, &model, 5, 2, &mut make_rng())
            .unwrap()
            .unwrap();
        for (v1, v2) in result1.iter().zip(result2.iter()) {
            assert_eq!(v1.location, v2.location);
            assert_eq!(v1.variant_type, v2.variant_type);
        }
    }
}
