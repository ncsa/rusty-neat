use log::{error, debug};
use common::rng::NeatRng;
use crate::{
    common::{
        models::mutation_model::MutationModel,
        structs::{
            sequence_block::SequenceBlock,
            variants::Variant
        }
    },
    gen_reads::errors::GenerateReadsError
};

pub fn generate_variants(
    sequence_block: &SequenceBlock,
    region_weights: &Vec<f64>,
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsError> {
    let block_len = sequence_block.get_len();
    if region_weights.len() != block_len {
        error!("Region weights did not match block length");
        return Err(GenerateReadsError::GenerateVariantsError);
    }

    let mut block_variants: Vec<Variant> = Vec::new();

    // Build compact segments of consecutive equal non-zero weights.
    // Avoids the 5-copy overhead of constructing a DiscreteDistribution over
    // the full chromosome-length weights vector.
    let mut segments: Vec<(usize, usize, f64)> = Vec::new();
    let mut i = 0;
    while i < block_len {
        let rate = region_weights[i];
        if rate > 0.0 {
            let start = i;
            while i < block_len && (region_weights[i] - rate).abs() < f64::EPSILON {
                i += 1;
            }
            segments.push((start, i, rate));
        } else {
            i += 1;
        }
    }

    if segments.is_empty() {
        return Ok(Some(block_variants));
    }

    // Cumulative weights for weighted segment selection via bisect.
    let mut cumulative: Vec<f64> = Vec::with_capacity(segments.len());
    let mut acc = 0.0_f64;
    for &(start, end, rate) in &segments {
        acc += (end - start) as f64 * rate;
        cumulative.push(acc);
    }
    let total_weight = acc;

    for _ in 0..num_mutations {
        let target = rng.random()? * total_weight;
        let seg_idx = cumulative.partition_point(|&c| c <= target)
            .min(segments.len() - 1);
        let (seg_start, seg_end, seg_rate) = segments[seg_idx];
        let weight_before = if seg_idx > 0 { cumulative[seg_idx - 1] } else { 0.0 };
        let weight_in_seg = (target - weight_before).max(0.0);
        let pos_in_seg = (weight_in_seg / seg_rate) as usize;
        let location = (seg_start + pos_in_seg).min(seg_end - 1);

        debug!("location: {}", location);
        block_variants.push(
            mutation_model.generate_mutation(
                &sequence_block.sequence,
                location,
                ploidy,
                rng,
            )?
        );
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
    use common::structs::nucleotides::Nucleotide;
    use common::rng::NeatRng;

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
        ]).unwrap()
    }

    #[test]
    fn test_generate_variants_count() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &weights, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn test_generate_variants_positions_in_range() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let variants = generate_variants(&block, &weights, &model, 10, 2, &mut rng)
            .unwrap()
            .unwrap();
        for v in variants {
            assert!(v.location < 200, "Variant location {} out of range", v.location);
        }
    }

    #[test]
    fn test_generate_variants_weight_mismatch_errors() {
        let block = make_block(200);
        let weights = vec![1.0f64; 100]; // wrong length
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        assert!(generate_variants(&block, &weights, &model, 5, 2, &mut rng).is_err());
    }

    #[test]
    fn test_generate_variants_deterministic() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let result1 = generate_variants(&block, &weights, &model, 5, 2, &mut make_rng())
            .unwrap().unwrap();
        let result2 = generate_variants(&block, &weights, &model, 5, 2, &mut make_rng())
            .unwrap().unwrap();
        for (v1, v2) in result1.iter().zip(result2.iter()) {
            assert_eq!(v1.location, v2.location);
            assert_eq!(v1.variant_type, v2.variant_type);
        }
    }
}