use crate::{
    common::{
        models::{mutation_model::MutationModel, snp_trinuc_model::TrinucFrame},
        structs::{sequence_block::SequenceBlock, variants::Variant},
    },
    gen_reads::errors::GenerateReadsError,
};
use common::rng::NeatRng;
use log::debug;
use std::collections::HashMap;

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
    if rate_segments.is_empty() {
        return Ok(Some(Vec::new()));
    }

    // #372: weight WHERE mutations land by the local trinucleotide's mutation
    // propensity w(ctx), so context-specific signatures (e.g. APOBEC's TpC hotspots)
    // reproduce in the output. When the model's context weights are uniform (default /
    // untrained model), placement is context-neutral and we take the fast
    // segment-arithmetic path — byte-identical to the pre-#372 behavior.
    let ctx_weights = mutation_model.context_weights()?;
    if context_is_uniform(&ctx_weights) {
        return generate_variants_uniform(
            sequence_block,
            rate_segments,
            mutation_model,
            num_mutations,
            ploidy,
            rng,
        );
    }

    let mut block_variants: Vec<Variant> = Vec::new();
    let seq = &sequence_block.sequence;
    let mean_w = ctx_weights.values().sum::<f64>() / (ctx_weights.len().max(1) as f64);

    // Per-position cumulative weight = rate × w(ctx). `seg_flat_start[i]` is the flat
    // index at which rate_segments[i] begins, so a sampled flat index maps back to a
    // genomic position without storing a parallel positions vector.
    let total_positions: usize = rate_segments.iter().map(|&(s, e, _)| e.saturating_sub(s)).sum();
    let mut cum: Vec<f64> = Vec::with_capacity(total_positions);
    let mut seg_flat_start: Vec<usize> = Vec::with_capacity(rate_segments.len());
    let mut acc = 0.0_f64;
    for &(start, end, rate) in rate_segments {
        seg_flat_start.push(cum.len());
        for p in start..end {
            // Interior positions have a trinucleotide context; contig-edge positions
            // fall back to the mean weight (no ±1 flank available).
            let w = if p >= 1 && p + 1 < seq.len() {
                let ctx = TrinucFrame::from(&[seq[p - 1], seq[p], seq[p + 1]]);
                *ctx_weights.get(&ctx).unwrap_or(&mean_w)
            } else {
                mean_w
            };
            acc += rate * w;
            cum.push(acc);
        }
    }
    let total_weight = acc;
    if total_weight <= 0.0 || cum.is_empty() {
        return Ok(Some(block_variants));
    }

    for _ in 0..num_mutations {
        let target = rng.random()? * total_weight;
        let flat_idx = cum.partition_point(|&c| c <= target).min(cum.len() - 1);
        let seg_idx = seg_flat_start
            .partition_point(|&s| s <= flat_idx)
            .saturating_sub(1);
        let seg_start = rate_segments[seg_idx].0;
        let location = seg_start + (flat_idx - seg_flat_start[seg_idx]);
        debug!("location (context-weighted): {}", location);
        block_variants.push(mutation_model.generate_mutation(seq, location, ploidy, rng)?);
    }

    Ok(Some(block_variants))
}

/// True when the per-context weights carry no signal (uniform / empty) — the default,
/// untrained-model case. Placement then stays context-neutral (pre-#372 behavior).
fn context_is_uniform(weights: &HashMap<TrinucFrame, f64>) -> bool {
    let mut vals = weights.values().copied();
    let Some(first) = vals.next() else {
        return true;
    };
    let (mut lo, mut hi) = (first, first);
    for v in vals {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    (hi - lo) <= 1e-9 * hi.max(1e-12)
}

/// Context-neutral placement: positions weighted only by the regional `rate_segments`,
/// uniform within a segment. This is the original (pre-#372) algorithm, used whenever
/// the model's context weights are uniform.
fn generate_variants_uniform(
    sequence_block: &SequenceBlock,
    rate_segments: &[(usize, usize, f64)],
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsError> {
    let mut block_variants: Vec<Variant> = Vec::new();

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
    fn context_weighted_placement_concentrates_at_upweighted_context() {
        use crate::common::models::snp_trinuc_model::TrinucFrame;
        use common::structs::nucleotides::Nucleotide::{A, C, G};
        use common::structs::transition_matrix::TransitionMatrix;
        use std::collections::HashMap;

        // make_block is ACGTACGT…, so C positions (p % 4 == 1) carry the context ACG.
        // Up-weight ACG heavily → #372 placement should concentrate mutations there;
        // context-neutral placement would put only ~1/4 at those positions.
        let mut trinuc_freq: HashMap<TrinucFrame, f64> = HashMap::new();
        trinuc_freq.insert(TrinucFrame::from((A, C, G)), 1000.0);
        let model = MutationModel::from_raw_data(
            0.01,
            0.0,
            vec![1.0, 0.0, 0.0], // all SNP
            HashMap::new(),
            trinuc_freq,
            HashMap::new(),
            vec![],
            vec![],
            vec![],
            vec![],
            Some(TransitionMatrix::default().unwrap()),
        )
        .unwrap();

        let block = make_block(400);
        let segments = vec![(0usize, 400usize, 1.0f64)];
        let mut rng = make_rng();
        let variants = generate_variants(&block, &segments, &model, 200, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(variants.len(), 200);
        let at_acg = variants.iter().filter(|v| v.location % 4 == 1).count();
        let frac = at_acg as f64 / variants.len() as f64;
        assert!(
            frac > 0.8,
            "expected >80% of placements at the up-weighted ACG context, got {:.2} ({}/{})",
            frac,
            at_acg,
            variants.len()
        );
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
