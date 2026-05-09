// This is the core functionality of NEAT. Generate reads turns a mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to give the contig the proper
// read coverage. Generate reads uses this to create a list of coordinates to take slices from
// the mutated fasta file. These will either be read-length fragments or fragment model length
// fragments.
use crate::{
    common::models::fragment_length::FragmentLengthModel, 
    gen_reads::errors::GenerateReadsError
};
use log::*;
use std::collections::VecDeque;
use common::models::gc_bias_model::GcBiasModel;
use common::rng::NeatRng;
use common::structs::sequence_block::SequenceBlock;
use common::structs::nucleotides::Nucleotide;

pub fn generate_fragments(
    sequence_length: usize,
    read_length: usize,
    max_del_len: usize,
    start: usize,
    coverage: usize,
    paired_ended: bool,
    fragment_model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Takes:
    // sequence_length: The length of the sequence to generate reads for.
    // read_length: the length of the reads for this run
    // coverage: the average depth of coverage for this run
    // paired_ended: when false, read_length is used as the fragment length so that
    //   coverage is computed in read-length steps rather than insert-size steps.
    // rng: the random number generator for the run
    // Returns:
    // A vector of (start, end) pairs representing fragment coordinates.
    let mut fragment_pool: Vec<usize> = Vec::new();
    // First thing we'll do is throw out regions where we can't get a full read.
    // Adding the 2 * max_del_len should ensure we are sufficiently padded.
    if sequence_length <= (read_length + (max_del_len * 2)) {
        debug!("Sequence length was too short, maybe because of a small bed region.");
        return Ok(Vec::new())
    }

    let num_frags = sequence_length.saturating_mul(coverage) / read_length;
    if paired_ended {
        // For paired-end reads the physical fragment length (insert size) determines spacing.
        for _ in 0..num_frags {
            let frag = fragment_model.generate_fragment(rng.random()?)?;
            if frag >= read_length {
                fragment_pool.push(frag);
            }
        }
    } else {
        // For single-end reads the fragment is exactly one read long; the fragment model
        // is irrelevant for spacing and would cause under-coverage if its mean >> read_length.
        fragment_pool = vec![read_length; num_frags];
    }
    if fragment_pool.is_empty() {
        debug!("Fragment pool is empty.");
        return Ok(Vec::new())
    }

    // Generate a vector of read positions
    debug!("Generating read coordinates.");
    let fragments: Vec<(usize, usize)> = cover_dataset(
        sequence_length,
        read_length,
        start,
        fragment_pool,
        coverage,
        rng
    )?;

    if fragments.is_empty() {
        debug!("No fragments generated!");
        return Ok(Vec::new())
    } else {
        Ok(fragments)
    }
}

/// Builds an unnormalized prefix-sum CDF over per-position GC bias weights for a region.
///
/// Each entry `prefix_sum[i+1] = prefix_sum[i] + weight(window starting at region_start + i)`.
/// The model's `window_size` is clamped to `region_len` so short regions are handled gracefully.
/// Returns `None` when the region is empty; the caller maps that to zero fragments.
///
/// Uses a sliding window to maintain a running GC and N count, avoiding per-position
/// Vec allocations that would dominate cost on large chromosomes.
fn build_gc_weight_prefix_sum(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    gc_bias_model: &GcBiasModel,
) -> Result<Option<Vec<f64>>, GenerateReadsError> {
    if region_end <= region_start {
        return Ok(None);
    }
    let region_len = region_end - region_start;
    let window = gc_bias_model.window_size().min(region_len);
    let n_positions = region_len - window + 1;

    let seq = &sequence_block.sequence[region_start..region_end];

    // Initialise counts for the first window.
    let mut gc_count: usize = seq[..window]
        .iter()
        .filter(|&&n| n == Nucleotide::G || n == Nucleotide::C)
        .count();
    let mut n_count: usize = seq[..window]
        .iter()
        .filter(|&&n| n == Nucleotide::N)
        .count();

    let weight_at = |gc: usize, n: usize| -> f64 {
        let called = window - n;
        if called == 0 {
            1.0  // all-N window: neutral, never filtered
        } else {
            gc_bias_model.weight_for_gc_fraction(gc as f64 / called as f64)
        }
    };

    let mut prefix_sum = vec![0.0f64; n_positions + 1];
    prefix_sum[1] = weight_at(gc_count, n_count);

    // Slide the window one base at a time, updating counts incrementally.
    for i in 1..n_positions {
        let outgoing = seq[i - 1];
        let incoming = seq[i + window - 1];

        match outgoing {
            Nucleotide::G | Nucleotide::C => gc_count -= 1,
            Nucleotide::N => n_count -= 1,
            _ => {}
        }
        match incoming {
            Nucleotide::G | Nucleotide::C => gc_count += 1,
            Nucleotide::N => n_count += 1,
            _ => {}
        }

        prefix_sum[i + 1] = prefix_sum[i] + weight_at(gc_count, n_count);
    }

    Ok(Some(prefix_sum))
}

/// Samples a 0-based position offset from the prefix-sum CDF using a single uniform draw.
fn sample_from_prefix_sum(
    prefix_sum: &[f64],
    total_weight: f64,
    rng: &mut NeatRng,
) -> Result<usize, GenerateReadsError> {
    let v = rng.random()? * total_weight;
    // Find the last index k where prefix_sum[k] <= v, i.e. the interval [k, k+1) containing v.
    let k = prefix_sum.partition_point(|&ps| ps <= v).saturating_sub(1);
    Ok(k)
}

/// Generates fragments for a region by sampling start positions directly from a per-position
/// GC bias weight CDF, with no rejection step.
///
/// Start positions are drawn with probability proportional to the model weight of the window
/// at each position, so high-weight positions are naturally over-sampled and low-weight
/// positions are under-sampled, exactly matching the sequencer's bias.
///
/// When `normalize` is true the fragment count is inflated by `1/mean_weight` (capped at
/// `100×`) so that the expected total coverage still approaches the target even in low-weight
/// regions. When `normalize` is false the fragment count equals what a uniform sequencer
/// would generate, and coverage will be proportionally lower in low-weight regions.
pub fn generate_weighted_fragments(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    read_length: usize,
    max_del_len: usize,
    coverage: usize,
    gc_bias_model: &GcBiasModel,
    fragment_model: &FragmentLengthModel,
    normalize: bool,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    const MAX_COVERAGE_MULTIPLIER: usize = 100;

    let region_len = region_end.saturating_sub(region_start);
    if region_len <= read_length + max_del_len * 2 {
        debug!("Region too short for weighted fragment generation.");
        return Ok(Vec::new());
    }

    let prefix_sum = match build_gc_weight_prefix_sum(sequence_block, region_start, region_end, gc_bias_model)? {
        Some(ps) => ps,
        None => return Ok(Vec::new()),
    };

    let total_weight = *prefix_sum.last().unwrap();
    if total_weight == 0.0 {
        debug!("All positions in region have zero GC bias weight; no fragments generated.");
        return Ok(Vec::new());
    }

    let n_valid_positions = prefix_sum.len() - 1;
    let mean_weight = total_weight / n_valid_positions as f64;

    let num_frags_base = region_len.saturating_mul(coverage) / read_length;
    let num_frags = if normalize {
        let inflated = (num_frags_base as f64 / mean_weight).round() as usize;
        let cap = num_frags_base.saturating_mul(MAX_COVERAGE_MULTIPLIER);
        if inflated > cap {
            warn!(
                "GC bias coverage inflation capped at {}x (mean weight {:.4}); \
                 target coverage may not be reached in this region",
                MAX_COVERAGE_MULTIPLIER, mean_weight
            );
        }
        inflated.min(cap)
    } else {
        num_frags_base
    };

    if num_frags == 0 {
        return Ok(Vec::new());
    }

    let mut fragments = Vec::with_capacity(num_frags);
    for _ in 0..num_frags {
        let offset = sample_from_prefix_sum(&prefix_sum, total_weight, rng)?;
        let start = region_start + offset;
        let frag_len = fragment_model.generate_fragment(rng.random()?)?;
        if frag_len < read_length {
            continue;
        }
        let end = (start + frag_len).min(region_end);
        if end - start < read_length {
            continue;
        }
        fragments.push((start, end));
    }

    fragments.sort_by_key(|&(s, _)| s);
    Ok(fragments)
}

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    read_start: usize,
    mut fragment_pool: Vec<usize>,
    coverage: usize,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Place exactly (span_length * coverage / read_length) fragments so that mean depth
    // equals the target regardless of inter-fragment jitter size.  Stopping on fragment
    // count rather than sweep count means the wildcard gap can be proportional to
    // read_length without reducing mean coverage.
    let target_count = span_length.saturating_mul(coverage) / read_length;

    rng.shuffle_in_place(&mut fragment_pool)?;
    let mut cover_fragment_pool = VecDeque::from(fragment_pool);

    let mut fragment_set: Vec<(usize, usize)> = Vec::with_capacity(target_count);
    let mut start = (rng.rand_int()? as usize) % (read_length / 4).max(1);

    while fragment_set.len() < target_count {
        let fragment_length = cover_fragment_pool.pop_front().unwrap();
        cover_fragment_pool.push_back(fragment_length);

        let temp_end = start + fragment_length;
        if temp_end > span_length {
            // Restart from the beginning of the contig for the next sweep.
            start = 0;
            continue;
        }

        fragment_set.push((start + read_start, temp_end + read_start));
        let wildcard = (rng.rand_u32()? % (read_length as u32 / 4).max(1)) as usize;
        start = temp_end + wildcard;
        if start >= span_length {
            start = 0;
        }
    }

    fragment_set.sort_by_key(|&(s, _)| s);
    Ok(fragment_set)
}


#[cfg(test)]
mod tests {
    use super::*;
    use common::rng::NeatRng;
    use common::models::gc_bias_model::GcBiasModel;
    use common::structs::sequence_block::{RegionType, SequenceBlock, SequenceMap};
    use common::structs::nucleotides::Nucleotide::{A, C, G, T};

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec!["gc-bias-test".to_string()]).unwrap()
    }

    fn make_sequence_block(sequence: Vec<common::structs::nucleotides::Nucleotide>) -> SequenceBlock {
        let len = sequence.len();
        SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, len)],
        }
    }

    fn weights_with_value(value: f64) -> Vec<f64> {
        vec![value; 101]
    }

    // ── generate_fragments (uniform path) ────────────────────────────────────

    #[test]
    fn test_cover_dataset() {
        let span_length = 100;
        let read_length = 10;
        let fragment_pool = vec![10];
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let cover = cover_dataset(
            span_length,
            read_length,
            0,
            fragment_pool,
            coverage,
            &mut rng
        ).unwrap();
        assert_eq!(cover[0], (0, 10))
    }

    #[test]
    fn test_gap_function() {
        let span_length = 100_000;
        let read_length = 100;
        let fragment_pool = vec![300];
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let cover = cover_dataset(
            span_length,
            read_length,
            0,
            fragment_pool,
            coverage,
            &mut rng
        ).unwrap();
        // target_count = span / read_length * coverage = 1000 fragments.
        assert_eq!(cover.len(), span_length / read_length * coverage);
        // After a wrap the cursor resets to 0, so the sorted first fragment starts there.
        assert_eq!(cover[0], (0, 300));
        // All fragments must fit within the span and have the pool fragment length.
        assert!(cover.iter().all(|&(s, e)| e <= span_length && e - s == 300));
    }

    #[test]
    fn test_generate_reads_single() {
        let read_length = 10;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_fragments(
            2000,
            read_length,
            0,
            0,
            coverage,
            false,
            &fragment_model,
            &mut rng,
        ).unwrap();
        // With paired_ended=false every fragment is exactly read_length wide.
        assert!(!reads.is_empty());
        assert!(reads.iter().all(|(s, e)| e - s == read_length));
    }

    #[test]
    fn test_seed_rng() {
        let sequence = vec![
            0, 0, 2, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 2, 2, 2, 4, 4, 4, 4
        ];
        let read_length = 10;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let run1 = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            false,
            &fragment_model,
            &mut rng,
        ).unwrap();

        // Reset to the same seed for reproducibility check.
        let mut rng2 = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let run2 = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            false,
            &fragment_model,
            &mut rng2,
        ).unwrap();

        assert_eq!(run1, run2)
    }

    #[test]
    fn test_generate_reads_paired() {
        let sequence: Vec<u8> = std::iter::repeat(0_u8).take(100_000).collect();
        let read_length = 100;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            true,
            &fragment_model,
            &mut rng,
        ).unwrap();
        assert!(!reads.is_empty())
    }

    // ── generate_weighted_fragments ──────────────────────────────────────────

    #[test]
    fn test_weighted_fragments_short_region_returns_empty() {
        // region_len=50, read_length=30, max_del_len=10 → guard: 50 <= 30 + 20
        let sequence_block = make_sequence_block(
            std::iter::repeat([A, C, G, T]).take(25).flatten().collect()
        );
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let result = generate_weighted_fragments(
            &sequence_block, 0, 50, 30, 10, 5,
            &GcBiasModel::default(), &fragment_model, false, &mut rng,
        ).unwrap();

        assert!(result.is_empty(), "Region too short should return no fragments");
    }

    #[test]
    fn test_weighted_fragments_zero_weight_region_returns_empty() {
        // All-A sequence (0% GC), model has weight only for 50% GC → all positions weight 0
        let sequence_block = make_sequence_block(vec![A; 500]);
        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;
        let model = GcBiasModel::from_weights(weights, 10).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let result = generate_weighted_fragments(
            &sequence_block, 0, 500, 10, 0, 5, &model, &fragment_model, false, &mut rng,
        ).unwrap();

        assert!(result.is_empty(), "All-zero-weight region should produce no fragments");
    }

    #[test]
    fn test_weighted_fragments_produces_fragments_for_valid_region() {
        let sequence_block = make_sequence_block(
            std::iter::repeat([A, C, G, T]).take(500).flatten().collect()
        );
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let fragments = generate_weighted_fragments(
            &sequence_block, 0, 2000, 100, 0, 10,
            &GcBiasModel::default(), &fragment_model, false, &mut rng,
        ).unwrap();

        assert!(!fragments.is_empty(), "Valid region with uniform model should produce fragments");
    }

    #[test]
    fn test_weighted_fragments_all_within_bounds() {
        let sequence_block = make_sequence_block(
            std::iter::repeat([A, C, G, T]).take(1000).flatten().collect()
        );
        let (region_start, region_end) = (200, 3800);
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let fragments = generate_weighted_fragments(
            &sequence_block, region_start, region_end, 100, 0, 5,
            &GcBiasModel::default(), &fragment_model, false, &mut rng,
        ).unwrap();

        for (start, end) in &fragments {
            assert!(*start >= region_start, "start {} < region_start {}", start, region_start);
            assert!(*end <= region_end, "end {} > region_end {}", end, region_end);
        }
    }

    #[test]
    fn test_weighted_fragments_deterministic() {
        let sequence_block = make_sequence_block(
            std::iter::repeat([A, C, G, T]).take(500).flatten().collect()
        );
        let fragment_model = FragmentLengthModel::default().unwrap();

        let mut rng1 = make_rng();
        let run1 = generate_weighted_fragments(
            &sequence_block, 0, 2000, 100, 0, 5,
            &GcBiasModel::default(), &fragment_model, false, &mut rng1,
        ).unwrap();

        let mut rng2 = make_rng();
        let run2 = generate_weighted_fragments(
            &sequence_block, 0, 2000, 100, 0, 5,
            &GcBiasModel::default(), &fragment_model, false, &mut rng2,
        ).unwrap();

        assert_eq!(run1, run2, "Same seed must produce identical fragment sets");
    }

    #[test]
    fn test_weighted_fragments_normalize_true_generates_more_than_false() {
        // 50% GC sequence; model mean weight ≈ 0.5 → normalize should roughly double count.
        let sequence_block = make_sequence_block(
            std::iter::repeat([A, C, G, T]).take(2500).flatten().collect()
        );
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;
        let model = GcBiasModel::from_weights(weights, 100).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();

        let mut rng = make_rng();
        let unnormalized = generate_weighted_fragments(
            &sequence_block, 0, 10000, 100, 0, 10, &model, &fragment_model, false, &mut rng,
        ).unwrap();

        let mut rng = make_rng();
        let normalized = generate_weighted_fragments(
            &sequence_block, 0, 10000, 100, 0, 10, &model, &fragment_model, true, &mut rng,
        ).unwrap();

        assert!(
            normalized.len() > unnormalized.len(),
            "normalize=true ({}) should produce more fragments than normalize=false ({})",
            normalized.len(), unnormalized.len()
        );
    }

    #[test]
    fn test_weighted_fragments_near_zero_mean_weight_is_capped() {
        // 100% GC sequence; model weights[100]=1e-6 → mean_weight ≈ 1e-6
        // → normalize would inflate by ~1e6, capped at 100x.
        let sequence_block = make_sequence_block(
            std::iter::repeat([C, G]).take(5000).flatten().collect()
        );
        let mut weights = weights_with_value(0.0);
        weights[100] = 1e-6;
        weights[50] = 1.0; // ensures model is valid (at least one positive weight)
        let model = GcBiasModel::from_weights(weights, 100).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let read_len = 100;
        let target_coverage = 10;

        let mut rng = make_rng();
        let result = generate_weighted_fragments(
            &sequence_block, 0, sequence_block.sequence.len(),
            read_len, 0, target_coverage, &model, &fragment_model, true, &mut rng,
        ).unwrap();

        let base_frags = (sequence_block.sequence.len() / read_len) * target_coverage;
        let cap = base_frags * 100;
        assert!(
            result.len() <= cap,
            "Fragment count {} should not exceed cap {}",
            result.len(), cap
        );
    }
}
