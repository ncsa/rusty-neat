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
use common::structs::fasta_map::{FastaMapError, SequenceBlock};

pub fn generate_fragments(
    sequence_length: usize,
    read_length: usize,
    max_del_len: usize,
    start: usize,
    coverage: usize,
    fragment_model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Takes:
    // sequence_length: The length of the sequence to generate reads for.
    // read_length: the length ef the reads for this run
    // coverage: the average depth of coverage for this run
    // rng: the random number generator for the run
    // Returns:
    // HashSet of vectors representing the read sequences, stored on the heap in box.
    //
    // This takes a mutated sequence and produces a set of reads based on the mutated sequence. For
    // paired ended reads, this will generate a set of reads from each end, by taking the reverse
    // complement in the output
    //
    let mut fragment_pool: Vec<usize> = Vec::new();
    // First thing we'll do is throw out regions where we can't get a full read.
    // Adding the 2 * max_del_len should ensure we are sufficiently padded.
    if sequence_length <= (read_length + (max_del_len * 2)) {
        debug!("Sequence length was too short, maybe because of a small bed region.");
        return Ok(Vec::new())
    }
    
    let num_frags = (sequence_length / read_length) * coverage;
    // add fragments to the fragment pool
    for _ in 0..num_frags {
        let frag = fragment_model.generate_fragment(rng.random()?)?;
        // filter fragments
        if frag < read_length {
            continue
        } else {
            fragment_pool.push(frag);
        }
    }
    if fragment_pool.is_empty() {
        // Not sure if this should be a recoverable error.
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

/// Probabilistically filter fragments using acceptance-rejection sampling against the GC bias
/// model. Each fragment is accepted with probability `weight(fragment) / max_weight`, so the
/// highest-weighted GC bin is always accepted (probability 1.0) and all others are retained at
/// a rate proportional to their relative weight. This preserves the shape of the model without
/// requiring an absolute probability scale.
///
/// Call `estimate_region_mean_gc_weight` first and inflate the fragment count accordingly so
/// that the surviving fragments still hit the target coverage depth.
pub fn apply_gc_bias_to_fragments(
    fragments: Vec<(usize, usize)>,
    sequence_block: &SequenceBlock,
    gc_bias_model: &GcBiasModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    let max_weight = gc_bias_model.max_weight();
    // from_weights guarantees at least one positive weight, so max_weight > 0 always holds
    debug_assert!(max_weight > 0.0, "GcBiasModel guarantees at least one positive weight");
    let mut retained = Vec::new();

    for (start, end) in fragments {
        let seq = sequence_block.get_subseq(start, end)?;
        let weight = gc_bias_model.weight_for_sequence(&seq);
        let accept_prob = weight / max_weight;

        if rng.random()? < accept_prob {
            retained.push((start, end));
        }
    }

    Ok(retained)
}

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    read_start: usize,
    mut fragment_pool: Vec<usize>,
    coverage: usize,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Takes:
    // span_length: Total number of bases in the sequence
    // read_length: The length of the reads for this run
    // fragment_pool: a vector of sizes for the fragments. If empty, it will instead be filled
    // by the read_length (single ended reads)
    // paired_ended: true or false if the run is paired ended mode or not.
    // coverage: The coverage depth for the reads
    // Returns:
    // A vector of tuples (usize, usize), denoting the start and end positions of the fragment of
    // DNA that was sequenced.
    //
    // This function selects the positions of the reads. It starts at the beginning and goes out
    // one read length, then picks a random jump between 0 and half the read length to move
    // And picks those coordinates for a second read. Once the tail of the read is past the end,
    // we start over again at 0.
    //
    // Reads that will be start and end of the fragment.
    let mut fragment_set: Vec<(usize, usize)> = vec![];

    let mut cover_fragment_pool: VecDeque<usize>;
    // shuffle the fragment pool
    rng.shuffle_in_place(&mut fragment_pool)?;
    cover_fragment_pool = VecDeque::from(fragment_pool);
    // Gap size to keep track of how many uncovered bases we have per layer, to help decide if we
    // need more layers
    let mut layer_count: usize = 0;
    // start this party off somewhere random.
    let mut start = ((rng.rand_int()?) as usize) % (read_length/4);
    // create coverage number of layers
    // TODO - make sure this doesn't get caught in an infinite loop.
    while layer_count < coverage {
        let fragment_length = cover_fragment_pool.pop_front().unwrap();
        let temp_end = start + fragment_length;
        cover_fragment_pool.push_back(fragment_length);
        if temp_end > span_length {
            // TODO some variation on this modulo idea will work for bacterial reads
            start = temp_end % span_length;
            layer_count += 1;
            continue
        }
        // Read start ensures our fragments are within the sequence block
        fragment_set.push((start+read_start, temp_end+read_start));
        // Picks a number between zero and a quarter of a read length
        let wildcard: usize = (rng.rand_u32().unwrap() % 10) as usize;
        // adds to the start to give it some spice
        start += temp_end + wildcard;
        // sanity check. If we are out of bounds, take the modulo
        if start >= span_length {
            // get us back in bounds
            start = start % span_length;
            layer_count += 1;
        }
    }
    fragment_set.sort_by(
        |a, b| a.0.cmp(&b.0)
    );
    Ok(fragment_set)
}

/// Estimates the mean GC bias weight across a genomic region by sliding the model's
/// stored window across the region and averaging the model weight at each position.
///
/// The result is used by the caller to inflate `effective_coverage` before fragment
/// generation so that rejection sampling in `apply_gc_bias_to_fragments` still hits
/// the target depth after filtering.
///
/// Returns at least `MIN_WEIGHT` (1e-6) to prevent division-by-zero in the caller
/// when every sampled window has zero weight.
pub fn estimate_region_mean_gc_weight(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    gc_bias_model: &GcBiasModel,
) -> Result<f64, FastaMapError> {
    // Stride of 25 bp gives a dense enough sample for typical read lengths (100–300 bp)
    // without scanning every base. When the window is smaller than 25 bp the stride
    // clamps to the window size so every position is still covered.
    const DEFAULT_STRIDE: usize = 25;
    const MIN_WEIGHT: f64 = 1e-6;

    if region_start >= region_end {
        return Err(FastaMapError::BadCoordinatesError);
    }

    let region_len = region_end - region_start;

    let window_size = gc_bias_model.window_size().min(region_len);

    if window_size == 0 {
        return Ok(MIN_WEIGHT);
    }

    let stride = DEFAULT_STRIDE.min(window_size).max(1);

    let mut total_weight = 0.0;
    let mut sampled_windows = 0usize;

    let last_start = region_end - window_size;
    let mut start = region_start;

    while start <= last_start {
        let end = start + window_size;
        let sequence = sequence_block.get_subseq(start, end)?;
        let weight = gc_bias_model.weight_for_sequence(&sequence);

        total_weight += weight;
        sampled_windows += 1;

        start += stride;
    }

    if sampled_windows == 0 {
        return Ok(MIN_WEIGHT);
    }

    let mean_weight = total_weight / sampled_windows as f64;

    Ok(mean_weight.max(MIN_WEIGHT))
}


#[cfg(test)]
mod tests {
    use super::*;
    use common::rng::NeatRng;
    use common::models::gc_bias_model::GcBiasModel;
    use common::structs::fasta_map::{RegionType, SequenceBlock, SequenceMap};
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

    #[test]
    fn test_apply_gc_bias_to_fragments_keeps_all_fragments_with_default_model() {
        let sequence_block = make_sequence_block(vec![
            A, C, G, T, A, C, G, T, A, C, G, T,
        ]);
        let fragments = vec![(0, 4), (2, 6), (4, 8), (8, 12)];
        let model = GcBiasModel::default();
        let mut rng = make_rng();

        let retained = apply_gc_bias_to_fragments(
            fragments.clone(),
            &sequence_block,
            &model,
            &mut rng,
        ).unwrap();

        assert_eq!(retained, fragments);
    }

    #[test]
    fn test_apply_gc_bias_to_fragments_rejects_zero_weight_fragments() {
        let sequence_block = make_sequence_block(vec![
            A, A, A, A, // 0% GC
            C, G, C, G, // 100% GC
            A, C, G, T, // 50% GC
        ]);

        let fragments = vec![
            (0, 4),  // 0% GC
            (4, 8),  // 100% GC
            (8, 12), // 50% GC
        ];

        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();
        let mut rng = make_rng();

        let retained = apply_gc_bias_to_fragments(
            fragments,
            &sequence_block,
            &model,
            &mut rng,
        ).unwrap();

        assert_eq!(retained, vec![(8, 12)]);
    }

    #[test]
    fn test_apply_gc_bias_to_fragments_uses_relative_weight_against_max_weight() {
        let sequence_block = make_sequence_block(vec![
            A, C, G, T, // 50% GC, weight 0.5
            C, G, C, G, // 100% GC, weight 1.0
        ]);

        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;
        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.max_weight() - 1.0).abs() < 1e-12);
        assert!((model.weight_for_sequence(&sequence_block.sequence[0..4]) - 0.5).abs() < 1e-12);
        assert!((model.weight_for_sequence(&sequence_block.sequence[4..8]) - 1.0).abs() < 1e-12);

        // Seed "gc-bias-test": first draw >= 0.5, so the 50% GC fragment is rejected.
        // Only the max-weight (100% GC) fragment survives.
        let mut rng = make_rng();
        let retained = apply_gc_bias_to_fragments(
            vec![(0, 4), (4, 8)],
            &sequence_block,
            &model,
            &mut rng,
        ).unwrap();
        assert_eq!(retained, vec![(4, 8)], "Seed gc-bias-test: expected only max-weight fragment");

        // Seed "Hello": first draw < 0.5, so the 50% GC fragment is also accepted.
        // Both fragments survive, confirming the accept path works too.
        let mut rng_accept = NeatRng::new_from_seed(&vec!["Hello".to_string()]).unwrap();
        let retained_accept = apply_gc_bias_to_fragments(
            vec![(0, 4), (4, 8)],
            &sequence_block,
            &model,
            &mut rng_accept,
        ).unwrap();
        assert_eq!(retained_accept, vec![(0, 4), (4, 8)], "Seed Hello: expected both fragments retained");
    }

    #[test]
    fn test_apply_gc_bias_to_fragments_returns_error_for_out_of_bounds_fragment() {
        let sequence_block = make_sequence_block(vec![A, C, G, T]);
        let fragments = vec![(0, 4), (2, 10)];
        let model = GcBiasModel::default();
        let mut rng = make_rng();

        let result = apply_gc_bias_to_fragments(
            fragments,
            &sequence_block,
            &model,
            &mut rng,
        );

        assert!(
            result.is_err(),
            "Expected out-of-bounds fragment to return an error"
        );
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_uniform_model() {
        let sequence_block = make_sequence_block(vec![
            A, A, A, A,
            A, C, G, T,
            C, G, C, G,
            A, C, G, T,
        ]);
        let model = GcBiasModel::default();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            sequence_block.sequence.len(),
            &model,
        ).unwrap();

        assert!((mean_weight - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_averages_sampled_windows() {
        let sequence_block = make_sequence_block(vec![
            A, A, A, A, // 0% GC, weight 0.2
            A, C, G, T, // 50% GC, weight 1.0
            C, G, C, G, // 100% GC, weight 0.4
        ]);

        let mut weights = weights_with_value(0.0);
        weights[0] = 0.2;
        weights[50] = 1.0;
        weights[100] = 0.4;

        let model = GcBiasModel::from_weights(weights, 4).unwrap();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            12,
            &model,
        ).unwrap();

        let expected = (0.2 + 1.0 + 0.4) / 3.0;
        assert!(
            (mean_weight - expected).abs() < 1e-12,
            "Expected {}, got {}",
            expected,
            mean_weight
        );
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_uses_whole_region_when_window_is_larger() {
        let sequence_block = make_sequence_block(vec![
            A, C, G, T, // 50% GC
        ]);

        let mut weights = weights_with_value(0.0);
        weights[50] = 0.75;

        let model = GcBiasModel::from_weights(weights, 100).unwrap();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            4,
            &model,
        ).unwrap();

        assert!((mean_weight - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_rejects_empty_region() {
        let sequence_block = make_sequence_block(vec![A, C, G, T]);
        let model = GcBiasModel::default();

        let result = estimate_region_mean_gc_weight(
            &sequence_block,
            2,
            2,
            &model,
        );

        assert!(
            result.is_err(),
            "Expected empty region to return an error"
        );
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_rejects_reversed_region() {
        let sequence_block = make_sequence_block(vec![A, C, G, T]);
        let model = GcBiasModel::default();

        let result = estimate_region_mean_gc_weight(
            &sequence_block,
            3,
            2,
            &model,
        );

        assert!(
            result.is_err(),
            "Expected reversed region to return an error"
        );
    }

    #[test]
    fn test_estimate_region_mean_gc_weight_clamps_zero_mean_to_minimum() {
        let sequence_block = make_sequence_block(vec![
            A, A, A, A,
            A, C, G, T,
            C, G, C, G,
        ]);

        let mut weights = weights_with_value(0.0);
        weights[0] = 0.0;
        weights[50] = 0.0;
        weights[100] = 1.0;

        let model = GcBiasModel::from_weights(weights, 4).unwrap();

        // Only sample the first two windows, both of which have zero weight.
        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            8,
            &model,
        ).unwrap();

        assert!(
            mean_weight > 0.0,
            "Expected zero mean to be clamped to a small positive value"
        );
        assert!(
            mean_weight <= 1e-5,
            "Expected clamped mean to be very small; got {}",
            mean_weight
        );
    }

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
        assert_eq!(cover[0], (15, 315))
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
            &fragment_model,
            &mut rng,
        ).unwrap();
        assert!(reads.contains(&(0, 382)));
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
            &fragment_model,
            &mut rng,
        ).unwrap();

        let fragment_model = FragmentLengthModel::default().unwrap();
        let run2 = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            &fragment_model,
            &mut rng,
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
            &fragment_model,
            &mut rng,
        ).unwrap();
        assert!(!reads.is_empty())
    }

    // Mirrors the runner's per-region logic: optionally inflate effective_coverage using
    // estimate_region_mean_gc_weight, generate fragments, then apply rejection sampling.
    // Uniform models short-circuit both steps, matching runner.rs behaviour exactly.
    fn fragment_count_after_gc_pipeline(
        sequence_block: &SequenceBlock,
        read_len: usize,
        target_coverage: usize,
        model: &GcBiasModel,
        normalize: bool,
        rng: &mut NeatRng,
    ) -> usize {
        let seq_len = sequence_block.sequence.len();
        let effective_coverage = if !model.is_uniform() && normalize {
            let mean_weight = estimate_region_mean_gc_weight(
                sequence_block, 0, seq_len, model,
            ).unwrap();
            (target_coverage as f64 / mean_weight).round() as usize
        } else {
            target_coverage
        };
        let fragment_model = FragmentLengthModel::default().unwrap();
        let fragments = generate_fragments(
            seq_len, read_len, 0, 0, effective_coverage, &fragment_model, rng,
        ).unwrap();
        if model.is_uniform() {
            fragments.len()
        } else {
            apply_gc_bias_to_fragments(fragments, sequence_block, model, rng).unwrap().len()
        }
    }

    #[test]
    fn test_coverage_normalization_compensates_for_gc_bias() {
        // 10,000 bp of repeating ACGT — exactly 50% GC throughout.
        let sequence: Vec<_> = std::iter::repeat([A, C, G, T])
            .take(2500)
            .flatten()
            .collect();
        let sequence_block = make_sequence_block(sequence);
        let read_len = 100;
        let target_coverage = 10;

        // Model: 50% GC has weight 0.5, max_weight = 1.0.
        // Normalization inflates coverage by 2x; rejection keeps 50% → net ≈ target.
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;
        let model = GcBiasModel::from_weights(weights, 100).unwrap();

        // Baseline: uniform model, no rejection — establishes the "target fragment count"
        // for this sequence length and coverage without any GC influence.
        let mut rng = make_rng();
        let baseline = fragment_count_after_gc_pipeline(
            &sequence_block, read_len, target_coverage, &GcBiasModel::default(), false, &mut rng,
        );

        // Normalized: inflation + rejection should land near the baseline.
        let mut rng = make_rng();
        let normalized = fragment_count_after_gc_pipeline(
            &sequence_block, read_len, target_coverage, &model, true, &mut rng,
        );

        let tolerance = baseline / 4; // ±25%
        assert!(
            normalized >= baseline.saturating_sub(tolerance) && normalized <= baseline + tolerance,
            "Normalized count ({}) should be within ±25% of baseline ({})",
            normalized, baseline
        );
    }

    #[test]
    fn test_without_normalization_coverage_reflects_gc_bias() {
        // Same sequence and model as the normalization test.
        let sequence: Vec<_> = std::iter::repeat([A, C, G, T])
            .take(2500)
            .flatten()
            .collect();
        let sequence_block = make_sequence_block(sequence);
        let read_len = 100;
        let target_coverage = 10;

        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;
        let model = GcBiasModel::from_weights(weights, 100).unwrap();

        // Without normalization: effective_coverage == target, then 50% rejection
        // → roughly half the fragment count compared to normalized.
        let mut rng = make_rng();
        let unnormalized = fragment_count_after_gc_pipeline(
            &sequence_block, read_len, target_coverage, &model, false, &mut rng,
        );

        let mut rng = make_rng();
        let normalized = fragment_count_after_gc_pipeline(
            &sequence_block, read_len, target_coverage, &model, true, &mut rng,
        );

        assert!(
            normalized > unnormalized * 3 / 2,
            "Normalized ({}) should be at least 1.5x more than unnormalized ({})",
            normalized, unnormalized
        );
    }

    #[test]
    fn test_uniform_model_coverage_unaffected_by_normalize_flag() {
        // For a uniform model, both paths short-circuit identically: no inflation,
        // no rejection sampling. The flag has no effect on the fragment count.
        let sequence: Vec<_> = std::iter::repeat([A, C, G, T])
            .take(2500)
            .flatten()
            .collect();
        let sequence_block = make_sequence_block(sequence);
        let model = GcBiasModel::default();

        let mut rng = make_rng();
        let with_flag = fragment_count_after_gc_pipeline(
            &sequence_block, 100, 10, &model, true, &mut rng,
        );
        let mut rng = make_rng();
        let without_flag = fragment_count_after_gc_pipeline(
            &sequence_block, 100, 10, &model, false, &mut rng,
        );

        assert_eq!(
            with_flag, without_flag,
            "Uniform model must produce identical counts regardless of the normalize flag"
        );
    }

    #[test]
    fn test_near_zero_weight_region_does_not_overflow() {
        // A 100% GC sequence with a model that gives zero weight to 100% GC bins.
        // estimate_region_mean_gc_weight clamps to MIN_WEIGHT (1e-6), which without a
        // cap would inflate effective_coverage to coverage / 1e-6 = 10_000_000.
        // This test verifies the pipeline completes in finite time and memory.
        let sequence: Vec<_> = std::iter::repeat([C, G]).take(5000).flatten().collect();
        let sequence_block = make_sequence_block(sequence);
        let read_len = 100;
        let target_coverage = 10;
        const MAX_COVERAGE_MULTIPLIER: usize = 100;

        // Only bin 50 has weight; every 100% GC window scores 0.0 → mean clamps to MIN_WEIGHT.
        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;
        let model = GcBiasModel::from_weights(weights, 100).unwrap();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block, 0, sequence_block.sequence.len(), &model,
        ).unwrap();

        let inflated = (target_coverage as f64 / mean_weight).round() as usize;
        let cap = target_coverage * MAX_COVERAGE_MULTIPLIER;
        let effective_coverage = inflated.min(cap);

        // Inflation would be ~10_000_000 without the cap; with it, max is 1_000.
        assert!(
            effective_coverage <= cap,
            "Effective coverage {} exceeds cap {}",
            effective_coverage, cap
        );

        // Pipeline must complete without OOM or panic.
        let mut rng = make_rng();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let fragments = generate_fragments(
            sequence_block.sequence.len(), read_len, 0, 0,
            effective_coverage, &fragment_model, &mut rng,
        ).unwrap();
        let _ = apply_gc_bias_to_fragments(fragments, &sequence_block, &model, &mut rng).unwrap();
    }
}
