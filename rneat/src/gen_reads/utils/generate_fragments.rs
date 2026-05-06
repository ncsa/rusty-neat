// This is the core functionality of NEAT. Generate reads turns a mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to give the contig the proper
// read coverage. Generate reads uses this to create a list of coordinates to take slices from
// the mutated fasta file. These will either be read-length fragments or fragment model length
// fragments.
use crate::{
    common::models::fragment_length::FragmentLengthModel, 
    gen_reads::errors::GenerateReadsErrors
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
) -> Result<Vec<(usize, usize)>, GenerateReadsErrors> {
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

pub fn apply_gc_bias_to_fragments(
    fragments: Vec<(usize, usize)>,
    sequence_block: &SequenceBlock,
    gc_bias_model: &GcBiasModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsErrors> {
    let max_weight = gc_bias_model.max_weight();
    let mut retained = Vec::new();

    for (start, end) in fragments {
        let seq = sequence_block.get_subseq(start, end)?;
        let weight = gc_bias_model.weight_for_sequence(&seq);
        let accept_prob = if max_weight > 0.0 {
            weight / max_weight
        } else {
            0.0
        };

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
) -> Result<Vec<(usize, usize)>, GenerateReadsErrors> {
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

pub fn estimate_region_mean_gc_weight(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    gc_window_size: usize,
    gc_bias_model: &GcBiasModel,
) -> Result<f64, FastaMapError> {
    const DEFAULT_STRIDE: usize = 25;
    const MIN_WEIGHT: f64 = 1e-6;

    if region_start >= region_end {
        return Err(FastaMapError::BadCoordinatesError);
    }

    let region_len = region_end - region_start;

    if region_len == 0 {
        return Err(FastaMapError::BadCoordinatesError);
    }

    let window_size = gc_window_size.min(region_len);

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

        let model = GcBiasModel::from_weights(weights).unwrap();
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

        let fragments = vec![
            (0, 4),
            (4, 8),
        ];

        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.max_weight() - 1.0).abs() < 1e-12);
        assert!((model.weight_for_sequence(&sequence_block.sequence[0..4]) - 0.5).abs() < 1e-12);
        assert!((model.weight_for_sequence(&sequence_block.sequence[4..8]) - 1.0).abs() < 1e-12);

        let mut rng = make_rng();

        let retained = apply_gc_bias_to_fragments(
            fragments,
            &sequence_block,
            &model,
            &mut rng,
        ).unwrap();

        // The 100% GC fragment has accept_prob 1.0 and must always be retained.
        assert!(
            retained.contains(&(4, 8)),
            "Expected max-weight fragment to be retained"
        );

        // The 50% GC fragment has accept_prob 0.5, so it may or may not be retained
        // depending on the deterministic RNG draw. We do not assert either outcome here.
        assert!(
            retained.iter().all(|fragment| *fragment == (0, 4) || *fragment == (4, 8)),
            "Unexpected retained fragment list: {:?}",
            retained
        );
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
            4,
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

        let model = GcBiasModel::from_weights(weights).unwrap();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            12,
            4,
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

        let model = GcBiasModel::from_weights(weights).unwrap();

        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            4,
            100,
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
            4,
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
            4,
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

        let model = GcBiasModel::from_weights(weights).unwrap();

        // Only sample the first two windows, both of which have zero weight.
        let mean_weight = estimate_region_mean_gc_weight(
            &sequence_block,
            0,
            8,
            4,
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
}
