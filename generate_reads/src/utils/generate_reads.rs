// This is the core functionality of NEAT. Generate reads turns a mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to give the contig the proper
// read coverage. Generate reads uses this to create a list of coordinates to take slices from
// the mutated fasta file. These will either be read-length fragments or fragment model length
// fragments.
use crate::common::{
    models::fragment_length::FragmentLengthModel,
};
use log::debug;
use std::collections::VecDeque;
use simple_rng::NeatRng;
use crate::errors::GenerateReadsErrors;

pub fn generate_reads(
    sequence_length: usize,
    read_length: usize,
    coverage: usize,
    paired_ended: bool,
    fragment_model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize, Option<usize>, usize)>, GenerateReadsErrors> {
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

    // todo We need to select for non-N areas for the reads, or at least have reads be more than 1/2
    //  true bases, I think.
    let mut fragment_pool: Vec<usize> = Vec::new();
    if paired_ended {
        let num_frags = (sequence_length / read_length) * coverage;
        // add fragments to the fragment pool
        for _ in 0..num_frags {
            let frag = fragment_model.generate_fragment(rng.random()?)?;
            fragment_pool.push(frag);
        }
    }
    // Generate a vector of read positions
    debug!("Generating read coordinates.");
    let fragments: Vec<(usize, usize)> = cover_dataset(
        sequence_length,
        read_length,
        fragment_pool,
        coverage,
        rng
    );

    // TODO Generate single or paired ended reads from the fragments
    let mut read_set: Vec<(usize, usize, Option<usize>, usize)> = Vec::new();
    if paired_ended {
        for fragment in fragments {
            // Push a paired ended read off the fragment
            read_set.push((
                fragment.0, 
                fragment.0 + read_length, 
                Some(fragment.1 - read_length), 
                fragment.1
            ))
        }
    } else {
        // Push a single ended read
        for fragment in fragments {
            read_set.push((
                fragment.0, 
                fragment.0 + read_length,
                None,
                fragment.1,
            ))
        }
    }

    debug!("Outputting read set");
    if read_set.is_empty() {
        Err(GenerateReadsErrors::GenerateReadsError)
    } else {
        Ok(read_set)
    }
}

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    mut fragment_pool: Vec<usize>,
    coverage: usize,
    rng: &mut NeatRng,
) -> Vec<(usize, usize)> {
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
    if fragment_pool.is_empty() {
        // set the shuffled fragment pool just equal to an instance of read_length
        cover_fragment_pool = VecDeque::from([read_length*2]);
    } else {
        // shuffle the fragment pool
        rng.shuffle_in_place(&mut fragment_pool)
            .expect("Generate reads failed while shuffling the fragment pool");
        cover_fragment_pool = VecDeque::from(fragment_pool)
    }
    // Gap size to keep track of how many uncovered bases we have per layer, to help decide if we
    // need more layers
    let mut layer_count: usize = 0;
    // start this party off at zero.
    let mut start: usize = 0;
    // create coverage number of layers
    while layer_count < coverage {
        let fragment_length = cover_fragment_pool.pop_front().unwrap();
        let temp_end = start + fragment_length;
        cover_fragment_pool.push_back(fragment_length.clone());
        if temp_end > span_length {
            // TODO some variation on this modulo idea will work for bacterial reads
            start = temp_end % span_length;
            layer_count += 1;
            continue
        }
        fragment_set.push((start, temp_end));
        // Picks a number between zero and a quarter of a read length
        let wildcard: usize = (rng.rand_u32().unwrap() % 10) as usize;
        // adds to the start to give it some spice
        start += temp_end + wildcard;
        // sanity check. If we are out of bounds, take the modulo
        if start >= span_length {
            // get us back in bounds
            start = start % span_length;
        }
    }
    fragment_set
}

#[cfg(test)]
mod tests {
    use super::*;
    use simple_rng::NeatRng;

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
        let cover = cover_dataset(span_length, read_length, fragment_pool, coverage, &mut rng);
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
        let cover = cover_dataset(span_length, read_length, fragment_pool, coverage, &mut rng);
        assert_eq!(cover[0], (0, 300))
    }

    #[test]
    fn test_generate_reads_single() {
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_reads(
            2000,
            read_length,
            coverage,
            paired_ended,
            &fragment_model,
            &mut rng,
        ).unwrap();
        assert!(reads.contains(&(0, 10, None, 20)));
    }

    #[test]
    fn test_seed_rng() {
        let sequnce = vec![
            0, 0, 2, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 2, 2, 2, 4, 4, 4, 4
        ];
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let run1 = generate_reads(
            sequnce.len(),
            read_length,
            coverage,
            paired_ended,
            &fragment_model,
            &mut rng,
        ).unwrap();

        let fragment_model = FragmentLengthModel::default().unwrap();
        let run2 = generate_reads(
            sequnce.len(),
            read_length,
            coverage,
            paired_ended,
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
        let paired_ended = true;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_reads(
            sequence.len(),
            read_length,
            coverage,
            paired_ended,
            &fragment_model,
            &mut rng,
        ).unwrap();
        assert!(!reads.is_empty())
    }
}
