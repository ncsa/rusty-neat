// This is the core functionality of NEAT. Generate reads turns a mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to give the contig the proper
// read coverage. Generate reads uses this to create a list of coordinates to take slices from
// the mutated fasta file. These will either be read-length fragments or fragment model length
// fragments.
use std::collections::{HashSet, VecDeque};
use rand::RngCore;
use rand::seq::SliceRandom;
use rand_distr::{Normal, Distribution};
use utils::neat_rng::NeatRng;
use utils::nucleotides::Nuc;
fn cover_dataset(
    span_length: usize,
    read_length: usize,
    mut fragment_pool: Vec<usize>,
    coverage: usize,
    mut rng: &mut NeatRng,
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
    let mut read_set: Vec<(usize, usize)> = vec![];

    let mut cover_fragment_pool: VecDeque<usize>;
    if fragment_pool.is_empty() {
        // set the shuffled fragment pool just equal to an instance of read_length
        cover_fragment_pool = VecDeque::from([read_length]);
    } else {
        // shuffle the fragment pool
        fragment_pool.shuffle(&mut rng);
        cover_fragment_pool = VecDeque::from(fragment_pool)
    }
    // Gap size to keep track of how many uncovered bases we have per layer, to help decide if we
    // need more layers
    let mut gap_size: usize = 0;
    let mut layer_count: usize = 0;
    // start this party off at zero.
    let mut start: usize = 0;
    // create coverage number of layers
    while layer_count <= coverage {
        let fragment_length = cover_fragment_pool[0];
        cover_fragment_pool.push_back(fragment_length);
        let temp_end = start+fragment_length;
        if temp_end > span_length {
            // TODO some variation on this modulo idea will work for bacterial reads
            start = temp_end % span_length;
            gap_size += start;
            //
            if gap_size >= span_length {
                // if we have accumulated enough gap, then we need to run the same layer again.
                // We'll reset gap size but not increment layer_count.
                gap_size = gap_size % span_length;
                continue;
            } else {
                layer_count += 1;
                continue
            }
        }
        read_set.push((start, temp_end));
        // insert size is the number of bases between reads in the fragment for paired ended reads
        // if these are singled ended reads, then the insert size will always be -read_length
        if fragment_length > (read_length * 2) {
            // if there's any insert size on paired ended reads, we'll add
            // that to the gap to ensure adequate coverage.
            gap_size += fragment_length - (read_length * 2)
        };
        // Picks a number between zero and a quarter of a read length
        let wildcard: usize = (rng.next_u32() % (read_length/4) as u32) as usize;
        // adds to the start to give it some spice
        start += temp_end + wildcard;
        // sanity check. If we are already out of bounds, take the modulo
        if start >= span_length {
            // get us back in bounds
            start = start % span_length;
            // add the gap
            gap_size += start;
        } else {
            // still in bounds, just add the gap
            gap_size += wildcard;
        }
    }
    read_set
}

pub fn generate_reads(
    mutated_sequence: &Vec<Nuc>,
    read_length: &usize,
    coverage: &usize,
    paired_ended: bool,
    mean: Option<f64>,
    st_dev: Option<f64>,
    mut rng: &mut NeatRng,
) -> Result<Box<HashSet<Vec<Nuc>>>, &'static str>{
    // Takes:
    // mutated_sequence: a vector of u8's representing the mutated sequence.
    // read_length: the length ef the reads for this run
    // coverage: the average depth of coverage for this run
    // rng: the random number generator for the run
    // Returns:
    // HashSet of vectors representing the read sequences, stored on the heap in box.
    //
    // This takes a mutated sequence and produces a set of reads based on the mutated sequence. For
    // paired ended reads, this will generate a set of reads from each end, by taking the reverse
    // complement int the output
    let mut fragment_pool: Vec<usize> = Vec::new();
    if paired_ended {
        let num_frags = (mutated_sequence.len() / read_length) * (coverage * 2);
        let fragment_distribution = Normal::new(mean.unwrap(), st_dev.unwrap()).unwrap();
        // add fragments to the fragment pool
        for _ in 0..num_frags {
            let frag = fragment_distribution.sample(&mut rng).round() as usize;
            fragment_pool.push(frag);
        }
    }
    // set up some defaults and storage
    let mut read_set: HashSet<Vec<Nuc>> = HashSet::new();
    // length of the mutated sequence
    let seq_len = mutated_sequence.len();
    // Generate a vector of read positions
    let read_positions: Vec<(usize, usize)> = cover_dataset(
        seq_len,
        *read_length,
        fragment_pool,
        *coverage,
        &mut rng,
    );
    // Generate the reads from the read positions.
    for (start, end) in read_positions {
        read_set.insert(mutated_sequence[start..end].into());
    }
    // puts the reads in the heap.
    if read_set.is_empty() {
        Err("No reads generated")
    } else {
        Ok(Box::new(read_set))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use utils::neat_rng::NeatRng;
    use utils::nucleotides::Nuc::*;

    #[test]
    fn test_cover_dataset() {
        let span_length = 100;
        let read_length = 10;
        let fragment_pool = vec![10];
        let coverage = 1;
        let mut rng = NeatRng::seed_from_u64(0);

        let cover = cover_dataset(
            span_length,
            read_length,
            fragment_pool,
            coverage,
            &mut rng,
        );
        assert_eq!(cover[0], (0,10))
    }

    #[test]
    fn test_gap_function() {
        let span_length = 100_000;
        let read_length = 100;
        let fragment_pool = vec![300];
        let coverage = 1;
        let mut rng = NeatRng::seed_from_u64(0);

        let cover = cover_dataset(
            span_length,
            read_length,
            fragment_pool,
            coverage,
            &mut rng,
        );
        assert_eq!(cover[0], (0, 300))
    }

    #[test]
    fn test_generate_reads_single() {
        let mutated_sequence = vec![
            A, A, C, A, T, T, T, T, A, A, A, A, A, G, G, G, N, N, N, N
        ];
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mean = None;
        let st_dev = None;
        let mut rng = NeatRng::seed_from_u64(0);
        let reads = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        ).unwrap();
        println!("{:?}", reads);
        assert!(reads.contains(&(vec![A, A, C, A, T, T, T, T, A, A])));
    }

    #[test]
    fn test_seed_rng() {
        let mutated_sequence = vec![
            A, A, C, A, T, T, T, T, A, A, A, A, A, G, G, G, N, N, N, N
        ];
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mean = None;
        let st_dev = None;
        let mut rng = NeatRng::seed_from_u64(0);
        let run1 = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        ).unwrap();

        let run2 = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        ).unwrap();

        assert_eq!(run1, run2)
    }

    #[test]
    fn test_generate_reads_paired() {
        let mutated_sequence: Vec<Nuc> = std::iter::repeat(A).take(100_000).collect();
        let read_length = 100;
        let coverage = 1;
        let paired_ended = true;
        let mean = Some(200.0);
        let st_dev = Some(1.0);
        let mut rng = NeatRng::seed_from_u64(0);
        let reads = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        );
        println!("{:?}", reads);
        assert!(!reads.unwrap().is_empty())
    }
}