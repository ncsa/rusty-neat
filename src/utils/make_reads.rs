use std::collections::HashSet;
use rand::RngCore;
use rand::rngs::ThreadRng;

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    coverage: usize,
    rng: &mut ThreadRng,
) -> Vec<(usize, usize)> {
    /*
    Takes:
        span_length: Total number of bases in the sequence
        read_length: The length of the reads for this run
        coverage: The coverage depth for the reads
    Returns:
        A vector of tuples (usize, usize), denoting the start and end positions of the reads

    This function selects the positions of the reads. It starts at the beginning and goes out
    one read length, then picks a random jump between 0 and half the read length to move
    And picks those coordinates for a second read. Once the tail of the read is past the end,
    we start over again at 0.
     */
    // todo currently the reads look a little light so we may need to improve this calculation.
    let mut read_set: Vec<(usize, usize)> = vec![];
    let mut start: usize = 0;
    'outer: for _ in 0..coverage {
        while start < span_length {
            let temp_end = start+read_length;
            if temp_end > span_length {
                start = 0;
                continue 'outer;
            }
            read_set.push((start, temp_end));
            // Picks a number between zero and half of a read length
            let wildcard: usize = (rng.next_u32() % (read_length/4) as u32) as usize;
            start += read_length + wildcard; // adds to the start to give it some spice
        }
    }
    read_set
}

fn complement(nucleotide: u8) -> u8 {
    /*
    0 = A, 1 = C, 2 = G, 3 = T,
    matches with the complement of each nucleotide.

    Todo: make this part of a struct to standardize across the program.
     */
    return match nucleotide {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

fn reverse_complement(sequence: &Vec<u8>) -> Vec<u8> {
    /*
    Returns the reverse complement of a vector of u8's representing a DNA sequence.
     */
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(complement(sequence[i]))
    }
    rev_comp
}

pub fn generate_reads(
    mutated_sequence: &Vec<u8>,
    read_length: &usize,
    coverage: &usize,
    paired_ended: bool,
    mut rng: &mut ThreadRng,
) -> Box<HashSet<Vec<u8>>> {
    /*
    Takes:
        mutated_sequence: a vector of u8's representing the mutated sequence.
        read_length: the length ef the reads for this run
        coverage: the average depth of coverage for this run
        rng: the random number generator for the run
    Returns:
        HashSet of vectors representing the read sequences, stored on the heap in box.

    This takes a mutated sequence and produces a set of reads based on the mutated sequence. For
    paired ended reads, this will generate a set of reads from each end, by taking the reverse
    complement int the output
     */

    if paired_ended {
        let rev_comp = reverse_complement(mutated_sequence);
        todo!();
    }

    // set up some defaults and storage
    let mut read_set: HashSet<Vec<u8>> = HashSet::new();

    let seq_len = mutated_sequence.len();
    // Generate a vector of read positions
    let read_positions: Vec<(usize, usize)> = cover_dataset(
        seq_len,
        *read_length,
        *coverage,
        &mut rng,
    );

    // Generate the reads from the read positions.
    for (start, end) in read_positions {
        read_set.insert(mutated_sequence[start..end].into());
    }
    // puts the reads in the heap.
    Box::new(read_set)
}