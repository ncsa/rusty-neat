use std::collections::HashSet;
use rand::RngCore;
use rand::rngs::ThreadRng;

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    coverage: usize,
    rng: &mut ThreadRng,
) -> Vec<(usize, usize)> {
    // This function selects the positions of the reads. It starts at the beginning and goes out
    // one read length, then picks a random jump between 0 and half the read length to move
    // And picks those coordinates for a second read. Once the tail of the read is past the end,
    // we start over again at 0.
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

pub fn generate_reads(
    mutated_sequences: &Vec<Vec<u8>>,
    read_length: &usize,
    coverage: &usize,
    rng: &mut ThreadRng,
) -> Box<HashSet<Vec<u8>>> {
    // Infer ploidy
    let ploidy = mutated_sequences.len();
    // Need X/ploidy reads per ploid, but probably fewer
    let coverage_per_ploid = coverage / ploidy + 1;

    // set up some defaults and storage
    let mut read_set: HashSet<Vec<u8>> = HashSet::new();

    for i in 0..ploidy {
        let seq_len = *(&mutated_sequences[i].len());
        // Generate a vector of read positions
        let read_positions: Vec<(usize, usize)> = cover_dataset(
            seq_len,
            *read_length,
            coverage_per_ploid,
            rng,
        );

        // Generate the reads from the read positions.
        for (start, end) in read_positions {
            read_set.insert(mutated_sequences[i][start..end].into());
        }
    }
    // puts the reads in the heap.
    Box::new(read_set)
}