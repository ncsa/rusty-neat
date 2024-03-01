use std::collections::HashSet;
use rand::RngCore;
use rand::rngs::ThreadRng;

fn cover_dataset(
    span_length: usize,
    read_length: &usize,
    strict_read_length: bool,
    coverage: &usize,
    rng: &mut ThreadRng,
) -> Vec<(usize, usize)> {
    let mut read_set: Vec<(usize, usize)> = vec![];
    let half_read = read_length/2;
    let mut start: usize = 0;
    'outer: for _ in 0..*coverage {
        while start < span_length {
            let mut length_wildcard: usize = 0;
            if strict_read_length == false {
                length_wildcard = (rng.next_u32() % ((read_length / 10) as u32)) as usize;
            }
            let temp_end = start+read_length+length_wildcard;
            if temp_end > span_length {
                start = 0;
                continue 'outer;
            }
            read_set.push((start, temp_end));
            let wildcard: usize = (rng.next_u32() % (half_read/2) as u32) as usize; // Picks a number between zero and quarter of a read length
            start += half_read + wildcard; // adds to the start to give it some spice
        }
    }
    read_set
}

pub fn generate_reads(
    mutated_sequence: &Vec<u8>,
    read_length: &usize,
    coverage: &usize,
    rng: &mut ThreadRng,
    strict_read_length: Option<bool>,
) -> Box<HashSet<Vec<u8>>> {
    strict_read_length.unwrap_or(true);

    // Generate a vector of read positions
    let read_positions: Vec<(usize, usize)> = cover_dataset(
        mutated_sequence.len(),
        read_length,
        strict_read_length.unwrap(),
        coverage,
        rng,
    );

    // set up some defaults and storage
    let mut read_set: HashSet<Vec<u8>> = HashSet::new();

    // Generate the reads from the read positions.
    for (start, end) in read_positions {
        read_set.insert(mutated_sequence[start..end].into());
    }
    Box::new(read_set)
}