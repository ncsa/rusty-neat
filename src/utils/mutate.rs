use std::cmp::max;
use std::collections::HashMap;
use rand::prelude::{ThreadRng, IndexedRandom};
use rand::seq::SliceRandom;
use rand::Rng;
use log::debug;
use itertools::izip;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
struct NucModel {
    // This simple nucleotide model simple tracks the base to mutate from
    // and the weights of each are turned into a vector. For example, if the weight count for "A"
    // is 8, then the vector will be [0, 0, 0, 0, 0, 0, 0, 0].
    // 0: A, 1: C, 2: G, 3: T
    base: u8,
    // vector of 0s, the length of which is the weight of the A in the vector.
    a: Vec<u8>,
    // vector of 1s, the length of which is the weight of the C in the vector.
    c: Vec<u8>,
    // vector of 2s, the length of which is the weight of the G in the vector.
    g: Vec<u8>,
    // vector of 3s, the length of which is the weight of the T in the vector.
    t: Vec<u8>,
}

impl NucModel {
    fn from(base: u8, weights: Vec<usize>) -> Self {
        Self {
            base,
            a: vec![0; weights[0]],
            c: vec![1; weights[1]],
            g: vec![2; weights[2]],
            t: vec![3; weights[3]],
        }
    }

    fn choose_new_nuc(&self, mut rng: &mut ThreadRng) -> u8 {
        // concatenate the vectors together
        let mut choices = [&self.a[..], &self.c[..], &self.g[..], &self.t[..]].concat();
        // shuffle the weighted list, based on whatever statistics we come up with.
        choices.shuffle(&mut rng);
        // start with the first choice
        let mut i = 0;
        // Look for the first nuc that isn't the old nuc.
        while choices[i] == self.base {
            i += 1;
        }
        choices[i]
    }
}

pub fn mutate_fasta(
    file_struct: &HashMap<String, Vec<u8>>,
    mut rng: &mut ThreadRng
) -> (Box<HashMap<String, Vec<u8>>>, Box<HashMap<String, Vec<(usize, u8, u8)>>>) {
    /*
    Takes:
        file_struct: a hashmap of contig names (keys) and a vector
            representing the reference sequence.
        ploidy: The number of copies of the genome within an organism's cells
        rng: random number generator for the run

    Returns:
        A tuple with pointers to:
          A hashmap with keys that are contig names and a vector with the mutated sequence
          A vector of tuples containing the location and alt of each variant

    This function performs a basic calculation (length x mutation rate +/- a random amount)
    and chooses that many positions along the sequence to mutate. It then builds a return
    string that represents the altered sequence and stores all the variants.
     */
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    let mut return_struct: HashMap<String, Vec<u8>> = HashMap::new(); // the mutated sequences

    // hashmap with keys of the contig names with a list of positions and alts under the contig.
    let mut all_variants: HashMap<String, Vec<(usize, u8, u8)>> = HashMap::new();

    for (name, sequence) in file_struct {
        // Mutations for this contig
        let contig_mutations: Vec<(usize, u8, u8)>;
        // The length of this sequence
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Clone the reference to create mutations
        let mut mutated_record: Vec<u8> = sequence.clone();
        // Calculate how many mutations to add
        let mut rough_num_positions: f64 = sequence_length as f64 * MUT_RATE;
        // Add or subtract a few extra positions.
        rough_num_positions += {
            // A random amount up to 10% of the reads
            let factor: f64 = rng.gen::<f64>() * 0.10;
            // 25% of the time subtract, otherwise we'll add.
            let sign: f64 = if rng.gen_bool(0.25) { -1.0 } else { 1.0 };
            // add or subtract up to 10% of the reads.
            rough_num_positions * (sign * factor)
        };
        // Round the number of positions to the nearest usize.
        // If negative or no reads, we still want at least 1 mutation per contig.
        let num_positions = max(1, rough_num_positions.round() as usize);

        (mutated_record, contig_mutations) = mutate_sequence(
            mutated_record, num_positions, &mut rng
        );
        // Add to the return struct and variants map.
        return_struct.entry(name.clone()).or_insert(mutated_record.clone());
        all_variants.entry(name.clone()).or_insert(contig_mutations);
    }

    (Box::new(return_struct), Box::new(all_variants))
}

fn mutate_sequence(
    sequence: Vec<u8>,
    num_positions: usize,
    mut rng: &mut ThreadRng
) -> (Vec<u8>, Vec<(usize, u8, u8)>) {
    /*
    Takes:
        sequence: A u8 vector representing a sequence of DNA
        num_positions: The number of mutations to add to this sequence
        rng: random number generator for the run

    returns a tuple with:
        Vec<u8> = the sequence itself
        Vec(usize, u8) = the position of the snp and the alt allele for that snp.

    Takes a vector of u8's and mutate a few positions at random. Returns the mutated sequence and
    a list of tuples with the position and the alts of the SNPs.
     */
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.clone();
    // Randomly select num_positions from positions, weighted by gc bias and whatever. For now
    // all he weights are just equal.
    let weights = vec![1; mutated_record.len()];
    // zip together weights and positions
    let mut weighted_positions: Vec<(usize, i32)> = Vec::new();
    // izip!() accepts iterators and/or values with IntoIterator.
    for (x, y) in izip!(&mut (0..mutated_record.len()), &weights) {
        weighted_positions.push((x, *y))
    }
    // now choose a random selection of num_positions without replacement
    let mutated_elements: Vec<&(usize, i32)> = weighted_positions
        .choose_multiple_weighted(&mut rng, num_positions, |x| x.1)
        .unwrap()
        .collect::<Vec<_>>();
    // Build mutation model
    /*
    Using NEAT's original mutation model of
    [[0.0, 0.15, 0.7, 0.15],
     [0.15, 0.0, 0.15, 0.7],
     [0.7, 0.15, 0.0, 0.15],
     [0.15, 0.7, 0.15, 0.0]]

     multiply each value by 20, which produces integer values
     */
    let mut_a = NucModel::from(0, vec![0, 3, 14, 3]);
    let mut_c = NucModel::from(1, vec![3, 0, 3, 14]);
    let mut_g = NucModel::from(2, vec![14, 3, 0, 3]);
    let mut_t = NucModel::from(3, vec![3, 14, 3, 0]);

    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<(usize, u8, u8)> = Vec::new();
    // for each index, picks a new base
    for (index, _weight) in mutated_elements {
        let reference_base = sequence[*index];
        mutated_record[*index] = match reference_base {
            0 => mut_a.choose_new_nuc(&mut rng),
            1 => mut_c.choose_new_nuc(&mut rng),
            2 => mut_g.choose_new_nuc(&mut rng),
            3 => mut_t.choose_new_nuc(&mut rng),
            _ => 4
        };
        // add the location and alt base for the variant
        sequence_variants.push((*index, mutated_record[*index], reference_base))
    }
    (mutated_record, sequence_variants)
}