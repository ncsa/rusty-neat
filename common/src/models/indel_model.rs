use simple_rng::{DiscreteDistribution, NeatRng};

// The following section are the models for each type of variant. In order to create the variant,
// we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
// Note that indels are actually two types, insertions and deletions, but they are usually classed
// together in the literature and are usually short. Insertions are often slips that cause sections
// to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
// deletions, since they are treated similarly by variant calling software.
//
// We have the variants separated out in the variants struct, into Insertion and Deletion, but they
// will share this model to avoid duplicating code, since code-wise they are closely linked.

#[derive(Debug, Clone)]
pub struct IndelModel {
    // Based what was in the original NEAT
    insertion_probability: f64,
    ins_lengths: Vec<u32>,
    ins_dist: DiscreteDistribution,
    del_lengths: Vec<u32>,
    del_dist: DiscreteDistribution,
}

impl IndelModel {
    pub fn default() -> Self {
        // Default Indel model from the original NEAT, scaled by the lowest value and rounded.
        let insertion_probability = 0.6;
        let ins_lengths: Vec<u32> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ins_dist: DiscreteDistribution = DiscreteDistribution::new(
            &vec![1.0, 1.0, 2.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        );
        let del_lengths: Vec<u32> = vec![1, 2, 3, 4, 5];
        let del_dist: DiscreteDistribution = DiscreteDistribution::new(
            &vec![3.0, 2.0, 2.0, 2.0, 1.0]
        );

        IndelModel {
            insertion_probability,
            ins_lengths,
            ins_dist,
            del_lengths,
            del_dist,
        }
    }

    pub fn from(ins_prob: f64, ins_lengths: Vec<u32>, ins_weights: Vec<f64>, del_lengths: Vec<u32>, del_weights: Vec<f64>) -> Self {
        let ins_dist = DiscreteDistribution::new(&ins_weights);
        let del_dist = DiscreteDistribution::new(&del_weights);
        IndelModel {
            insertion_probability: ins_prob,
            ins_lengths,
            ins_dist,
            del_lengths,
            del_dist,
        }
    }

    pub fn is_insertion(&self, rng: &mut NeatRng) -> bool {
        rng.random() < self.insertion_probability
    }

    pub fn new_insert_length(&self, mut rng: &mut NeatRng) -> u32 {
        // This function uses the insertion weights to choose an insertion length from the list.
        self.ins_lengths[self.ins_dist.sample(&mut rng)]
    }

    pub fn new_delete_length(&self, mut rng: &mut NeatRng) -> u32 {
        // This function is the same as above, but uses the deletion lengths instead.
        self.del_lengths[self.del_dist.sample(&mut rng)]
    }
}

pub fn generate_random_insertion(length: u32, rng: &mut NeatRng) -> Vec<u8> {
    // We could refine this with a nucleotide bias matrix. Maybe it would make a difference,
    // but probably not, since the presence of the insertion is more important than it's content,
    // for this use. If there was a call for it, maybe.
    let mut insertion_vec = Vec::new();
    for _ in 0..length {
        // since our range is restricted, we are covered
        insertion_vec.push(rng.range_i64(0, 4).unwrap() as u8)
    }
    insertion_vec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_indel_model() {
        let indel_model = IndelModel {
            insertion_probability: 0.75,
            ins_lengths: vec![1, 2],
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0]),
            del_lengths: vec![3, 4],
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0]),
        };
        let mut rng = NeatRng::new_from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let length = indel_model.new_insert_length(&mut rng);
        if length > 0 {
            assert!((length == 1) || (length == 2));
        }
    }

    #[test]
    fn test_generate_insertion () {
        let mut rng = NeatRng::new_from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);

        let insertion = generate_random_insertion(10, &mut rng);
        assert_eq!(insertion.len(), 10);
        assert!(!insertion.contains(&4));

        let insertion2 = generate_random_insertion(12, &mut rng);
        assert_eq!(insertion2.len(), 12);
        assert!(!insertion2.contains(&4));

        println!("{:?} and {:?}", insertion, insertion2);
    }
}
