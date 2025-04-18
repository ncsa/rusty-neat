// Throughout this program, we are standardizing the use of a u8 representation of the nucleotides
//     A = 0
//     C = 1
//     G = 2
//     T = 3
//     N = 4
// This is intended to make it easier to store them. We thought about using the u8 representation
// of the character as built into Rust, but we'd then have to figure out the translations and keep
// track of extra numbers. So this is intended to simplify everything
use simple_rng::{DiscreteDistribution, Rng};

pub fn base_to_u8(char_of_interest: char) -> u8 {
    // This defines the relationship between the 4 possible nucleotides in DNA and
    // a simple u8 numbering system. Everything that isn't a recognized base is a 4.
    // Note that NEAT ignores soft masking.
    return match char_of_interest {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}

pub fn u8_to_base(nuc_num: u8) -> char {
    // Canonical conversion from base u8 representation back into the character.
    // We're returning a string instead of a char to facilitate. No attempt to preserve or display
    // any soft masking.
    return match nuc_num {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

#[derive(Debug, Clone)]
pub struct NucModel {
    // Simple mutation model. The letter in the model represents
    // the base we are mutating and the vector are the weights for mutating
    // to another base (in the same a, c, g, t order)
    //
    // By definition, the model is a 4x4 matrix with zeros along the diagonal
    // because, e.g., A can't "mutate" to A.
    a: Vec<u32>,
    c: Vec<u32>,
    g: Vec<u32>,
    t: Vec<u32>,
}

impl NucModel {
    pub fn new() -> Self {
        // Default mutation model based on the original from NEAT 2.0
        Self {
            a: vec![0, 17, 69, 14],
            c: vec![16, 0, 17, 67],
            g: vec![67, 17, 0, 16],
            t: vec![14, 69, 16, 0],
        }
    }

    #[allow(dead_code)]
    // todo, once we have numbers we can implement this.
    pub fn from(weights: Vec<Vec<u32>>) -> Self {
        // Supply a vector of 4 vectors that define the mutation chance
        // from the given base to the other 4 bases.

        // First some safety checks. This should be a 4x4 matrix defining mutation from
        // ACGT (top -> down) to ACGT (left -> right)
        if weights.len() != 4 {
            panic!("Weights supplied to NucModel is wrong size");
        }
        for weight_vec in &weights {
            if weight_vec.len() != 4 {
                panic!("Weights supplied to NucModel is wrong size");
            }
        }
        Self {
            a: weights[0].clone(),
            c: weights[1].clone(),
            g: weights[2].clone(),
            t: weights[3].clone(),
        }
    }

    pub fn choose_new_nuc(&self, base: u8, rng: &mut Rng) -> u8 {

        // the canonical choices for DNA, as defined above
        let choices: [u8; 4] = [0, 1, 2, 3];
        // Pick the weights list for the base that was input
        let weights: Vec<u32> = match base {
            0 => self.a.clone(),
            1 => self.c.clone(),
            2 => self.g.clone(),
            3 => self.t.clone(),
            // anything else we return the N value of 4
            _ => { return 4; },
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = DiscreteDistribution::new(&weights, false);
        choices[dist.sample(rng)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nuc_model_build() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];

        let model = NucModel {
            a: a_weights.clone(),
            c: c_weights.clone(),
            g: g_weights.clone(),
            t: t_weights.clone(),
        };

        let str = format!("{:?}", model);
        let str_repr = String::from("NucModel { a: [0, 20, 1, 20], c: [20, 0, 1, 1], g: [1, 1, 0, 20], t: [20, 1, 20, 0] }");
        assert_eq!(str, str_repr);
        assert_eq!(model.a, a_weights);
    }

    #[test]
    fn test_nuc_model_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let test_model = NucModel::from(
            vec![a_weights, c_weights, g_weights, t_weights]
        );
        // It actually mutates the base
        assert_ne!(test_model.choose_new_nuc(0, &mut rng), 0);
        assert_ne!(test_model.choose_new_nuc(1, &mut rng), 1);
        assert_ne!(test_model.choose_new_nuc(2, &mut rng), 2);
        assert_ne!(test_model.choose_new_nuc(3, &mut rng), 3);
        // It gives back N when you give it N
        assert_eq!(test_model.choose_new_nuc(4, &mut rng), 4);
    }

    #[test]
    #[should_panic]
    fn test_nuc_model_too_many_vecs() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let u_weights = vec![20, 1, 20, 0];
        NucModel::from(vec![a_weights, c_weights, g_weights, t_weights, u_weights]);
    }

    #[test]
    #[should_panic]
    fn test_nuc_model_too_many_bases() {
        let a_weights = vec![0, 20, 1, 20, 1];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        NucModel::from(vec![a_weights, c_weights, g_weights, t_weights]);
    }
}