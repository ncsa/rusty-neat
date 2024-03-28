use utils::nucleotides::Nuc;
use utils::neat_rng::NeatRng;
use log::debug;
use rand::distributions::WeightedIndex;
use rand::prelude::Distribution;
use utils::transition_matrix::TransitionMatrix;

pub struct SequencingErrorModel {
    transition_matrix: TransitionMatrix,
}

impl SequencingErrorModel {
    fn generate_snp_error(
        &self,
        base: Nuc,
        rng: &mut NeatRng
    ) -> Nuc {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let weights: &Vec<usize> = match base {
            Nuc::A => &self.transition_matrix.a_weights,
            Nuc::C => &self.transition_matrix.c_weights,
            Nuc::G => &self.transition_matrix.g_weights,
            Nuc::T => &self.transition_matrix.t_weights,
            // return the N value for N with no further computation.
            Nuc::N => { return Nuc::N; },
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = WeightedIndex::new(weights).unwrap();
        match dist.sample(rng) {
            0 => Nuc::A,
            1 => Nuc::C,
            2 => Nuc::G,
            3 => Nuc::T,
            _ => Nuc::N,
        }
    }
}