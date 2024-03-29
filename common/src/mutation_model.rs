use std::collections::HashMap;
use rand::distributions::Distribution;
use crate::neat_rng::NeatRng;
use crate::nucleotides::Nuc;
use crate::vars::variants::VariantType;
use crate::transition_matrix::TransitionMatrix;
use crate::vars::variants::*;

pub struct MutationModel {
    // This is the model for mutations, the same construct used by the python version, basically.
    //
    // mutation rate is the average rate of mutation for the dataset. It is used to calculate how
    // many variants to add to the dataset.
    mutation_rate: f64,
    // Homozygous frequency is the fraction of mutations that are homozygous, meaning the alternate
    // was inherited from both parents (in the case of humans). The definition of "homozygous"
    // is ambiguous in polyploid organisms. We'll take it to mean "on all ploids"
    homozygous_frequency: f64,
    // If a variant occurs, this is the probability that it will be a SNP, the most common type
    // of variant. There are only 2 types at the moment, but this will get expanded out in time.
    variant_probs: HashMap<VariantType, usize>,
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
}
impl MutationModel {
    pub fn new() -> Self {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. To make it easier to sample, we
        // will use weights. This would have been 0.95 and 0.05, so we divided both by 0.05 to get
        // 19 and 1 as our weights.
        let variant_probs = HashMap::from([
            (VariantType::SNP, 19),
            (VariantType::Indel, 1),
        ]);

        let statistical_models = StatisticalModels::new();

        MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_probs,
            statistical_models,
        }
    }
    // We may need several ways to create the mutation model. Here is one for testing.
    pub fn from_transition_matrix(transition_matrix: TransitionMatrix) -> Self {
        let mut model = MutationModel::new();
        model.statistical_models.transition_matrix = transition_matrix;
        model
    }

    pub fn generate_mutation(
        &self,
        variant_type: VariantType,
        variant_location: usize,
        input_sequence: &Vec<Nuc>,
        rng: &mut NeatRng
    ) -> Vec<Nuc> {
        self.statistical_models.get_variant(variant_type, variant_location, input_sequence, rng)
    }
}

#[cfg(test)]
mod tests {
    use rand_core::SeedableRng;
    use lib::mutation_model::{MutationModel, TransitionMatrix};
    use utils::neat_rng::NeatRng;
    use lib::nucleotides::Nuc::*;
    use lib::variants::VariantType;

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let mut rng = NeatRng::seed_from_u64(0);
        let matrix = TransitionMatrix::from(
            vec![a_weights, c_weights, g_weights, t_weights]
        );
        let test_model = MutationModel::from_transition_matrix(matrix);
        // It actually mutates the base
        assert_ne!(test_model.generate_mutation(
            VariantType::SNP, 0, &vec![A], &mut rng), vec![A]
        );
        assert_ne!(test_model.generate_mutation(
            VariantType::SNP, 0, &vec![C], &mut rng), vec![C]
        );
        assert_ne!(test_model.generate_mutation(
            VariantType::SNP, 0, &vec![G], &mut rng), vec![G]
        );
        assert_ne!(test_model.generate_mutation(
            VariantType::SNP, 0, &vec![T], &mut rng), vec![T]
        );
        // It gives back N when you give it N
        assert_eq!(test_model.mutate(
            VariantType::SNP, 0, &vec![N], &mut rng), vec![N]
        );
        // test alternate mutation methods
        let mutation = test_model.generate_mutation(
            VariantType::SNP, 1, &vec![A, C, G], &mut rng
        );
        assert_eq!(mutation[0], A);
        assert_eq!(mutation[2], G);
        assert_ne!(mutation[3], C);
        let mutation = test_model.generate_mutation(
            VariantType::Indel,
            0,
            &vec![A, C, C, G, T, T, A, C, G],
            &mut rng
        );
        // todo something with this mutation variable

    }
}