use std::collections::HashMap;
use rand::distributions::{WeightedIndex, Distribution};
use rand::Rng;

use structs::nucleotides::Nuc;
use structs::transition_matrix::TransitionMatrix;
use structs::variants::{Variant, VariantType};
use models::indel_model::IndelModel;
use models::quality_scores::QualityScoreModel;
use models::sequencing_error_model::SequencingErrorModel;
use models::snp_model::SnpModel;
use neat_rng::NeatRng;

#[derive(Clone)]
pub struct MutationModel {
    // This is the model for mutations, the same construct used by the python version, basically.
    // The idea is to store all the data we need to create mutations in one construct.
    //
    // mutation rate is the average rate of mutation for the dataset. It is used to calculate how
    // many variants to add to the dataset.
    pub mutation_rate: f64,
    // Homozygous frequency is the fraction of mutations that are homozygous, meaning the alternate
    // was inherited from both parents (in the case of humans). The definition of "homozygous"
    // is ambiguous in polyploid organisms. We'll take it to mean "on all ploids"
    pub homozygous_frequency: f64,
    // If a variant occurs, this is how frequently each type occurs. There are only 2 types at the
    // moment, but this will get expanded out in time.
    pub variant_weights: [(VariantType, usize); 2],
    // These hold all the statistical model data we need to apply the mutations with this model
    pub statistical_models: StatisticalModels,
    // This is the rng for the run, which may not be needed depending on how this struct is
    // being used.
    pub rng: Option<NeatRng>,
}
impl MutationModel {
    pub fn new(rng: NeatRng) -> Self {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. To make it easier to sample, we
        // will use weights. This would have been 0.95 and 0.05, so we divided both by 0.05 to get
        // 19 and 1 as our weights. This means there is roughly 1 indel every 20 mutations.
        let variant_weights = [(VariantType::SNP, 19), (VariantType::Indel, 1)];
        let statistical_models = StatisticalModels::new();

        MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_weights,
            statistical_models,
            rng: Some(rng),
        }
    }

    #[allow(dead_code)]
    // Todo I think I need this to generate the mutation model from data
    pub fn default_no_rng() -> Self {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. To make it easier to sample, we
        // will use weights. This would have been 0.95 and 0.05, so we divided both by 0.05 to get
        // 19 and 1 as our weights. This means there is roughly 1 indel every 20 mutations.
        let variant_weights = [(VariantType::SNP, 19), (VariantType::Indel, 1)];
        let statistical_models = StatisticalModels::new();

        MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_weights,
            statistical_models,
            rng: None,
        }
    }

    pub fn add_rng(
        &mut self,
        rng: &mut NeatRng
    ) {
        // Cloning the rng is my only choice here. I'm curious to see how the results turn out.
        self.rng = Some(rng.clone());
    }

    pub fn get_mut_rng(&self) -> &mut NeatRng { &mut self.rng.clone().unwrap() }

    fn generate_genotype(&self, ploidy: usize) -> (Vec<u8>, bool) {
        // "Homozygous" is ambiguous for polyploid organisms, so we'll just take "heterozygous" to
        // mean roughly half the reads will have the variant, to keep it simple
        let is_homozygous = self.get_mut_rng().gen_bool(self.homozygous_frequency);
        if is_homozygous {
            (vec![1; ploidy], true)
        } else {
            let num_ploids = ploidy/2;
            ([vec![1; num_ploids], vec![0; ploidy-num_ploids]].concat(), false)
        }
    }

    pub fn generate_mutation(
        &self,
        reference_sequence: &Vec<Nuc>,
        variant_location: usize,
        ploidy: usize,
    ) -> Variant {
        // Select a genotype for the variant
        let (genotype, is_homozygous) = self.generate_genotype(ploidy);
        // Select a type of mutation.
        let mut dist = WeightedIndex::new(
            self.variant_weights
                .iter()
                .map(|item| item.0)
        ).unwrap();
        let variant_type = self.variant_weights[dist.sample(self.get_mut_rng())].0;
        let (reference, alternate) = match variant_type {
            VariantType::SNP => {
                let trinuc_reference = reference_sequence.get(
                    variant_location-1..=variant_location+1
                ).unwrap();

                let alternate_base = self.statistical_models.snp_model.generate_snp(
                    &trinuc_reference,
                    self.get_mut_rng(),
                );

                let reference = vec![trinuc_reference[1]];
                let alternate = vec![alternate_base];

                (reference, alternate)
            },
            VariantType::Indel => {
                let length = self.statistical_models.indel_model
                    .generate_new_indel_length(self.get_mut_rng());
                if length > 0 { // insertion
                    let insertion_vec = self.statistical_models.indel_model
                        .generate_random_insertion(length.abs() as usize, self.get_mut_rng())
                        .clone();
                    let reference = vec![
                        reference_sequence.get(variant_location)
                            .unwrap()
                            .clone()
                    ];
                    let alternate = [reference.clone(), insertion_vec.clone()].concat();
                    (reference, alternate)
                } else { // deletion
                    // plus one because the first base is not deleted, but serves as a reference
                    // in the vcf
                    let reference: Vec<Nuc> = vec![
                        *(reference_sequence.get(
                        &(variant_location..(variant_location + length.abs() as usize + 1))
                    ).unwrap())];
                    let alternate = vec![reference[0].clone()];

                    (reference, alternate)
                }
            },
        };
        Variant::new(variant_type, &reference, &alternate, genotype, is_homozygous)
    }
}

#[derive(Clone)]
struct StatisticalModels {
    // This struct links the Mutation model to the other statistical models for modeling different
    // variant types. If new variant types are added, then we will need to expand this struct to
    // include them.
    transition_matrix: TransitionMatrix,
    quality_score_model: QualityScoreModel,
    sequencing_error_model: SequencingErrorModel,
    indel_model: IndelModel,
    snp_model: SnpModel,
}

impl StatisticalModels {
    pub fn new() -> Self {
        // use the default transition matrix, snp model and indel model
        let transition_matrix = TransitionMatrix::new();
        let snp_model = SnpModel::new();
        let indel_model = IndelModel::new();
        let quality_score_model = QualityScoreModel::new();
        let sequencing_error_model = SequencingErrorModel::new();

        StatisticalModels {
            transition_matrix,
            quality_score_model,
            snp_model,
            indel_model,
            sequencing_error_model,
        }
    }

    pub fn from(
        transition_matrix: TransitionMatrix,
        quality_score_model: QualityScoreModel,
        snp_model: SnpModel,
        indel_model: IndelModel,
        sequencing_error_model: SequencingErrorModel,
    ) -> Self {
        StatisticalModels {
            transition_matrix,
            quality_score_model,
            snp_model,
            indel_model,
            sequencing_error_model,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nuc::*;
    use rand_core::SeedableRng;
    use neat_rng::NeatRng;

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let rng = NeatRng::seed_from_u64(0);
        let matrix = TransitionMatrix::from(vec![a_weights, c_weights, g_weights, t_weights]);
        // It actually mutates the base
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![A],
        //     &mut rng
        // ), vec![A]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![C],
        //     &mut rng
        // ), vec![C]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![G],
        //     &mut rng
        // ), vec![G]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![T],
        //     &mut rng
        // ), vec![T]);
        // // It gives back N when you give it N
        // assert_eq!(test_model.mutate(
        //     VariantType::SNP, 0, &vec![N], &mut rng), vec![N]
        // );
        // // test alternate mutation methods
        // let mutation = test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     1,
        //     &vec![A, C, G],
        //     &mut rng
        // );
        // assert_eq!(mutation[0], A);
        // assert_eq!(mutation[2], G);
        // assert_ne!(mutation[3], C);
        // let mutation = test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::Indel,
        //     0,
        //     &vec![A, C, C, G, T, T, A, C, G],
        //     &mut rng
        // );
        // // todo something with this mutation variable
    }
}
