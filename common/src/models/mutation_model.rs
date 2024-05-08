use rand::distributions::{WeightedIndex, Distribution};
use rand::Rng;
use std::borrow::BorrowMut;
use rand_chacha::ChaCha20Rng;

use structs::nucleotides::Nuc;
use structs::transition_matrix::TransitionMatrix;
use structs::variants::{Variant, VariantType};
use models::indel_model::{IndelModel, generate_random_insertion};
use models::quality_scores::QualityScoreModel;
use models::sequencing_error_model::SequencingErrorModel;
use models::snp_model::SnpModel;
use structs::fasta_map::{SequenceBlock, insert_variant};

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
    pub variant_weights: [(VariantType, usize); 3],
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
}
impl MutationModel {
    pub fn new() -> Self {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. To make it easier to sample, we
        // will use weights. This would have been 0.95 and 0.05, combine with an insertion
        // frequency (relative to their being an indel variant) was 0.6.
        // This gives us a probability of insertion of 0.6 * 0.05 = 0.03 * 100 = 3
        //                 probability of deletion of  0.4 * 0.05 = 0.02 * 100 = 2
        //                 probability of SNP of       0.95              * 100 = 95
        // which works out to roughly 1 indel every 20 mutations (as originally envisioned).
        let variant_weights = [
            (VariantType::SNP, 95), (VariantType::Insertion, 3), (VariantType::Deletion, 2)
        ];
        let statistical_models = StatisticalModels::new();

        MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_weights,
            statistical_models,
        }
    }

    fn generate_genotype(&mut self, ploidy: usize, mut rng: ChaCha20Rng) -> Vec<u8> {
        // "Homozygous" is ambiguous for polyploid organisms, so we'll just take "heterozygous" to
        // mean roughly half the reads will have the variant, to keep it simple
        let is_homozygous = rng.borrow_mut().gen_bool(self.homozygous_frequency);
        if is_homozygous {
            vec![1; ploidy]
        } else {
            let num_ploids = ploidy/2;
            [vec![1; num_ploids], vec![0; ploidy-num_ploids]].concat()
        }
    }

    pub fn generate_mutation(
        &mut self,
        reference_sequence: &Vec<Nuc>,
        fasta_blocks: &Vec<SequenceBlock>,
        variant_location: usize,
        ploidy: usize,
        mut rng: ChaCha20Rng,
    ) -> Vec<SequenceBlock> {
        // Select a genotype for the variant
        let genotype= self.generate_genotype(ploidy, rng.clone());
        // Select a type of mutation.
        let dist = WeightedIndex::new(
            self.variant_weights
                .iter()
                .map(|item| item.1)
        ).unwrap();
        let variant_type = self.variant_weights[dist.sample(rng.borrow_mut())].0;
        // todo figure out which block to mutate
        //   figure out the mutation to add
        //   Add a function in FastaMap to add a block to a vector
        let (reference, alternate) = match variant_type {
            VariantType::SNP => {
                let trinuc_reference = &reference_sequence.get(
                    variant_location-1..=variant_location+1
                ).unwrap();

                let alternate_base = self.statistical_models.snp_model.generate_snp(
                    &trinuc_reference,
                    &mut rng,
                );

                let reference = vec![trinuc_reference[1]];
                let alternate = vec![alternate_base];

                (reference, alternate)
            },
            VariantType::Insertion => {
                let length = self.statistical_models.indel_model.new_insert_length(&mut rng);
                let insertion_vec = generate_random_insertion(length, rng.borrow_mut());
                let reference = vec![
                    reference_sequence.get(variant_location)
                        .unwrap()
                        .clone()
                ];
                let alternate = [reference.clone(), insertion_vec.clone()].concat();
                (reference, alternate)
            },
            VariantType::Deletion => {
                // todo Deletions are a bitch. I need to think about them.
                let length = self.statistical_models.indel_model.new_delete_length(&mut rng);
                // +1 is so that we grab a base for the reference in the VCF. This is similar
                // to how we appended bases to the reference in the insertion model.
                let reference: Vec<Nuc> = reference_sequence
                    .get(variant_location..(variant_location + length + 1))
                    .unwrap()
                    .iter()
                    .map(|item| *item)
                    .collect();
                let alternate: Vec<Nuc> = vec![reference[0].clone()];

                (reference, alternate)
            },
        };
        let variant_to_insert = Variant::new(
            variant_type, &reference, &alternate, genotype
        );
        let return_blocks: Vec<SequenceBlock> = insert_variant(
            fasta_blocks, variant_to_insert, variant_location
        );
        return_blocks
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
    use rand_core::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        TransitionMatrix::from(vec![a_weights, c_weights, g_weights, t_weights]);
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
