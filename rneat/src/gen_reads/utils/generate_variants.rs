use log::{error, debug};
use simple_rng::NeatRng;
use crate::{
    common::{
        models::mutation_model::MutationModel,
        structs::{
            distributions::DiscreteDistribution, 
            fasta_map::SequenceBlock, 
            variants::Variant
        }
    }, 
    gen_reads::errors::GenerateReadsErrors
};

pub fn generate_variants(
    sequence_block: &SequenceBlock,
    region_weights: &Vec<f64>,
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsErrors>{
    // Variants for this contig, organized by location. Each location may have only one variant.
    // And we may need to implement a further check that there aren't other types of overlaps.
    // This will simply overwrite any existing value, so we just capture the last variant at each
    // location.
    // We take as input the number of mutations to add to this block. The number to add will be
    // calculated on the contig level.
    // Region weights are calculated from the contig map (they essentially just mask N-regions for
    // now, but we'll add more machine bias and trinuc bias nuance in future releases)
    let block_len = sequence_block.get_len();
    if region_weights.len() != block_len {
        error!("Region weights did not match block length");
        return Err(GenerateReadsErrors::GenerateVariantsError)
    }
    // The keys are the position of each variant added.
    let mut block_variants: Vec<Variant> = Vec::new();
    let sequence_dist = DiscreteDistribution::new(
        &region_weights[0..block_len].to_vec(),
        &(0..block_len).collect()
    )?;
    for _ in 0..num_mutations {
        let location: usize = sequence_dist.sample(rng.random()?)?;
        debug!("location: {}", location);
        block_variants.push(
            mutation_model.generate_mutation(
                &sequence_block.get_seq_clone()?,
                location,
                ploidy,
                rng,
            )?
        );
    }
    Ok(Some(block_variants))
}


#[cfg(test)]
mod tests {
    use common::structs::fasta_map::RegionType::NonNRegion;
    use common::structs::fasta_map::SequenceMap;
    use common::structs::nucleotides::Nucleotide::*;
    use super::*;

    fn setup() -> (SequenceBlock, NeatRng) {
        let sequence_map = SequenceMap{
            region_type: NonNRegion,
            start: 0,
            end: 10
        };
        let sequence_block = SequenceBlock{
            contig: "chr1".to_string(),
            ref_start: 10,
            ref_end: 20,
            sequence: vec![A, A, A, G, T, A, G, C, C, T],
            sequence_map: vec![sequence_map]
        };
        let rng = NeatRng::new_from_seed(
            &vec![
                "Goodbye".to_string(),
                "Cruel".to_string(),
                "World".to_string(),
            ]
        ).unwrap();
        (sequence_block, rng)
    }
    // Now the tests
    #[test]
    fn test_generate_variants_success() {
        let (sequence_block, mut rng) = setup();
        let region_weights = vec![1.0; 10];
        let mutation_model = MutationModel::default().unwrap();
        let num_mutations = 3;
        let ploidy = 2;

        let result = generate_variants(
            &sequence_block,
            &region_weights,
            &mutation_model,
            num_mutations,
            ploidy,
            &mut rng,
        );

        assert!(result.is_ok());
        let variants = result.unwrap().unwrap();
        assert_eq!(variants.len(), num_mutations);
        for variant in variants {
            assert!(variant.location < sequence_block.get_len());
        }
    }

    #[test]
    fn test_generate_variants_region_weights_length_mismatch() {
        let (sequence_block, mut rng) = setup();
        let region_weights = vec![1.0; 7]; // length mismatch
        let mutation_model = MutationModel::default().unwrap();
        let num_mutations = 1;
        let ploidy = 2;

        let result = generate_variants(
            &sequence_block,
            &region_weights,
            &mutation_model,
            num_mutations,
            ploidy,
            &mut rng,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_generate_variants_zero_mutations() {
        let (sequence_block, mut rng) = setup();
        let region_weights = vec![1.0; 10];
        let mutation_model = MutationModel::default().unwrap();
        let num_mutations = 0;
        let ploidy = 2;

        let result = generate_variants(
            &sequence_block,
            &region_weights,
            &mutation_model,
            num_mutations,
            ploidy,
            &mut rng,
        );

        assert!(result.is_ok());
        let variants = result.unwrap();
        assert!(variants.is_some());
        assert_eq!(variants.unwrap().len(), 0);
    }
}