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
    use super::*;
    use crate::common::{
        models::mutation_model::MutationModel,
        structs::fasta_map::{RegionType, SequenceBlock, SequenceMap},
    };
    use common::structs::nucleotides::Nucleotide;
    use simple_rng::NeatRng;

    fn make_block(len: usize) -> SequenceBlock {
        let sequence: Vec<Nucleotide> = (0..len)
            .map(|i| match i % 4 {
                0 => Nucleotide::A,
                1 => Nucleotide::C,
                2 => Nucleotide::G,
                _ => Nucleotide::T,
            })
            .collect();
        SequenceBlock {
            contig: "test".to_string(),
            ref_start: 0,
            ref_end: len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, len)],
        }
    }

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap()
    }

    #[test]
    fn test_generate_variants_count() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &weights, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn test_generate_variants_positions_in_range() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let variants = generate_variants(&block, &weights, &model, 10, 2, &mut rng)
            .unwrap()
            .unwrap();
        for v in variants {
            assert!(v.location < 200, "Variant location {} out of range", v.location);
        }
    }

    #[test]
    fn test_generate_variants_weight_mismatch_errors() {
        let block = make_block(200);
        let weights = vec![1.0f64; 100]; // wrong length
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        assert!(generate_variants(&block, &weights, &model, 5, 2, &mut rng).is_err());
    }

    #[test]
    fn test_generate_variants_deterministic() {
        let block = make_block(200);
        let weights = vec![1.0f64; block.get_len()];
        let model = MutationModel::default().unwrap();
        let result1 = generate_variants(&block, &weights, &model, 5, 2, &mut make_rng())
            .unwrap().unwrap();
        let result2 = generate_variants(&block, &weights, &model, 5, 2, &mut make_rng())
            .unwrap().unwrap();
        for (v1, v2) in result1.iter().zip(result2.iter()) {
            assert_eq!(v1.location, v2.location);
            assert_eq!(v1.variant_type, v2.variant_type);
        }
    }
}