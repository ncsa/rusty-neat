use log::error;
use simple_rng::NeatRng;
use crate::common::{
    models::mutation_model::MutationModel,
    structs::{
        variants::Variant,
        distributions::DiscreteDistribution,
        fasta_map::SequenceBlock,
    }
};
use crate::errors::GenerateReadsErrors;

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
    if region_weights.len() != sequence_block.get_len() {
        error!("Region weights did not match block length");
        return Err(GenerateReadsErrors::GenerateVariantsError)
    }
    // The keys are the position of each variant added.
    let mut block_variants: Vec<Variant> = Vec::new();
    let sequence_dist = DiscreteDistribution::new(
        &region_weights,
        &(0..sequence_block.get_len()).collect()
    )?;
    for _ in 0..num_mutations {
        let location: usize = sequence_dist.sample(rng.random()?)?;
        block_variants.push(
            mutation_model.generate_mutation(
                &sequence_block.get_seq_clone().unwrap(),
                location,
                ploidy,
                rng,
            )?
        );
    }
    Ok(Some(block_variants))
}