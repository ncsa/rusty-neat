use std::collections::HashMap;
use itertools::Itertools;
use simple_rng::{DiscreteDistribution, NeatRng};
use common::models::mutation_model::MutationModel;
use common::structs::variants::{Variant, VariantType};
use common::structs::fasta_map::SequenceBlock;
use common::structs::fasta_map::RegionType;
use common::structs::variants::VariantType::SNP;

pub fn generate_variants(
    sequence_block: SequenceBlock,
    region_weights: Vec<f64>,
    mutation_model: &mut MutationModel,
    num_mutations: usize,
    ploidy: usize,
    mut rng: NeatRng,
) -> Result<Option<HashMap<usize, Variant>>, Err>{
    // Variants for this contig, organized by location. Each location may have only one variant.
    // And we may need to implement a further check that there aren't other types of overlaps.
    // This will simply overwrite any existing value, so we just capture the last variant at each
    // location.
    // We take as input the number of mutations to add to this block. The number to add will be
    // calculated on the contig level.
    // Region weights are calculated from the contig map (they essentially just mask N-regions for
    // now, but we'll add more machine bias and trinuc bias nuance in future releases)
    if region_weights.len() != sequence_block.len() {
        return Err("Region weights did not match block length")
    }
    if num_mutations <=0 {
        return Ok(None)
    }
    let mut block_variants: HashMap<usize, Variant> = HashMap::new();
    let i = 0;
    let dist = DiscreteDistribution::new(&region_weights);
    while i < num_mutations {
        let location: usize = dist.sample(&mut rng);
        block_variants.insert(location, mutation_model.generate_mutation(
            &sequence_block.retrieve_all().unwrap(),
            location,
            ploidy,
            &mut rng,
        ));

    }
    Ok(Some(block_variants))
}