use std::{collections::HashMap, path::PathBuf};
use itertools::Itertools;
use log::*;

use common::{
    models::snp_trinuc_model::TrinucFrame, 
    structs::{
        bed_record::BedRecord, 
        fasta_map::{FastaMap, RegionType, SequenceBlock}, 
        nucleotides::Nucleotide, 
        variants::{Variant, VariantType}
    }
};

use crate::gen_mut_model::errors::GenMutationModelError;

pub fn runner(
    fasta_map: Box<FastaMap>,
    filtered_mutations: HashMap<String, Vec<Variant>>,
    bed_table: HashMap<String, Vec<BedRecord>>,
    output_file: &PathBuf,
) -> Result<(), GenMutationModelError> {
    // This will generate a MutationModel and write it to the output file 
    // provided.
    // params:
    //   fasta_map: The database for the fasta file, with links back to the 
    //      original sequences
    //   filtered_mutations: The mutations from the original VCF, after any 
    //      appropriate bed filtering has been applied
    //   bed_table: The hash map representing thi input bed
    //   output_file: Where to write the output model
    //
    // This will store the count of trinucleotide combo A transitioning to combo B
    // This is the raw data for the trinuc model.
    let mut trinuc_transition_count: HashMap<(TrinucFrame, TrinucFrame), usize> = HashMap::new();
    let mut snp_count = 0;
    let mut snp_transition_count: HashMap<(Nucleotide, Nucleotide), usize> = HashMap::new();
    let mut insertion_count: HashMap<usize, usize> = HashMap::new();
    let mut deletion_count: HashMap<usize, usize> = HashMap::new();
    let mut homozygous_count = 0;
    let mut common_variants: Vec<(String, Variant, f64)> = Vec::new();
    let trinuc_count = count_trinculeotides (
        &fasta_map,
        bed_table,
    )?;
    if trinuc_count.is_empty() {
        error!("Trinucleotide counts were empty");
        return Err(
            GenMutationModelError::TrinucCountError(
                "Trinuc counts are empty. Unkown error".to_string()
            )
        )
    }
    let mut total_reflen = 0;
    for contig in &fasta_map.contigs{
        let mut non_n_region_count: usize = 0;
        let countable_regions = contig.get_non_n_regions()?;
        for region in countable_regions {
            // I think the coordinates are forward on all regions
            total_reflen += region.end - region.start;
        }
        let matching_variants = &filtered_mutations[&contig.name];
        if matching_variants.is_empty() {
            debug!("No variants found for {}", contig.name);
            continue
        }
        for variant in matching_variants {
            match variant.variant_type {
                VariantType::SNP => {
                    snp_count += 1;
                    // Retrieve the trinucleotide sequence from the sequence block
                    // If this is slow we may need to load the entire block and comb it.
                    let trinuc_to_analyze = contig.get_block_subseq(variant.location-1, variant.location+2)?;
                    // I'm curious to see what happens if the above is 4 bases long.
                    // does impl From implement some error check under the hood?
                    let ref_array: [&Nucleotide; 3] = trinuc_to_analyze.iter().collect_array().unwrap(); 
                    let ref_frame = TrinucFrame::from(ref_array);
                    let alt_array = [ref_array[0], &variant.alternate[0], ref_array[2]];
                    let alt_frame = TrinucFrame::from(alt_array);
                    if *ref_array[1] != variant.reference[0] {
                        return Err(GenMutationModelError::BaseMismatch(format!("{:?}", variant)))
                    }
                    *trinuc_transition_count.entry((ref_frame, alt_frame)).or_default() += 1;
                    *snp_transition_count.entry((variant.reference[0], variant.alternate[0])).or_default() += 1;
                },
                VariantType::Insertion => {
                    let variant_len = variant.alternate.len() - variant.reference.len();
                    *insertion_count.entry(variant_len).or_default() += 1;
                },
                VariantType::Deletion => {
                    let variant_len = variant.reference.len() - variant.alternate.len();
                    *deletion_count.entry(variant_len).or_default() += 1;
                },
                _ => debug!("Unknown variant type, skipping for this analyis."),
            }
            let var_pop_freq = find_caf(&variant.alternate);

        }
    }
    Ok(())
}

fn find_caf(alt: &Vec<Nucleotide>) -> Result<f64, GenMutationModelError> {
    // Python:
    // info_split = [a.split('=') for a in candidate_field.split(';')]
    // for item in info_split:
    //     if item[0].upper() == 'CAF':
    //         if ',' in item[1]:
    //             return float(item[1].split(',')[1])
    // return VCF_DEFAULT_POP_FREQ
    todo!()
}

fn count_trinculeotides(
    fasta_map: &FastaMap,
    bed_table: HashMap<String, Vec<BedRecord>>,
) -> Result<HashMap<TrinucFrame, usize>, GenMutationModelError> {
    // Counts the frequency of the various trinucleotide combinations in the dataset
    // fasta_map - Database for reading back the sequenced data.
    // bed_table - A table for the bed records

    // HashMap of trinuce to count, e.g., "AAA": 12
    let mut trinuc_count: HashMap<TrinucFrame, usize> = HashMap::new();
    if !bed_table.is_empty() {
        // There is actual bed input
        info!("Counting trinucleotide combinations in bed regions");
        for (chrom, bed_vec) in bed_table {
            let contig_block = fasta_map.retrieve_contig(chrom.clone())
                .expect(&format!("Error retrieving contig {}", chrom));
            let regions_of_interest = contig_block.get_non_n_regions()
                .expect("Error retrieving non-n regions!");
            for region in regions_of_interest {
                // iterate over the regions, grab trinucleotides, profit
                // The first trincu is one base in. The last is one base from the end
                for i in (region.start+1)..(region.end-1) {
                    let trinuc = contig_block.get_block_subseq(i-1, i+1)?;
                    let frame = TrinucFrame::from((trinuc[0], trinuc[1], trinuc[2]));
                    *trinuc_count.entry(frame.to_owned()).or_default() += 1;
                }
            }
        }
    } else {
        // No bed input
    }
    Ok(trinuc_count)
}