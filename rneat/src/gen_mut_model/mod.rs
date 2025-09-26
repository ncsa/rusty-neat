/// This file will create a mutation model based on input data. 
///
pub mod errors;
pub mod utils;

use log::*;
use std::{
    collections::HashMap, 
    path::PathBuf
};
use common::{
    file_tools::fasta_reader::read_fasta, 
    structs::variants::Variant,
};
use crate::gen_mut_model::{
    errors::GenMutationModelError, 
    utils::config::RunConfiguration,
    utils::runner::runner,
};

pub fn main(config_file: &PathBuf) -> Result<(), GenMutationModelError> {
    let run_config = RunConfiguration::from(config_file)?;
    // filter the variants with the bed_table, if applicable
    let mut filtered_mutations: HashMap<String, Vec<Variant>> = HashMap::new();
    if !run_config.bed_table.is_empty() {
        for chrom in run_config.mutations.keys() {
            if run_config.bed_table.contains_key(chrom) {
                let search_region = &run_config.bed_table[chrom];
                let contig_mutations = &run_config.mutations[chrom];
                for variant in contig_mutations {
                    for record in search_region {
                        if record.contains(chrom, variant.location) {
                            if !filtered_mutations.contains_key(chrom) {
                                filtered_mutations.insert(chrom.clone(), Vec::new());
                            }
                            filtered_mutations
                                .get_mut(chrom)
                                .unwrap()
                                .push(variant.clone());
                            break;
                        }
                    }
                }
            } else {
                // Throw out all mutations in the vcf with this chrom
                debug!("Filtered out contig: {chrom}")
            }
        }
    } else {
        filtered_mutations = run_config.mutations
    }
    // this will store our fasta file pieces
    let temp_dir = tempfile::tempdir().unwrap();
    // We want to read this fasta in as raw as possible, with no overlap. 
    // We don't need to generate reads, so we only want the sequence.
    let fasta_map = read_fasta(
        &run_config.reference,
        None,
        0,
        &temp_dir,
        None,
    ).expect("Error reading fasta file");

    let result = runner(
        fasta_map,
        filtered_mutations,
        run_config.bed_table,
        &run_config.output_file,
    ).unwrap();

    Ok(())
}