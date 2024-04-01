use crate::utils;
use common;

use std::thread;
use std::time::Duration;

use log::info;
use rand::prelude::SliceRandom;
use std::collections::{HashMap, HashSet};
use common::models::mutation_model::MutationModel;
use utils::config::RunConfiguration;
use utils::fasta_tools::{read_fasta, write_fasta};
use utils::fastq_tools::write_fastq;
use utils::make_reads::generate_reads;
use utils::mutate::mutate_fasta;
use common::neat_rng::NeatRng;
use common::structs::nucleotides::Nuc;
use common::structs::variants::Variant;
use utils::read_models::read_quality_score_model_json;
use utils::vcf_tools::write_vcf;
use crate::utils::generate_variants::generate_variants;

#[derive(Debug)]
pub enum RunNeatError {
    GeneralBullshit,
}

pub fn run_neat(config: Box<RunConfiguration>, rng: NeatRng) -> Result<(), RunNeatError> {
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir.display(), config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    // Trying to think of how to avoid copying this. We could read the file, get the order, but
    // skip the reading in for now, let the threads handle that. Or we could create an index of byte
    // positions of the sequences, like a FAI, but that's a different type of problem to solve.
    // Alternative to that, require an index file and then use that to find the read positions?
    // That might be the easiest way to go.
    let (fasta_map, mut fasta_order) = read_fasta(&config.reference)
        .unwrap();

    let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
    for contig in &fasta_order {
        fasta_lengths.insert(contig.to_string(), fasta_map.get(&contig).unwrap().len())
    }

    // Load models that will be used for the runs.
    // For now, we will use the one supplied, pulled directly from NEAT2.0's original model.
    let default_quality_score_model_file = "model_data/neat_quality_score_model.json";
    let quality_score_model = read_quality_score_model_json(
        default_quality_score_model_file
    );

    // Todo load all models and set up the run.
    // load mutation model
    let mutation_model = MutationModel::new(rng);

    // Mutating the reference and recording the variant locations.
    info!("Generating variants");
    // I think looping like this will allow us to multithread easier
    let mut all_variants: Box<HashMap<String, HashMap<usize, Variant>>> = Box::new(HashMap::new());
    let mut all_reads: Box<HashMap<String, (usize, usize)>>;
    let handle = thread::spawn(move || {
        for contig in &fasta_order {
            info!("Generating variants for {chromosome}");
            let contig_sequence = fasta_map.get(contig).unwrap();
            let contig_variants = generate_variants(
                contig_sequence, &mutation_model, config.ploidy,
            );
            info!("Finished generating variants for {chromosome}");
            all_variants.insert(contig.to_owned(), contig_variants);

            let contig_reads = generate_reads(
                contig_sequence.len(),
                &contig_variants,
                config.read_len,
                config.coverage,
                config.paired_ended,
                config.fragment_mean,
                config.fragment_st_dev,
                mutation_model.get_mut_rng(),
            );

            // This all adds up to if this run is fasta only. If it is, then we don't
            // need to generate reads
            if config.produce_fasta &&
                !config.produce_fastq &&
                !config.produce_vcf &&
                !config.produce_bam {
                all_reads.insert(contig.to_owned(), *contig_reads.unwrap()[0])
            }
        }
    });

    handle.join().unwrap();

    if config.produce_fasta {
        write_fasta(
            &fasta_map,
            &all_variants,
            &fasta_order,
            config.overwrite_output,
            &output_file,
        ).expect("Error writing fasta file!")
    }

    // If we produce any more files, we will need to know the genotype. It will be used in the
    // output vcf, bam, and in the fastq.

    if config.produce_fastq {
        write_fastq(
            &fasta_map,
            &all_reads,
            config.overwrite_output,
            &output_file,
            config.paired_ended,
            quality_score_model,
        ).expect("Error writing fastq file(s)!")
    }

    // todo for the vcf and bam, we need to know things like genotype
    if config.produce_vcf {
        write_vcf(
            &all_variants,
            &fasta_order,
            &fasta_lengths,
            &config.reference,
            config.overwrite_output,
            &output_file,
        ).expect("Error writing vcf file!")
    }

    // let (mutated_map, variant_locations) =
    //     mutate_fasta(&fasta_map, config.minimum_mutations, &mut rng);
    //
    // if config.produce_fasta {
    //     info!("Outputting fasta file");
    //     write_fasta(
    //         &mutated_map,
    //         &fasta_order,
    //         config.overwrite_output,
    //         &output_file,
    //     )
    //     .unwrap();
    // }
    //
    // if config.produce_vcf {
    //     info!("Writing vcf file");
    //     write_vcf(
    //         &variant_locations,
    //         &fasta_order,
    //         config.ploidy,
    //         &config.reference,
    //         config.overwrite_output,
    //         &output_file,
    //         &mut rng,
    //     )
    //     .unwrap();
    // }
    //
    // let mut read_sets: HashSet<Vec<Nuc>> = HashSet::new();
    // info!("Generating reads");
    // for (_name, sequence) in mutated_map.iter() {
    //     // defined as a set of read sequences that should cover
    //     // the mutated sequence `coverage` number of times
    //     let data_set = generate_reads(
    //         sequence,
    //         &config.read_len,
    //         &config.coverage,
    //         config.paired_ended,
    //         config.fragment_mean,
    //         config.fragment_st_dev,
    //         &mut rng,
    //     )
    //     .unwrap();
    //
    //     read_sets.extend(*data_set);
    // }
    //
    // if config.produce_fastq {
    //     info!("Shuffling output fastq data");
    //     let mut outsets: Box<Vec<&Vec<Nuc>>> = Box::new(read_sets.iter().collect());
    //     outsets.shuffle(&mut rng);
    //
    //     info!("Writing fastq");
    //     write_fastq(
    //         &output_file,
    //         config.overwrite_output,
    //         config.paired_ended,
    //         *outsets,
    //         quality_score_model,
    //         rng,
    //     )
    //     .unwrap();
    //     info!("Processing complete")
    // }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::SeedableRng;
    use std::fs;
    use std::path::PathBuf;
    use utils::config::ConfigBuilder;

    #[test]
    fn test_runner() {
        let mut config = RunConfiguration::build();
        config.reference = Some("test_data/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test");
        fs::create_dir("test").unwrap();
        let config = config.build();
        let _ = run_neat(Box::new(config), &mut NeatRng::seed_from_u64(0)).unwrap();
        fs::remove_dir_all("test").unwrap();
    }

    #[test]
    fn test_runner_files_messages() {
        let mut config = ConfigBuilder::new();
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.produce_vcf = true;
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("output");
        fs::create_dir("output").unwrap();
        let config = config.build();
        let _ = run_neat(Box::new(config), &mut NeatRng::seed_from_u64(0)).unwrap();
        fs::remove_dir_all("output").unwrap();
    }
}
