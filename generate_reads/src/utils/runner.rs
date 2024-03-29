use std::collections::HashSet;
use log::info;
use rand::prelude::SliceRandom;
use utils::config::RunConfiguration;
use utils::fasta_tools::{read_fasta, write_fasta};
use utils::fastq_tools::write_fastq;
use utils::make_reads::generate_reads;
use utils::mutate::mutate_fasta;
use utils::neat_rng::NeatRng;
use utils::vcf_tools::write_vcf;
use utils::read_models::read_quality_score_model_json;

pub fn run_neat(config: Box<RunConfiguration>, mut rng: &mut NeatRng) -> Result<(), &'static str>{
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir.display(), config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    let (fasta_map, fasta_order) = read_fasta(&config.reference)
        .unwrap();

    // Load models that will be used for the runs.
    // For now we will use the one supplied, pulled directly from NEAT2.0's original model.
    let default_quality_score_model_file = "model_data/neat_quality_score_model.json";
    let quality_score_model = read_quality_score_model_json(
        default_quality_score_model_file
    );

    // Mutating the reference and recording the variant locations.
    info!("Mutating reference.");
    let (mutated_map, variant_locations) = mutate_fasta(
        &fasta_map,
        config.minimum_mutations,
        &mut rng
    );

    if config.produce_fasta {
        info!("Outputting fasta file");
        write_fasta(
            &mutated_map,
            &fasta_order,
            config.overwrite_output,
            &output_file,
        ).unwrap();
    }

    if config.produce_vcf {
        info!("Writing vcf file");
        write_vcf(
            &variant_locations,
            &fasta_order,
            config.ploidy,
            &config.reference,
            config.overwrite_output,
            &output_file,
            &mut rng
        ).unwrap();
    }

    let mut read_sets: HashSet<Vec<u8>> = HashSet::new();
    info!("Generating reads");
    for (_name, sequence) in mutated_map.iter() {
        // defined as a set of read sequences that should cover
        // the mutated sequence `coverage` number of times
        let data_set = generate_reads(
            sequence,
            &config.read_len,
            &config.coverage,
            config.paired_ended,
            config.fragment_mean,
            config.fragment_st_dev,
            &mut rng
        ).unwrap();

        read_sets.extend(*data_set);
    }

    if config.produce_fastq {
        info!("Shuffling output fastq data");
        let mut outsets: Box<Vec<&Vec<u8>>> = Box::new(read_sets.iter().collect());
        outsets.shuffle(&mut rng);

        info!("Writing fastq");
        write_fastq(
            &output_file,
            config.overwrite_output,
            config.paired_ended,
            *outsets,
            quality_score_model,
            rng,
        ).unwrap();
        info!("Processing complete")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use rand_core::SeedableRng;
    use utils::config::ConfigBuilder;
    use super::*;

    #[test]
    fn test_runner() {
        let mut config = RunConfiguration::build();
        config.reference = Some("data/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test");
        fs::create_dir("test").unwrap();
        let config = config.build();
        let _ = run_neat(
            Box::new(config),
            &mut NeatRng::seed_from_u64(0),
        ).unwrap();
        fs::remove_dir_all("test").unwrap();
    }

    #[test]
    fn test_runner_files_messagse() {
        let mut config = ConfigBuilder::new();
        config.reference = Some("data/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.produce_vcf = true;
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("output");
        fs::create_dir("output").unwrap();
        let config = config.build();
        let _ = run_neat(
            Box::new(config),
            &mut NeatRng::seed_from_u64(0),
        ).unwrap();
        fs::remove_dir_all("output").unwrap();
    }
}