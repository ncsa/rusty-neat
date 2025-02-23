use std::collections::HashSet;
use log::info;
use simple_rng::Rng;
use super::config::RunConfiguration;
use super::fasta_tools::{read_fasta, write_fasta};
use super::fastq_tools::write_fastq;
use super::make_reads::generate_reads;
use super::mutate::mutate_fasta;
use super::vcf_tools::write_vcf;
use super::read_models::read_quality_score_model_json;

pub fn run_neat(config: Box<RunConfiguration>, mut rng: &mut Rng) -> Result<(), &'static str>{
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir.display(), config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    let (fasta_map, fasta_order) = read_fasta(&config.reference)
        .unwrap();

    // Load models that will be used for the runs.
    // For now we will use the one supplied, pulled directly from NEAT2.0's original model.
    let default_quality_score_model_file = "models/neat_quality_score_model.json";
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
        let outsets: Box<Vec<&Vec<u8>>> = Box::new(read_sets.iter().collect());
        let mut outsets_order: Vec<usize> = (0..outsets.len()).collect();
        rng.shuffle(&mut outsets_order);

        info!("Writing fastq");
        write_fastq(
            &output_file,
            config.overwrite_output,
            config.paired_ended,
            *outsets,
            outsets_order,
            quality_score_model,
            rng,
        ).unwrap();
        info!("Processing complete")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use super::super::config::ConfigBuilder;

    #[test]
    fn test_runner() {
        let mut config = RunConfiguration::build();
        config.reference = Some("test_data/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test");
        fs::create_dir("test").unwrap();
        let config = config.build();
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let _ = run_neat(
            Box::new(config),
            &mut rng,
        ).unwrap();
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
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let _ = run_neat(
            Box::new(config),
            &mut rng,
        ).unwrap();
        fs::remove_dir_all("output").unwrap();
    }
}