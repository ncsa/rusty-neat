use crate::utils;
use crate::data;
use common;

use std::thread;
use std::time;

use log::info;
use std::collections::{HashMap, VecDeque};
use std::sync::{Arc, Mutex};
use rand_chacha::ChaCha20Rng;

use utils::config::RunConfiguration;
use utils::fasta_tools::{read_fasta, write_fasta};
use utils::generate_reads::generate_reads;
use utils::read_models::{read_quality_score_model_file, read_quality_score_raw_data};
use utils::generate_variants::generate_variants;
use common::structs::nucleotides::Nuc;
use common::structs::variants::Variant;
use common::models::mutation_model::MutationModel;
use data::quality_score_data::RawQualityScoreData;

#[derive(Debug)]
pub enum RunNeatError {
    GeneralBullshit,
}

pub fn run_neat(config: Box<RunConfiguration>, rng: ChaCha20Rng) -> Result<(), RunNeatError> {
    let now = time::Instant::now();
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir.display(), config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    // Trying to think of how to avoid copying this. We could read the file, get the order, but
    // skip the reading in for now, let the threads handle that. Or we could create an index of byte
    // positions of the sequences, like a FAI, but that's a different type of problem to solve.
    // Alternative to that, require an index file and then use that to find the read positions?
    // That might be the easiest way to go.
    let (fasta_map, fasta_order) =
        read_fasta(&config.reference)
            .unwrap();

    let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
    for contig in &fasta_order {
        fasta_lengths.insert(contig.to_string(), fasta_map.get(contig).unwrap().len());
    }

    // Load models that will be used for the runs.
    // For now, we will use the one supplied, pulled directly from NEAT2.0's original model.
    // Later these variables will take the config input.
    let input_quality_score_model = false;
    let input_quality_model = String::new();

    let quality_score_model = {
        if input_quality_score_model {
            read_quality_score_model_file(
                &input_quality_model
            )
        } else {
            read_quality_score_raw_data(
                RawQualityScoreData::new()
            )
        }
    };

    // Todo load all models and set up the run.
    // load mutation model
    let mutation_model = MutationModel::new();

    // Mutating the reference and recording the variant locations.
    info!("Generating variants");

    // initialize all variants with the names of the contig, so it's not empty later
    let mut all_variants: Box<HashMap<String, HashMap<usize, Variant>>> =
        Box::new(HashMap::new());
    for contig in &fasta_order {
        all_variants.entry(contig.clone()).or_insert(HashMap::new());
    }
    let mut all_reads:Box<HashMap<String, Vec<(usize, usize)>>> = Box::new(HashMap::new());

    for contig in &fasta_order {
        info!("Generating variants for {}", contig);
        let contig_sequence = fasta_map.get(contig).unwrap();
        let contig_variants = generate_variants(
            contig_sequence,
            &mut mutation_model.clone(),
            config.ploidy,
            rng.clone()
        );

        info!("Finished generating variants for {}", contig);
        let _ = {
            all_variants.insert(contig.to_owned(), contig_variants)
        };

        // This all means if this run is fasta only, we can skip generating reads
        if config.produce_fastq ||
            config.produce_vcf ||
            config.produce_bam {
            let contig_sequence_len = fasta_map[contig].len().clone();
            let contig_reads = generate_reads(
                contig_sequence_len,
                config.read_len,
                config.coverage,
                config.paired_ended,
                config.fragment_mean,
                config.fragment_st_dev,
                rng.clone(),
            ).expect("Error generating reads");

            all_reads.insert(contig.to_owned(), contig_reads);
        };
    }

    if config.produce_fasta {
        info!("Producing fasta file");
        write_fasta(
            &fasta_map,
            &all_variants,
            &fasta_order,
            config.overwrite_output,
            &output_file,
        ).expect("Error writing fasta file!")
    }
    //
    // // If we produce any more files, we will need to know the genotype. It will be used in the
    // // output vcf, bam, and in the fastq.
    //
    // if config.produce_fastq {
    //     write_fastq(
    //         &fasta_map,
    //         &all_reads,
    //         &all_variants,
    //         config.read_len,
    //         config.overwrite_output,
    //         &output_file,
    //         config.paired_ended,
    //         quality_score_model,
    //         rng.clone()
    //     ).expect("Error writing fastq file(s)!")
    // }
    //
    // // todo for the vcf and bam, we need to know things like genotype
    // if config.produce_vcf {
    //     write_vcf(
    //         &all_variants,
    //         &fasta_order,
    //         &fasta_lengths,
    //         &config.reference,
    //         config.overwrite_output,
    //         &output_file,
    //     ).expect("Error writing vcf file!")
    // }

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
    let elapsed_time = now.elapsed();
    info!("Total run time {} seconds.", elapsed_time.as_secs());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::SeedableRng;
    use std::fs;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use utils::config::ConfigBuilder;
    use simplelog;
    use simplelog::*;

    #[test]
    fn test_runner() {
        let mut config = RunConfiguration::build();
        config.reference = Some("test_data/references/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test_runner");
        fs::create_dir("test_runner").unwrap();
        let config = config.build();
        let _ = run_neat(Box::new(config), ChaCha20Rng::seed_from_u64(0)).unwrap();
        fs::remove_dir_all("test_runner").unwrap();
    }

    #[test]
    fn test_runner_files_message() {
        let mut config = ConfigBuilder::new();
        config.reference = Some("test_data/references/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.produce_vcf = true;
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test_run_output");

        TermLogger::init(
            LevelFilter::Trace,
            Config::default(),
            TerminalMode::Stdout,
            ColorChoice::Auto,
        )
            .unwrap();

        fs::create_dir("test_run_output").unwrap();
        let config = Box::new(config.build());
        let rng = ChaCha20Rng::seed_from_u64(0);
        run_neat(config, rng.clone()).unwrap();
        let file_path = "test_run_output/neat_out.fasta";
        let input = File::open(file_path).unwrap();
        let buffered = BufReader::new(input);
        let line_count = buffered.lines().count();
        assert!(line_count > 0);
        fs::remove_dir_all("test_run_output").unwrap();
    }
}
