use crate::errors::GenerateReadsErrors;
use tempfile;
use std::time;
use log::{info, debug};
use std::collections::HashMap;
use std::path::PathBuf;
use simple_rng::NeatRng;

use crate::common::{
    file_tools::{
        fasta_reader::read_fasta,
        fastq_tools::write_fastq,
        vcf_tools::write_vcf,
    },
    structs::{
        variants::Variant,
        fasta_map::SequenceBlock,
        mutated_map::MutatedMap,
    },
    models::{
        mutation_model::MutationModel,
        fragment_length::FragmentLengthModel,
        quality_scores::QualityScoreModel,
        sequencing_error_model::{
            SequencingErrorModel
        }
    },
};
use crate::utils::{
    config::RunConfiguration,
    generate_reads::generate_reads,
    generate_variants::generate_variants,
};

pub fn run_neat(config: &Box<RunConfiguration>, rng: &mut NeatRng) -> Result<(), GenerateReadsErrors> {
    let start = time::Instant::now();
    info!("////////////// Welcome to rusty-neat read generator!");
    info!("Processing started!");
    // We'll need a temp dir to store file fragments
    let working_dir = tempfile::tempdir().unwrap();
    info!("Created temp dir at {:?}", working_dir);
    // Load models that will be used for the runs.
    //
    // Quality score model
    let quality_score_model = {
        match &config.quality_score_model {
            Some(filename) => {
                 QualityScoreModel::from_file(&filename)?
            },
            None => {
                QualityScoreModel::default()?
            }
        }
    };

    // load mutation model
    let mutation_model = {
        match &config.mutation_model {
            Some(filename) => {
                MutationModel::from_file(&filename)?
            },
            None => {
                MutationModel::default()?
            }
        }
    };

    // Fragment Length Model
    let fragment_length_model: FragmentLengthModel = {
        // FragmentLengthModel is an enum that allows us to choose one of two models.
        match &config.fragment_model {
            Some(filename) => {
                FragmentLengthModel::discrete_from_file(&filename)?.into()
            },
            None => {
                match config.fragment_mean {
                    Some(mean) => {
                        FragmentLengthModel::new_normal(
                            mean,
                            // Config should already have caught issues with missing data
                            config.fragment_st_dev.unwrap()
                        )?.into()
                    },
                    None => {
                        FragmentLengthModel::default()?.into()
                    },
                }
            }
        }
    };

    // Load Sequencing Error Model
    let seq_error_model: SequencingErrorModel = {
        match &config.sequence_error_model {
            Some(filename) => {
                SequencingErrorModel::from_file(&filename)?
            },
            None => {
                SequencingErrorModel::default()?
            }
        }
    };

    // The FastaMap struct consists of a vector of Contig objects, each describing a
    // block of Sequence, broken up into chunks in the SequenceBlock objects. It also
    // holds a name_map hashmap that links tho contig name from the Fasta with the shortname
    let fasta_map = read_fasta(
        &config.reference,
        config.read_len,
        &working_dir,
    ).unwrap();

    // Mutating the reference and recording the variant locations.
    info!("Generating simulated dataset");

    // all variants will be a hashmap of the contig name indexing a variants.
    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    // all reads is a hashmap of contigs and read tuples
    // Todo think about if we want to make a read container struct. 
    //    We have sequece retrieval covered by fasta map, so I'm feeling like it's
    //    just extra work at this point.
    let mut all_reads: Vec<(String, usize, usize, Option<usize>, usize)> = Vec::new();
    // Count the number of blocks

    // iterate over contigs
    for contig in &fasta_map.contigs {
        // Iterate over blocks within the contig
        // This is probably where we want to parallelize
        let contig_blocks = &contig.blocks;
        let contig_name = &contig.name;
        let mut contig_reads: Vec<(usize, usize, Option<usize>, usize)> = Vec::new();
        let mut contig_maps: Vec<MutatedMap> = Vec::new();
        debug!("Processing {}", contig_name);
        for block_filename in contig_blocks {
            let mut block_variants: Vec<Variant> = Vec::new();
            debug!("Processing block {:?}", block_filename);
            // Set up current block
            let current_block = SequenceBlock::from(&block_filename)?;
            // First we have to create the region weights data based on the fasta 
            // map and maybe gc-bias later at some point
            // filter down to minimal regions to search
            debug!("    > Generating bias map.");
            let regions_of_interest = current_block.get_non_n_regions()?;
            if regions_of_interest.is_empty() {
                // This block is all N's so we can skip
                continue
            }
            let mut bias_map: Vec<f64> = Vec::with_capacity(current_block.get_len());
            for _ in 0..current_block.get_len() {
                bias_map.push(0.0);
            }
            for region in regions_of_interest {
                // all of these have RegionType::NonNRegion
                // This is where we might apply bias
                for j in region.1..region.2 {
                    bias_map[j] = 1.0;
                }
            }
            // determine the number of mutations for this segment
            let num_mutations = (current_block.get_len() as f64 * config.mutation_rate)
                .trunc()
                as usize;
            debug!("Adding {} mutations to block {:?}", num_mutations, block_filename);
            // Generate any relevant mutations
            if num_mutations > 0 {
                // first we generate variants for the block
                let variants: Vec<Variant> = generate_variants(
                    &current_block, 
                    &bias_map, 
                    &mutation_model, 
                    num_mutations, 
                    config.ploidy, 
                    rng
                )?.unwrap();
                // Add variants 
                block_variants.extend(variants);
            }
            
            let block_reads: Vec<(usize, usize, Option<usize>, usize)> = generate_reads(
                current_block.get_len(),
                config.read_len,
                config.coverage,
                config.paired_ended,
                &fragment_length_model,
                rng,
            )?;

            contig_reads.extend(block_reads);

            let mutated_map = MutatedMap::new(
                block_filename.to_path_buf(),
                block_variants,
            )?;

            contig_maps.push(mutated_map);
        }
        mutated_maps.insert(contig_name.clone(), contig_maps.clone());
        // Add reads to all_reads
        for read in contig_reads {
            all_reads.push((contig_name.clone(), read.0, read.1, read.2, read.3))
        }
    }

    if config.produce_fastq {
        if let Some(filename_r1) = &config.output_fastq_1 {
            if let Some(filename_r2) = &config.output_fastq_2 {
                info!("Producing paired-ended fastq files");
                write_fastq(
                    &mut all_reads,
                    &mutated_maps,
                    config.read_len,
                    config.paired_ended,
                    (PathBuf::from(filename_r1), Some(PathBuf::from(filename_r2))),
                    &quality_score_model,
                    &seq_error_model,
                    rng,
                ).expect("Error writing fastq file(s)!")
            } else {
                info!("Producing single-ended fastq file");
                write_fastq(
                    &mut all_reads,
                    &mutated_maps,
                    config.read_len,
                    config.paired_ended,
                    (PathBuf::from(filename_r1), None),
                    &quality_score_model,
                    &seq_error_model,
                    rng,
                ).expect("Error writing fastq file(s)!")
            }
        } else {
            panic!("Error produce fastq set to true but no fastq filename set")
        }
    }

    match &config.output_vcf {
        Some(filename) => {
            info!("Writing output vcf file");
            // Maps contig to a total contig size, a required entry for a valid vcf file.
            let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
            let contigs = fasta_map.contigs.clone();
            for contig in contigs {
                fasta_lengths.insert(contig.name.clone(), contig.contig_len);
            }
            write_vcf(
                &mutated_maps,
                &fasta_map.contig_order,
                &fasta_lengths,
                &config.reference,
                config.overwrite_output,
                &filename,
            ).expect("Error writing vcf file!")
        },
        None => {
            info!("Skipping vcf creation")
        }
    }

    let elapsed_time = time::Instant::now() - start;
    info!("Processing finished in {} milliseconds", elapsed_time.as_millis());
    Ok(())
}

#[cfg(test)]
mod tests {
    use log::LevelFilter;
    use simplelog::{ColorChoice, Config, TermLogger, TerminalMode};

    use super::*;
    use std::{fs::{self, File}, io::{BufRead, BufReader}};
    use crate::utils::config::ConfigBuilder;

    #[test]
    fn test_runner() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test_runner");
        fs::create_dir("test_runner").unwrap();
        config.update_and_log().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let _ = run_neat(&Box::new(config), &mut rng);
        fs::remove_dir_all("test_runner").unwrap();
    }

    #[test]
    fn test_runner_files_message() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
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
        config.update_and_log().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        run_neat(&Box::new(config), &mut rng).unwrap();
        let file_path = "test_run_output/neat_out.fasta";
        let input = File::open(file_path).unwrap();
        let buffered = BufReader::new(input);
        let line_count = buffered.lines().count();
        assert!(line_count > 0);
        fs::remove_dir_all("test_run_output").unwrap();
    }
}
