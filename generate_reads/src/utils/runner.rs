use crate::errors::GenerateReadsErrors;
use tempfile;
use std::time;
use log::{info, debug};
use std::collections::HashMap;
use simple_rng::NeatRng;

use crate::common::{
    file_tools::fasta_reader::read_fasta,
    structs::{
        variants::Variant,
        fasta_map::SequenceBlock
    },
    models::{
        mutation_model::MutationModel,
        fragment_length::{
            DiscreteFragmentLengthModel,
            FragmentLengthModel,
            NormalFragmentLengthModel,
        },
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
        match config.quality_score_model {
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
        match config.mutation_model {
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
        match config.fragment_model {
            Some(filename) => {
                DiscreteFragmentLengthModel::discrete_from_file(&filename)?.into()
            },
            None => {
                match config.fragment_mean {
                    Some(mean) => {
                        NormalFragmentLengthModel::new_from_mean(
                            mean,
                            // Config should already have caught issues with missing data
                            config.fragment_st_dev.unwrap()
                        )?.into()
                    },
                    None => {
                        NormalFragmentLengthModel::default()?.into()
                    },
                }
            }
        }
    };

    // Load Sequencing Error Model
    let seq_error_model: SequencingErrorModel = {
        match config.sequence_error_model {
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
    let all_variants: HashMap<String, Vec<Variant>> = HashMap::new();
    // all reads is a hashmap of contigs and read tuples
    // Todo think about if we want to make a read container struct. 
    //    We have sequece retrieval covered by fasta map, so I'm feeling like it's
    //    just extra work at this point.
    let all_reads: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    // iterate over contigs
    for (i, contig) in fasta_map.contigs.iter().enumerate() {
        // Iterate over blocks within the contig
        // This is probably where we want to parallelize
        let contig_blocks = contig.blocks;
        debug!("Processing {}", contig.name);
        for (j, block_filename) in contig.blocks.iter().enumerate() {
            let mut block_variants: Vec<Variant> = Vec::new();
            let mut block_reads: Vec<(usize, usize)> = Vec::new();
            debug!("Processing block {:?}", block_filename);
            // Set up current block
            let current_block = SequenceBlock::from(&block_filename)?;
            // First we have to create the region weights data based on the fasta 
            // map and maybe gc-bias later at some point
            let block_start = current_block.ref_start;
            let block_end = current_block.ref_end;
            // filter down to minimal regions to search
            debug!("    > Generating bias map.");
            let regions_of_interest = current_block.get_non_n_regions()?;
            if regions_of_interest.is_empty() {
                // This block is all N's so we can skip
                continue
            }
            let mut bias_map: Vec<f64> = Vec::with_capacity(current_block.get_len());
            for i in 0..current_block.get_len() {
                bias_map.push(0.0);
            }
            for region in regions_of_interest {
                // all of these have RegionType::NonNRegion
                // This is where we might apply bias
                for i in region.1..region.2 {
                    bias_map[i] = 1.0;
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
                let variants = generate_variants(
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
            
            let contig_reads = generate_reads(
                current_block.get_len(),
                config.read_len,
                config.coverage,
                config.paired_ended,
                fragment_length_model,
                rng,
            );

            /// TODO update to write block fastq file
            if config.produce_fastq {
                write_fastq(
                    &current_block,
                    &mut all_reads,
                    config.read_len,
                    config.overwrite_output,
                    &,
                    config.paired_ended,
                    quality_score_model,
                    rng.clone()
                ).expect("Error writing fastq file(s)!")
            }

            /// TODO update to write block vcf flie
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
        }
    }
    let elapsed_time = time::Instant::now() - start;
    info!("Processing finished in {} milliseconds", elapsed_time.as_millis());
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
