use tempfile::{self, tempdir};
use std::time;
use log::{info, debug, error};
use std::collections::HashMap;
use std::path::PathBuf;
use simple_rng::NeatRng;
use flate2::{ 
    Compression,
    write::GzEncoder
};

use crate::errors::GenerateReadsErrors;
use crate::common::{
    file_tools::{
        fasta_reader::read_fasta,
        vcf_tools::write_vcf,
        fastq_tools::{write_block_fastq, combine_temp_fastqs},
        file_io::{append_to_file, VectorBuffer}
    },
    structs::{
        variants::Variant,
        fasta_map::SequenceBlock,
        mutated_map::MutatedMap,
        nucleotides::NucleotideSelector
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
    generate_fragments::generate_fragments,
    generate_variants::generate_variants,
};

pub fn run_neat(config: &Box<RunConfiguration>, rng: &mut NeatRng) -> Result<(), GenerateReadsErrors> {
    let start = time::Instant::now();
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

    // This model provides an alternate base for N based on a simple equal weight distribution.
    let nuc_sub_model: NucleotideSelector = {
        NucleotideSelector::new()
    };

    // The FastaMap struct consists of a vector of Contig objects, each describing a
    // block of Sequence, broken up into chunks in the SequenceBlock objects. It also
    // holds a name_map hashmap that links tho contig name from the Fasta with the shortname
    let fasta_map = read_fasta(
        &config.reference,
        nuc_sub_model,
        config.read_len,
        &working_dir,
        rng,
    )?;

    // Mutating the reference and recording the variant locations.
    info!("Generating simulated dataset");

    // all variants will be a hashmap of the contig name indexing a variants.
    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    // All files written indexed by contig_name
    let mut all_fastq_files: HashMap<String, Vec<PathBuf>> = HashMap::new(); 
    // iterate over contigs
    for contig in &fasta_map.contigs {
        // Iterate over blocks within the contig
        // This is probably where we want to parallelize
        let contig_blocks = &contig.blocks;
        let contig_name = &contig.name;
        let mut contig_files: Vec<PathBuf> = Vec::new();
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
                for j in region.start..region.end {
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
                let result = generate_variants(
                    &current_block, 
                    &bias_map, 
                    &mutation_model, 
                    num_mutations, 
                    config.ploidy, 
                    rng
                )?;
                match result {
                    Some(vec) => {
                        // Add variants 
                        block_variants.extend(vec);
                    },
                    None => {
                        // No variants for this block
                    }
                }
            }
            
            // blocks get returned sorted
            let block_fragments: Vec<(usize, usize)> = generate_fragments(
                current_block.get_len(),
                config.read_len,
                config.coverage,
                config.paired_ended,
                &fragment_length_model,
                rng,
            )?;

            let mutated_map = MutatedMap::new(
                block_filename.to_path_buf(),
                block_variants,
            )?;

            if config.produce_fastq {
                match &config.output_fastq_1 {
                    Some(_) => {
                        // Need a file name that helps us ID this fastq later
                        let mut file_to_write_1 = PathBuf::from(working_dir.path());
                        file_to_write_1
                            .push(&config.output_filename);
                        file_to_write_1
                            .push(
                                format!(
                                    "_{}_{}_{}_r1_tmp.fastq.gz",
                                    contig_name,
                                    current_block.get_start()?, 
                                    current_block.get_end()?,
                                )
                            );
                        let file1 = append_to_file(&file_to_write_1)?;
                        let mut buffer1 = GzEncoder::new(
                            file1, Compression::default()
                        );
                        let read_name_prefix = format!(
                            "neat_generated_{}_{}_{}_",
                            current_block.contig,
                            current_block.get_start()?,
                            current_block.get_end()?,
                        );
                        if config.paired_ended {
                            match &config.output_fastq_2 {
                                Some(_) => {
                                    let mut file_to_write_2 = PathBuf::from(working_dir.path());
                                    file_to_write_2
                                        .push(&config.output_filename);
                                    file_to_write_2
                                        .push(
                                            format!(
                                                "_{}_{}_{}_r2_tmp.fastq.gz",
                                                contig_name,
                                                current_block.get_start()?, 
                                                current_block.get_end()?
                                            )
                                        );
                                    let file2 = append_to_file(&file_to_write_2)?;
                                    let mut buffer2 = GzEncoder::new(
                                        file2, Compression::default()
                                    );
                                    debug!("Writing paired-ended block fastq files");
                                    write_block_fastq (
                                        block_fragments,
                                        &mutated_map,
                                        true,
                                        &mut buffer1,
                                        &mut buffer2,
                                        config.read_len,
                                        &quality_score_model,
                                        &seq_error_model,
                                        rng,
                                    )?;
                                    contig_files
                                        .push(file_to_write_1);
                                    contig_files
                                        .push(file_to_write_2);
                                },
                                None =>{
                                    error!("Produce fastq true and paired-ended true, but output_fastq_2 was empty.");
                                    return Err(GenerateReadsErrors::ConfigError)
                                }
                            }
                        } else {
                            debug!("Writing single-ended block fastq file");
                            // Single-ended mode. We create a fake vector to hold the data we won't write.
                            // dummy will be destroyed at the end of this block without being used.
                            let dummy_data: VectorBuffer = VectorBuffer::new();
                            let mut buffer2 = GzEncoder::new(
                                dummy_data, Compression::default()
                            );
                            write_block_fastq (
                                block_fragments,
                                &mutated_map,
                                false,
                                &mut buffer1,
                                &mut buffer2,
                                config.read_len,
                                &quality_score_model,
                                &seq_error_model,
                                rng,
                            )?;
                            contig_files
                                .push(file_to_write_1)
                        }
                    },
                    None => {
                        error!("Produce fastq true, but output_fastq_1 was empty.");
                        return Err(GenerateReadsErrors::ConfigError)
                    },
                }
                
            } 
            contig_maps.push(mutated_map);
        }
        mutated_maps.insert(contig_name.clone(), contig_maps.clone());
        all_fastq_files.insert(contig_name.clone(), contig_files);
    }

    if config.produce_fastq {
        // First we need to shuffle the output order using a HashMap that maps
        // The original filename: read number to a final read number. Will need
        // to keep pairs together during this
        //
        // We know it's max 2 files for fastqs
        let mut final_fastqs = Vec::with_capacity(2);

        // Shuffle all reads to randomize fastq output
        
        if let Some(filename_r1) = &config.output_fastq_1 {
            if let Some(filename_r2) = &config.output_fastq_2 {
                info!("Producing paired-ended fastq files");
                info!("Writing fastq 1");
                combine_temp_fastqs( todo!());
                info!("Writing fastq 2");
                combine_temp_fastqs(todo!());
                final_fastqs.push(filename_r1.clone());
                final_fastqs.push(filename_r2.clone());
            } else {
                info!("Producing single-ended fastq file");
                combine_temp_fastqs(todo!());
                final_fastqs.push(filename_r1.clone());
            }
        } else {
            panic!("Error produce fastq set to true but no fastq filename set")
        };
        info!("Fastq(s) successfully written: {:?}!", final_fastqs);
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
