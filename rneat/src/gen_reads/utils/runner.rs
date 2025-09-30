use common::file_tools::file_io::create_output_file;
use common::structs::variants::VariantType;
use tempfile;
use log::{info, debug, error};
use simple_rng::NeatRng;
use std::collections::HashMap;
use std::path::PathBuf;
use std::io::BufWriter;

use indicatif::ProgressBar;
use flate2::{ 
    Compression,
    write::GzEncoder
};

use crate::{
    common::{
        file_tools::{
            fasta_reader::read_fasta, 
            fastq_tools::{
                combine_temp_fastqs, 
                write_block_fastq
            }, 
            file_io::{append_to_file, VectorBuffer}, 
            vcf_tools::write_vcf
        }, models::{
            fragment_length::FragmentLengthModel, 
            mutation_model::MutationModel, 
            sequencing_error_model::SequencingErrorModel
        }, structs::{
            fasta_map::SequenceBlock, 
            mutated_map::MutatedMap, 
            nucleotides::NucleotideSelector, 
            variants::Variant
        }
    }, 
    gen_reads::{
        errors::GenerateReadsErrors, 
        utils::{
            config::RunConfiguration, 
            generate_variants::generate_variants,
            generate_fragments::generate_fragments,
        }
    }
};


pub fn run_neat(config: &Box<RunConfiguration>, rng: &mut NeatRng) -> Result<Vec<PathBuf>, GenerateReadsErrors> {
    // We'll need a temp dir to store file fragments
    let working_dir = tempfile::tempdir().unwrap();
    info!("Created temp dir at {:?}", working_dir);
    // Load models that will be used for the runs
    //
    // load mutation model
    info!("Generate mutation model");
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
    info!("Generate fragment length model");
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
    info!("Generate sequencing error model");
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
    info!("Initialize Nucleotide selector");
    let nuc_sub_model: NucleotideSelector = {
        NucleotideSelector::new()
    };

    // The FastaMap struct consists of a vector of Contig objects, each describing a
    // block of Sequence, broken up into chunks in the SequenceBlock objects. It also
    // holds a name_map hashmap that links tho contig name from the Fasta with the shortname
    info!("Reading fasta file: {}", &config.reference.display());    
    let fasta_map = read_fasta(
        &config.reference,
        Some(&nuc_sub_model),
        config.read_len,
        &working_dir,
        Some(rng),
    )?;

    // Mutating the reference and recording the variant locations.
    info!("Generating simulated dataset");

    // all variants will be a hashmap of the contig name indexing a variants.
    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    // All files written indexed by contig_name, separated by read1/read2
    let mut all_fastq_files: HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)> = HashMap::new();
    let bar: ProgressBar = ProgressBar::new(fasta_map.contigs.len() as u64);
    // to display bar
    bar.tick();
    // iterate over contigs. 
    for contig in &fasta_map.contigs {
        // Iterate over blocks within the contig
        // This is probably where we want to parallelize
        let contig_blocks = &contig.blocks;
        let contig_name = &contig.name;
        let mut contig_files_r1: Vec<PathBuf> = Vec::new();
        let mut contig_files_r2: Vec<PathBuf> = Vec::new();
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
            let mut bias_map: Vec<f64> = vec![0.0; current_block.get_len()];
            let mut num_mutable_area = 0.0;
            for region in &regions_of_interest {
                // all of these have RegionType::NonNRegion
                // This is where we might apply bias
                for j in region.start..region.end {
                    bias_map[j] = 1.0;
                    num_mutable_area += 1.0;
                }
            }
            // determine the number of mutations for this segment
            let num_mutations = (num_mutable_area * config.mutation_rate)
                .trunc()
                as usize;
            debug!("Adding {} mutations to block {:?}", num_mutations, block_filename);
            let mut max_del_len = 0;
            // Generate any relevant mutations
            if num_mutations > 0 {
                // first we generate variants for the block
                let result: Option<Vec<Variant>> = generate_variants(
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
                        for variant in vec {
                            if variant.variant_type == VariantType::Deletion {
                                // -1 because we're ignoring the reference base
                                if variant.reference.len() - 1 > max_del_len {
                                    // This will get us the maximum del length.
                                    max_del_len = variant.reference.len() - 1
                                }
                            }
                            block_variants.push(variant);
                        }
                    },
                    None => {
                        // No variants for this block
                    }
                }
            }

            let block_fragments: Vec<(usize, usize)> = {
                let mut temp_regions: Vec<(usize, usize)> = Vec::new();
                // find matching region.
                // if next region is within a read, maybe we'll use the whole thing.
                for region in regions_of_interest {
                    temp_regions.push((region.start, region.end));
                }
                let mut block_frags = Vec::new();
                for (frag_start, frag_end) in temp_regions {
                    let result = generate_fragments(
                        frag_end-frag_start,
                        config.read_len,
                        max_del_len,
                        frag_start,
                        config.coverage,
                        &fragment_length_model,
                        rng,
                    );
                    match result {
                        Ok(frags) => {
                            if !frags.is_empty() {
                                block_frags.extend_from_slice(&frags);
                            }
                        },
                        Err(error) => return Err(error),
                    }
                    
                }
                block_frags.clone()
            };

            let mutated_map = MutatedMap::new(
                block_filename.to_path_buf(),
                block_variants,
            )?;

            if config.produce_fastq {
                // Need a file name that helps us ID this fastq later
                let mut file_to_write_1 = PathBuf::from(working_dir.path());
                file_to_write_1
                    .push(
                        format!(
                            "temp_{}_{:010}_{:010}_r1_tmp.fastq.gz",
                            contig_name,
                            current_block.get_start()?, 
                            current_block.get_end()?,
                        )
                    );
                let file1 = append_to_file(&file_to_write_1)?;
                let writer1 = BufWriter::new(&file1);
                let mut buffer1 = GzEncoder::new(
                    writer1, Compression::default()
                );
                let read_name_prefix = format!(
                    "neat_generated_{}",
                    current_block.contig,
                );
                if config.paired_ended {
                    let mut file_to_write_2 = PathBuf::from(working_dir.path());
                    file_to_write_2
                        .push(
                            format!(
                                "temp_{}_{:010}_{:010}_r2_tmp.fastq.gz",
                                contig_name,
                                current_block.get_start()?, 
                                current_block.get_end()?
                            )
                        );
                    let file2 = append_to_file(&file_to_write_2)?;
                    let writer2 = BufWriter::new(&file2);
                    let mut buffer2 = GzEncoder::new(
                        writer2, Compression::default()
                    );
                    debug!("Writing paired-ended block fastq files");
                    write_block_fastq (
                        block_fragments,
                        &mutated_map,
                        true,
                        &mut buffer1,
                        &mut buffer2,
                        config.read_len,
                        &read_name_prefix,
                        &seq_error_model,
                        rng,
                    )?;
                    contig_files_r1
                        .push(file_to_write_1);
                    contig_files_r2
                        .push(file_to_write_2);
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
                        &read_name_prefix,
                        &seq_error_model,
                        rng,
                    )?;
                    contig_files_r1
                        .push(file_to_write_1)
                }
            } 
            contig_maps.push(mutated_map);
        }
        mutated_maps.insert(contig_name.clone(), contig_maps.clone());
        all_fastq_files.insert(contig_name.clone(), (contig_files_r1, contig_files_r2));
        bar.inc(1);
    }

    bar.finish_and_clear();
    
    if config.produce_fastq {
        // First we need to shuffle the output order using a HashMap that maps
        // The original filename: read number to a final read number. Will need
        // to keep pairs together during this
        //
        info!("Producing final fastq(s) file(s)");
        let bar: ProgressBar = ProgressBar::new(all_fastq_files.len() as u64);
        match &config.output_fastq_1 {
            Some(filename1) => {
                // If this is already a file, this will clear it.
                create_output_file(filename1, config.overwrite_output)?;
                if config.paired_ended {
                    match &config.output_fastq_2 {
                        Some(filename2) => {
                            // If this is already a file, this will clear it.
                            create_output_file(filename2, config.overwrite_output)?;
                            for (_contig, temp_fastqs) in all_fastq_files {
                                combine_temp_fastqs(
                                    temp_fastqs.0, 
                                    &filename1, 
                                    false
                                )?;
                                combine_temp_fastqs(
                                    temp_fastqs.1,
                                    &filename2,
                                    false,
                                )?;
                                bar.inc(1);
                            }
                        },
                        None => {
                            error!("Produce fastq true and paired-ended true, but output_fastq_2 was missing.");
                            return Err(GenerateReadsErrors::ConfigError)
                        }
                    }
                } else {
                    for (_contig, temp_fastqs) in all_fastq_files {
                        combine_temp_fastqs(
                            temp_fastqs.0, 
                            &filename1, 
                            false
                        )?;
                        bar.inc(1);
                    }
                }
            },
            None => {
                error!("Produce fastq true but output_fastq_1 was missing.");
                return Err(GenerateReadsErrors::ConfigError)
            },
        } 
    }

    bar.finish_and_clear();

    let mut files_written = Vec::new();
    if config.paired_ended {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", &filename1);
            files_written.push(filename1.clone());
            if let Some(filename2) = &config.output_fastq_2 {
                info!("Successfully wrote fastq file: {:?}", &filename2);
                files_written.push(filename2.clone());
            }
        }
    } else {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", &filename1);
            files_written.push(filename1.clone());
        }
    }

    if let Some(filename) = &config.output_vcf {
        info!("Writing output vcf file");
        // Maps contig to a total contig size, a required entry for a valid vcf file.
        let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
        let contigs = fasta_map.contigs.clone();
        for contig in contigs {
            fasta_lengths.insert(contig.name.clone(), contig.contig_len);
        }
        let result = write_vcf(
            &mutated_maps,
            &fasta_map.contig_order,
            &fasta_lengths,
            &config.reference,
            config.overwrite_output,
            &filename,
        );
        match result {
            Ok(()) => {
                info!("Successfully wrote vcf file: {:?}", filename);
                files_written.push(filename.clone());
            },
            Err(error) => {
                error!("Error writing vcf file!");
                return Err(GenerateReadsErrors::IoError(error))
            },
        }
    }
    Ok(files_written.clone())
}
