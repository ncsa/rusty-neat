pub mod errors;
pub mod utils;
use log::*;
use errors::GenerateReadsErrors;

use clap::Parser;
use common::file_tools::bed_reader::read_bed;
use common::structs::bed_record::BedRecord;
use std::env;
use simplelog::*;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use crate::{
    filter_reads::utils::runner::{
        filter_fastq, 
        filter_paird_fastqs, 
        filter_vcf,
    }, 
    gen_reads::utils::{
        config::RunConfiguration, 
        runner::run_neat
    }
};
use simple_rng::NeatRng;
use std::{thread, time};


pub fn main(config: &PathBuf) -> Result<(), GenerateReadsErrors> {  
    if !config.is_file() {
        return Err(
            GenerateReadsErrors::CliError(
                "Must supply either a configuration file or a reference file to run NEAT!".to_string()
            )
        )
    }
    
    info!("////////////// Welcome to rusty-neat read generator! \\\\\\\\\\\\\\\\\\\\\\\\");
    thread::sleep(time::Duration::from_millis(1000));
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if &config.display().to_string() != "" {
        info!("Using Configuration file input: {:?}", &config);
        RunConfiguration::from_yaml_file(config)
            .expect("Error generating run configuration from input yaml")
    } else {
        panic!("Failed to supply config file!");
    };

    // Decide whether to run bed filtering or nah
    let mut bed_filter_enabled = false;
    let bed_filename = {
        if let Some(bed_file) = &config.filter_output {
            bed_filter_enabled = true;
            Some(PathBuf::from(bed_file))
        } else {
            None
        }
    };
    
    let mut files_created: Vec<String> = Vec::new();
    // Compile files to write
    if let Some(filename) = &config.output_fastq_1 {
        files_created.push(filename.display().to_string());
    };
    if let Some(filename) = &config.output_fastq_2 {
        files_created.push(filename.display().to_string());
    }
    if let Some(filename) = &config.output_vcf {
        files_created.push(filename.display().to_string());
    }

    info!("////////////// Configuration successuful! Ready to run! \\\\\\\\\\\\\\\\\\\\\\\\\\");
    thread::sleep(time::Duration::from_millis(100));
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let mut rng: NeatRng = NeatRng::new_from_seed(&config.seed_vec)
        .expect("Neat failed during rng creation!");
    // run the generate reads main script
    let result = run_neat(&Box::new(config.clone()), &mut rng);
    match result {
        Ok(()) => {
            // Continue on for bed filtering
            info!("Successfully produced unfiltered output file: {:?}", files_created)
        },
        Err(error) => {
            panic!("Main returned an error {:?}", error);
        },
    }

    if bed_filter_enabled {
        if let Some(bed_file) = bed_filename {
            info!("Filtering output based on input bed: {:?}", bed_file);
            let bed_table: HashMap<String, Vec<BedRecord>> = read_bed(&bed_file)?;
            if config.produce_fastq {
                if config.paired_ended {
                    // filter both
                    filter_paird_fastqs(
                        &bed_table, 
                        &PathBuf::from(files_created[0].clone()), 
                        &PathBuf::from(files_created[1].clone()),
                    ).expect("Error filtering fastq files");
                    info!("Successfully filtered fastq files: {:?}, {:?}", files_created[0], files_created[1]);
                } else {
                    // filter only one
                    filter_fastq(&bed_table, &PathBuf::from(files_created[0].clone()))
                        .expect("Error filtering fastq file");
                    info!("Successfully filtered fastq file: {:?}", files_created[0]);
                }
            }

            if config.produce_vcf {
                // filter vcf
                let last_file = files_created.len() - 1;
                filter_vcf(&bed_table, &PathBuf::from(files_created[last_file].clone()))
                    .expect("Error filtering vcf file");
                info!("Successfully filtered fastq file: {:?}", files_created[last_file]);
            }

        } else {
            panic!("Bed filtering was enabled, but no bed file name specified")
        }
    }
    Ok(())
}