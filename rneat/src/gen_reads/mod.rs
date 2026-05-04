pub mod errors;
pub mod utils;
use log::*;
use errors::GenerateReadsErrors;

use std::path::PathBuf;
use crate::{
    gen_reads::utils::{
        config::RunConfiguration, 
        runner::run_neat
    }
};
use simple_rng::NeatRng;

/// gen-reads is the primary read generation function of rneat. It reads a fasta file and generates a set of fastqs and/or a set of variants. It can now also filter reads by bed file.
pub fn main(config: &PathBuf) -> Result<(), GenerateReadsErrors> {   
    info!("////////////// Welcome to rusty-neat read generator! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if &config.display().to_string() != "" {
        info!("Using Configuration file input: {:?}", &config);
        RunConfiguration::from_yaml_file(config)
            .expect("Error generating run configuration from input yaml")
    } else {
        panic!("Failed to supply config file!");
    };

    // Check that we are clear to write outputs
    if !config.overwrite_output {
        if let Some(filename) = &config.output_fastq_1 {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
        if let Some(filename) = &config.output_fastq_2 {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
        if let Some(filename) = &config.output_vcf {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
    }
    
    info!("////////////// Configuration successuful! Ready to run! \\\\\\\\\\\\\\\\\\\\\\\\\\");
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let mut rng: NeatRng = NeatRng::new_from_seed(&config.seed_vec)
        .expect("Neat failed during rng creation!");
    // run the generate reads main script
    let result = run_neat(&Box::new(config.clone()), &mut rng);
    match result {
        Ok(_) => {
            // Continue on for bed filtering
            Ok(())
        },
        Err(error) => {
            error!("runner returned an error {:?}", error);
            Err(GenerateReadsErrors::RunnerError)
        },
    }
}