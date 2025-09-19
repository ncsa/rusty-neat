pub mod errors;
pub mod utils;
use log::*;
use errors::GenerateReadsErrors;

use common::file_tools::bed_reader::read_bed;
use common::structs::bed_record::BedRecord;
use std::collections::HashMap;
use std::path::PathBuf;
use crate::{
    filter_reads::utils::filter_lib::{
        filter_fastq, 
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
            for file_name in files_created {
                let unfiltered_file = PathBuf::from(&file_name);
                let mut outfile = PathBuf::from(file_name.clone());
                // This gets us to <PREFIX>.vcf
                outfile = outfile.with_extension("");
                let extension = outfile.extension();
                let mut is_fastq = true;
                match extension {
                    Some(ext) =>{
                        if ext == "vcf" {
                            is_fastq = false;
                        }
                    },
                    None => { 
                        unreachable!("Somehow a file made it to filtering with messed up extensions: {file_name}")
                    },
                }
                // This gets us to <PREFIX>
                outfile = outfile.with_extension("");
                let mut out_file_str = outfile.display().to_string();
                out_file_str.push_str("_filtered");
                if is_fastq {
                    out_file_str.push_str(".fastq.gz");
                    info!("Filtering output: {}", out_file_str);
                    filter_fastq(
                        &bed_table,
                        &unfiltered_file,
                        true,
                        &PathBuf::from(out_file_str.clone()),
                    ).expect(
                        &format!("Error filtering {file_name}")
                    );
                } else {
                    out_file_str.push_str(".vcf.gz");
                    filter_vcf(
                        &bed_table,
                        &unfiltered_file,
                        true,
                        &PathBuf::from(out_file_str.clone()),
                    ).expect(
                        &format!("Error filtering {file_name}")
                    );
                }
        
                info!("Filtered {file_name}, wrote: {}", out_file_str)
            }
        }
        info!("All files successfully filtered")
    };
    Ok(())
}