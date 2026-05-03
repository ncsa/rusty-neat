// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use serde_yml::Value;
use std::collections::HashMap;
use std::path::PathBuf;
use std::string::String;
use std::fs;
use common::{
    file_tools::{
        bed_reader::read_bed, 
        vcf_tools::read_vcf
    }, 
    structs::{
        bed_record::BedRecord, 
        variants::Variant
    }
};
use crate::gen_mut_model::errors::GenMutationModelError;

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this filter-reads run. It is built from user-supplied input
    // in the form of a configuration yaml file
    //
    // bed_file: The path to the bed_file for the run.
    // files_to_filter: The list of files to filter.
    // filter_key: The key to add to the filtered file names so you know they have been filtered.
    pub reference: PathBuf,
    pub mutations: HashMap<String, Vec<Variant>>,
    pub bed_table: HashMap<String, Vec<BedRecord>>,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenMutationModelError> {
        // Reads an input configuration file from yaml using the serde package. Then sets the
        // parameters based on the inputs. A "." value means to use the default value.
        //
        // Opens file for reading
        let f = fs::File::open(&yml_file);
        let file = match f {
            Ok(l) => l,
            Err(error) => panic!(
                "Problem reading the config file: {:?}. System error: {}",
                &yml_file,
                error,
            ),
        };
        // Uses serde_yaml to read the file into a HashMap
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)
            .expect("Error reading yaml file!");
        // Fill in the bed_file first, all hinges on that
        let reference = PathBuf::from(scrape_config["reference"].as_str().unwrap());
        if !reference.is_file() {
            panic!("Invalid reference file {:?}", reference)
        }
        let vcf_file = PathBuf::from(scrape_config["vcf_file"].as_str().unwrap());
        if !vcf_file.is_file() {
            panic!("Invalid bed file {:?}", vcf_file)
        }
        let bed_file_raw = scrape_config["bed_file"].as_str().unwrap_or(".");
        let bed_table = if bed_file_raw == "." {
            HashMap::new()
        } else {
            let bed_file = PathBuf::from(bed_file_raw);
            if !bed_file.is_file() {
                panic!("Invalid bed file {:?}", bed_file)
            }
            read_bed(&bed_file).expect("Error reading bed file!")
        };
        let overwrite_output = scrape_config["overwrite_output"].as_bool().unwrap_or(false);
        let output_file = PathBuf::from(scrape_config["output_file"].as_str().unwrap());
        if !overwrite_output {
            if output_file.is_file() {
                panic!("Attempting to overwrite an existing file {:?}", output_file)
            }
        }
        let mutations = read_vcf(vcf_file)?;

        Ok(RunConfiguration {
            reference,
            mutations,
            bed_table,
            output_file,
            overwrite_output,
        })
    }
}

pub fn create_map_item(
    split_string: Vec<&str>, 
    raw_string: &str,
    length: usize,
    filter_key: &str,
) -> (PathBuf, PathBuf) {
    let old_element = split_string[length-3];
    let new_element = format!("{old_element}{filter_key}");
    let mut output_name = String::new();
    // stop one short of the end
    for i in 0..length-1 {
        if i == length-3 {
            output_name.push_str(&new_element);
            output_name.push('.');
        } else {
            output_name.push_str(split_string[i]);
            output_name.push('.');
        }
    }
    // End with the last extension with no traling dot.
    output_name.push_str(split_string[length-1]);
    let temp_path = PathBuf::from(raw_string);
    if !temp_path.is_file() {
        panic!("Input file not found! {:?}", temp_path)
    }
    (temp_path, PathBuf::from(output_name))
}

#[cfg(test)]
mod tests {
    // use super::*;

    #[test]
    fn test_run_configuration() {
        // TODO create test
        assert!(true)
    }
}
