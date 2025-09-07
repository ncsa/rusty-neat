// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use chrono::Utc;

use log::{info, warn};
use serde_yml::Value;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::string::String;
use std::{env, fs};
use crate::{
    common::file_tools::folder_tools::check_create_dir,
    utils::cli::Cli
};

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this particular run. It is derived from input either
    // from a configuration file or from command line inputs. This is not built directly in the code,
    // but is constructed by a builder to enable default values.
    //
    // reference: The path to the reference for the run.
    // read_len: The length of reads in the output fastq.
    // coverage: The average depth of coverage for the output fastq file.
    // mutation_rate: The rate of mutation for the file.
    // ploidy: The number of copies of each chromosome in the target organism. Mutation process will
    // be replicated this number of times.
    // paired_ended: If the run will be in paired-ended mode.
    // fragment_mean: Mean size of the fragments in paired-ended mode.
    // fragment_st_dev: Standard deviation of the fragment mean describing the sample set to sample
    // fragments from.
    // produce_fastq: True or false on whether to produce an output fastq file.
    // produce_vcf: True or false on whether to produce an output VCF file, with genotyped variants.
    // produce_bam: True or false on whether to produce an output BAM file, which will be aligned to
    // the reference.
    // overwrite_output: if true, will overwrite output. If false will error and exit you attempt to
    // overwrite files with the same name.
    // output_dir: The directory, relative or absolute, path to the directory to place output.
    // output_prefix: The name to use for the output files.
    pub reference: String,
    pub read_len: usize,
    pub coverage: usize,
    pub mutation_rate: f64,
    pub ploidy: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    pub produce_fastq: bool,
    pub produce_vcf: bool,
    pub produce_bam: bool,
    pub rng_seed: Option<String>,
    pub seed_vec: Vec<String>,
    pub overwrite_output: bool,
    pub minimum_mutations: usize,
    pub output_dir: PathBuf,
    pub output_prefix: String,
    pub output_fastq_1: Option<String>,
    pub output_fastq_2: Option<String>,
    pub output_vcf: Option<String>,
    pub output_bam: Option<String>, 
    // model input
    pub quality_score_model: Option<PathBuf>,
    pub mutation_model: Option<PathBuf>,
    pub fragment_model: Option<PathBuf>,
    pub sequence_error_model: Option<PathBuf>,
}

impl RunConfiguration {
    pub (crate) fn default() -> RunConfiguration {
        let mut config = RunConfiguration { 
            reference: String::new(),
            read_len: 151, 
            coverage: 10, 
            mutation_rate: 0.001, 
            ploidy: 2, 
            paired_ended: false, 
            fragment_mean: None, 
            fragment_st_dev: None, 
            produce_fastq: true, 
            produce_vcf: false, 
            produce_bam: false, 
            rng_seed: None,
            seed_vec: Vec::new(),
            overwrite_output: false, 
            minimum_mutations: 0,
            output_dir: env::current_dir().unwrap(), 
            output_prefix: "neat_out".to_string(),
            output_fastq_1: None,
            output_fastq_2: None,
            output_vcf: None,
            output_bam: None, 
            quality_score_model: None,
            mutation_model: None,
            fragment_model: None,
            sequence_error_model: None,
        };
        config.update_and_log();
        config
    }

    pub fn from_yaml_file(yaml_file: String) -> RunConfiguration {
        // Reads an input configuration file from yaml using the serde package. Then sets the parameters
        // based on the inputs. A "." value means to use the default value.

        // Opens file for reading
        let f = fs::File::open(&yaml_file);
        let file = match f {
            Ok(l) => l,
            Err(error) => panic!(
                "Problem reading the config file: {}. System error: {}",
                &yaml_file,
                error,
            ),
        };
        // Uses serde_yaml to read the file into a HashMap
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)
            .expect("Could not read values");
        // create a default and update it
        let mut configuration = RunConfiguration::default();

        for (key, value) in scrape_config {
            match key.as_str() {
                // The reference is required and can't be skipped.
                "reference" => {
                    let reference_path = Path::new(value.as_str().unwrap());
                    if !reference_path.is_file() {
                        panic!("Reference file not found: {}", value.as_str().unwrap())
                    } else {
                        // Okay, serde Value is weird for strings
                        configuration.reference = reference_path.display().to_string();
                    }
                }
                _ => {
                    match &value.as_str() {
                        // Any item with a . for a value we keep the default
                        Some(".") => continue,
                        _ => match key.as_str() {
                            "read_len" => {
                                configuration.read_len = value
                                    .as_u64()
                                    .expect(&generate_error(&key, "integer", &value))
                                    as usize
                            }
                            "coverage" => {
                                configuration.coverage = value
                                    .as_u64()
                                    .expect(&generate_error(&key, "integer", &value))
                                    as usize
                            }
                            "mutation_rate" => {
                                configuration.mutation_rate = value
                                    .as_f64()
                                    .expect(&generate_error(&key, "float", &value))
                            }
                            "ploidy" => {
                                configuration.ploidy = value
                                    .as_u64()
                                    .expect(&generate_error(&key, "integer", &value))
                                    as usize
                            }
                            "paired_ended" => {
                                configuration.paired_ended = value
                                    .as_bool()
                                    .expect(&generate_error(&key, "boolean", &value))
                            }
                            "fragment_mean" => {
                                configuration.fragment_mean = value
                                    .as_f64()
                                    .expect(&generate_error(&key, "float", &value))
                                    .into() // to make it an option
                            }
                            "fragment_st_dev" => {
                                configuration.fragment_st_dev = value
                                    .as_f64()
                                    .expect(&generate_error(&key, "float", &value))
                                    .into() // to make it an option
                            }
                            "produce_fastq" => {
                                configuration.produce_fastq = value
                                    .as_bool()
                                    .expect(&generate_error(&key, "boolean", &value))
                            }
                            "produce_vcf" => {
                                configuration.produce_vcf = value
                                    .as_bool()
                                    .expect(&generate_error(&key, "boolean", &value))
                            }
                            "produce_bam" => {
                                configuration.produce_bam = value
                                    .as_bool()
                                    .expect(&generate_error(&key, "boolean", &value))
                            }
                            "rng_seed" => {
                                configuration.rng_seed = Some(value
                                    .as_str()
                                    .expect(&generate_error(&key, "string", &value))
                                    .to_string()
                                )
                            }
                            "overwrite_output" => {
                                configuration.overwrite_output = value
                                    .as_bool()
                                    .expect(&generate_error(&key, "boolean", &value))
                            }
                            "minimum_mutations" => {
                                configuration.minimum_mutations =
                                    value
                                        .as_u64()
                                        .expect(&generate_error(&key, "Valid integer", &value))
                                        .try_into()
                                        .unwrap();
                            }
                            "output_dir" => {
                                let output_path = value.as_str().unwrap().to_string();
                                if !Path::new(&output_path).is_dir() {
                                    generate_error(&key, "Path", &value);
                                }
                                configuration.output_dir = PathBuf::from(output_path);
                            }
                            "output_prefix" => {
                                configuration.output_prefix = value
                                    .as_str()
                                    .expect(&generate_error(&key, "String", &value))
                                    .to_string()
                            }
                            "quality_score_model" => {
                                let output_path = value
                                    .as_str()
                                    .expect("Value for quality_score_model must be a string")
                                    .to_string();
                                if !Path::new(&output_path).is_file() {
                                    generate_error(&key, "Path", &value);
                                }
                                configuration.quality_score_model = Some(PathBuf::from(output_path))
                            }
                            "mutation_model" => {
                                let output_path = value
                                    .as_str()
                                    .expect("Value for mutation model must be a string")
                                    .to_string();
                                if !Path::new(&output_path).is_file() {
                                    generate_error(&key, "Path", &value);
                                }
                                configuration.mutation_model = Some(PathBuf::from(output_path))
                            }
                            "fragment_model_datafile" => {
                                let output_path = value
                                    .as_str()
                                    .expect("Value for fragment lengt model must be a string")
                                    .to_string();
                                if !Path::new(&output_path).is_file() {
                                    generate_error(&key, "Path", &value);
                                }
                                configuration.mutation_model = Some(PathBuf::from(output_path))
                            },
                            "sequence_error_model" => {
                                let output_path = value
                                    .as_str()
                                    .expect("Value for mutation model must be a string")
                                    .to_string();
                                if !Path::new(&output_path).is_file() {
                                    generate_error(&key, "Path", &value);
                                }
                                configuration.mutation_model = Some(PathBuf::from(output_path))
                            },
                            _ => continue,
                        },
                    }
                }
            }
        }
        configuration.update_and_log();
        configuration
    }

    pub fn from_args(args: Cli) -> RunConfiguration {
        // Takes in a bunch of args from a clap CLI and builds a config based on that. More CLI options
        // will need additional items entered here. To add them to the config, so they can be implemented.

        // Create the ConfigBuilder object with default values
        let mut configuration = RunConfiguration::default();
        // Can't do a run without a reference
        if !args.reference.is_empty() {
            configuration.reference = args.reference.into();
        } else {
            panic!("No reference specified");
        }
        // The default value works directly for the config builder and CLI handles the type checking
        configuration.read_len = args.read_length;
        configuration.coverage = args.coverage;
        // default is empty string, in which case the config builder controls the default
        if args.output_dir == "" {
            configuration.output_dir = env::current_dir()
                .expect("Error finding current directory. Please specify --output-dir (-o) option.")
        } else {
            let output_path = Path::new(&args.output_dir);
            check_create_dir(output_path);
            configuration.output_dir = PathBuf::from(output_path);
        };
        // If this is unset, sets the default value of "neat_out" by CLI
        configuration.output_prefix = args.output_file_prefix;
        // To set a minimum mutation rate, such as for debugging, or for small datasets, use this
        if !args.minimum_mutations.is_none() {
            let input_min_muts: usize = args.minimum_mutations.unwrap() as usize;
            configuration.minimum_mutations = input_min_muts;
        }
        // Wraps things in a Box to move this object to the heap
        configuration.update_and_log();
        configuration
    }

    pub fn update_and_log(&mut self) {
        // This does a final check of the configuration for valid items. It will print info
        // message of the items, to work as a record and to assist in debugging any issues that
        // come up.
        if self.reference.is_empty() {
            panic!("No reference was specified.")
        }
        info!(
            "Running rusty-neat to generate reads on {} with...",
            &self.reference
        );
        info!("  >read length: {}", &self.read_len);
        info!("  >coverage: {}", &self.coverage);
        info!("  >mutation rate: {}", &self.mutation_rate);
        info!("  >ploidy: {}", &self.ploidy);
        info!("  >paired ended: {}", &self.paired_ended);
        if self.overwrite_output {
            warn!("Overwriting any existing files.")
        }
        if self.minimum_mutations > 0 {
            info!(
                "  >minimum mutations per contig: {}",
                &self.minimum_mutations
            )
        }
        let output_path = &self.output_dir;
        // This check may be overkill, but here it is. Let's make sure we ended up with something
        if !output_path.as_path().is_dir() {
            info!("Creating output directory: {:?}", self.output_dir);
            check_create_dir(output_path);
        }
        let file_prefix = format!("{:?}/{}", self.output_dir, self.output_prefix);

        // No point in running if we aren't producing files
        if !(self.produce_fastq | self.produce_vcf | self.produce_bam) {
            panic!("All file types set to false, no files would be produced.");
        }

        if self.paired_ended {
            if self.fragment_mean.is_none() | self.fragment_st_dev.is_none() {
                panic!(
                    "Paired ended is set to true, but fragment mean \
                    and standard deviation were not set."
                );
            }
            if self.produce_fastq {
                info!("\t> fragment mean: {}", self.fragment_mean.unwrap());
                info!(
                    "\t> fragment standard deviation: {}",
                    self.fragment_st_dev.unwrap()
                );
                info!("Producing fastq files:");

                let fastq_1 = format!("{}_r1.fastq", file_prefix);
                info!("\t> {}", &fastq_1);
                self.output_fastq_1 = Some(fastq_1);
                let fastq_2 = format!("{}_r2.fastq", file_prefix);
                info!("\t> {}", fastq_2);
                self.output_fastq_2 = Some(fastq_2)
            }
        } else {
            info!("Producing fastq file:");
            let fastq_1 = format!("{}_r1.fastq", file_prefix);
            info!("\t> {}", fastq_1);
            self.output_fastq_1 = Some(fastq_1);
        }
        if self.produce_vcf {
            let vcf: String = format!("{}.vcf", file_prefix);
            info!("Producing vcf file: {}.vcf", vcf);
            self.output_vcf = Some(vcf);
        }
        if self.produce_bam {
            let bam: String = format!("{}.bam", file_prefix);
            info!("Produce bam file: {}", bam);
            self.output_bam = Some(bam);
        }

        match &self.rng_seed {
            Some(seed) => {
                // User supplied seed
                // Convert the string to a string vec split at the white space.
                // Seeds can be any whitespace-separated series of strings.
                for seed_term in seed.split_whitespace() {
                    self.seed_vec.push(seed_term.to_string());
                }
                info!("Seed string to regenerate these exact results: {}", &seed);
            },

            None => {
                // Since no seed was provided, we'll use a datetime stamp with nanoseconds
                // The seed can be any space separated or tab separated series of strings
                // e.g., "Every good boy does Fine"
                // seeds are case-sensitive
                let raw_string = Utc::now().format("%Y %m %d %H %M %S %f").to_string();
                for item in raw_string.split_whitespace() {
                    self.seed_vec.push(item.to_string());
                }
                info!("Seed string to regenerate these exact results: {}", &raw_string);
            },
        }
    }
}

fn generate_error(key: &str, key_type: &str, value: &Value) -> String {
    format!(
        "Input {} could not be converted to {}: {:?}",
        key, key_type, value
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_run_configuration() {
        let test_configuration = RunConfiguration {
            reference: String::from("Hello.world"),
            read_len: 100,
            coverage: 22,
            mutation_rate: 0.09,
            ploidy: 3,
            paired_ended: true,
            fragment_mean: Option::from(333.0),
            fragment_st_dev: Option::from(33.0),
            produce_fastq: false,
            produce_bam: true,
            produce_vcf: true,
            rng_seed: None,
            seed_vec: Vec::new(),
            overwrite_output: true,
            minimum_mutations: 0,
            output_dir: PathBuf::from("/my/my"),
            output_prefix: String::from("Hey.hey"),
            output_fastq_1: None,
            output_fastq_2: None,
            output_vcf: None,
            output_bam: None, 
            quality_score_model: None,
            mutation_model: None,
            fragment_model: None,
            sequence_error_model: None,
        };

        println!("{:?}", test_configuration);
        assert_eq!(test_configuration.reference, "Hello.world".to_string());
        assert_eq!(test_configuration.read_len, 100);
        assert_eq!(test_configuration.coverage, 22);
        assert_eq!(test_configuration.mutation_rate, 0.09);
        assert_eq!(test_configuration.ploidy, 3);
        assert_eq!(test_configuration.paired_ended, true);
        assert_eq!(test_configuration.fragment_mean.unwrap(), 333.0);
        assert_eq!(test_configuration.fragment_st_dev.unwrap(), 33.0);
        assert_eq!(test_configuration.produce_fastq, false);
        assert_eq!(test_configuration.produce_vcf, true);
        assert_eq!(test_configuration.produce_bam, true);
        assert_eq!(test_configuration.rng_seed, None);
        assert_eq!(test_configuration.overwrite_output, true);
        assert_eq!(test_configuration.output_dir, PathBuf::from("/my/my"));
        assert_eq!(test_configuration.output_prefix, "Hey.hey".to_string());
    }

    #[test]
    fn test_build() {
        use super::*;
        let x = RunConfiguration::default();
        assert!(x.reference.is_empty())
    }

    #[test]
    fn test_read_config_yaml() {
        let yaml = String::from("test_data/configs/neat_test.yml");
        let test_config = RunConfiguration::from_yaml_file(yaml);
        assert_eq!(test_config.reference, "test_data/references/ecoli.fa".to_string());
        assert_eq!(test_config.coverage, 3);
    }

    #[test]
    #[should_panic]
    fn test_bad_yaml() {
        let yaml = String::from("test_data/fake_file.yml");
        RunConfiguration::from_yaml_file(yaml);
    }

    #[test]
    #[should_panic]
    fn test_missing_ref() {
        let yaml = String::from("test_data/configs/simple_template.yml");
        RunConfiguration::from_yaml_file(yaml);
    }

    #[test]
    fn test_creates_out_dir() {
        let yaml = String::from("test_data/configs/neat_test_bad.yml");
        RunConfiguration::from_yaml_file(yaml);
        assert!(Path::new("fake").is_dir());
        fs::remove_dir("fake").unwrap()
    }

    #[test]
    fn test_command_line_inputs() {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/references/ecoli.fa"),
            output_dir: String::from("data"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        let test_config = RunConfiguration::from_args(args);
        assert_eq!(test_config.reference, "test_data/references/ecoli.fa".to_string());
        fs::remove_dir("data").unwrap();
    }

    #[test]
    #[should_panic]
    fn test_bad_config_builder() {
        let mut config = RunConfiguration::default();
        config.update_and_log();
    }

    #[test]
    fn test_creat_nonexisting_out() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.output_dir = PathBuf::from("contig/");
        config.update_and_log();
        fs::remove_dir("contig").unwrap()
    }

    #[test]
    #[should_panic]
    fn test_cl_missing_ref() {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from(""),
            output_dir: String::from("test_dir"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        RunConfiguration::from_args(args);
        fs::remove_dir("test_dir").unwrap()
    }

    #[test]
    fn test_overwrite_warn() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.overwrite_output = true;
        config.update_and_log();
    }

    #[test]
    fn test_produce_fastq_messages() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.fragment_st_dev = Some(10.0);
        // tests first branch of if statement for paired_ended & produce_fastq = true
        config.update_and_log();
        // Checks the alternative pe = true, produce_fastq = false
        config.produce_fastq = false;
        // need to produce at least one file or check will panic
        config.update_and_log();
    }

    #[test]
    #[should_panic]
    fn test_no_files() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.produce_fastq = false;
        config.update_and_log();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean_or_stdev() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        // paired end set to true, by default, fragment mean and st dev are None
        config.paired_ended = true;
        config.update_and_log();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.paired_ended = true;
        config.fragment_st_dev = Some(10.0);
        config.update_and_log();
    }

    #[test]
    #[should_panic]
    fn test_no_stdev() {
        let mut config = RunConfiguration::default();
        config.reference = "test_data/references/H1N1.fa".to_string();
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.update_and_log();
    }

    #[test]
    fn no_output_dir_given() {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/references/H1N1.fa"),
            output_dir: String::new(),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        let config = RunConfiguration::from_args(args);
        assert_eq!(env::current_dir().unwrap().as_path(), config.output_dir);
    }

    #[test]
    fn test_minimum_mutations_and_others() {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/references/H1N1.fa"),
            output_dir: String::new(),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: Some(10),
            read_length: 120,
            coverage: 13,
        };

        let config = RunConfiguration::from_args(args);
        assert_eq!(10, config.minimum_mutations);
        assert_eq!(120, config.read_len);
        assert_eq!(13, config.coverage);
    }
}
