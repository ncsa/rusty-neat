// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use chrono::Utc;

use log::{info, warn, error};
use serde_yml::Value;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::string::String;
use std::{env, fs};
use crate::{
    common::file_tools::folder_tools::check_create_dir,
    utils::cli::Cli
};
use crate::errors::GenerateReadsErrors;

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
    pub output_filename: String,
    pub output_fastq_1: Option<PathBuf>,
    pub output_fastq_2: Option<PathBuf>,
    pub output_vcf: Option<PathBuf>,
    pub output_bam: Option<PathBuf>, 
    // model input
    pub quality_score_model: Option<PathBuf>,
    pub mutation_model: Option<PathBuf>,
    pub fragment_model: Option<PathBuf>,
    pub sequence_error_model: Option<PathBuf>,
}

// The config builder allows us to construct a config in multiple different ways, depending
// on the input.
pub struct ConfigBuilder {
    // The point of this is so we don't have to check for the reference to exist in the rest of the code
    // We set it as an option in this, then run some checks, for yaml v command line inputs, then when
    // we have the reference, we build the real config object, which requires a reference. Any other arguments
    // can be made required this way. It gives us some flexibility in the CLI v config. Probably this was more
    // ambitious than it needed to be, though.
    pub(crate) reference: Option<String>,
}

impl ConfigBuilder {
    pub fn new() -> ConfigBuilder {
        ConfigBuilder {
            reference: None,
        }
    }

    // Function to build the actual configuration.
    pub fn build(self) -> RunConfiguration {
        match self.reference {
            Some(reference) => {
                info!("Running rusty-neat to generate reads on {:?} with...", &reference);

                RunConfiguration {
                    reference: reference,
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
                    output_filename: "neat_out".to_string(),
                    output_fastq_1: None,
                    output_fastq_2: None,
                    output_vcf: None,
                    output_bam: None, 
                    quality_score_model: None,
                    mutation_model: None,
                    fragment_model: None,
                    sequence_error_model: None,
                }
            },
            None => {
                panic!("Must specify a reference with -r or -c with a valid config file")
            },
        }
    }
}

impl RunConfiguration {

    #[allow(dead_code)]
    // The purpose of this function is to redirect you to the ConfigBuilder
    pub fn build() -> ConfigBuilder {
        ConfigBuilder::new()
    }

    pub fn from_yaml_file(yaml_file: String) -> Result<RunConfiguration, GenerateReadsErrors> {
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
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)?;
        // create a default and update it
        let mut config_builder = ConfigBuilder::new();
        // Fill in the reference first, all hinges on that
        let reference = scrape_config["reference"].clone();
        let ref_str = reference.as_str();
        match ref_str {
            Some(str) => {
                let ref_path = Path::new(str);
                if ref_path.is_file() {
                    config_builder.reference = Some(ref_path.display().to_string());
                } else {
                    return Err(GenerateReadsErrors::FileNotFound(str.to_owned()))
                }
            },
            None => {
                return Err(GenerateReadsErrors::MissingReferenceError)
            }
        }
        // Creates a default configuration with the reference value, which we can then update.
        let mut configuration = config_builder.build();
        for (key, value) in scrape_config {
            match &value.as_str() {
                // Any item with a . for a value we keep the default
                Some(".") => continue,
                _ => match key.as_str() {
                    "read_len" => {
                        let rl_option = value.as_u64();
                        match rl_option {
                            Some(rl_opt) => {
                                configuration.read_len = rl_opt as usize;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("read_len".to_string(), "integer".to_string()))
                            }
                        }
                    },
                    "coverage" => {
                        let cov_option = value.as_u64();
                        match cov_option {
                            Some(cov_opt) => {
                                configuration.coverage = cov_opt as usize;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("coverage".to_string(), "integer".to_string()))
                            }
                        }
                    },
                    "mutation_rate" => {
                        let mr_option = value.as_f64();
                        match mr_option {
                            Some(mr_opt) => {
                                configuration.mutation_rate = mr_opt;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("mutation_rate".to_string(), "float".to_string()))
                            }
                        }
                    },
                    "ploidy" => {
                        let pl_option = value.as_u64();
                        match pl_option {
                            Some(pl_opt) => {
                                configuration.ploidy = pl_opt as usize;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("ploidy".to_string(), "integer".to_string()))
                            }
                        }
                    },
                    "paired_ended" => {
                        let pe_option = value.as_bool();
                        match pe_option {
                            Some(pe_opt) => {
                                configuration.paired_ended = pe_opt;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("paired_ended".to_string(), "boolean".to_string()))
                            }
                        }
                    },
                    "fragment_mean" => {
                        let fm_option = value.as_f64();
                        match fm_option {
                            Some(fm_opt) => {
                                configuration.fragment_mean = Some(fm_opt);
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("fragment_mean".to_string(), "float".to_string()))
                            }
                        }
                    },
                    "fragment_st_dev" => {
                        let fsd_option = value.as_f64();
                        match fsd_option {
                            Some(fsd_opt) => {
                                configuration.fragment_st_dev = Some(fsd_opt);
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("frament_st_dev".to_string(), "float".to_string()))
                            }
                        }
                    },
                    "produce_fastq" => {
                        let pfq_option = value.as_bool();
                        match pfq_option {
                            Some(pfq_opt) => {
                                configuration.produce_fastq = pfq_opt;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("produce_fastq".to_string(), "boolean".to_string()))
                            }
                        }
                    },
                    "produce_vcf" => {
                        let pvc_option = value.as_bool();
                        match pvc_option {
                            Some(pvc_opt) => {
                                configuration.produce_vcf = pvc_opt;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("produce_vcf".to_string(), "boolean".to_string()))
                            }
                        }
                    },
                    "produce_bam" => {
                        let pb_option = value.as_bool();
                        match pb_option {
                            Some(pb_opt) => {
                                configuration.produce_bam = pb_opt;
                                warn!("Produce BAM is currently not functional, but will be implemented soon!")
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("produce_bam".to_string(), "boolean".to_string()))
                            }
                        }
                    },
                    "rng_seed" => {
                        let rs_option = value.as_str();
                        match rs_option {
                            Some(rs_opt) => {
                                configuration.rng_seed = Some(rs_opt.to_string());
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("rng_seed".to_string(), "String".to_string()))
                            }
                        }
                    },
                    "overwrite_output" => {
                        let ow_option = value.as_bool();
                        match ow_option {
                            Some(ow_opt) => {
                                configuration.overwrite_output = ow_opt;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("overwrite_output".to_string(), "boolean".to_string()))
                            }
                        }
                    },
                    "minimum_mutations" => {
                        let mm_option = value.as_u64();
                        match mm_option {
                            Some(mm_opt) => {
                                configuration.minimum_mutations = mm_opt as usize;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("minimum_mutations".to_string(), "Integer".to_string()))
                            }
                        }
                    },
                    "output_dir" => {
                        let output_option = value.as_str();
                        match output_option {
                            Some(output_path) => {
                                let path = PathBuf::from(output_path);
                                configuration.output_dir = path;
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("output_dir".to_string(), "String".to_string()))
                            }
                        }
                    },
                    // output_prefix is for backward compatability
                    "output_filename" | "output_prefix "=> {
                        let pref_val = value.as_str();
                        match pref_val {
                            Some(name) => {
                                configuration.output_filename = name.to_string();
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("output_filename".to_string(), "String".to_string()))
                            }
                        }
                    },
                    "quality_score_model" => {
                        let qs_val = value.as_str();
                        match qs_val {
                            Some(name) => {
                                let filename = PathBuf::from(name);
                                if filename.is_file() {
                                    configuration.quality_score_model = Some(filename);
                                } else {
                                    return Err(GenerateReadsErrors::FileNotFound(name.to_string()))
                                }
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("quality_score_model".to_string(), "path to file".to_string()))
                            }
                        }
                    },
                    "mutation_model" => {
                        let mm_val = value.as_str();
                        match mm_val {
                            Some(name) => {
                                let filename = PathBuf::from(name);
                                if filename.is_file() {
                                    configuration.mutation_model = Some(filename);
                                } else {
                                    return Err(GenerateReadsErrors::FileNotFound(name.to_string()))
                                }
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("mutation_model".to_string(), "path to file".to_string()))
                            }
                        }
                    },
                    "fragment_model" => {
                        let fr_val = value.as_str();
                        match fr_val {
                            Some(name) => {
                                let filename = PathBuf::from(name);
                                if filename.is_file() {
                                    configuration.fragment_model = Some(filename);
                                } else {
                                    return Err(GenerateReadsErrors::FileNotFound(name.to_string()))
                                }
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("fragment_model".to_string(), "path to file".to_string()))
                            }
                        }
                    },
                    "sequence_error_model" => {
                        let se_val = value.as_str();
                        match se_val {
                            Some(name) => {
                                let filename = PathBuf::from(name);
                                if filename.is_file() {
                                    configuration.sequence_error_model = Some(filename);
                                } else {
                                    return Err(GenerateReadsErrors::FileNotFound(name.to_string()))
                                }
                            },
                            None => {
                                return Err(GenerateReadsErrors::ConfigReadError("sequence_error_model".to_string(), "path to file".to_string()))
                            }
                        }
                    },
                    _ => continue,
                },
            }
        }
        configuration.update_and_log()?;
        Ok(configuration)
    }

    pub fn from_args(args: Cli) -> Result<RunConfiguration, GenerateReadsErrors> {
        // Takes in a bunch of args from a clap CLI and builds a config based on that. More CLI options
        // will need additional items entered here. To add them to the config, so they can be implemented.

        // Create the ConfigBuilder object with default values
        let mut config_builder = ConfigBuilder::new();
        // Can't do a run without a reference
        if !args.reference.is_empty() {
            let reference: String = args.reference.into();
            let ref_path = Path::new(&reference);
            if ref_path.is_file() {
                config_builder.reference = Some(ref_path.display().to_string());
            } else {
                return Err(GenerateReadsErrors::FileNotFound(reference))
            }
        } else {
            return Err(GenerateReadsErrors::MissingReferenceError);
        }
        // We confirmed the reference is a valid file, proceed with the rest of the CLI
        let mut configuration = config_builder.build();

        // The default value works directly for the config builder and CLI handles the type checking
        configuration.read_len = args.read_length;
        configuration.coverage = args.coverage;
        // default is empty string, in which case the config builder controls the default
        if args.output_dir == "" {
            configuration.output_dir = env::current_dir()?
        } else {
            let output_path = Path::new(&args.output_dir);
            check_create_dir(output_path);
            configuration.output_dir = PathBuf::from(output_path);
        };
        // If this is unset, sets the default value of "neat_out" by CLI
        configuration.output_filename = args.output_filename;
        // Wraps things in a Box to move this object to the heap
        configuration.update_and_log()?;
        Ok(configuration)
    }

    pub fn update_and_log(&mut self) -> Result<(), GenerateReadsErrors> {
        // This does a final check of the configuration for valid items. It will print info
        // message of the items, to work as a record and to assist in debugging any issues that
        // come up.
        if self.reference.is_empty() {
            return Err(GenerateReadsErrors::MissingReferenceError)
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

        // No point in running if we aren't producing files
        if !(self.produce_fastq | self.produce_vcf | self.produce_bam) {
            error!("All file types set to false, no files would be produced.");
            return Err(GenerateReadsErrors::MissingReferenceError)
        }

        if self.paired_ended {
            if self.fragment_mean.is_none() | self.fragment_st_dev.is_none() {
                error!(
                    "Paired ended is set to true, but fragment mean \
                    and standard deviation were not set."
                );
                return Err(GenerateReadsErrors::MissingReferenceError)
            }
            if self.produce_fastq {
                info!("\t> fragment mean: {}", self.fragment_mean.unwrap());
                info!(
                    "\t> fragment standard deviation: {}",
                    self.fragment_st_dev.unwrap()
                );
                info!("Producing fastq files:");

                let mut fastq_1 = self.output_dir.clone();
                fastq_1.push(format!("{}_r1.fastq.gz", self.output_filename));
                info!("\t> {:?}", &fastq_1);
                self.output_fastq_1 = Some(fastq_1);
                let mut fastq_2 = self.output_dir.clone();
                fastq_2.push(format!("{}_r2.fastq.gz", self.output_filename));
                info!("\t> {:?}", fastq_2);
                self.output_fastq_2 = Some(fastq_2);
            }
        } else {
            info!("Producing fastq file:");
            let mut fastq_1 = self.output_dir.clone();
            fastq_1.push(format!("{}_r1.fastq.gz", self.output_filename));
            info!("\t> {:?}", &fastq_1);
            self.output_fastq_1 = Some(fastq_1);
        }
        if self.produce_vcf {
            let mut vcf = self.output_dir.clone();
            vcf.push(format!("{}.vcf", self.output_filename));
            info!("\t> {:?}", &vcf);
            self.output_vcf = Some(vcf);
        }
        if self.produce_bam {
            let mut bam = self.output_dir.clone();
            bam.push(format!("{}.bam", self.output_filename));
            info!("\t> {:?}", &bam);
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
                Ok(())
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
                Ok(())
            },
        }
    }
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
            output_filename: String::from("Hey.hey"),
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
        assert_eq!(test_configuration.output_filename, "Hey.hey".to_string());
    }

    #[test]
    fn test_build() {
        use super::*;
        let mut x = ConfigBuilder::new();
        x.reference = Some("test_data/H1N1.fa".to_string());
        let config = x.build();
        assert_eq!(config.reference, "test_data/H1N1.fa".to_string())
    }

    #[test]
    fn test_read_config_yaml() {
        let yaml = String::from("test_data/configs/neat_test.yml");
        let test_config = RunConfiguration::from_yaml_file(yaml).unwrap();
        assert_eq!(test_config.reference, "test_data/references/ecoli.fa".to_string());
        assert_eq!(test_config.coverage, 3);
    }

    #[test]
    #[should_panic]
    fn test_bad_yaml() {
        let yaml = String::from("test_data/fake_file.yml");
        RunConfiguration::from_yaml_file(yaml).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_missing_ref() {
        let yaml = String::from("test_data/configs/simple_template.yml");
        RunConfiguration::from_yaml_file(yaml).unwrap();
    }

    #[test]
    fn test_creates_out_dir() {
        let yaml = String::from("test_data/configs/neat_test_bad.yml");
        RunConfiguration::from_yaml_file(yaml).unwrap();
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
            output_filename: String::from("test"),
            read_length: 150,
            coverage: 10,
        };

        let test_config = RunConfiguration::from_args(args).unwrap();
        assert_eq!(test_config.reference, "test_data/references/ecoli.fa".to_string());
        fs::remove_dir("data").unwrap();
    }

    #[test]
    fn test_no_config_builder() {
        let config = RunConfiguration::build();
        assert_eq!(config.reference, None)
    }

    #[test]
    fn test_creat_nonexisting_out() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.output_dir = PathBuf::from("contig/");
        config.update_and_log().unwrap();
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
            output_filename: String::from("test"),
            read_length: 150,
            coverage: 10,
        };

        RunConfiguration::from_args(args).unwrap();
        fs::remove_dir("test_dir").unwrap()
    }

    #[test]
    fn test_overwrite_warn() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.overwrite_output = true;
        config.update_and_log().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_produce_fastq_messages() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.fragment_st_dev = Some(10.0);
        // tests first branch of if statement for paired_ended & produce_fastq = true
        config.update_and_log().unwrap();
        // Checks the alternative pe = true, produce_fastq = false
        config.produce_fastq = false;
        // need to produce at least one file or check will panic
        config.update_and_log().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_no_files() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.produce_fastq = false;
        config.update_and_log().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean_or_stdev() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        // paired end set to true, by default, fragment mean and st dev are None
        config.paired_ended = true;
        config.update_and_log().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.paired_ended = true;
        config.fragment_st_dev = Some(10.0);
        config.update_and_log().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_no_stdev() {
        let mut config_builder = ConfigBuilder::new();
        config_builder.reference = Some("test_data/references/H1N1.fa".to_string());
        let mut config = config_builder.build();
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.update_and_log().unwrap();
    }

    #[test]
    fn no_output_dir_given() {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/references/H1N1.fa"),
            output_dir: String::new(),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_filename: String::from("test"),
            read_length: 150,
            coverage: 10,
        };

        let config = RunConfiguration::from_args(args).unwrap();
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
            output_filename: String::from("test"),
            read_length: 120,
            coverage: 13,
        };

        let config = RunConfiguration::from_args(args).unwrap();
        assert_eq!(120, config.read_len);
        assert_eq!(13, config.coverage);
    }
}
