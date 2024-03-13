extern crate log;

use std::collections::{HashMap};
use std::env;
use std::string::String;
use crate::utils::cli::Cli;
use log::{debug, error};
use serde_yaml::Error;

#[derive(Debug)]
pub struct RunConfiguration {
    /*
    This struct holds all the parameters for this particular run. It is derived from input either
    from a configuration file or from command line inputs. This is not built directly in the code,
    but is constructed by a builder to enable default values.

    reference: The path to the reference for the run.
    read_len: The length of reads in the output fastq.
    coverage: The average depth of coverage for the output fastq file.
    mutation_rate: The rate of mutation for the file.
    ploidy: The number of copies of each chromosome in the target organism. Mutation process will
        be replicated this number of times.
    paired_ended: If the run will be in paired-ended mode.
    fragment_mean: Mean size of the fragments in paired-ended mode.
    fragment_st_dev: Standard deviation of the fragment mean describing the sample set to sample
        fragments from.
    produce_fastq: True or false on whether to produce an output fastq file.
    produce_fasta: True or false on whether to produce an output fasta file, 1 per ploid.
    produce_vcf: True or false on whether to produce an output VCF file, with genotyped variants.
    produce_bam: True or false on whether to produce an output BAM file, which will be aligned to
        the reference.
    output_dir: The directory, relative or absolute, path to the directory to place output.
    output_prefix: The name to use for the output files.
     */
    pub reference: String,
    pub read_len: usize,
    pub coverage: usize,
    pub mutation_rate: f64,
    pub ploidy: usize,
    pub paired_ended: bool,
    pub fragment_mean: f64,
    pub fragment_st_dev: f64,
    pub produce_fastq: bool,
    pub produce_fasta: bool,
    pub produce_vcf:  bool,
    pub produce_bam: bool,
    pub output_dir: String,
    pub output_prefix: String,
}

#[warn(dead_code)]
impl RunConfiguration {
    // The purpose of this function is to redirect you to the ConfigBuilder
    pub fn build() -> ConfigBuilder {
        ConfigBuilder::new()
    }
}

// The config builder allows us to construct a config in multiple different ways, depending
// on the input.
#[derive(Default)]
pub struct ConfigBuilder {
    reference: String,
    read_len: usize,
    coverage: usize,
    mutation_rate: f64,
    ploidy: usize,
    paired_ended: bool,
    fragment_mean: f64,
    fragment_st_dev: f64,
    produce_fastq: bool,
    produce_fasta: bool,
    produce_vcf:  bool,
    produce_bam: bool,
    output_dir: String,
    output_prefix: String,
}

#[derive(Debug)]
pub enum ConfigError {
    FileReadError,
    ConfigurationError,
}

impl ConfigBuilder {
    pub fn new() -> ConfigBuilder {
        let current_dir = env::current_dir()
            .expect("Could not find current working directory")
            .display()
            .to_string();
        ConfigBuilder {
            // Setting default values
            reference: String::from("data/H1N1.fa"),
            read_len: 150,
            coverage: 10,
            mutation_rate: 0.001,
            ploidy: 2,
            paired_ended: false,
            fragment_mean: 300.0,
            fragment_st_dev: 30.0,
            produce_fastq: true,
            produce_fasta: true,
            produce_vcf: false,
            produce_bam: false,
            output_dir: current_dir,
            output_prefix: String::from("neat_out"),
        }
    }

    // Basically a bunch of setters
    pub fn set_reference(mut self, reference: String) -> ConfigBuilder {
        self.reference = reference;
        self
    }

    pub fn set_read_len(mut self, read_len: usize) -> ConfigBuilder {
        self.read_len = read_len;
        self
    }

    pub fn set_coverage(mut self, coverage: usize) -> ConfigBuilder {
        self.coverage = coverage;
        self
    }

    pub fn set_mutation_rate(mut self, mutation_rate: f64) -> ConfigBuilder {
        self.mutation_rate = mutation_rate;
        self
    }

    pub fn set_ploidy(mut self, ploidy: usize) -> ConfigBuilder {
        self.ploidy = ploidy;
        self
    }

    pub fn set_paired_ended(mut self, paired_ended: bool) -> ConfigBuilder {
        self.paired_ended = paired_ended;
        self
    }

    pub fn set_fragment_mean(mut self, fragment_mean: f64) -> ConfigBuilder {
        self.fragment_mean = fragment_mean;
        self
    }

    pub fn set_fragment_st_dev(mut self, fragment_st_dev: f64) -> ConfigBuilder {
        self.fragment_st_dev = fragment_st_dev;
        self
    }

    pub fn set_produce_fastq(mut self, produce_fastq: bool) -> ConfigBuilder {
        self.produce_fastq = produce_fastq;
        self
    }

    pub fn set_produce_fasta(mut self, produce_fasta: bool) -> ConfigBuilder {
        self.produce_fasta = produce_fasta;
        self
    }

    pub fn set_produce_vcf(mut self, produce_vcf: bool) -> ConfigBuilder {
        self.produce_vcf = produce_vcf;
        self
    }

    pub fn set_produce_bam(mut self, produce_bam: bool) -> ConfigBuilder {
        self.produce_bam = produce_bam;
        self
    }

    pub fn set_output_dir(mut self, output_dir: String) -> ConfigBuilder {
        self.output_dir = String::from(output_dir);
        self
    }

    pub fn set_output_prefix(mut self, output_prefix: String) -> ConfigBuilder {
        self.output_prefix = String::from(output_prefix);
        self
    }

    pub fn check_and_print_config(&self) -> Result<(), ConfigError> {
        debug!("Running rusty-neat to generate reads on {} with...", self.reference);
        debug!("\t> read length: {}", self.read_len);
        debug!("\t> coverage: {}", self.coverage);
        debug!("\t> mutation rate: {}", self.mutation_rate);
        debug!("\t> ploidy: {}", self.ploidy);
        debug!("\t> paired_ended: {}", self.paired_ended);
        let file_prefix = format!("{}/{}", self.output_dir, self.output_prefix);
        if !(self.produce_fastq | self.produce_fasta | self.produce_vcf | self.produce_bam) {
            error!("All file types set to false, no files would be produced.");
            return Err(ConfigError::ConfigurationError);
        }
        if self.paired_ended {
            if self.fragment_mean.is_nan() | self.fragment_st_dev.is_nan() {
                error!(
                    "Paired ended is set to true, but fragment mean \
                    and standard deviation were not set."
                );
                return Err(ConfigError::ConfigurationError);
            }
            if self.produce_fastq {
                debug!("\t> fragment mean: {}", self.fragment_mean);
                debug!("\t> fragment standard deviation: {}", self.fragment_st_dev);
                debug!("Producing fastq files:\n\t> {}_r1.fastq\n\t {}_r2.fastq",
                    file_prefix, file_prefix
                )
            } else {
                debug!("Producing fastq file:\n\t> {}_r1.fastq", file_prefix)
            }
        }
        if self.produce_fasta {
            debug!("Producing fasta file: {}.fasta", file_prefix);
        }
        if self.produce_vcf {
            debug!("Producing vcf file: {}.vcf", file_prefix)
        }
        if self.produce_bam {
            debug!("Produce bam file: {}.bam", file_prefix)
        }
        Ok(())
    }

    // Function to build the actual configuration.
    pub fn build(self) -> RunConfiguration {
        RunConfiguration {
            reference: self.reference,
            read_len: self.read_len,
            coverage: self.coverage,
            mutation_rate: self.mutation_rate,
            ploidy: self.ploidy,
            paired_ended: self.paired_ended,
            fragment_mean: self.fragment_mean,
            fragment_st_dev: self.fragment_st_dev,
            produce_fastq: self.produce_fastq,
            produce_fasta: self.produce_fasta,
            produce_vcf: self.produce_vcf,
            produce_bam: self.produce_bam,
            output_dir: self.output_dir,
            output_prefix: self.output_prefix,
        }
    }
}

pub fn read_config_yaml(yaml: String) -> Result<Box<RunConfiguration>, &'static str> {
    /*
    Reads an input configuration file from yaml using the serde package. Then sets the parameters
    based on the inputs. A "." value means to use the default value.
     */

    // Opens file for reading
    let f = std::fs::File::open(yaml);
    let file = match f {
        Ok(l) => l,
        Err(_) => { return Err("problem reading configuration file"); },
    };
    // Uses serde_yaml to read the file into a HashMap
    let scrape_config: HashMap<String, String> = serde_yaml::from_reader(file)
        .expect("Could not read values");
    debug!("{:?}", scrape_config);
    // Create the config builder then update any items from the configuration file in the
    // configuration object and returns it.
    let mut config_builder = ConfigBuilder::new();
    for (key, value) in scrape_config {
        match key.as_str() {
            "reference" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_reference(value)
                }
            },
            "read_len" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_read_len(value.parse::<usize>().unwrap())
                }
            },
            "coverage" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_coverage(value.parse::<usize>().unwrap())
                }
            },
            "mutation_rate" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_mutation_rate(value.parse::<f64>().unwrap())
                }
            },
            "ploidy" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_ploidy(value.parse::<usize>().unwrap())
                }
            },
            "paired_ended" => {
                // since it is false by default, all we need to handle is the true case.
                if value.to_lowercase() == "true" {
                    config_builder = config_builder.set_paired_ended(true)
                }
            },
            "fragment_mean" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_fragment_mean(value.parse::<f64>().unwrap())
                }
            },
            "fragment_st_dev" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_fragment_st_dev(value.parse::<f64>().unwrap())
                }
            },
            "produce_fastq" => {
                // since produce_fastq is true by default, we only need to handle the false case
                if value.to_lowercase() == "false" {
                    config_builder = config_builder.set_produce_fastq(false)
                }
            },
            "produce_fasta" => {
                // since produce_fasta is false by default, we only need to handle the true case
                if value.to_lowercase() == "true" {
                    config_builder = config_builder.set_produce_fasta(true)
                }
            },
            "produce_vcf" => {
                // since produce_vcf is false by default, we only need to handle the true case
                if value.to_lowercase() == "true" {
                    config_builder = config_builder.set_produce_vcf(true)
                }
            },
            "produce_bam" => {
                // since produce_bam is false by default, we only need to handle the true case
                if value.to_lowercase() == "true" {
                    config_builder = config_builder.set_produce_bam(true)
                }
            },
            "output_dir" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_output_dir(value)
                }
            },
            "output_prefix" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_output_prefix(value)
                }
            },
            _ => continue,
        }
    }
    let _ = &config_builder.check_and_print_config().unwrap();
    Ok(Box::new(config_builder.build()))
}

pub fn build_config_from_args(args: Cli) -> Result<Box<RunConfiguration>, &'static str> {
    /*
    Takes in a bunch of args from a clap CLI and builds a config based on that. More CLI options
    will need additional items entered here. To add them to the config, so they can be implemented.
     */

    // Create the ConfigBuilder object with default values
    let mut config_builder = ConfigBuilder::new();
    // The default assigned by CLI will keep it as data/H1N1.fa
    if args.reference != "" {
        config_builder = config_builder.set_reference(args.reference)
    }
    // The default value works directly for the config builder
    config_builder = config_builder.set_read_len(args.read_length);
    config_builder = config_builder.set_coverage(args.coverage);
    // default is empty string, in which case the config builder controls the default
    if args.output_dir != "" {
        config_builder = config_builder.set_output_dir(args.output_dir)
    };
    // If this is unset, sets the default value of neat_out
    if args.output_file_prefix != "" {
        config_builder = config_builder.set_output_prefix(args.output_file_prefix)
    };
    // Wraps things in a Box to move this object to the heap
    let _ = &config_builder.check_and_print_config().unwrap();
    Ok(Box::new(config_builder.build()))
}

#[test]
fn test_read_config_yaml() {
    let yaml = String::from("config/neat_test.yml");
    let test_config = read_config_yaml(yaml).unwrap();
    assert_eq!(test_config.reference, "data/ecoli.fa".to_string());
    assert_eq!(test_config.coverage, 10);
}

#[test]
fn test_command_line_inputs() {
    let args: Cli = Cli{
        config: String::new(),
        reference: String::from("data/ecoli.fa"),
        output_dir: String::from("test_dir"),
        output_file_prefix: String::from("test"),
        read_length: 150,
        coverage: 10,
    };

    let test_config = build_config_from_args(args).unwrap();
    assert_eq!(test_config.reference, "data/ecoli.fa".to_string())
}