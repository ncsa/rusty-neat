use serde::{Serialize, Deserialize};
use std::collections::{HashMap};
use std::env;
use std::string::String;
use crate::utils::cli::Cli;

#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    pub reference: String,
    pub read_len: u32,
    pub coverage: u32,
    pub mutation_rate: f64,
    pub ploidy: u32,
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

impl Config {
    pub fn build() -> ConfigBuilder {
        ConfigBuilder::default()
    }
}

#[derive(Default)]
struct ConfigBuilder {
    reference: String,
    read_len: u32,
    coverage: u32,
    mutation_rate: f64,
    ploidy: u32,
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

impl ConfigBuilder {
    pub fn new() -> ConfigBuilder {
        let current_dir = env::current_dir()
            .expect("Could not find current working directory")
            .display()
            .to_string();
        ConfigBuilder {
            reference: String::from("data/H1N1.fa"),
            read_len: 150,
            coverage: 10,
            mutation_rate: 0.001,
            ploidy: 2,
            paired_ended: false,
            fragment_mean: 300.0,
            fragment_st_dev: 30.0,
            produce_fastq: true,
            produce_fasta: false,
            produce_vcf: false,
            produce_bam: false,
            output_dir: current_dir,
            output_prefix: String::from("neat_out"),
        }
    }

    pub fn set_ref(mut self, reference: String) -> ConfigBuilder {
        self.reference = reference;
        self
    }

    pub fn set_rl(mut self, read_len: u32) -> ConfigBuilder {
        self.read_len = read_len;
        self
    }

    pub fn set_cov(mut self, coverage: u32) -> ConfigBuilder {
        self.coverage = coverage;
        self
    }

    pub fn set_outdir(mut self, output_dir: String) -> ConfigBuilder {
        self.output_dir = output_dir;
        self
    }

    pub fn set_outpref(mut self, output_prefix: String) -> ConfigBuilder {
        self.output_prefix = output_prefix;
        self
    }

    pub fn build(self) -> Config {
        Config {
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

pub fn read_config_yaml(yaml: String) -> Result<Config, serde_yaml::Error> {
    let f = std::fs::File::open(yaml).expect("Could not open config file.");
    let scrape_config: HashMap<String, String> = serde_yaml::from_reader(f).expect("Could not read values");
    println!("{:?}", scrape_config);
    let mut config_builder = ConfigBuilder::new();
    for (key, value) in scrape_config {
        match key.as_str() {
            "reference" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_ref(value)
                } else {
                    continue
                }
            },
            "read_len" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_rl(value.parse::<u32>().unwrap())
                } else {
                    continue
                }
            },
            "coverage" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_cov(value.parse::<u32>().unwrap())
                } else {
                    continue
                }
            },
            "output_dir" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_outdir(value)
                } else {
                    continue
                }
            },
            "output_prefix" => {
                if value != ".".to_string() {
                    config_builder = config_builder.set_outpref(value)
                } else {
                    continue
                }
            },
            _ => continue,
        }
    }
    Ok(config_builder.build())
}

pub fn build_config_from_args(args: Cli) -> Result<Config, serde_yaml::Error> {
    let mut config_builder = ConfigBuilder::new();
    if args.reference != "" { config_builder = config_builder.set_ref(args.reference) }
    config_builder = config_builder.set_rl(args.read_length);
    config_builder = config_builder.set_cov(args.coverage);
    if args.output_dir != "" { config_builder = config_builder.set_outdir(args.output_dir); };
    if args.output_file_prefix != "" { config_builder = config_builder.set_outpref(args.output_file_prefix) };
    Ok(config_builder.build())
}

#[test]
fn test_read_config_yaml() {
    let yaml = String::from("config/neat_test.yml");
    let test_config = read_config_yaml(yaml).unwrap();
    assert_eq!(test_config.reference, "data/ecoli.fa".to_string());
    assert_eq!(test_config.coverage, 10);
}