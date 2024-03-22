// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.

use std::collections::{HashMap};
use std::string::String;
use utils::cli::Cli;
use log::{debug, warn};
use std::env;
use std::path::Path;

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
    overwrite_output: if true, will overwrite output. If false will error and exit you attempt to
        overwrite files with the same name.
    output_dir: The directory, relative or absolute, path to the directory to place output.
    output_prefix: The name to use for the output files.
     */
    pub reference: String,
    pub read_len: usize,
    pub coverage: usize,
    pub mutation_rate: f64,
    pub ploidy: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    pub produce_fastq: bool,
    pub produce_fasta: bool,
    pub produce_vcf:  bool,
    pub produce_bam: bool,
    pub rng_seed: Option<u64>,
    pub overwrite_output: bool,
    pub output_dir: String,
    pub output_prefix: String,
}
#[allow(dead_code)]
impl RunConfiguration {
    // The purpose of this function is to redirect you to the ConfigBuilder
    pub fn build() -> ConfigBuilder {
        ConfigBuilder::new()
    }
}

// The config builder allows us to construct a config in multiple different ways, depending
// on the input.
pub struct ConfigBuilder {
    reference: Option<String>,
    read_len: usize,
    coverage: usize,
    mutation_rate: f64,
    ploidy: usize,
    paired_ended: bool,
    fragment_mean: Option<f64>,
    fragment_st_dev: Option<f64>,
    produce_fastq: bool,
    produce_fasta: bool,
    produce_vcf:  bool,
    produce_bam: bool,
    rng_seed: Option<u64>,
    overwrite_output: bool,
    output_dir: String,
    output_prefix: String,
}

impl ConfigBuilder {
    pub fn new() -> ConfigBuilder {
        ConfigBuilder {
            // Setting default values
            reference: None,
            read_len: 150,
            coverage: 10,
            mutation_rate: 0.001,
            ploidy: 2,
            paired_ended: false,
            fragment_mean: None,
            fragment_st_dev: None,
            produce_fastq: true,
            produce_fasta: false,
            produce_vcf: false,
            produce_bam: false,
            rng_seed: None,
            overwrite_output: false,
            output_dir: String::from(env::current_dir().unwrap().to_str().unwrap()),
            output_prefix: String::from("neat_out"),
        }
    }

    // Basically a bunch of setters
    pub fn set_reference(mut self, reference: String) -> ConfigBuilder {
        self.reference = Option::from(reference);
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
        self.fragment_mean = Option::from(fragment_mean);
        self
    }

    pub fn set_fragment_st_dev(mut self, fragment_st_dev: f64) -> ConfigBuilder {
        self.fragment_st_dev = Option::from(fragment_st_dev);
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

    pub fn set_rng_seed(mut self, rng_seed: u64) -> ConfigBuilder {
        self.rng_seed = Option::from(rng_seed);
        self
    }

    pub fn set_overwrite_output(mut self) -> ConfigBuilder {
        self.overwrite_output = true;
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

    pub fn check_and_print_config(&self) {
        if self.reference.is_none() {
            panic!("No reference was specified.")
        }
        debug!("Running rusty-neat to generate reads on {} with...", self.reference.clone().unwrap());
        debug!("  >read length: {}", self.read_len);
        debug!("  >coverage: {}", self.coverage);
        debug!("  >mutation rate: {}", self.mutation_rate);
        debug!("  >ploidy: {}", self.ploidy);
        debug!("  >paired_ended: {}", self.paired_ended);
        if self.overwrite_output {
            warn!("Overwriting any existing files.")
        }
        let file_prefix = format!("{}/{}", self.output_dir, self.output_prefix);
        if !(self.produce_fastq | self.produce_fasta | self.produce_vcf | self.produce_bam) {
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
                debug!("\t> fragment mean: {}", self.fragment_mean.unwrap());
                debug!("\t> fragment standard deviation: {}", self.fragment_st_dev.unwrap());
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
        if self.rng_seed.is_some() {
            debug!("Using rng seed: {}", self.rng_seed.unwrap())
        }
    }

    // Function to build the actual configuration.
    pub fn build(self) -> RunConfiguration {
        RunConfiguration {
            reference: self.reference.unwrap(),
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
            rng_seed: self.rng_seed,
            overwrite_output: self.overwrite_output,
            output_dir: self.output_dir,
            output_prefix: self.output_prefix,
        }
    }
}

pub fn read_config_yaml(yaml: String) -> Box<RunConfiguration> {
    /*
    Reads an input configuration file from yaml using the serde package. Then sets the parameters
    based on the inputs. A "." value means to use the default value.
     */

    // Opens file for reading
    let f = std::fs::File::open(yaml);
    let file = match f {
        Ok(l) => l,
        Err(error) => panic!("Problem reading the cofig file: {}", error),
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
                if value != ".".to_string() && value != "".to_string() {
                    let reference_path = Path::new(&value);
                    if !reference_path.exists() {
                        panic!("Reference file not found: {}", value)
                    } else {
                        config_builder = config_builder.set_reference(value);
                    }
                } else {
                    panic!("Reference was not specified in config.");
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
            "rng_seed" => {
                // Optional rng seed to replicate results
                if value != ".".to_string() {
                    config_builder = config_builder.set_rng_seed(value.parse::<u64>().unwrap())
                }
            },
            "overwrite_output" => {
                // overwrite_output is false by default, to prevent data loss, setting this flag
                // will instead enable rusty-neat to overwrite the output
                if value.to_lowercase() == "true" {
                    config_builder = config_builder.set_overwrite_output()
                }
            }
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
    let _ = &config_builder.check_and_print_config();
    Box::new(config_builder.build())
}

pub fn build_config_from_args(args: Cli) -> Result<Box<RunConfiguration>, &'static str> {
    /*
    Takes in a bunch of args from a clap CLI and builds a config based on that. More CLI options
    will need additional items entered here. To add them to the config, so they can be implemented.
     */

    // Create the ConfigBuilder object with default values
    let mut config_builder = ConfigBuilder::new();
    // Can't do a run without a reference
    if args.reference != "" {
        config_builder = config_builder.set_reference(args.reference)
    } else {
        panic!("No reference specified");
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
    let _ = &config_builder.check_and_print_config();
    Ok(Box::new(config_builder.build()))
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
            produce_fasta: true,
            produce_vcf: true,
            rng_seed: None,
            overwrite_output: true,
            output_dir: String::from("/my/my"),
            output_prefix: String::from("Hey.hey")
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
        assert_eq!(test_configuration.produce_fasta, true);
        assert_eq!(test_configuration.rng_seed, None);
        assert_eq!(test_configuration.overwrite_output, true);
        assert_eq!(test_configuration.output_dir, "/my/my".to_string());
        assert_eq!(test_configuration.output_prefix, "Hey.hey".to_string());
    }

    #[test]
    fn test_build() {
        use super::*;
        let x = RunConfiguration::build();
        assert!(x.reference.is_none())
    }

    #[test]
    fn test_read_config_yaml() {
        let yaml = String::from("config/neat_test.yml");
        let test_config = read_config_yaml(yaml);
        assert_eq!(test_config.reference, "data/ecoli.fa".to_string());
        assert_eq!(test_config.coverage, 3);
    }

    #[test]
    #[should_panic]
    fn test_bad_yaml() {
        let yaml = String::from("fake_file.yml");
        read_config_yaml(yaml);
    }

    #[test]
    #[should_panic]
    fn test_missing_ref() {
        let yaml = String::from("config/simple_template.yml");
        read_config_yaml(yaml);
    }

    #[test]
    fn test_setters() {
        let mut config = ConfigBuilder::new();
        assert_eq!(config.reference, None);
        assert_eq!(config.mutation_rate, 0.001);
        assert_eq!(config.ploidy, 2);
        assert_eq!(config.paired_ended, false);
        assert_eq!(config.fragment_mean, None);
        assert_eq!(config.fragment_st_dev, None);
        assert_eq!(config.produce_fasta, false);
        assert_eq!(config.produce_fastq, true);
        assert_eq!(config.produce_vcf, false);
        assert_eq!(config.produce_bam, false);
        assert_eq!(config.rng_seed, None);
        config = config.set_reference("data/H1N1.fa".to_string())
            .set_mutation_rate(0.111)
            .set_ploidy(3)
            .set_paired_ended(true)
            .set_fragment_mean(111.0)
            .set_fragment_st_dev(0.011)
            .set_produce_fastq(false)
            .set_produce_fasta(true)
            .set_produce_vcf(true)
            .set_produce_bam(true)
            .set_rng_seed(11393)
            .set_overwrite_output();
        assert_eq!(config.reference, Some("data/H1N1.fa".to_string()));
        assert_eq!(config.mutation_rate, 0.111);
        assert_eq!(config.ploidy, 3);
        assert_eq!(config.paired_ended, true);
        assert_eq!(config.fragment_mean, Some(111.0));
        assert_eq!(config.fragment_st_dev, Some(0.011));
        assert_eq!(config.produce_fasta, true);
        assert_eq!(config.produce_fastq, false);
        assert_eq!(config.produce_vcf, true);
        assert_eq!(config.produce_bam, true);
        assert_eq!(config.rng_seed, Some(11393));
        config.check_and_print_config()
    }

    #[test]
    fn test_command_line_inputs() {
        let args: Cli = Cli{
            config: String::new(),
            reference: String::from("data/ecoli.fa"),
            output_dir: String::from("test_dir"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            read_length: 150,
            coverage: 10,
        };

        let test_config = build_config_from_args(args).unwrap();
        assert_eq!(test_config.reference, "data/ecoli.fa".to_string())
    }

    #[test]
    #[should_panic]
    fn test_cl_missing_ref() {
        let args: Cli = Cli{
            config: String::new(),
            reference: String::from(""),
            output_dir: String::from("test_dir"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            read_length: 150,
            coverage: 10,
        };

        build_config_from_args(args).unwrap();
    }

    #[test]
    fn test_overwrite_warn() {
        let mut config = ConfigBuilder::new();
        config = config.set_reference("data/H1N1.fa".to_string());
        config = config.set_overwrite_output();
        config.check_and_print_config();
    }

    #[test]
    fn test_produce_fastq_messages() {
        let mut config = ConfigBuilder::new();
        config = config.set_reference("data/H1N1.fa".to_string());
        config = config.set_paired_ended(true);
        config = config.set_fragment_mean(100.0);
        config = config.set_fragment_st_dev(10.0);
        // tests first branch of if statement for paired_ended & produce_fastq = true
        config.check_and_print_config();
        // Checks the alternative pe = true, produce_fastq = false
        config = config.set_produce_fastq(false);
        // need to produce at least one file or check will panic
        config = config.set_produce_fasta(true);
        config.check_and_print_config();
    }

    #[test]
    fn test_produce_fasta_messages() {
        let mut config = ConfigBuilder::new();
        config = config.set_reference("data/H1N1.fa".to_string());
        config = config.set_produce_fasta(true);
        config.check_and_print_config();
        config = config.set_produce_vcf(true);
        config.check_and_print_config();
        config = config.set_produce_bam(true);
        config.check_and_print_config();
        // If it passes all the checks, we're good.
    }

    #[test]
    #[should_panic]
    fn test_no_files() {
        let mut config = ConfigBuilder::new();
        config = config.set_produce_fastq(false);
        config.check_and_print_config();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean_or_stdev() {
        let mut config = ConfigBuilder::new();
        // paired end set to true, by default, fragment mean and st dev are None
        config = config.set_paired_ended(true);
        config.check_and_print_config();
    }

    #[test]
    #[should_panic]
    fn test_no_frag_mean() {
        let mut config = ConfigBuilder::new();
        config = config.set_paired_ended(true);
        config = config.set_fragment_st_dev(10.0);
        config.check_and_print_config();
    }

    #[test]
    #[should_panic]
    fn test_no_stdev() {
        let mut config = ConfigBuilder::new();
        config = config.set_paired_ended(true);
        config = config.set_fragment_mean(100.0);
        config.check_and_print_config();
    }

}