// This is a pretty basic implementation of Clap CLI
// The idea of this interface and supporting code is that the user can enter an optional
// config file that will take the place of the other command line options, except the logging
// features, which are handled separately. Either way, these options are read into a configuration
// struct that holds the variables for the run. Logging, meanwhile, is handled separately,
// outside run configuration parsing.

extern crate clap;

use clap::Parser;
use std::env;

#[derive(Parser, Debug)]
pub struct Cli {
    /*
    Command line interface for neat. The current configuration items allowed are listed below,
    please add to the list below as new things are added (then update this message when the options
    are stable

    config <String> = the full path to a configuration yaml file. If entered, it will override
        all other command line inputs. No default.

    reference <String> = The relative path or full path to the reference file. Must be in fasta
        file format. Default "data/H1N1.fa"
    output_dir <String> = The directory where output files will be written. If nothing is entered,
        it will write output files to the current working directory.
    output_file_prefix <String> = output files will start with this name. Default = neat_out
    read_length <usize> = the length of the reads in the output fastq file. Default = 150
    coverage <usize> = The average depth per read of the fastq files. Default = 10

    The following commands are independent of the config and not affected by it one way or another:
    log_level <String> = Set a log level for the run. Everything at and above the level chosen will
        be displayed in both logs. See simplelog docs for more info:
        https://docs.rs/simplelog/latest/simplelog/enum.Level.html
    log_dest <String> = Full path filename where to write the log. The default is current working
        dir, filename "neat_out.log," which is set during config parsing.
     */
    #[arg(short='C', long="configuration_yaml", default_value_t=String::new(),
    help="Enter a full path and filename to a configuration file. \
    This will override most other options")]
    pub config: String,

    // All of these arguments are overridden by the config file
    #[arg(short='r', long="reference", default_value_t=String::from("data/H1N1.fa"))]
    pub reference: String,
    #[arg(short='o', long="output_dir", default_value_t=String::new())]
    pub output_dir: String,
    #[arg(short='f', long="output_file_prefix", default_value_t=String::new())]
    pub output_file_prefix: String,
    #[arg(short='l', long="read_len", default_value_t = 150)]
    pub read_length: usize,
    #[arg(short='c', long="coverage", default_value_t = 10)]
    pub coverage: usize,

    // These options relate to the logging features and are not overridden by a config
    #[arg(long="log-level", default_value_t=String::from("Trace"), help="Enter one of Trace, Debug, Info, Warn, Error, Off")]
    pub log_level: String,
    #[arg(long="log-dest", default_value_t=env::current_dir().unwrap().display().to_string() + "neat_out.log", help="Full path and name to log file")]
    pub log_dest: String,
}

// Tests are handled in other places.