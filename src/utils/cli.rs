extern crate clap;

use clap::Parser;

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
     */
    #[arg(short='C', long="configuration_yaml", default_value_t=String::new())]
    pub config: String,
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
}