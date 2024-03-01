use std::ffi::CString;
use clap::Parser;

#[derive(Parser, Debug)]
pub struct Cli {
    #[arg(short='y', long="configuration_yaml", default_value_t=String::new())]
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