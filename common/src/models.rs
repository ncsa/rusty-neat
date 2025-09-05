//! Here are various models from the original NEAT project, written by Zach Stephens in Python2
//! These are rust implementations, and contain a lot of statistical data extracted from the original
//! compressed model data. Some of this may need to be converted to Json ultimately.

pub mod indel_model;
pub mod mutation_model;
pub mod quality_scores;
pub mod sequencing_error_model;
pub mod snp_trinuc_model;
pub mod fragment_length;
pub mod mutated_map;
mod lib;