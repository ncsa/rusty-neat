use crate::utils;

use std::io::Write;
use std::fs::File;
use std::io;

use utils::fasta_tools::sequence_array_to_string;

pub fn write_fastq(
    fastq_filename: &str,
    dataset: Vec<&Vec<u8>>,
) -> io::Result<()> {
    let name_prefix = "neat_generated_".to_string();
    let fastq_filename = String::from(fastq_filename) + ".fastq";
    let mut outfile = File::options().create(true).append(true).open(fastq_filename)?;
    let mut index = 1;
    for sequence in dataset {
        let sequence_len = sequence.len();
        // sequence name
        writeln!(&mut outfile, "@{}", name_prefix.clone() + &index.to_string())?;
        // Array as a string
        writeln!(&mut outfile, "{}", sequence_array_to_string(&sequence))?;
        // The stupid plus sign
        writeln!(&mut outfile, "+")?;
        // Qual score of all F's for the whole thing.
        writeln!(&mut outfile, "{}", std::iter::repeat("F").take(sequence_len).collect::<String>())?;
        index += 1;
    };
    Ok(())
}