use std::io::Write;
use std::fs::File;
use std::io;

use utils::fasta_tools::sequence_array_to_string;

pub fn write_fastq(
    fastq_filename: &str,
    paired_ended: bool,
    dataset: Vec<&Vec<u8>>,
) -> io::Result<()> {
    /*
    Writes fastq files. At the moment, it only writes out single r1 file, but will eventually write
    out r1 and r2 files.
     */

    // name_prefix is for the prefix for the read names. Reads are numbered in output order
    // (Although this feature is currently untested and unknown).
    // (May need sorting.)
    let name_prefix = "neat_generated_".to_string();
    // vector to hold the fastq filenames, since we don't know ahead of time if we need one or two
    let mut fastq_files: Vec<String> = Vec::new();
    fastq_files.push(String::from(fastq_filename) + "_r1.fastq");
    if paired_ended {
        fastq_files.push(String::from(fastq_filename) + "_r1.fastq");
    }
    // open the filename
    let mut outfile = File::options().create_new(true).append(true).open(fastq_filename)?;
    // write out sequence by sequence using index to track the numbering
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