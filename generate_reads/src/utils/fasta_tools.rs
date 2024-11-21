// This library contains tools needed to process fasta files as input and output.
use common;

use std::{io, str, usize};
use std::io::Write;
use std::*;
use std::collections::{HashMap, VecDeque};
use fasta_reader::read_fasta;



pub fn write_fasta(
    mutated_fasta: &Box<HashMap<String, Vec<u8>>>,
    fasta_order: &VecDeque<String>,
    overwrite_output: bool,
    output_file: &str,
) -> io::Result<()> {
    // // reference_fasta: the hashmap of mutated sequences with contig names
    // // variants: the location, type, and size of the variants that have been inserted in the fasta.
    // // fasta_lengths: A vector with the proper order for the fasta elements and lengths of the
    // //      sequences.
    // // overwrite_output: Boolean that controls whether an existing file is overwritten.
    // // output_file: the prefix for the output file name
    // const LINE_LENGTH: usize = 70; // fasta is usually writes with 70 characters per line
    // let mut output_fasta = format!("{}.fasta", output_file);
    // let mut outfile = open_file(&mut output_fasta, overwrite_output)
    //     .expect(&format!("Error opening {}", output_fasta));
    // for contig in fasta_order {
    //     let sequence = &mutated_fasta[contig];
    //     let contig_length = sequence.len();
    //     // Write contig name
    //     writeln!(&mut outfile, ">{}", contig)?;
    //     // write sequences[ploid] to this_fasta
    //     let mut i = 0;
    //     while i < contig_length {
    //         let mut outlines = String::new();
    //         // Generate a few lines before writing
    //         'num_lines: for _ in 0..50 {
    //             for j in 0..LINE_LENGTH {
    //
    //                 // Check to make sure we are not out of bounds, which will happen if we
    //                 // reach the end and there are not exactly 70 bases left (basically, every time)
    //                 if (i + j) >= contig_length {
    //                     // append the remaining bases
    //                     for base in sequence[i..].iter() {
    //                         outlines += &base.to_string();
    //                     }
    //                     outlines += "\n";
    //                     // ensure outer loop breaks after writing final line, break inner loop.
    //                     i += j;
    //                     // todo this is missing a line break after the last full 70 chars.
    //                     break 'num_lines
    //                 }
    //                 outlines += &sequence[i+j].to_string();
    //                 i += 1
    //             }
    //             outlines += "\n";
    //         }
    //         writeln!(&mut outfile, "{}", outlines)?;
    //     }
    // }
    // Ok(())
    todo!();
}

#[cfg(test)]
mod tests {
    use std::path::Path;
    use super::*;

    #[test]
    fn test_write_fasta() -> Result<(), Box<dyn error::Error>> {
        // let mutated_seq1 = vec![A, G, T, A, C, T, C, A, G, T, G, T, T, C, C, T];
        // let mutated_seq2 = vec![G, C, G, C, G, C, T, T, T, T, G, G, C, A, C, G, T, A, A];
        // let fasta_map = Box::new(HashMap::from([
        //     ("chr1".to_string(), mutated_seq1.clone()),
        //     ("chr2".to_string(), mutated_seq2.clone()),
        // ]));
        //
        // let fasta_order = VecDeque::from([
        //     String::from("chr1"),
        //     String::from("chr2"),
        // ]);
        // let output_file = "test";
        // let test_write = write_fasta(
        //     &fasta_map,
        //     &fasta_order,
        //     false,
        //     output_file
        // ).unwrap();
        // let file_name = "test.fasta";
        // assert_eq!(test_write, ());
        // let attr = fs::metadata(file_name).unwrap();
        // assert!(attr.len() > 0);
        // fs::remove_file(file_name)?;
        // Ok(())
    }
}
