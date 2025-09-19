use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Lines, Write};
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;
use common::structs::bed_record::BedRecord;
use flate2::read::GzDecoder;
use crate::filter_reads::FilterReadsError;

use common::file_tools::file_io::{read_gzip_lines, read_lines, create_output_file};
/// The purpose of this module is to filter a rneat-generated fastq and vcf file, so that they only show regions
/// of overlap with an input bed file.
enum ReaderType {
    GzipReader(Lines<BufReader<GzDecoder<File>>>),
    NonZipReader(Lines<BufReader<File>>),
}

impl Iterator for ReaderType {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::GzipReader(iterator) => {
                let next = iterator.next();
                match next {
                    Some(result) => {
                        Some(result.expect("Error reading gzipped file."))
                    },
                    None => None,
                }
            },
            Self::NonZipReader(iterator) => {
                let next = iterator.next();
                match next {
                    Some(result) => {
                        Some(result.expect("Error reading non-zipped file."))
                    },
                    None => None,
                }
            },
        }
    }
}

fn prep_files_for_filtering(file_in: &PathBuf, is_gzip: bool, file_out: &PathBuf) -> Result<(ReaderType, GzEncoder<File>), FilterReadsError> {
    // Open file file for reading.
    let lines: ReaderType = {
        if is_gzip {
            open_gzipped(&file_in)
        } else {
            open_unzipped(&file_in)
        }
    };
    let out_file = create_output_file(file_out, true)?;
    let buffer = GzEncoder::new(out_file, Compression::default());
    Ok((lines, buffer))
}

fn open_gzipped(filename: &PathBuf) -> ReaderType {
    ReaderType::GzipReader(read_gzip_lines(filename).expect("Error decompressing gzipped fastq file."))
}

fn open_unzipped(filename: &PathBuf) -> ReaderType {
    ReaderType::NonZipReader(read_lines(filename).expect("Error reading fastq file while filtering."))
}

fn parse_fastq_record_name(line: &str) -> Result<(&str, usize, usize), FilterReadsError> {
    // The rneat-generated record will have the pattern 
    // "@neat_generated_<contig_name>_<10-digit first coord>_<10-digit second coord>_<read number>:<1 for Forward, 2 for Reverse>"
    // Our goal is to extract contig_name and the two coordinates
    let split_line: Vec<&str> = line.split('_').collect();
    // [0: @neat, 1: generated, 2: contig_name, etc]
    let contig_name = split_line[2];
    let start: usize = split_line[3].parse()?;
    let end: usize = split_line[4].parse()?;
    Ok((contig_name, start, end))
}

pub fn filter_fastq(
    bed_table: &HashMap<String, Vec<BedRecord>>, 
    fastq_in: &PathBuf, 
    is_gzip: bool,
    fastq_out: &PathBuf,
)-> Result<(), FilterReadsError> {
    // bed_file: path to bed file to use for filtering
    // fastq_in: fastq ready for filtering, can be gzipped or not
    // is_gzip: whether the input file is gzipped or not
    // output_file: will be gzipped.
    // Run prep for filtering to open fastq
    let (in_file, mut out_file) = prep_files_for_filtering(fastq_in, is_gzip, fastq_out)?;
    // For each record in fastq, find name, parse, and check against bed records
    let mut include = false;

    for line in in_file {
        // Technically a quality score can start with @
        if line.starts_with("@neat") {
            // Record name found
            let fragment_region = parse_fastq_record_name(&line)?;
            if bed_table.contains_key(fragment_region.0) {
                let mut found = false;
                for record in &bed_table[fragment_region.0] {
                    if record.overlaps(fragment_region) {
                        // Found a match
                        found = true;
                        out_file.write(format!("{}\n", line).as_bytes())?;
                        break
                    }
                }
                if found {
                    include = true;
                } else {
                    include = false;
                }
            } else {
                // not in the bed at all, skip
                include = false;
            }
        } else {
            // this indicates an actual record line. Fastq has a format of
            //     @name
            //     ACCCATTAC
            //     +
            //     ?????????
            // where the name is unique and in our case contains coordinates of the read
            // The second line is a sequence 
            // the third line is a literal plus sign spacer
            // the foruth line is an encoded ascii string of the quality scores (usually 1-42 once translated)
            // So once we find a name we want, we need to read in the next three lines or we messed up the read
            // Once we hit the next name, we'll turn this off if that one is not included
            if include {
                out_file.write(format!("{}\n", line).as_bytes())?;
            }
        }
    }
    Ok(())
}

pub fn filter_vcf(
    bed_table: &HashMap<String, Vec<BedRecord>>, 
    vcf_in: &PathBuf, 
    is_gzip: bool,
    vcf_out: &PathBuf, 
) -> Result<(), FilterReadsError> {
    // Run prep for filtering to open vcf
    let (infile, mut outfile) = prep_files_for_filtering(vcf_in, is_gzip, vcf_out)?;
    for line in infile {
        if line.starts_with("#") {
            // header line, skip
            continue;
        } else {
            // break the line into a vector by tabs
            let line_vec: Vec<&str> = line.split('\t').collect();
            // A vcf line has the format ["CHROM", "POS", "ID", "REF", "ALT", etc]
            // The only ones we care about are CHROM [0] and POS [1]
            let chrom = line_vec[0];
            let pos: usize = line_vec[1].parse()?;
            if bed_table.contains_key(line_vec[0]) {
                for record in &bed_table[line_vec[0]] {
                    if record.contains(chrom, pos) {
                        // Found a match
                        outfile.write(format!("{}\n", line).as_bytes())?;
                        break
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_prep_file_for_filtering() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();
        let mut temp_file: PathBuf = PathBuf::from(temp_dir.path());
        temp_file.push("Test_data.txt\n");
        let mut file_obj = create_output_file(&temp_file, true).unwrap();
        file_obj.write("Test data!\nThis is only a test\nPlease disregard all information within\nthe end".as_bytes()).unwrap();
        file_obj.flush().unwrap();
        let mut outfile = PathBuf::from(temp_dir.path());
        outfile.push("test_data_filtered.txt");
        let result = prep_files_for_filtering(
            &temp_file, 
            false, 
            &outfile);
        match result {
            Ok(_tuple) => assert!(true),
            Err(_error) => assert!(false),
        }
    }

    #[test]
    fn test_parse_record_name() {
        let record_name = "@neat_generated_chrom1_001000_002000_223:1".to_string();
        assert_eq!(
            parse_fastq_record_name(&record_name).unwrap(),
            ("chrom1", 1000, 2000),
        );
    }
}