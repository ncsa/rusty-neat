use std::collections::HashMap;
use log::*;
use std::fs::File;
use std::io::{BufReader, Lines, Write};
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;
use common::structs::bed_record::BedRecord;
use flate2::read::GzDecoder;
use crate::filter_reads::FilterReadsError;
use thiserror::Error;

use common::file_tools::file_io::{read_gzip_lines, read_lines, create_output_file};
/// The purpose of this module is to filter a rneat-generated fastq and vcf file, so that they only show regions
/// of overlap with an input bed file.
enum ReaderType {
    GzipReader(Lines<BufReader<GzDecoder<File>>>),
    NonZipReader(Lines<BufReader<File>>),
}


#[derive(Debug, Error)]
enum FilterLibError {
    #[error("Possible Q score")]
    InvalidReadName,
    #[error("Line parsing was not understood: {0}")]
    MalformedLine(String)
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

fn parse_fastq_record_name(line: &str) -> Result<(String, usize, usize), FilterLibError> {
    // The rneat-generated record will have the pattern 
    // "@neat_generated_<contig_name>_<10-digit first coord>_<10-digit second coord>_<read number>:<1 for Forward, 2 for Reverse>"
    // Our goal is to extract contig_name and the two coordinates
    //
    // Start by checking if this is even a read
    let result = line.find("@neat_generated_");
    match result {
        Some(index) => {
            if index == 0 {
                // valid read name
            } else {
                // We don't know what this is, calling it an error.
                return Err(FilterLibError::MalformedLine(line.to_string()))
            }
        },
        None => {
            // This is probably just a quality score
            return Err(FilterLibError::InvalidReadName)
        }
    }
    let end = line.len();
    // We can remove the strand number
    let mut trimmed_line = line[..(end - 2)].to_string();
    trimmed_line = trimmed_line[16..].to_string();
    // What we have left should be <name_parts>_<start>_<end>_<read_number>
    let split_line: Vec<&str> = trimmed_line.split('_').collect();
    let length = split_line.len();
    // start is 2 from the end
    let start: usize = split_line[length-3].parse().unwrap();
    // end is 1 from the end
    let end: usize = split_line[length-2].parse().unwrap();
    // However "contig_name" might contain an underscore, so we make sure we have the right coordinates
    // Grab everything up to the start
    let mut contig = String::new();
    for i in 0..(length-3) {
        // No underscore in the last loop
        if i == (length-4) {
            contig.push_str(split_line[i])
        } else {
            contig.push_str(split_line[i]);
            // replace what "split" removed
            contig.push('_');
        }
    }
    Ok((contig, start, end))

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
        if line.starts_with("@") {
            // Record name candidate found
            let candidate_result = parse_fastq_record_name(&line);
            // It's possible that this is not a read name, but rather a quality score starting with @
            let (contig, start, end) = match candidate_result {
                Ok(region) => {
                    Some(region)
                },
                Err(error) => {
                    // If this is a quality score and we need to write it, we do, then we run the next loop
                    debug!("Had a weird record: {line}");
                    match error {
                        FilterLibError::MalformedLine(line) => {
                            return Err(FilterReadsError::FastqFilterError(format!("Filter error with line: {}", line)))
                        },
                        FilterLibError::InvalidReadName => {
                            // In this case this is a q-score and, if included, we should write it
                            if include {
                                out_file.write(format!("{}\n", line).as_bytes())?;
                            };
                        },
                    }
                    // That line is processed and we can move on
                    continue
                },
            }.unwrap(); // This unwrap will either continue us or be a region, so this is good to go.
            if bed_table.contains_key(&contig) {
                let mut found = false;
                for record in &bed_table[&contig] {
                    if record.overlaps(&contig, start, end) {
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
            // header line, write out with no changes
            outfile.write(format!("{}\n", line).as_bytes())?;;
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
        let record_name = "@neat_generated_chrom_1_0000001000_0000002000_223:1".to_string();
        assert_eq!(
            parse_fastq_record_name(&record_name).unwrap(),
            ("chrom_1".to_string(), 1000, 2000),
        );
    }
}