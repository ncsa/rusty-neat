use std::collections::HashMap;
use log::*;
use std::fs::File;
use std::io::{BufReader, Lines, Write};
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;
use common::structs::bed_record::BedRecord;
use flate2::read::GzDecoder;
use crate::filter_reads::errors::FilterReadsError;
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
    // Open file for reading.
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
    // gen_reads produces: "@RNEAT_generated_<contig>_<10-digit start>_<10-digit end>/N"
    if !line.starts_with("@RNEAT_generated_") {
        return if line.contains("@RNEAT_generated_") {
            Err(FilterLibError::MalformedLine(line.to_string()))
        } else {
            Err(FilterLibError::InvalidReadName)
        };
    }
    // Slice off the "@RNEAT_generated_" prefix (17 chars) and "/N" strand suffix (2 chars).
    let trimmed = &line[17..line.len() - 2];
    let split_line: Vec<&str> = trimmed.split('_').collect();
    let length = split_line.len();
    let start: usize = split_line[length - 2].parse().unwrap();
    let end: usize = split_line[length - 1].parse().unwrap();
    let contig = split_line[..length - 2].join("_");
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
                            if include {
                                out_file.write_all(line.as_bytes())?;
                                out_file.write_all(b"\n")?;
                            }
                        },
                    }
                    // That line is processed and we can move on
                    continue
                },
            }.unwrap(); // This unwrap will either continue us or be a region, so this is good to go.
            include = if let Some(records) = bed_table.get(&contig) {
                let found = records.iter().any(|r| r.overlaps(&contig, start, end));
                if found {
                    out_file.write_all(line.as_bytes())?;
                    out_file.write_all(b"\n")?;
                }
                found
            } else {
                false
            };
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
                out_file.write_all(line.as_bytes())?;
                out_file.write_all(b"\n")?;
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
            outfile.write_all(line.as_bytes())?;
            outfile.write_all(b"\n")?;
        } else {
            let line_vec: Vec<&str> = line.split('\t').collect();
            let chrom = line_vec[0];
            let pos: usize = line_vec[1].parse()?;
            if let Some(records) = bed_table.get(line_vec[0]) {
                if records.iter().any(|r| r.contains(chrom, pos)) {
                    outfile.write_all(line.as_bytes())?;
                    outfile.write_all(b"\n")?;
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
        let record_name = "@RNEAT_generated_chrom_1_0000001000_0000002000/1".to_string();
        assert_eq!(
            parse_fastq_record_name(&record_name).unwrap(),
            ("chrom_1".to_string(), 1000, 2000),
        );
    }

    #[test]
    fn test_filter_fastq() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();

        // Two reads: first inside BED chr1:0-2000, second outside
        let fastq_content = concat!(
            "@RNEAT_generated_chr1_0000001000_0000002000/1\n",
            "ACGTACGT\n",
            "+\n",
            "IIIIIIII\n",
            "@RNEAT_generated_chr1_0000005000_0000006000/1\n",
            "TTTTTTTT\n",
            "+\n",
            "IIIIIIII\n",
        );
        let input = temp_dir.path().join("input.fastq");
        std::fs::write(&input, fastq_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("output.fastq.gz");
        filter_fastq(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(lines.len(), 4, "Only the in-range read (4 lines) should be in output");
        assert!(lines[0].contains("0000001000"));
    }

    #[test]
    fn test_filter_vcf() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();

        // Two variants: pos 1500 inside BED chr1:0-2000, pos 5000 outside
        let vcf_content = concat!(
            "##fileformat=VCFv4.1\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            "chr1\t1500\t.\tA\tG\t37\tPASS\t.\tGT\t0/1\n",
            "chr1\t5000\t.\tT\tC\t37\tPASS\t.\tGT\t0/1\n",
        );
        let input = temp_dir.path().join("input.vcf");
        std::fs::write(&input, vcf_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("output.vcf.gz");
        filter_vcf(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let variant_lines: Vec<&String> = lines.iter().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(variant_lines.len(), 1, "Only the in-range variant should appear");
        assert!(variant_lines[0].contains("1500"));
    }
}