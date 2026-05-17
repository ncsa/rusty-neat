use log::*;
use thiserror::Error;
use std::{ 
    collections::HashMap,
    io::{self, BufReader, Lines, Read}, 
    num::ParseIntError, 
    path::PathBuf
};
use crate::{
    file_tools::file_io::{read_gzip_lines, read_lines, is_gzipped_file},
    structs::bed_record::{BedErrors, BedRecord}
};

#[derive(Debug, Error)]
pub enum BedReaderError {
    #[error("Error reading the bed file. Check file format: {0}")]
    MalformedBed(String),
    #[error("IO error from the bed reader: {0}")]
    IoError(#[from] io::Error),
    #[error("Error creating the bed record: {0}")]
    BedRecordError(#[from] BedErrors),
    #[error("Error parsing ints from the bed file: {0}")]
    ParserError(#[from] ParseIntError),
}

pub fn read_bed(
    bed_path: &PathBuf,
    mut_bed: bool,
) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    if !bed_path.is_file() {
        return Err(BedReaderError::MalformedBed(format!(
            "Input bed file does not exist: {:?}", 
            bed_path
        )))
    } 
    if is_gzipped_file(bed_path)? {
        process_gzip_bed(bed_path, mut_bed)
    } else {
        process_bed(bed_path, mut_bed)
    }
}

fn process_gzip_bed(filename: &PathBuf, mut_bed: bool) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let reader = read_gzip_lines(filename)?;
    read_open_bed(reader, mut_bed)
}

fn process_bed(filename: &PathBuf, mut_bed: bool) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let reader = read_lines(filename)?;
    read_open_bed(reader, mut_bed)
}
    
fn read_open_bed<P: Read> (reader: Lines<BufReader<P>>, mut_bed: bool) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let mut file_records: HashMap<String, Vec<BedRecord>> = HashMap::new();
    for result in reader {
        let rec_str: String = result?;
        let mut fields: Vec<&str> = rec_str.split_whitespace().collect();
        let contig: String = fields[0].to_string();
        let start: isize = fields[1].parse()?;
        let end: isize = fields[2].parse()?;
        if start < 0 || end < 0{
            // Skip this record
            warn!("Skipped bed region with negative coordinates: {:?}", rec_str);
            continue
        }
        let record = if mut_bed {
            // We only care about the "other" part in a mut_regions_bed
            let mut other = String::new();
            for remainder in &mut fields[3..] {
                other = other + remainder + "\t";
            }
            BedRecord::new_mut_region_record(
                contig.clone(),
                start as usize,
                end as usize,
                &other,
            )?
        } else {
            BedRecord::new_bed_record(
                contig.clone(),
                start as usize,
                end as usize,
            )?
        };
        if !file_records.contains_key(&contig) {
            file_records.insert(contig.clone(), Vec::new());
        }
        file_records.get_mut(&contig)
            .expect(&format!("BUG: The contig was not found in file records {}", contig))
            .push(record)
    }
    Ok(file_records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_process_bed() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bed_path = temp_dir.path().join("test.bed");
        let mut f = std::fs::File::create(&bed_path).unwrap();
        writeln!(f, "chr1\t100\t200\textra_field").unwrap();
        writeln!(f, "chr1\t500\t600").unwrap();
        writeln!(f, "chr2\t300\t400").unwrap();
        drop(f);

        let result = read_bed(&bed_path, false).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result["chr1"].len(), 2);
        assert_eq!(result["chr2"].len(), 1);
        assert_eq!(result["chr1"][0].start, 100);
        assert_eq!(result["chr1"][0].end, 200);
        assert_eq!(result["chr2"][0].start, 300);
        assert_eq!(result["chr2"][0].end, 400);
    }

    #[test]
    fn test_process_mut_bed() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bed_path = temp_dir.path().join("test_mut.bed");
        let mut f = std::fs::File::create(&bed_path).unwrap();
        writeln!(f, "chr1\t100\t200\tmut_rate=0.01").unwrap();
        writeln!(f, "chr2\t300\t400\tname=region2;mut_rate=0.005").unwrap();
        drop(f);

        let result = read_bed(&bed_path, true).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result["chr1"][0].mut_rate, Some(0.01));
        assert_eq!(result["chr2"][0].mut_rate, Some(0.005));
    }

    #[test]
    fn test_process_bed_no_extension_succeeds() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bad_path = temp_dir.path().join("no_extension");
        let mut f = std::fs::File::create(&bad_path).unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        drop(f);
        let result = read_bed(&bad_path, false).unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result["chr1"].len(), 1);
        assert_eq!(result["chr1"][0].start, 100);
        assert_eq!(result["chr1"][0].end, 200);
    }

    #[test]
    fn test_read_bed_file_not_found() {
        let bed_path = PathBuf::from("non_existent_file.bed");
        let result = read_bed(&bed_path, false);
        assert!(matches!(result, Err(BedReaderError::MalformedBed(_))));
    }

    #[test]
    fn test_read_bed_negative_coordinates() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bed_path = temp_dir.path().join("negative.bed");
        let mut f = std::fs::File::create(&bed_path).unwrap();
        writeln!(f, "chr1\t-100\t200").unwrap();
        writeln!(f, "chr1\t100\t-200").unwrap();
        writeln!(f, "chr2\t300\t400").unwrap();
        drop(f);

        let result = read_bed(&bed_path, false).unwrap();
        // chr1 records should be skipped
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("chr2"));
        assert!(!result.contains_key("chr1"));
    }

    #[test]
    fn test_read_bed_malformed_int() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bed_path = temp_dir.path().join("malformed.bed");
        let mut f = std::fs::File::create(&bed_path).unwrap();
        writeln!(f, "chr1\tabc\t200").unwrap();
        drop(f);

        let result = read_bed(&bed_path, false);
        assert!(matches!(result, Err(BedReaderError::ParserError(_))));
    }
}
