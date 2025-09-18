use log::*;
use thiserror::Error;
use std::{ 
    collections::HashMap,
    io::{self, BufReader, Lines, Read}, 
    num::ParseIntError, 
    path::PathBuf
};
use crate::{
    file_tools::file_io::{read_gzip_lines, read_lines}, 
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
    #[error("File has an unknown or missing file extension: {0}")]
    FileExtensionUnknown(String)
}

pub fn read_bed(
    bed_path: &PathBuf,
) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    if !bed_path.is_file() {
        return Err(BedReaderError::MalformedBed(format!(
            "Input bed file does not exist: {:?}", 
            bed_path
        )))
    } 
    let ext = bed_path.extension();
    match ext {
        Some(ext) => {
            if ext == "gz" {
                return process_gzip_bed(&bed_path)
            } else {
                return process_bed(&bed_path)
            }
        },
        None => {
            return Err(BedReaderError::FileExtensionUnknown(
                format!("file: {:?}", bed_path)
            ))
        }
    }
}

fn process_gzip_bed(filename: &PathBuf) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let reader = read_gzip_lines(filename)?;
    return Ok(read_open_bed(reader).expect("Problem processing lines of gzipped bed file"))
}

fn process_bed(filename: &PathBuf) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let reader = read_lines(filename)?;
    return Ok(read_open_bed(reader).expect("Error processing bed file"))
}
    
fn read_open_bed<P: Read> (reader: Lines<BufReader<P>>) -> Result<HashMap<String, Vec<BedRecord>>, BedReaderError> {
    let mut file_records: HashMap<String, Vec<BedRecord>> = HashMap::new();
    for result in reader {
        let rec_str: String = result?;
        let mut fields: Vec<&str> = rec_str.split_whitespace().collect();
        let contig: String = fields[0].to_string();
        let start: isize = fields[1].parse()?;
        let end: isize = fields[2].parse()?;
        if start < 0 || end < 0{
            // Skip this record
            warn!("Not sure what to do with negative coordinates yet: {:?}", rec_str);
            continue
        }
        let mut other = Vec::new();
        for remainder in &mut fields[3..] {
            other.push(remainder.to_string());
        }
        let record = BedRecord::new(
            contig.clone(), 
            start as usize, 
            end as usize, 
            other,
        )?;
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

    #[test]
    fn test_process_bed() {
        // TODO add some tests
        assert!(true)
    }
}