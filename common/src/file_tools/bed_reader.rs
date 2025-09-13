use log::*;
use thiserror::Error;
use std::{ 
    io, num::ParseIntError, path::PathBuf
};
use crate::{
    file_tools::file_io::{read_gzip_lines, read_lines}, 
    structs::bed::{BedErrors, BedRecord}
};

#[derive(Debug, Error)]
pub enum BedReaderError {
    #[error("Error reading the bed file. Check file format: {0}")]
    MalformedBed(String),
    #[error("IO error from the bed reader: {0}")]
    IoError(#[from] io::Error),
    #[error("Error creating the bed record: {0}")]
    BedRecordError(#[from] BedErrors),
    #[error("Error deserializing the CSV: {0}")]
    CsvError(#[from] csv::Error),
    #[error("Error parsing ints from the bed file: {0}")]
    ParserError(#[from] ParseIntError),
    #[error("File has an unknown or missing file extension: {0}")]
    FileExtensionUnknown(String)
}

pub fn read_bed(
    bed_path: &PathBuf,
) -> Result<Vec<BedRecord>, BedReaderError> {
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

fn process_gzip_bed (filename: &PathBuf) -> Result<Vec<BedRecord>, BedReaderError> {
    let reader = read_gzip_lines(filename)?;
    let mut file_records = Vec::new();
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
            contig, 
            start as usize, 
            end as usize, 
            other,
        )?;
        println!("{:?}", record);
        file_records.push(record)
    }
    Ok(file_records)
}

fn process_bed (filename: &PathBuf) -> Result<Vec<BedRecord>, BedReaderError> {
    let open_file = read_lines(filename)?;
    let mut file_records = Vec::new();
    for result in open_file {
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
            contig, 
            start as usize, 
            end as usize, 
            other,
        )?;
        println!("{:?}", record);
        file_records.push(record)
    }
    Ok(file_records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_bed() {
        let file = "/home/joshfactorial/data/Escherichia_coli_K12_strain_K-12_substrain_MG1655_id4242.bed.gz";
        let pathbuf = PathBuf::from(file);
        let records = read_bed(&pathbuf).unwrap();
        for record in records {
            println!("{:?}", record)
        }
    }
}