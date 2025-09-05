//! This library contains some data and functions that are useful in the other models
use std;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use flate2::Compression;
use flate2::write::GzEncoder;
use crate::file_tools::file_io::create_output_file;

pub fn model_writer<T: Serialize>(model: T, filename: &PathBuf) -> std::io::Result<()> 
{
    // This will take any serializable model and write it to file.
    let data = serde_json::to_vec(&model).unwrap();
    let fileout = create_output_file(&filename, false)
        .expect(&format!("Error creating output {:?}. If file alread exists, remove.", &filename));
    let writer = BufWriter::new(fileout);
    let mut encoder = GzEncoder::new(writer, Compression::default());
    encoder.write_all(&data)?;
    Ok(())
}

pub fn model_reader<T>(filename: &PathBuf) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein = File::open(filename)
        .expect(&format!("Error opening file {:?}", &filename));
    let reader = GzDecoder::new(BufReader::new(filein));
    let value: T = serde_json::from_reader(reader).unwrap();
    Ok(value)
}

#[allow(unused)]
pub fn model_unzipped_reader<T>(filename: &PathBuf) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein = File::open(filename)
        .expect(&format!("Error opening file {:?}", &filename));
    let reader = BufReader::new(filein);
    let value: T = serde_json::from_reader(reader).unwrap();
    Ok(value)
}