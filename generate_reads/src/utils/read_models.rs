use serde_json;
use std::fs;
use common::models::quality_scores::QualityScoreModel;
use crate::data::quality_score_data::RawQualityScoreData;

pub fn read_quality_score_model_file(filename: &str) -> QualityScoreModel {
    let f = fs::File::open(filename);
    let file = match f {
        Ok(l) => l,
        Err(error) => panic!("Problem reading the quality json file: {}", error),
    };
    serde_json::from_reader(file).expect("Problem with json file format.")
}

pub fn read_quality_score_raw_data(data: RawQualityScoreData) -> QualityScoreModel {
    serde_json::from_str(&data.data).unwrap()
}
