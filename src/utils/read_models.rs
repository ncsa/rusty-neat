use utils::quality_scores::QualityScoreModel;
use std::fs;
use serde_json;
pub fn read_quality_score_model_json(filename: &str) -> QualityScoreModel {
    let f = fs::File::open(filename);
    let file = match f {
        Ok(l) => l,
        Err(error) => panic!("Problem reading the quality json file: {}", error),
    };
    serde_json::from_reader(file).expect("Problem with json file format.")
}