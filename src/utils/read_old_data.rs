use std::fs;
use serde::{Deserialize, Serialize};
use serde_json::*;
use utils::quality_scores::QualityScoreModel;

#[derive(Serialize, Deserialize)]
pub struct QualityRaw {
    seed_stats: Vec<f64>,
    stats_from_one: Vec<Vec<Vec<f64>>>
}

impl QualityRaw {
    fn convert_to_weights(&self) -> QualityScoreModel {
        let mut seed_weights: Vec<usize> = Vec::new();
        for item in &self.seed_stats {
            let weight = (item * 10e5).round() as usize;
            seed_weights.push(weight);
        }
        let mut weights_from_one = Vec::new();
        for i in 0..self.stats_from_one.len() {
            let mut position_arrays = Vec::new();
            for j in 0..self.stats_from_one[i].len() {
                let mut score_arrays = Vec::new();
                for k in 0..self.stats_from_one[i][j].len() {
                    let weight = (self.stats_from_one[i][j][k] * 10e5).round() as usize;
                    score_arrays.push(weight);
                }
                position_arrays.push(score_arrays.clone())
            }
            weights_from_one.push(position_arrays.clone())
        }
        let assumed_read_length = seed_weights.len();
        let quality_score_options: Vec<usize> = (0..42).collect();
        let binned_scores = false;

        QualityScoreModel {
            quality_score_options,
            binned_scores,
            assumed_read_length,
            seed_weights,
            weights_from_one,
        }
    }
}

pub fn parse_neat_quality_scores(filename: &str) -> QualityScoreModel {
    let f = fs::File::open(filename);
    let file = match f {
        Ok(l) => l,
        Err(error) => panic!("Problem reading the quality json file: {}", error),
    };
    let quality_raw: QualityRaw = from_reader(file).expect("Problem with json file format.");
    let quality_model: QualityScoreModel = quality_raw.convert_to_weights();
    quality_model
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_yaml() {
        let filename = "test_data/test_model.json".to_string();
        let qual_data = parse_neat_quality_scores(&filename);
        assert_eq!(qual_data.seed_weights.len(), 42);
        assert_eq!(qual_data.weights_from_one.len(), 100);
        assert_eq!(qual_data.weights_from_one[0].len(), 42);
        assert_eq!(qual_data.weights_from_one[0][0].len(), 42);
        let mut out_file = "test_data/test.json".to_string();
        qual_data.write_out_quality_model(&mut out_file).unwrap();
        fs::remove_file(out_file).expect("Failed in removing test json file");
    }
}