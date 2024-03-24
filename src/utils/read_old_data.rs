use std::collections::HashMap;
use std::fs;
use serde::{Deserialize, Serialize};
use serde_json::*;

#[derive(Debug, Serialize, Deserialize)]
pub struct QualityData {
    seed_weights: Vec<usize>,
    weights_from_one: Vec<Vec<usize>>,
}

pub fn parse_neat_quality_scores(filename: &str) -> QualityData {
    let f = fs::File::open(filename);
    let file = match f {
        Ok(l) => l,
        Err(error) => panic!("Problem reading the quality yaml file: {}", error),
    };
    // Uses serde_json to read the file into a HashMap
    let scrape_config: HashMap<String, serde_yaml::Value> = serde_json::from_reader(file)
        .expect("Could not read values");
    // Need to unwrap the inner values
    let mut seed_weights: Vec<usize> = Vec::new();
    for weight in scrape_config.get("seed_weights").unwrap().as_sequence().unwrap() {
        let unwrapped: usize = (weight.as_f64().unwrap() * 10000.0).round() as usize;
        seed_weights.push(unwrapped)
    }
    let mut weights_from_one: Vec<Vec<usize>> = Vec::new();
    for stat_array in scrape_config.get("stats_from_one").unwrap().as_sequence().unwrap() {
        let mut weight_array: Vec<usize> = Vec::new();
        for Some(weight) in stat_array {
            let unwrapped: usize = (weight.as_f64().unwrap() * 10000.0).round() as usize;
            weight_array.push(unwrapped)
        }
        weights_from_one.push(weight_array)
    }

    QualityData {
        seed_weights,
        weights_from_one
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_yaml() {
        let filename = "/home/joshfactorial/code/neat2/models/quality_score_model.json".to_string();
        let qual_data = parse_neat_quality_scores(&filename);
        println!("{:?}", qual_data);
    }
}