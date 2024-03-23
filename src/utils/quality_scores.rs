#[allow(dead_code)]
pub struct QualityScores {
    score_options: Vec<usize>,
    score_weights: Vec<Vec<usize>>,
    read_length: usize,
}