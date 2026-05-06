use std::path::PathBuf;
use serde_derive::{Serialize, Deserialize};
use thiserror::Error;
use crate::models::lib::{model_reader, model_writer};
use crate::structs::nucleotides::Nucleotide;

#[derive(Error, Debug)]
pub enum GcBiasModelError {
    #[error("GC bias model must contain exactly 101 weights")]
    InvalidBinCount,
    #[error("GC bias weights must be finite and non-negative")]
    InvalidWeight,
    #[error("GC bias model must contain at least one positive weight")]
    EmptyModel,
    #[error("GC bias model IO error: {0}")]
    IoError(#[from] std::io::Error),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcBiasModel {
    weights_by_percent_gc: Vec<f64>,
    is_uniform: bool,
}

impl GcBiasModel {
    pub fn default() -> Self {
        let weights = vec![1.0; 101];
        Self {
            weights_by_percent_gc: weights,
            is_uniform: true,
        }
    }

    pub fn from_file(path: &PathBuf) -> Result<Self, GcBiasModelError> {
        let data: GcBiasModel = model_reader(path)?;
        Ok(data)
    }

    pub fn write_to_file(&self, path: &PathBuf) -> Result<(), GcBiasModelError> {
        model_writer(self, path)?;
        Ok(())
    }

    pub fn from_weights(weights_by_percent_gc: Vec<f64>) -> Result<Self, GcBiasModelError> {
        // This function creates a GcBias model from a vector of weights. The length must
        // be exactly 101, where each index represents a percentage of GC
        // in the read (0-100,inclusive). Presumably the only thing different the
        // data-generated model will be that the list length is already vetted and
        // it tells us it is uniform or not.
        //
        // We have to enforce the length
        if weights_by_percent_gc.len() != 101 {
            return Err(GcBiasModelError::InvalidBinCount)
        }
        // check if uniform and for invalid values
        let first_item = weights_by_percent_gc[0];
        // fail right away in this case
        if first_item.is_infinite() || first_item.is_nan() || first_item < 0.0 {
            return Err(GcBiasModelError::InvalidWeight)
        }
        // assume uniform and look for an exception
        let mut is_uniform = true;
        for item in &weights_by_percent_gc[1..] {
            if item.is_nan() || (*item < 0.0) || item.is_infinite() {
                return Err(GcBiasModelError::InvalidWeight)
            } else if (*item != first_item) && (is_uniform == true) {
                is_uniform = false;
            }
        }
        if (is_uniform == true) && (first_item == 0.0) {
            // uniformly zero is not allowed
            return Err(GcBiasModelError::EmptyModel)
        }
        Ok(Self {
            weights_by_percent_gc,
            is_uniform,
        })
    }

    pub fn is_uniform(&self) -> bool {
        self.is_uniform
    }

    pub fn weight_for_gc_fraction(&self, gc_fraction: f64) -> f64 {
        let index: usize = if gc_fraction < 0.0 {
            0
        } else if gc_fraction > 1.0 {
            100
        } else {
            (gc_fraction * 100.0).round() as usize
        };
        self.weights_by_percent_gc[index]
    }

    pub fn weight_for_sequence(&self, sequence: &[Nucleotide]) -> f64 {
        let gc_count = sequence
            .iter()
            .filter(|&nuc| *nuc == Nucleotide::G || *nuc == Nucleotide::C)
            .count() as f64;
        let sequence_count = sequence
            .iter()
            .filter(|&base| *base != Nucleotide::N)
            .count() as f64;
        if sequence_count == 0.0 {
            return 1.0
        }
        let gc_fraction = gc_count / sequence_count;
        self.weight_for_gc_fraction(gc_fraction)
    }

    pub fn max_weight(&self) -> f64 {
        *self.weights_by_percent_gc.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
    }

    pub fn mean_weight(&self) -> f64 {
        self.weights_by_percent_gc.iter().sum::<f64>() / self.weights_by_percent_gc.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, C, G, N, T};

    fn weights_with_value(value: f64) -> Vec<f64> {
        vec![value; 101]
    }

    #[test]
    fn test_default_model_has_101_unit_weights() {
        let model = GcBiasModel::default();

        assert_eq!(model.max_weight(), 1.0);
        assert!((model.mean_weight() - 1.0).abs() < 1e-12);

        for percent_gc in 0..=100 {
            let gc_fraction = percent_gc as f64 / 100.0;
            assert!(
                (model.weight_for_gc_fraction(gc_fraction) - 1.0).abs() < 1e-12,
                "Expected uniform weight 1.0 for {}% GC",
                percent_gc
            );
        }
    }

    #[test]
    fn test_from_weights_rejects_too_few_bins() {
        let result = GcBiasModel::from_weights(vec![1.0; 100]);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidBinCount)),
            "Expected InvalidBinCount for 100 weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_too_many_bins() {
        let result = GcBiasModel::from_weights(vec![1.0; 102]);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidBinCount)),
            "Expected InvalidBinCount for 102 weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_negative_weight() {
        let mut weights = weights_with_value(1.0);
        weights[30] = -0.1;

        let result = GcBiasModel::from_weights(weights);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for negative weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_nan_weight() {
        let mut weights = weights_with_value(1.0);
        weights[50] = f64::NAN;

        let result = GcBiasModel::from_weights(weights);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for NaN weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_infinite_weight() {
        let mut weights = weights_with_value(1.0);
        weights[50] = f64::INFINITY;

        let result = GcBiasModel::from_weights(weights);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for infinite weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_all_zero_weights() {
        let result = GcBiasModel::from_weights(weights_with_value(0.0));

        assert!(
            matches!(result, Err(GcBiasModelError::EmptyModel)),
            "Expected EmptyModel for all-zero weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_weight_for_gc_fraction_uses_percent_bins() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.1;
        weights[30] = 0.3;
        weights[50] = 1.0;
        weights[100] = 0.2;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.weight_for_gc_fraction(0.00) - 0.1).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.30) - 0.3).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.50) - 1.0).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(1.00) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_rounds_to_nearest_percent() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[51] = 0.51;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.weight_for_gc_fraction(0.504) - 0.5).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.505) - 0.51).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_clamps_below_zero() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.25;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.weight_for_gc_fraction(-0.10) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_clamps_above_one() {
        let mut weights = weights_with_value(0.0);
        weights[100] = 0.75;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.weight_for_gc_fraction(1.10) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_zero_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.2;

        let model = GcBiasModel::from_weights(weights).unwrap();
        let sequence = vec![A, A, T, T, A, T];

        assert!((model.weight_for_sequence(&sequence) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_fifty_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;

        let model = GcBiasModel::from_weights(weights).unwrap();
        let sequence = vec![A, C, G, T];

        assert!((model.weight_for_sequence(&sequence) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_one_hundred_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[100] = 0.4;

        let model = GcBiasModel::from_weights(weights).unwrap();
        let sequence = vec![C, G, G, C, C, G];

        assert!((model.weight_for_sequence(&sequence) - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_ignores_n_bases_in_gc_fraction() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.8;

        let model = GcBiasModel::from_weights(weights).unwrap();

        // Ignoring N bases leaves A, C, G, T => 2 GC / 4 called bases = 50% GC.
        let sequence = vec![A, C, G, T, N, N];

        assert!((model.weight_for_sequence(&sequence) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_max_weight() {
        let mut weights = weights_with_value(0.1);
        weights[30] = 0.5;
        weights[50] = 1.25;
        weights[70] = 0.75;

        let model = GcBiasModel::from_weights(weights).unwrap();

        assert!((model.max_weight() - 1.25).abs() < 1e-12);
    }

    #[test]
    fn test_mean_weight() {
        let mut weights = weights_with_value(1.0);
        weights[0] = 0.0;
        weights[100] = 2.0;

        let model = GcBiasModel::from_weights(weights).unwrap();

        // Sum is still 101.0 because one bin changed 1 -> 0 and one changed 1 -> 2.
        assert!((model.mean_weight() - 1.0).abs() < 1e-12);
    }
}