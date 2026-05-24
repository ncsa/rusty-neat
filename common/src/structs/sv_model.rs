//! Statistical model for generating symbolic / structural variants
//! (`<DEL>` / `<DUP>` / `<CNV>`) de novo from a learned distribution.
//!
//! Learned by `gen-mut-model` from a real-world SV-rich VCF and consumed
//! by `gen-reads` to emit symbolic records that drive coverage modulation
//! through the depth-modulation path added in v1.9.0. Wrapped in
//! `Option<SvModel>` on [`crate::models::mutation_model::MutationModel`]
//! so v1.9 model JSON files load unchanged (default = `None`).
//!
//! v1 scope: `<DEL>` / `<DUP>` / `<CNV>` only. `<INV>` and breakends
//! observed in training data are counted-and-dropped (the trainer logs
//! how many were excluded); they may be added later as the read-content
//! modeling for inversions and junctions matures.

use crate::structs::variants::SvType;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Minimum SV length (in bases) considered an SV rather than an indel.
/// Matches the conventional cutoff used by Manta, DELLY, Lumpy, GRIDSS,
/// and the rest of the SV-calling ecosystem. Enforced on both sides:
///   - Training: input records with span `< MIN_SV_LENGTH_BP` are
///     excluded from the log-normal fit (logged at warn level with a
///     count).
///   - Generation: sampled lengths below the cutoff are rejected and
///     re-drawn, with a retry cap.
pub const MIN_SV_LENGTH_BP: usize = 50;

/// Learned statistical model for de novo symbolic-SV generation.
///
/// Constructed by `gen-mut-model::runner` after the trinucleotide /
/// SNP / indel fit completes. `Option<SvModel>` on `MutationModel` is
/// `None` when the training VCF had zero usable SV observations — see
/// [`SvModel::is_usable`] for the "do we have enough to generate?" check.
#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
pub struct SvModel {
    /// Per-base SV rate across the training corpus
    /// (= total observed SVs / total bases of training contigs seen).
    /// Used as the Poisson rate parameter when sampling SV count per
    /// generation contig: `lambda = per_base_rate * contig_len *
    /// sv_rate_scale`.
    pub per_base_rate: f64,

    /// Probability distribution over SV types. v1 keys are
    /// `SvType::Del`, `SvType::Dup`, and `SvType::Cnv`. Probabilities
    /// sum to ~1.0 across present keys; the generator uses a weighted
    /// pick.
    pub type_probabilities: HashMap<SvType, f64>,

    /// Length distribution per SV type — log-normal parameters
    /// `(mu, sigma)` in log-space. Sampling: `round(exp(N(mu, sigma)))`,
    /// rejecting any draw below [`MIN_SV_LENGTH_BP`] or above a
    /// caller-supplied upper bound (typically `contig_len / 4` for
    /// sanity).
    ///
    /// Conditioning length on type is the cheapest concession to
    /// biological realism: DELs trend smaller than DUPs and CNVs, so a
    /// single global length distribution would mis-bias against
    /// whichever type the training corpus had less of.
    pub length_log_normal: HashMap<SvType, (f64, f64)>,

    /// Copy-number distribution for `<CNV>` records: CN value →
    /// probability. CN=2 (the diploid reference) is excluded at training
    /// time — a CNV with CN=2 isn't a variant. Probabilities sum to
    /// ~1.0 across present keys.
    pub cnv_copy_number_distribution: HashMap<u32, f64>,

    /// Probability that a generated SV is homozygous (vs heterozygous)
    /// across all SV types. Estimated as the fraction of training
    /// observations whose genotype was homozygous.
    pub homozygous_frequency: f64,
}

impl SvModel {
    /// Returns `true` when this model has enough learned signal to drive
    /// generation — at least one type probability and a non-zero
    /// per-base rate. The trainer should return `None` instead of an
    /// unusable `Some(SvModel)`; this check exists as a guard for
    /// hand-constructed or partially-deserialized models.
    pub fn is_usable(&self) -> bool {
        self.per_base_rate > 0.0 && !self.type_probabilities.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sv_model_is_not_usable() {
        // `SvModel::default()` has per_base_rate=0 and empty maps — it's
        // a zero value, not a generatable model. `is_usable()` must say
        // so to protect downstream samplers from a NaN / divide-by-zero
        // path.
        let m = SvModel::default();
        assert!(!m.is_usable());
        assert_eq!(m.per_base_rate, 0.0);
        assert!(m.type_probabilities.is_empty());
        assert!(m.length_log_normal.is_empty());
        assert!(m.cnv_copy_number_distribution.is_empty());
        assert_eq!(m.homozygous_frequency, 0.0);
    }

    #[test]
    fn min_sv_length_is_the_industry_50bp_cutoff() {
        // Lock in the constant — the value itself isn't load-bearing
        // here, but downstream code (trainer + generator) keys off it,
        // and a silent change from 50 to anything else would shift the
        // entire SV/indel boundary across the pipeline.
        assert_eq!(MIN_SV_LENGTH_BP, 50);
    }

    #[test]
    fn populated_sv_model_round_trips_through_serde_json() {
        // The SvModel will sit inside MutationModel's serialized JSON;
        // serde-json round-trip with a known shape protects against
        // accidental enum-tag changes on SvType (which is the HashMap
        // key) breaking model files.
        let mut m = SvModel {
            per_base_rate: 2.3e-7,
            homozygous_frequency: 0.31,
            ..Default::default()
        };
        m.type_probabilities.insert(SvType::Del, 0.6);
        m.type_probabilities.insert(SvType::Dup, 0.3);
        m.type_probabilities.insert(SvType::Cnv, 0.1);
        m.length_log_normal.insert(SvType::Del, (7.5, 1.8));
        m.length_log_normal.insert(SvType::Dup, (8.1, 1.5));
        m.cnv_copy_number_distribution.insert(0, 0.05);
        m.cnv_copy_number_distribution.insert(1, 0.45);
        m.cnv_copy_number_distribution.insert(3, 0.30);
        m.cnv_copy_number_distribution.insert(4, 0.20);

        let json = serde_json::to_string(&m).unwrap();
        let back: SvModel = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);
        assert!(back.is_usable());
    }

    #[test]
    fn is_usable_requires_both_rate_and_type_probabilities() {
        // Either alone is insufficient — guards against partially
        // populated models that would crash a weighted pick.
        let mut only_rate = SvModel::default();
        only_rate.per_base_rate = 1e-6;
        assert!(!only_rate.is_usable());

        let mut only_types = SvModel::default();
        only_types.type_probabilities.insert(SvType::Del, 1.0);
        assert!(!only_types.is_usable());
    }
}
