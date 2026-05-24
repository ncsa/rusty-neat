//! Bundled default [`SvModel`] used by [`crate::models::mutation_model::MutationModel::default`].
//!
//! The parameters here are **rough approximations** of the gnomAD-SV v2.1
//! callset summary statistics (Collins et al., *Nature* 581, 444–451,
//! 2020), eyeballed from the paper's published distributions so users
//! can opt into de novo SV generation without first running
//! `gen-mut-model` against a real corpus. They are NOT the result of
//! refitting on the actual gnomAD VCFs — anyone doing serious
//! benchmarking should train their own model with
//! `gen-mut-model -c <config>.yml` pointed at the SV-rich VCF of their
//! choice and override the embedded default by passing
//! `mutation_model: <trained>.json.gz` in their gen-reads YAML.
//!
//! Provenance of each parameter is documented at the call site; the
//! intent is that a future regen against real gnomAD data is a search-
//! and-replace in this one file rather than a hunt across the codebase.

use crate::structs::sv_model::SvModel;
use crate::structs::variants::SvType;
use std::collections::HashMap;

/// Build the bundled default `SvModel`.
///
/// Returns a fully-populated, [`SvModel::is_usable`]-passing model with
/// `<DEL>`, `<DUP>`, and `<CNV>` entries. `gen-reads` consults this only
/// when `sv_rate_scale > 0.0`, so generation remains opt-in even though
/// the model itself is always loaded.
pub fn default_sv_model() -> SvModel {
    // Per-base rate. gnomAD-SV reports ~9000 PASS SVs per genome
    // averaged across ~14k individuals on a ~3.1Gb reference, so
    // ~3e-6 SVs/base is a reasonable composite. Real per-base rates
    // vary by chromosome (chr Y, mitochondrial) and population; this
    // value is a "typical autosome" approximation.
    let per_base_rate = 3.0e-6;

    // Type breakdown. gnomAD-SV's published per-type counts collapse
    // long-tail types into ~60% DEL / ~30% DUP / ~10% CNV after
    // dropping <INV> / breakends (which v1.10's sampler doesn't
    // generate). Probabilities sum to 1.0 within the supported set.
    let mut type_probabilities: HashMap<SvType, f64> = HashMap::new();
    type_probabilities.insert(SvType::Del, 0.60);
    type_probabilities.insert(SvType::Dup, 0.30);
    type_probabilities.insert(SvType::Cnv, 0.10);

    // Length log-normals fit by eye to the published length histograms
    // (Collins et al., Fig 2A). Means / sigmas chosen so the resulting
    // log-normal's median lands near the published per-type median:
    //   DEL: median ~270bp → mu = ln(270) ≈ 5.6, sigma ≈ 1.5
    //   DUP: median ~3kb   → mu = ln(3000) ≈ 8.0, sigma ≈ 1.7
    //   CNV: median ~50kb  → mu = ln(50000) ≈ 10.8, sigma ≈ 1.2
    // Sigmas are rough — the published distributions are heavier-
    // tailed than a pure log-normal, so we underestimate tail
    // probability slightly. Adequate for "looks-like-gnomAD" coverage
    // simulation but not for tail-event modeling.
    let mut length_log_normal: HashMap<SvType, (f64, f64)> = HashMap::new();
    length_log_normal.insert(SvType::Del, (5.6, 1.5));
    length_log_normal.insert(SvType::Dup, (8.0, 1.7));
    length_log_normal.insert(SvType::Cnv, (10.8, 1.2));

    // Copy-number distribution for <CNV> records. gnomAD's published
    // multi-allelic CNV calls cluster around CN={0,1,3,4} with a long
    // tail; rough proportions: 0 → 15%, 1 → 35%, 3 → 30%, 4 → 15%,
    // 5+ → 5%. Compressing the tail into CN=5 keeps the table small
    // without changing the simulated coverage signal much.
    let mut cnv_copy_number_distribution: HashMap<u32, f64> = HashMap::new();
    cnv_copy_number_distribution.insert(0, 0.15);
    cnv_copy_number_distribution.insert(1, 0.35);
    cnv_copy_number_distribution.insert(3, 0.30);
    cnv_copy_number_distribution.insert(4, 0.15);
    cnv_copy_number_distribution.insert(5, 0.05);

    // Homozygous fraction. gnomAD-SV reports het:hom roughly 4:1
    // averaged across the callset → homozygous_frequency ≈ 0.20.
    let homozygous_frequency = 0.20;

    SvModel {
        per_base_rate,
        type_probabilities,
        length_log_normal,
        cnv_copy_number_distribution,
        homozygous_frequency,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sv_model_is_usable() {
        // Smoke check: the bundled default must pass the same usability
        // bar the sampler keys off, otherwise `MutationModel::default()`
        // returns a model whose `sv_model: Some(...)` would silently
        // produce zero SVs even with `sv_rate_scale: 1.0`.
        assert!(default_sv_model().is_usable());
    }

    #[test]
    fn default_sv_model_carries_all_three_supported_types() {
        // Sampler skips types missing from `type_probabilities`. The
        // default must include DEL, DUP, and CNV so users see all three
        // when they enable generation.
        let m = default_sv_model();
        for t in [SvType::Del, SvType::Dup, SvType::Cnv] {
            assert!(
                m.type_probabilities.contains_key(&t),
                "missing type probability for {:?}",
                t
            );
            assert!(
                m.length_log_normal.contains_key(&t),
                "missing length distribution for {:?}",
                t
            );
        }
    }

    #[test]
    fn default_sv_model_type_probabilities_sum_to_one() {
        // Within float tolerance — if the hand-rolled probabilities
        // ever drift, fail loudly rather than silently sample from a
        // mis-normalized distribution.
        let total: f64 = default_sv_model().type_probabilities.values().sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "type_probabilities sum to {total}, expected 1.0"
        );
    }

    #[test]
    fn default_sv_model_cn_distribution_sums_to_one() {
        let total: f64 = default_sv_model()
            .cnv_copy_number_distribution
            .values()
            .sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "cnv_copy_number_distribution sums to {total}, expected 1.0"
        );
    }

    #[test]
    fn default_sv_model_parameters_in_plausible_ranges() {
        // Cheap sanity check — catches an accidental decimal-place
        // typo that would otherwise sail through type-checking.
        let m = default_sv_model();
        assert!(m.per_base_rate > 0.0 && m.per_base_rate < 1e-3);
        assert!(m.homozygous_frequency >= 0.0 && m.homozygous_frequency <= 1.0);
        for (_, (mu, sigma)) in m.length_log_normal.iter() {
            // mu is ln(median bp); medians range from ~tens to ~tens-of-thousands.
            assert!(*mu > 2.0 && *mu < 14.0, "mu={mu} out of plausible range");
            assert!(*sigma > 0.0 && *sigma < 4.0, "sigma={sigma} out of plausible range");
        }
        for cn in m.cnv_copy_number_distribution.keys() {
            // Excluding CN=2 (which is reference for diploid) is a
            // training-time invariant — the default should respect it.
            assert_ne!(*cn, 2, "default CN distribution must not include CN=2");
        }
    }
}
