//! Bundled default [`SvModel`] used by [`crate::models::mutation_model::MutationModel::default`].
//!
//! Most parameters are **refit** from the v1.10 validation run against
//! the gnomAD-SV v4.1 sites VCF on chr22 (~27k surviving SVs after the
//! standard INV/CPX/CTX/BND/INS:ME filter; `tools/validate_with_gnomad.sh
//! smoke` produces the same numbers). The two parameters that can't be
//! sourced from a sites-only VCF — `homozygous_frequency` and
//! `cnv_copy_number_distribution` — are kept at their original
//! literature-derived approximations (gnomAD-SV per-sample CN lives in
//! FORMAT/SAMPLE, not INFO/CN; het:hom ≈ 4:1 is taken from the
//! Collins et al. paper, *Nature* 581, 444–451, 2020).
//!
//! chr22 is gene-dense and may not exactly match genome-wide statistics
//! — a `full`-mode validation run can supersede the per_base_rate and
//! type-probability values here. Anyone doing serious benchmarking
//! should train their own model with `gen-mut-model` against the
//! SV-rich VCF of their choice and override the embedded default by
//! passing `mutation_model: <trained>.json.gz` in their gen-reads YAML.

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
    // Per-base rate: fitted from chr22 gnomAD-SV v4.1 (27,192 surviving
    // SVs over 50,818,468 bp = 5.35e-4 per chr22 bp — but the trainer
    // reports per_base_rate as `n / total_reflen_seen`, which uses the
    // full chr22 length not the simulated denominator. The 8.917e-6
    // value below is what the trainer wrote, and is what gen-reads
    // multiplies by `contig_len * sv_rate_scale`. chr22 is gene-dense
    // and may run hotter than the genome-wide average — a `full`-mode
    // run can refine this.
    let per_base_rate = 8.917e-6;

    // Type breakdown: fitted directly from chr22. gnomAD-SV chr22 is
    // overwhelmingly deletions, with CNV essentially a rounding error
    // because most multi-allelic CNV records on chr22 collapse into a
    // small number of sites. v1's hand-rolled 60/30/10 split was off
    // by ~20 percentage points on DEL.
    let mut type_probabilities: HashMap<SvType, f64> = HashMap::new();
    type_probabilities.insert(SvType::Del, 0.807);
    type_probabilities.insert(SvType::Dup, 0.191);
    type_probabilities.insert(SvType::Cnv, 0.002);

    // Length log-normals: fitted from chr22 via MLE on ln(span). Quick
    // gut-check on the medians (= exp(mu)):
    //   DEL: median ≈ 583bp (vs hand-rolled 270bp — too small)
    //   DUP: median ≈ 735bp (vs hand-rolled 3kb — way too large)
    //   CNV: median ≈ 18kb  (vs hand-rolled 50kb — too large)
    // Sigmas pulled in / out accordingly.
    let mut length_log_normal: HashMap<SvType, (f64, f64)> = HashMap::new();
    length_log_normal.insert(SvType::Del, (6.370, 1.433));
    length_log_normal.insert(SvType::Dup, (6.600, 2.375));
    length_log_normal.insert(SvType::Cnv, (9.806, 1.112));

    // Copy-number distribution for `<CNV>` records. gnomAD-SV's sites
    // VCF leaves `INFO/CN` unset (CN varies per-sample, in
    // FORMAT/SAMPLE), so the v1.10 chr22 validation run produced an
    // empty cnv_copy_number_distribution from training. We keep the
    // literature-derived distribution here so a user loading
    // MutationModel::default() with `sv_rate_scale > 0` still gets
    // CNVs with sensible CN values. SvModel::sample_variants also has
    // a `FALLBACK_CN_DISTRIBUTION` constant that mirrors this for the
    // case where the *user-trained* model has an empty CN dist —
    // keep both in sync.
    //
    // Rough per-CN proportions from Collins et al.: 0 → 15%, 1 → 35%,
    // 3 → 30%, 4 → 15%, 5+ → 5%. Compressing the tail into CN=5 keeps
    // the table small without changing the simulated coverage signal
    // much.
    let mut cnv_copy_number_distribution: HashMap<u32, f64> = HashMap::new();
    cnv_copy_number_distribution.insert(0, 0.15);
    cnv_copy_number_distribution.insert(1, 0.35);
    cnv_copy_number_distribution.insert(3, 0.30);
    cnv_copy_number_distribution.insert(4, 0.15);
    cnv_copy_number_distribution.insert(5, 0.05);

    // Homozygous fraction: kept at 0.20 (Collins et al.; gnomAD-SV
    // het:hom ≈ 4:1). The v1.10 chr22 validation run could not fit
    // this from data because gnomAD-SV is sites-only with no GT
    // column — the reader defaults to heterozygous and the trainer
    // reports `homozygous_frequency = 0`. A per-sample callset could
    // refit this; for now the literature value is the best estimate.
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
