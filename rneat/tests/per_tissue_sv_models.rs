//! Regression test for the bundled per-tissue cancer models (#202):
//! `tools/cosmic_per_tissue_{BRCA,skin,lung}.json.gz`.
//!
//! Each combines a **per-tissue COSMIC SNP/indel spectrum** (the #202 SNP/indel
//! half — trained from that tissue's COSMIC GenomeScreensMutant subset) with a
//! **per-tissue `sv_model`** (the #237 SV half, fitted from that tissue's PCAWG
//! donors: BRCA / skin-melanoma / lung). The test asserts both halves:
//!   1. the `sv_model` is structurally valid and cancer-shaped (the #218
//!      `cosmic_bundle.rs` canaries), so gen-reads can load it as `--tumor-model`;
//!   2. the SNP/indel half is genuinely per-tissue — the three models do not all
//!      share one pan-cancer base.

use std::path::PathBuf;

use ::common::models::mutation_model::MutationModel;
use ::common::structs::variants::SvType;

fn tissue_model_path(tissue: &str) -> PathBuf {
    PathBuf::from(format!(
        "{}/../tools/cosmic_per_tissue_{tissue}.json.gz",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

#[test]
fn bundled_per_tissue_models_carry_valid_sv_components() {
    for tissue in ["BRCA", "skin", "lung"] {
        let path = tissue_model_path(tissue);
        assert!(
            path.exists(),
            "bundled {tissue} model missing at {}",
            path.display()
        );

        let model = MutationModel::from_file(&path)
            .unwrap_or_else(|e| panic!("deserialize {tissue} model: {e:?}"));
        let sv = model
            .sv_model
            .as_ref()
            .unwrap_or_else(|| panic!("{tissue} model must carry a non-null sv_model"));

        // All six SV types present in both maps (otherwise type sampling silently
        // skips whichever is missing).
        for t in [
            SvType::Del,
            SvType::Dup,
            SvType::Cnv,
            SvType::Bnd,
            SvType::Inv,
            SvType::Ins,
        ] {
            assert!(
                sv.type_probabilities.contains_key(&t),
                "{tissue}: missing type_probabilities entry for {t:?}"
            );
            assert!(
                sv.length_log_normal.contains_key(&t),
                "{tissue}: missing length_log_normal entry for {t:?}"
            );
        }

        let prob_sum: f64 = sv.type_probabilities.values().sum();
        assert!(
            (prob_sum - 1.0).abs() < 1e-6,
            "{tissue}: type_probabilities sum to {prob_sum}, expected 1.0"
        );

        let cn_sum: f64 = sv.cnv_copy_number_distribution.values().sum();
        assert!(
            (cn_sum - 1.0).abs() < 1e-6,
            "{tissue}: cnv_copy_number_distribution sums to {cn_sum}, expected 1.0"
        );

        // Per-tumor (not corpus-aggregated) cancer rate. Per-tissue rates run a
        // bit higher than the pan-cancer 3.5e-8 (e.g. BRCA ~7e-8 from its higher
        // SV burden) but stay well under the germline-rate canary.
        assert!(
            sv.per_base_rate < 1e-6,
            "{tissue}: per_base_rate {} looks germline; cancer should be ≲ 1e-7",
            sv.per_base_rate
        );

        // BND stays a major class (the cosmic_bundle.rs ≥0.20 canary).
        let bnd = sv.type_probabilities[&SvType::Bnd];
        assert!(
            bnd >= 0.20,
            "{tissue}: BND proportion {bnd} too low for cancer"
        );
    }
}

#[test]
fn bundled_per_tissue_models_have_distinct_snp_indel_spectra() {
    // The #202 SNP/indel half: each model's SNP/indel spectrum is trained from
    // that tissue's COSMIC subset, so the three are not a shared pan-cancer base.
    // We fingerprint each by (mutation_rate, variant_dist weights) and assert the
    // three fingerprints are pairwise distinct.
    let fingerprints: Vec<(String, f64, Vec<f64>)> = ["BRCA", "skin", "lung"]
        .iter()
        .map(|tissue| {
            let model = MutationModel::from_file(&tissue_model_path(tissue))
                .unwrap_or_else(|e| panic!("deserialize {tissue} model: {e:?}"));
            let weights = model.variant_dist.weights().expect("variant_dist weights");
            (tissue.to_string(), model.mutation_rate, weights)
        })
        .collect();

    for i in 0..fingerprints.len() {
        for j in (i + 1)..fingerprints.len() {
            let (ti, ri, wi) = &fingerprints[i];
            let (tj, rj, wj) = &fingerprints[j];
            let same_rate = (ri - rj).abs() < 1e-12;
            let same_weights =
                wi.len() == wj.len() && wi.iter().zip(wj).all(|(a, b)| (a - b).abs() < 1e-12);
            assert!(
                !(same_rate && same_weights),
                "{ti} and {tj} have identical SNP/indel spectra — \
                 per-tissue training did not take effect"
            );
        }
    }
}
