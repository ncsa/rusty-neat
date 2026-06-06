//! Regression test for the bundled per-tissue cancer SV models (#202):
//! `tools/cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz`.
//!
//! Each is the pan-cancer COSMIC SNP/indel base with a tissue-specific `sv_model`
//! fitted from that tissue's PCAWG donors (BRCA / skin-melanoma / lung). They
//! must carry the same structurally-valid `sv_model` shape the pan-cancer bundle
//! does (issue #218's `cosmic_bundle.rs` checks), so gen-reads can load them as
//! `--tumor-model`.

use std::path::PathBuf;

use ::common::models::mutation_model::MutationModel;
use ::common::structs::variants::SvType;

fn tissue_model_path(tissue: &str) -> PathBuf {
    PathBuf::from(format!(
        "{}/../tools/cosmic_pancancer_sv_{tissue}.json.gz",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

#[test]
fn bundled_per_tissue_models_carry_valid_sv_components() {
    for tissue in ["BRCA", "skin", "lung"] {
        let path = tissue_model_path(tissue);
        assert!(path.exists(), "bundled {tissue} model missing at {}", path.display());

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
        assert!(bnd >= 0.20, "{tissue}: BND proportion {bnd} too low for cancer");
    }
}
