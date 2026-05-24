//! Statistical model for generating symbolic / structural variants
//! (`<DEL>` / `<DUP>` / `<CNV>`) de novo from a learned distribution.
//!
//! `SvModel` is fit by [`SvModel::fit_from_observations`] (called from
//! `gen-mut-model` after the SNP/indel pass completes) and lives on
//! [`crate::models::mutation_model::MutationModel`] as
//! `Option<SvModel>` — `None` whenever the training VCF didn't carry
//! enough usable SV observations, or whenever the model JSON was
//! written by a pre-v1.10 build (the field defaults to `None` via
//! `#[serde(default)]`).
//!
//! Supported types in this release: `<DEL>` / `<DUP>` / `<CNV>`. `<INV>`
//! and breakend records are dropped at fit time (logged with a count).
//! `gen-reads` does not yet consume the model — that wiring lives in a
//! later patch.

use crate::structs::variants::{Genotype, SvType, Variant};
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Minimum SV length (in bases) treated as an SV rather than an indel.
/// Matches the conventional cutoff used by Manta, DELLY, Lumpy, GRIDSS,
/// and the rest of the SV-calling ecosystem.
///
/// [`SvModel::fit_from_observations`] drops training records below this
/// span — they're indels, not SVs.
pub const MIN_SV_LENGTH_BP: usize = 50;

/// Per-type sample count required to fit a length distribution.
/// Log-normal MLE on a single observation yields `sigma = 0` — a
/// degenerate distribution. Types with fewer observations than this are
/// dropped from the trained model rather than carrying a fake `sigma`.
const MIN_OBS_FOR_LENGTH_FIT: usize = 2;

/// Learned statistical model for symbolic-SV generation.
///
/// The five fields together specify everything a sampler needs: how
/// many SVs per base (`per_base_rate`), which type each one is
/// (`type_probabilities`), how long it is (`length_log_normal`, keyed
/// by type), what copy number a `<CNV>` carries
/// (`cnv_copy_number_distribution`), and whether it's homozygous
/// (`homozygous_frequency`).
#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
pub struct SvModel {
    /// Per-base SV rate across the training corpus
    /// (= count of fit-eligible SVs / sum of training contig lengths
    /// the trainer iterated over).
    pub per_base_rate: f64,

    /// Probability that a generated SV is of a given type. Keys are
    /// limited to `SvType::Del` / `SvType::Dup` / `SvType::Cnv`;
    /// probabilities sum to ~1.0 across present keys. A type with
    /// fewer than [`MIN_OBS_FOR_LENGTH_FIT`] observations is dropped
    /// from this map (its length distribution couldn't be fit), and
    /// the remaining probabilities are renormalized.
    pub type_probabilities: HashMap<SvType, f64>,

    /// Length distribution per SV type — `(mu, sigma)` of a log-normal
    /// in log-space, fit by maximum likelihood on `ln(span)`. Same
    /// key set as `type_probabilities`.
    ///
    /// Conditioning length on type is the cheapest concession to
    /// biological realism: DELs trend smaller than DUPs / CNVs, and a
    /// single global length distribution would mis-bias whichever type
    /// the training corpus had less of.
    pub length_log_normal: HashMap<SvType, (f64, f64)>,

    /// Empirical copy-number distribution observed on `<CNV>` records
    /// (CN value → probability). Probabilities sum to ~1.0 if any
    /// CNVs in the training corpus carried `INFO/CN`; empty otherwise.
    pub cnv_copy_number_distribution: HashMap<u32, f64>,

    /// Probability that a generated SV is homozygous (vs heterozygous),
    /// estimated as the homozygous fraction across all fit-eligible
    /// observations of all types.
    pub homozygous_frequency: f64,
}

impl SvModel {
    /// Returns `true` when the model has enough learned signal to drive
    /// generation — at least one type probability and a positive
    /// per-base rate. Guards against partially-constructed or
    /// hand-rolled models reaching a sampler that would divide by zero
    /// or call a weighted pick on an empty distribution.
    pub fn is_usable(&self) -> bool {
        self.per_base_rate > 0.0 && !self.type_probabilities.is_empty()
    }

    /// Fit an `SvModel` from a set of observed SV records.
    ///
    /// Returns `None` when there isn't enough signal to build a usable
    /// model: zero in-scope observations, `total_reflen == 0`, or every
    /// type fell below [`MIN_OBS_FOR_LENGTH_FIT`]. Callers should treat
    /// `None` as "this corpus doesn't support de novo SV generation"
    /// and leave `MutationModel.sv_model` unset.
    ///
    /// Filtering applied (in order):
    ///   - records whose ALT isn't `AlternateType::Symbolic` are
    ///     ignored;
    ///   - `<INV>` / breakends / unknown tags are dropped (this release
    ///     only generates `<DEL>` / `<DUP>` / `<CNV>`);
    ///   - records with no derivable span (no `END`, no `SVLEN`) are
    ///     dropped;
    ///   - records below [`MIN_SV_LENGTH_BP`] are dropped (they're
    ///     indels, not SVs);
    ///   - types with fewer than [`MIN_OBS_FOR_LENGTH_FIT`] surviving
    ///     observations are dropped from the model entirely, so the
    ///     resulting `type_probabilities` and `length_log_normal` have
    ///     the same key set.
    ///
    /// Each filter stage logs a count so a "why is my SV model None?"
    /// debug is grep-able in the trainer output.
    pub fn fit_from_observations(
        observations: &[Variant],
        total_reflen: usize,
    ) -> Option<SvModel> {
        if total_reflen == 0 {
            warn!("SvModel: total_reflen is zero, cannot estimate per-base rate");
            return None;
        }

        let mut dropped_non_symbolic = 0usize;
        let mut dropped_unsupported_type = 0usize;
        let mut dropped_no_span = 0usize;
        let mut dropped_below_min_length = 0usize;

        // Per-type accumulators of (span, &Variant) so we can pull CN
        // and genotype back out without re-walking the source.
        let mut by_type: HashMap<SvType, Vec<(usize, &Variant)>> = HashMap::new();

        for v in observations {
            let sv = match v.alternate.as_symbolic() {
                Some(sv) => sv,
                None => {
                    dropped_non_symbolic += 1;
                    continue;
                }
            };
            match sv.sv_type {
                SvType::Del | SvType::Dup | SvType::Cnv => {}
                _ => {
                    dropped_unsupported_type += 1;
                    continue;
                }
            }
            let span = match sv.span(v.location) {
                Some(s) if s > 0 => s,
                _ => {
                    dropped_no_span += 1;
                    continue;
                }
            };
            if span < MIN_SV_LENGTH_BP {
                dropped_below_min_length += 1;
                continue;
            }
            by_type.entry(sv.sv_type).or_default().push((span, v));
        }

        // Drop types that don't have enough observations for a real
        // length fit. Renormalization happens below over what's left.
        let mut dropped_thin_types = 0usize;
        by_type.retain(|_, obs| {
            if obs.len() < MIN_OBS_FOR_LENGTH_FIT {
                dropped_thin_types += obs.len();
                false
            } else {
                true
            }
        });

        if dropped_non_symbolic > 0 {
            // Trainers should hand us symbolic-only records, but the
            // fitter is tolerant — a stray SNP just gets ignored.
            info!(
                "SvModel: ignored {} non-symbolic observation(s) (literal-base ALTs)",
                dropped_non_symbolic
            );
        }
        if dropped_unsupported_type > 0 {
            info!(
                "SvModel: dropped {} <INV>/breakend/unknown record(s) — not yet generatable",
                dropped_unsupported_type
            );
        }
        if dropped_no_span > 0 {
            info!(
                "SvModel: dropped {} symbolic SV(s) with no INFO/END or INFO/SVLEN",
                dropped_no_span
            );
        }
        if dropped_below_min_length > 0 {
            info!(
                "SvModel: dropped {} symbolic SV(s) below {}bp (treated as indels, not SVs)",
                dropped_below_min_length, MIN_SV_LENGTH_BP
            );
        }
        if dropped_thin_types > 0 {
            info!(
                "SvModel: dropped {} observation(s) from types with <{} samples \
                 (can't fit a length distribution)",
                dropped_thin_types, MIN_OBS_FOR_LENGTH_FIT
            );
        }

        if by_type.is_empty() {
            warn!(
                "SvModel: no SV type cleared the {}-observation fit bar; \
                 returning None (MutationModel.sv_model stays unset)",
                MIN_OBS_FOR_LENGTH_FIT
            );
            return None;
        }

        let total_usable: usize = by_type.values().map(|v| v.len()).sum();
        let n = total_usable as f64;

        let type_probabilities: HashMap<SvType, f64> = by_type
            .iter()
            .map(|(k, obs)| (*k, obs.len() as f64 / n))
            .collect();

        // Log-normal MLE: μ = mean(ln x), σ² = mean((ln x − μ)²).
        let length_log_normal: HashMap<SvType, (f64, f64)> = by_type
            .iter()
            .map(|(sv_type, obs)| {
                let log_lens: Vec<f64> =
                    obs.iter().map(|(s, _)| (*s as f64).ln()).collect();
                let mu = log_lens.iter().sum::<f64>() / log_lens.len() as f64;
                let variance = log_lens.iter().map(|x| (x - mu).powi(2)).sum::<f64>()
                    / log_lens.len() as f64;
                (*sv_type, (mu, variance.sqrt()))
            })
            .collect();

        // CN distribution from <CNV> records that carry INFO/CN.
        let mut cn_counts: HashMap<u32, usize> = HashMap::new();
        if let Some(cnv_obs) = by_type.get(&SvType::Cnv) {
            for (_, v) in cnv_obs {
                if let Some(sv) = v.alternate.as_symbolic()
                    && let Some(cn) = sv.copy_number
                {
                    *cn_counts.entry(cn).or_default() += 1;
                }
            }
        }
        let cn_total: usize = cn_counts.values().sum();
        let cnv_copy_number_distribution: HashMap<u32, f64> = if cn_total > 0 {
            cn_counts
                .iter()
                .map(|(k, c)| (*k, *c as f64 / cn_total as f64))
                .collect()
        } else {
            HashMap::new()
        };

        let hom_count = by_type
            .values()
            .flat_map(|obs| obs.iter())
            .filter(|(_, v)| matches!(v.genotype, Genotype::Homozygous))
            .count();
        let homozygous_frequency = hom_count as f64 / n;

        let per_base_rate = n / total_reflen as f64;

        info!(
            "SvModel: fit from {} observation(s) across {} type(s); \
             per_base_rate = {:.3e}, homozygous_frequency = {:.3}",
            total_usable,
            by_type.len(),
            per_base_rate,
            homozygous_frequency
        );

        Some(SvModel {
            per_base_rate,
            type_probabilities,
            length_log_normal,
            cnv_copy_number_distribution,
            homozygous_frequency,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide;
    use crate::structs::variants::{AlternateType, SvData, VariantType};

    #[test]
    fn default_sv_model_is_not_usable() {
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
        // Lock in the constant — downstream samplers and the trainer
        // both key off it; a silent change would shift the entire
        // SV/indel boundary across the pipeline.
        assert_eq!(MIN_SV_LENGTH_BP, 50);
    }

    #[test]
    fn populated_sv_model_round_trips_through_serde_json() {
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
        let mut only_rate = SvModel::default();
        only_rate.per_base_rate = 1e-6;
        assert!(!only_rate.is_usable());

        let mut only_types = SvModel::default();
        only_types.type_probabilities.insert(SvType::Del, 1.0);
        assert!(!only_types.is_usable());
    }

    /// Build a synthetic symbolic-SV `Variant` with a 1-based POS,
    /// a 1-based-inclusive END (or `None` to skip), the given type,
    /// genotype, and optional copy_number. Returns a value owned by
    /// the caller — tests construct a `Vec<Variant>` and pass a slice
    /// to the fitter.
    fn sv(
        pos_1based: usize,
        end_1based: Option<usize>,
        sv_type: SvType,
        genotype: Genotype,
        copy_number: Option<u32>,
    ) -> Variant {
        let raw_alt = match sv_type {
            SvType::Del => "<DEL>",
            SvType::Dup => "<DUP>",
            SvType::Cnv => "<CNV>",
            SvType::Ins => "<INS>",
            SvType::Inv => "<INV>",
            SvType::Bnd => "B[chr1:1000[",
            SvType::Unknown => "<NOVEL>",
        };
        let mut sv_data = SvData::new(raw_alt, sv_type);
        sv_data.end = end_1based;
        sv_data.copy_number = copy_number;
        let gt_str = match genotype {
            Genotype::Homozygous => "1/1",
            Genotype::Heterozygous => "0/1",
        };
        Variant {
            variant_type: VariantType::Complex,
            location: pos_1based,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(sv_data),
            genotype_str: gt_str.to_string(),
            genotype,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
        }
    }

    #[test]
    fn fit_returns_none_for_empty_observations() {
        let m = SvModel::fit_from_observations(&[], 10_000);
        assert!(m.is_none());
    }

    #[test]
    fn fit_returns_none_when_total_reflen_zero() {
        // Per-base rate would be 0/0 — `None` keeps the contract that a
        // `Some(SvModel)` is always usable.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None),
        ];
        assert!(SvModel::fit_from_observations(&obs, 0).is_none());
    }

    #[test]
    fn fit_drops_sub_50bp_records_and_unsupported_types() {
        // 5 input records, only the two ≥50bp DELs should survive into
        // the fit. The INV and BND are dropped as unsupported; the
        // 30bp DEL is dropped as an indel-disguised-as-SV.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None), // 101bp DEL — keep
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None), // 201bp DEL — keep
            sv(800, Some(829), SvType::Del, Genotype::Heterozygous, None), // 30bp — drop (sub-50)
            sv(1000, Some(2000), SvType::Inv, Genotype::Homozygous, None), // INV — drop
            sv(3000, None, SvType::Bnd, Genotype::Heterozygous, None),     // BND — drop
        ];
        let m = SvModel::fit_from_observations(&obs, 100_000).unwrap();
        assert!(m.is_usable());
        assert_eq!(m.type_probabilities.len(), 1);
        assert_eq!(m.type_probabilities.get(&SvType::Del), Some(&1.0));
        // 2 surviving SVs / 100kb = 2e-5 per base.
        assert!((m.per_base_rate - 2e-5).abs() < 1e-12);
        // 1 homozygous of 2 surviving → 0.5
        assert!((m.homozygous_frequency - 0.5).abs() < 1e-12);
    }

    #[test]
    fn fit_drops_types_with_one_observation() {
        // DEL has 2 observations → kept. DUP has 1 → dropped (can't fit σ
        // from a single sample). Probabilities renormalize over what's
        // left, so DEL becomes 1.0.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Homozygous, None),
            sv(1000, Some(1500), SvType::Dup, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 10_000).unwrap();
        assert_eq!(m.type_probabilities.len(), 1);
        assert!(m.type_probabilities.contains_key(&SvType::Del));
        assert!(!m.type_probabilities.contains_key(&SvType::Dup));
        assert!(m.length_log_normal.contains_key(&SvType::Del));
        assert!(!m.length_log_normal.contains_key(&SvType::Dup));
    }

    #[test]
    fn fit_recovers_log_normal_parameters_within_tolerance() {
        // For 4 deletions of length e^7, e^7, e^8, e^8 (i.e.
        // ln-lengths = 7, 7, 8, 8), MLE gives μ=7.5, σ²=0.25, σ=0.5.
        let lens = [
            (7.0_f64.exp() as usize), // ~1097
            (7.0_f64.exp() as usize),
            (8.0_f64.exp() as usize), // ~2981
            (8.0_f64.exp() as usize),
        ];
        let mut pos = 100usize;
        let mut obs = Vec::new();
        for &len in &lens {
            obs.push(sv(pos, Some(pos + len - 1), SvType::Del, Genotype::Heterozygous, None));
            pos += len + 1000; // separate them in coords (doesn't affect fit)
        }
        let m = SvModel::fit_from_observations(&obs, 1_000_000).unwrap();
        let (mu, sigma) = m.length_log_normal[&SvType::Del];
        // Rounding through usize introduces tiny error vs the exact
        // analytic value, but it should easily land within 5%.
        assert!((mu - 7.5).abs() < 0.05, "mu was {mu}, expected ~7.5");
        assert!((sigma - 0.5).abs() < 0.05, "sigma was {sigma}, expected ~0.5");
    }

    #[test]
    fn fit_populates_cnv_copy_number_distribution() {
        // Three CNVs: CN=1, CN=1, CN=3 → P(1)=2/3, P(3)=1/3.
        // Two CNVs without CN — they're kept but don't contribute to
        // the CN distribution (cn_counts is empty for those).
        let obs = vec![
            sv(100, Some(500), SvType::Cnv, Genotype::Heterozygous, Some(1)),
            sv(1000, Some(1500), SvType::Cnv, Genotype::Heterozygous, Some(1)),
            sv(2000, Some(2500), SvType::Cnv, Genotype::Heterozygous, Some(3)),
            sv(3000, Some(3500), SvType::Cnv, Genotype::Heterozygous, None),
            sv(4000, Some(4500), SvType::Cnv, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 100_000).unwrap();
        let cn = &m.cnv_copy_number_distribution;
        assert_eq!(cn.len(), 2);
        assert!((cn[&1] - 2.0 / 3.0).abs() < 1e-12);
        assert!((cn[&3] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn fit_returns_none_when_every_type_thin() {
        // Two records, two different types — neither type clears the
        // 2-observation bar, so the fit returns None.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(1500), SvType::Dup, Genotype::Heterozygous, None),
        ];
        assert!(SvModel::fit_from_observations(&obs, 10_000).is_none());
    }

    #[test]
    fn fit_ignores_non_symbolic_alts() {
        // A literal SNP slipped in shouldn't blow up the fitter — it's
        // not a symbolic record, so it's silently skipped.
        use crate::structs::variants::Variant as V;
        let snp = V::new(
            VariantType::SNP,
            42,
            &vec![Nucleotide::A],
            &vec![Nucleotide::C],
            &mut vec![1, 0],
        )
        .unwrap();
        let obs = vec![
            snp,
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 10_000).unwrap();
        assert_eq!(m.type_probabilities.len(), 1);
        assert!(m.is_usable());
    }
}
