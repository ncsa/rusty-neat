//! Equivalence detection: catch FPs and FNs that are different denotations of
//! the same edit.
//!
//! Port of NEAT 2.1's `vcf_compare.py` algorithm:
//!
//! For each raw FP, take a ±`window_bp` byte window around its 0-based
//! reference position. Find every other raw FP and raw FN whose position
//! falls inside that window. Apply each subset to a copy of the window
//! (left-to-right, tracking the cumulative length delta from earlier indels).
//! If the FP-applied window and FN-applied window are byte-identical, the
//! two subsets are alternative spellings of the same edit — promote every
//! consumed FN to TP and remove the consumed FPs from FP, and ditto for FNs.
//!
//! Why per-FP windows: left/right-aligned indels and equivalent multi-base
//! substitutions usually anchor on a single called-but-mis-denoted variant.
//! Centering each window on an FP catches the canonical case without an
//! O(n²) all-pairs scan.
use eidolon_core::structs::variants::Variant;
use std::collections::HashSet;

/// Indices of FP and FN entries that were "consumed" by the sweep — i.e.,
/// found to be denotation-different alternate spellings of golden variants.
pub struct SweepResult {
    /// Number of FNs that got promoted to TP. Each consumed FN counts once.
    pub equivalents_promoted: usize,
    /// Indices into the input `fps` slice to delete.
    pub fp_consumed: Vec<usize>,
    /// Indices into the input `fns` slice to delete.
    pub fn_consumed: Vec<usize>,
}

/// Run the equivalence sweep for one contig.
///
/// `reference` is the contig's raw byte sequence (case-preserved as read from
/// the FASTA). `fps` and `fns` are the raw, contig-restricted variant lists.
/// `window_bp` is the half-width of the equivalence window (e.g. 50 → 101 bp
/// total).
///
/// The comparison is case-folded to ASCII uppercase so that a reference with
/// soft-masked bases compares cleanly against ALT alleles emitted in
/// canonical uppercase.
pub fn sweep(fps: &[Variant], fns: &[Variant], reference: &[u8], window_bp: usize) -> SweepResult {
    let mut fp_consumed: HashSet<usize> = HashSet::new();
    let mut fn_consumed: HashSet<usize> = HashSet::new();
    let mut promoted = 0usize;

    for fp_idx in 0..fps.len() {
        if fp_consumed.contains(&fp_idx) {
            continue;
        }
        let center = fps[fp_idx].location;
        // Convert VCF 1-based pos to 0-based, then expand by `window_bp` on
        // each side. End is exclusive, capped at reference length.
        let center0 = center.saturating_sub(1);
        let start = center0.saturating_sub(window_bp);
        let end = (center0 + window_bp + 1).min(reference.len());
        if start >= end {
            continue;
        }
        let window = &reference[start..end];

        let fps_in_window = variants_in_window(fps, &fp_consumed, start, end);
        let fns_in_window = variants_in_window(fns, &fn_consumed, start, end);
        if fns_in_window.is_empty() {
            // Without a candidate FN partner, this window cannot be equivalent
            // (applying the FP alone shifts bytes; nothing on the FN side matches).
            continue;
        }

        let altered_fp = apply_variants(window, start, &fps_in_window, fps);
        let altered_fn = apply_variants(window, start, &fns_in_window, fns);
        if eq_ignore_case(&altered_fp, &altered_fn) {
            for i in &fps_in_window {
                fp_consumed.insert(*i);
            }
            promoted += fns_in_window.len();
            for i in &fns_in_window {
                fn_consumed.insert(*i);
            }
        }
    }

    SweepResult {
        equivalents_promoted: promoted,
        fp_consumed: fp_consumed.into_iter().collect(),
        fn_consumed: fn_consumed.into_iter().collect(),
    }
}

/// Indices of variants whose 0-based position falls inside `[start, end)`.
fn variants_in_window(
    variants: &[Variant],
    consumed: &HashSet<usize>,
    start: usize,
    end: usize,
) -> Vec<usize> {
    let mut out = Vec::new();
    for (i, v) in variants.iter().enumerate() {
        if consumed.contains(&i) {
            continue;
        }
        let pos0 = v.location.saturating_sub(1);
        if pos0 >= start && pos0 < end {
            out.push(i);
        }
    }
    out
}

/// Apply the variants at `indices` (resolved in `variants`) to a copy of
/// `window` (whose 0-based reference offset is `window_start`). Returns the
/// resulting byte sequence.
fn apply_variants(
    window: &[u8],
    window_start: usize,
    indices: &[usize],
    variants: &[Variant],
) -> Vec<u8> {
    let mut ordered: Vec<&Variant> = indices.iter().map(|&i| &variants[i]).collect();
    ordered.sort_by_key(|v| v.location);

    let mut result: Vec<u8> = window.to_vec();
    let mut adj: isize = 0;
    for v in ordered {
        let offset_signed = v.location as isize - 1 - window_start as isize + adj;
        if offset_signed < 0 {
            continue;
        }
        let offset = offset_signed as usize;
        let ref_len = v.reference.len();
        if offset + ref_len > result.len() {
            continue;
        }
        // filter_vcf strips symbolic ALTs into the `symbolic` skip bucket before
        // they reach the equivalence sweep, so anything here must be literal.
        debug_assert!(
            v.alternate.is_literal(),
            "symbolic ALT reached apply_variants at location {}",
            v.location
        );
        let alt_bytes: Vec<u8> = v
            .alternate
            .as_literal()
            .unwrap()
            .iter()
            .map(|n| Into::<char>::into(*n) as u8)
            .collect();
        let alt_len = alt_bytes.len();
        result.splice(offset..offset + ref_len, alt_bytes);
        adj += alt_len as isize - ref_len as isize;
    }
    result
}

fn eq_ignore_case(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|(x, y)| x.eq_ignore_ascii_case(y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::structs::variants::{AlternateType, Provenance};
    use eidolon_core::structs::{
        nucleotides::Nucleotide,
        variants::{Genotype, VariantType},
    };

    fn variant(
        loc: usize,
        ref_seq: &[Nucleotide],
        alt_seq: &[Nucleotide],
        vt: VariantType,
    ) -> Variant {
        Variant {
            variant_type: vt,
            location: loc,
            reference: ref_seq.to_vec(),
            alternate: AlternateType::Literal(alt_seq.to_vec()),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: Some("PASS".to_string()),
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        }
    }

    /// Reference: ACGTACGTACGT. A deletion at pos 5 (1-based) of "CGTA" can
    /// equivalently be denoted as a deletion at pos 1 of "ACGT" — same alt
    /// sequence either way. The sweep should catch this.
    #[test]
    fn detects_left_vs_right_aligned_deletion() {
        let reference = b"ACGTACGTACGT";
        // FN (golden): delete CGTA starting at 1-based pos 5 (REF=ACGTA, ALT=A)
        let fn_v = variant(
            5,
            &[
                Nucleotide::A,
                Nucleotide::C,
                Nucleotide::G,
                Nucleotide::T,
                Nucleotide::A,
            ],
            &[Nucleotide::A],
            VariantType::Deletion,
        );
        // FP (called): same deletion left-aligned at pos 1 (REF=ACGTA, ALT=A)
        let fp_v = variant(
            1,
            &[
                Nucleotide::A,
                Nucleotide::C,
                Nucleotide::G,
                Nucleotide::T,
                Nucleotide::A,
            ],
            &[Nucleotide::A],
            VariantType::Deletion,
        );
        let res = sweep(&[fp_v], &[fn_v], reference, 50);
        assert_eq!(res.equivalents_promoted, 1);
        assert_eq!(res.fp_consumed.len(), 1);
        assert_eq!(res.fn_consumed.len(), 1);
    }

    /// SNPs at the same position with the same ALT but different REF cases
    /// (uppercase vs soft-masked lowercase in the reference) should still
    /// equivalence-match because the case-fold compare absorbs the case
    /// mismatch.
    #[test]
    fn handles_case_mismatch_on_reference() {
        let reference = b"ACgTACGT"; // pos 3 is lowercase 'g'
        let fn_v = variant(3, &[Nucleotide::G], &[Nucleotide::A], VariantType::SNP);
        let fp_v = variant(3, &[Nucleotide::G], &[Nucleotide::A], VariantType::SNP);
        let res = sweep(&[fp_v], &[fn_v], reference, 50);
        assert_eq!(res.equivalents_promoted, 1);
    }

    /// Two FPs at adjacent positions that together produce the same alt
    /// sequence as a single FN multi-base substitution should
    /// equivalence-match. The Variant struct collapses MNP-style records
    /// into the `Complex` type.
    #[test]
    fn detects_compound_snp_vs_mnp() {
        let reference = b"ACGTACGT";
        // FN: pos 3 REF=GT → ALT=AA (a complex multi-base substitution)
        let fn_v = variant(
            3,
            &[Nucleotide::G, Nucleotide::T],
            &[Nucleotide::A, Nucleotide::A],
            VariantType::Complex,
        );
        // FP1: pos 3 G→A, FP2: pos 4 T→A
        let fp1 = variant(3, &[Nucleotide::G], &[Nucleotide::A], VariantType::SNP);
        let fp2 = variant(4, &[Nucleotide::T], &[Nucleotide::A], VariantType::SNP);
        let res = sweep(&[fp1, fp2], &[fn_v], reference, 50);
        assert_eq!(res.equivalents_promoted, 1);
        assert_eq!(res.fp_consumed.len(), 2);
        assert_eq!(res.fn_consumed.len(), 1);
    }

    /// Genuinely different variants should not equivalence-match.
    #[test]
    fn unrelated_variants_are_not_equivalent() {
        let reference = b"ACGTACGTACGT";
        let fn_v = variant(3, &[Nucleotide::G], &[Nucleotide::A], VariantType::SNP);
        let fp_v = variant(8, &[Nucleotide::T], &[Nucleotide::C], VariantType::SNP);
        let res = sweep(&[fp_v], &[fn_v], reference, 50);
        assert_eq!(res.equivalents_promoted, 0);
        assert!(res.fp_consumed.is_empty());
        assert!(res.fn_consumed.is_empty());
    }

    /// With no FNs in the window, a lonely FP cannot be equivalent.
    #[test]
    fn lonely_fp_is_not_consumed() {
        let reference = b"ACGTACGTACGT";
        let fp_v = variant(5, &[Nucleotide::A], &[Nucleotide::T], VariantType::SNP);
        let res = sweep(&[fp_v], &[], reference, 50);
        assert_eq!(res.equivalents_promoted, 0);
        assert!(res.fp_consumed.is_empty());
    }

    /// Variant near the contig start: window should clamp at 0 and the
    /// algorithm should still run.
    #[test]
    fn handles_window_clamp_at_start() {
        let reference = b"ACGTACGT";
        let fn_v = variant(1, &[Nucleotide::A], &[Nucleotide::G], VariantType::SNP);
        let fp_v = variant(1, &[Nucleotide::A], &[Nucleotide::G], VariantType::SNP);
        let res = sweep(&[fp_v], &[fn_v], reference, 50);
        assert_eq!(res.equivalents_promoted, 1);
    }
}
