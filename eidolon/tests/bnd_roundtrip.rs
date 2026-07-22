mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_with_bnd_variant_round_trips_to_vcf() {
    // End-to-end contract: a BND variant in the input VCF must be preserved
    // verbatim in the output golden VCF. BNDs currently have a coverage
    // multiplier of 1.0, so they shouldn't affect read depth.

    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_bnd.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        )
        .unwrap();
        // BND variant at H1N1_HA:500
        writeln!(
            f,
            "H1N1_HA\t500\t.\tG\tG]H1N1_HA:1500]\t60\tPASS\t.\tGT\t0/1"
        )
        .unwrap();
        // Another form of BND
        writeln!(
            f,
            "H1N1_HA\t600\t.\tC\t[H1N1_HA:1600[C\t60\tPASS\t.\tGT\t0/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "bnd_run");
    config.coverage = 10;
    config.produce_vcf = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_vcf_gz = work.join("bnd_run.vcf.gz");
    assert!(
        out_vcf_gz.exists(),
        "output VCF not produced: {out_vcf_gz:?}"
    );

    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_vcf_gz).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };

    // Check first BND
    assert!(
        vcf_lines
            .iter()
            .any(|l| l.contains("H1N1_HA\t500\t") && l.contains("G]H1N1_HA:1500]")),
        "expected first BND record in output VCF; got: {vcf_lines:?}"
    );

    // Check second BND
    assert!(
        vcf_lines
            .iter()
            .any(|l| l.contains("H1N1_HA\t600\t") && l.contains("[H1N1_HA:1600[C")),
        "expected second BND record in output VCF; got: {vcf_lines:?}"
    );
}
