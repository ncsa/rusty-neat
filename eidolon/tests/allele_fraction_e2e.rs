mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

/// End-to-end contract for #398: a literal SNP carrying `INFO/AF` in the input
/// VCF must drive the alt-allele fraction of the emitted reads, so the golden
/// VCF's measured `AF` tracks the requested value — not the `{0.5, 1.0}` the
/// Genotype-based default can produce. Two sites at 0.25 and 0.75 (both well
/// clear of the 0.5 het default) prove the fraction is honored per variant.
#[test]
fn input_vcf_allele_fraction_drives_golden_af() {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_af.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
        )
        .unwrap();
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
        // Two heterozygous SNPs at distinct target fractions. If AF were ignored,
        // both would land near 0.5 (the het coin flip).
        // H1N1_HA is 1701 bp, so keep both positions well inside it.
        writeln!(f, "H1N1_HA\t500\t.\tA\tT\t60\tPASS\tAF=0.25\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t1200\t.\tA\tT\t60\tPASS\tAF=0.75\tGT\t0/1").unwrap();
    }

    // High coverage so per-site AF estimation noise is small (~0.025 sd at 0.25),
    // and mutation_rate=0 so de-novo SNPs don't confound the measured sites.
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "af_run");
    config.coverage = 300;
    config.mutation_rate = Some(0.0);
    config.produce_vcf = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_vcf_gz = work.join("af_run.vcf.gz");
    assert!(
        out_vcf_gz.exists(),
        "output VCF not produced: {out_vcf_gz:?}"
    );

    let lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_vcf_gz).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };

    let measured_af = |pos: &str| -> f64 {
        let line = lines
            .iter()
            .find(|l| {
                let mut cols = l.split('\t');
                cols.next() == Some("H1N1_HA") && cols.next() == Some(pos)
            })
            .unwrap_or_else(|| panic!("no golden VCF record at H1N1_HA:{pos}\nlines: {lines:?}"));
        // FORMAT is GT:AD:DP:AF; the sample column is last, AF is its final field.
        let sample = line.split('\t').next_back().unwrap();
        let af_str = sample.split(':').next_back().unwrap();
        af_str
            .parse::<f64>()
            .unwrap_or_else(|_| panic!("unparseable AF {af_str:?} at H1N1_HA:{pos}"))
    };

    let af_low = measured_af("500");
    let af_high = measured_af("1200");

    assert!(
        (0.15..0.35).contains(&af_low),
        "site with AF=0.25 measured {af_low:.4} (expected ~0.25, clear of the 0.5 default)"
    );
    assert!(
        (0.65..0.85).contains(&af_high),
        "site with AF=0.75 measured {af_high:.4} (expected ~0.75, clear of the 0.5 default)"
    );
}
