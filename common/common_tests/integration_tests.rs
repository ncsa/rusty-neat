use std::path::PathBuf;
use walkdir::WalkDir;
use common::file_tools::fasta_reader::{count_fasta, read_fasta};
use common::structs::nucleotides::NucleotideSelector;
use simple_rng::NeatRng;

#[test]
fn check_data() {
    let test_fasta = PathBuf::from("test_data/H1N1.fa");

    // Verify count_fasta reports the correct number of contigs
    let contig_count = count_fasta(&test_fasta).unwrap();
    assert_eq!(contig_count, 8, "H1N1.fa should have 8 contigs");

    // read_fasta should produce one block file per contig in the temp dir
    let temp_dir = tempfile::tempdir().unwrap();
    let mut rng = NeatRng::new_from_seed(&vec![
        "integration".to_string(),
        "test".to_string(),
    ]).unwrap();
    let fasta_map = read_fasta(
        &test_fasta,
        Some(&NucleotideSelector::new()),
        350,
        &temp_dir,
        Some(&mut rng),
    ).unwrap();

    assert_eq!(fasta_map.contigs.len(), contig_count);
    assert_eq!(fasta_map.contig_order.len(), contig_count);

    let file_count = WalkDir::new(&temp_dir)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().is_file())
        .count();
    assert!(file_count >= contig_count, "Expected at least one block file per contig");
}