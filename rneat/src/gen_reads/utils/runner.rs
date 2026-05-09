use common::file_tools::file_io::create_output_file;
use common::structs::variants::VariantType;
use tempfile;
use log::{info, debug, warn, error};
use common::rng::NeatRng;
use std::collections::HashMap;
use std::path::PathBuf;
use std::io::BufWriter;

use flate2::{
    Compression,
    write::GzEncoder
};
use common::models::gc_bias_model::GcBiasModel;
use crate::{
    common::{
        file_tools::{
            bam_writer::BamWriter,
            bed_reader::read_bed,
            fasta_reader::{map_buffer, apply_n_substitution},
            fasta_stream::{FastaStream, scan_fasta_lengths},
            fastq_tools::{
                combine_temp_fastqs,
                write_block_fastq
            },
            file_io::{append_to_file, VectorBuffer},
            vcf_tools::{read_vcf, write_vcf}
        }, models::{
            fragment_length::FragmentLengthModel,
            mutation_model::MutationModel,
            quality_scores::QualityScoreModel,
            sequencing_error_model::SequencingErrorModel
        }, structs::{
            bed_record::BedRecord,
            fasta_map::{SequenceBlock, SequenceMap, RegionType},
            mutated_map::MutatedMap,
            nucleotides::NucleotideSelector,
            variants::Variant
        }
    },
    gen_reads::{
        errors::{
            GenerateReadsError,
        },
        utils::{
            config::RunConfiguration,
            generate_variants::generate_variants,
            generate_fragments::{
                generate_fragments,
                generate_weighted_fragments,
            },
        }
    }
};

pub fn run_neat(config: &Box<RunConfiguration>, rng: &mut NeatRng) -> Result<Vec<PathBuf>, GenerateReadsError> {
    // We'll need a temp dir to store file fragments
    let working_dir = tempfile::tempdir().unwrap();
    info!("Created temp dir at {:?}", working_dir);
    // Load models that will be used for the runs
    //
    // load mutation model
    info!("Generate mutation model");
    let mutation_model = {
        match &config.mutation_model {
            Some(filename) => {
                MutationModel::from_file(&filename)?
            },
            None => {
                MutationModel::default()?
            }
        }
    };
    // Putting code here to determine whether to use the model's mutation_rate, or the custom
    // mutation rate from the config.
    let run_mutation_rate = match config.mutation_rate {
        Some(rate) => rate,
        None => mutation_model.mutation_rate
    };

    // Fragment Length Model
    info!("Generate fragment length model");
    let fragment_length_model: FragmentLengthModel = {
        // FragmentLengthModel is an enum that allows us to choose one of two models.
        match &config.fragment_model {
            Some(filename) => {
                FragmentLengthModel::discrete_from_file(&filename)?.into()
            },
            None => {
                match config.fragment_mean {
                    Some(mean) => {
                        FragmentLengthModel::new_normal(
                            mean,
                            // Config should already have caught issues with missing data
                            config.fragment_st_dev.unwrap()
                        )?.into()
                    },
                    None => {
                        FragmentLengthModel::default()?.into()
                    },
                }
            }
        }
    };

    // Load Sequencing Error Model
    info!("Generate sequencing error model");
    let seq_error_model: SequencingErrorModel = {
        match &config.sequence_error_model {
            Some(filename) => {
                SequencingErrorModel::from_file(&filename)?
            },
            None => {
                SequencingErrorModel::default()?
            }
        }
    };

    // Load Quality Score Model
    info!("Generate quality score model");
    let quality_score_model: QualityScoreModel = {
        match &config.quality_score_model {
            Some(filename) => QualityScoreModel::from_file(filename)?,
            None => QualityScoreModel::default()?,
        }
    };

    // load GC-bias model, if provided.
    let gc_bias_model = match &config.gc_bias_model {
        Some(path) => {
            info!("Loading GC Bias model: {}", path.display());
            GcBiasModel::from_file(path)?
        },
        None => GcBiasModel::default(),
    };

    // This model provides an alternate base for N based on a simple equal weight distribution.
    info!("Initialize Nucleotide selector");
    let nuc_sub_model: NucleotideSelector = {
        NucleotideSelector::new()
    };

    // Load the target BED if provided. Regions outside the BED will be skipped entirely
    // during read-generation, so no reads or variants are produced for those areas.
    let target_bed = match &config.target_bed {
        Some(path) => {
            info!("Loading target BED: {}", path.display());
            Some(read_bed(path)?)
        },
        None => None,
    };

    // Load input VCF variants if provided. Complex variants (multi-base ref AND alt that are
    // not a simple indel) are skipped with a warning — they have no application code yet.
    let input_variants: Option<HashMap<String, Vec<Variant>>> = match &config.input_vcf {
        Some(path) => {
            info!("Loading input VCF: {}", path.display());
            let raw = read_vcf(path.to_path_buf())?;
            Some(filter_input_vcf(raw))
        },
        None => None,
    };

    // Pre-scan FASTA for contig lengths needed by the BAM header. This is a fast
    // single-pass read that does not store sequences.
    info!("Reading fasta file: {}", &config.reference.display());
    let mut bam_writer: Option<BamWriter> = if config.produce_bam {
        let contigs = scan_fasta_lengths(&config.reference)?;
        let bam_path = config.output_bam.as_ref().expect("output_bam must be set when produce_bam is true");
        Some(BamWriter::new(bam_path, &contigs)?)
    } else {
        None
    };

    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    let mut all_fastq_files: HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)> = HashMap::new();
    let mut contig_order: Vec<String> = Vec::new();
    let mut fasta_lengths: HashMap<String, usize> = HashMap::new();

    info!("Generating simulated dataset");

    let fasta = FastaStream::open(&config.reference)?;
    for result in fasta {
        let (contig_name, mut sequence) = result?;
        let contig_len = sequence.len();
        // Record all contigs for VCF header and output ordering.
        contig_order.push(contig_name.clone());
        fasta_lengths.insert(contig_name.clone(), contig_len);
        debug!("Processing {}", contig_name);

        // BED check before any expensive work on this contig.
        if let Some(bed) = &target_bed {
            if !bed.contains_key(&contig_name) {
                debug!("Skipping {} — not in target BED", contig_name);
                continue;
            }
        }

        if sequence.is_empty() {
            warn!("Contig {} has empty sequence, skipping", contig_name);
            continue;
        }

        // Apply N-substitution in-place, then build an in-memory SequenceBlock.
        apply_n_substitution(&mut sequence, &nuc_sub_model, rng)?;
        let sequence_map = map_buffer(&sequence)?;
        let current_block = SequenceBlock {
            contig: contig_name.clone(),
            ref_start: 0,
            ref_end: contig_len,
            sequence,
            sequence_map,
        };

        debug!("    > Generating bias map.");
        let raw_regions = current_block.get_non_n_regions()?;
        let regions_of_interest: Vec<SequenceMap> = if let Some(bed) = &target_bed {
            let contig_beds = bed.get(&contig_name).map(|v| v.as_slice()).unwrap_or(&[]);
            intersect_with_bed(&raw_regions, contig_beds, 0)
        } else {
            raw_regions.into_iter().cloned().collect()
        };
        if regions_of_interest.is_empty() {
            continue;
        }
        let mut bias_map: Vec<f64> = vec![0.0; current_block.get_len()];
        let mut num_mutable_area = 0.0_f64;
        for region in &regions_of_interest {
            for j in region.start..region.end {
                bias_map[j] = 1.0;
                num_mutable_area += 1.0;
            }
        }

        // VCF POS is 1-based; with ref_start=0 the local position is simply location-1.
        let block_end = contig_len;
        let mut block_variants: Vec<Variant> = if let Some(iv) = &input_variants {
            iv.get(&contig_name)
                .map(|vs| {
                    vs.iter()
                        .filter(|v| {
                            let pos0 = v.location.saturating_sub(1);
                            pos0 < block_end
                        })
                        .filter_map(|v| {
                            let local_pos = v.location - 1;
                            if local_pos >= bias_map.len() {
                                return None;
                            }
                            if bias_map[local_pos] > 0.0 {
                                bias_map[local_pos] = 0.0;
                                num_mutable_area -= 1.0;
                            }
                            let mut v2 = v.clone();
                            v2.location = local_pos;
                            Some(v2)
                        })
                        .collect()
                })
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        debug!("Seeded {} user variant(s) into contig {}", block_variants.len(), contig_name);

        let num_mutations = (num_mutable_area * run_mutation_rate).trunc() as usize;
        debug!("Adding {} mutations to contig {}", num_mutations, contig_name);
        let mut max_del_len = 0;
        for v in &block_variants {
            if v.variant_type == VariantType::Deletion && v.reference.len() > 1 {
                max_del_len = max_del_len.max(v.reference.len() - 1);
            }
        }
        if num_mutations > 0 {
            let result: Option<Vec<Variant>> = generate_variants(
                &current_block,
                &bias_map,
                &mutation_model,
                num_mutations,
                config.ploidy,
                rng
            )?;
            match result {
                Some(vec) => {
                    for variant in vec {
                        if variant.variant_type == VariantType::Deletion {
                            if variant.reference.len() - 1 > max_del_len {
                                max_del_len = variant.reference.len() - 1
                            }
                        }
                        block_variants.push(variant);
                    }
                },
                None => {}
            }
        }

        let block_fragments: Vec<(usize, usize)> = {
            let mut block_frags = Vec::new();
            for (frag_start, frag_end) in regions_of_interest.into_iter().map(|r| (r.start, r.end)) {
                let frags = if gc_bias_model.is_uniform() {
                    generate_fragments(
                        frag_end - frag_start,
                        config.read_len,
                        max_del_len,
                        frag_start,
                        config.coverage,
                        &fragment_length_model,
                        rng,
                    )?
                } else {
                    generate_weighted_fragments(
                        &current_block,
                        frag_start,
                        frag_end,
                        config.read_len,
                        max_del_len,
                        config.coverage,
                        &gc_bias_model,
                        &fragment_length_model,
                        config.gc_bias_normalize_coverage,
                        rng,
                    )?
                };
                if !frags.is_empty() {
                    block_frags.extend_from_slice(&frags);
                }
            }
            block_frags
        };

        let mutated_map = MutatedMap::from_interval(0, block_end, block_variants)?;

        let mut contig_files_r1: Vec<PathBuf> = Vec::new();
        let mut contig_files_r2: Vec<PathBuf> = Vec::new();

        if config.produce_fastq {
            let mut file_to_write_1 = PathBuf::from(working_dir.path());
            file_to_write_1.push(format!(
                "temp_{}_{:010}_{:010}_r1_tmp.fastq.gz",
                contig_name,
                current_block.get_start()?,
                current_block.get_end()?,
            ));
            let file1 = append_to_file(&file_to_write_1)?;
            let writer1 = BufWriter::new(&file1);
            let mut buffer1 = GzEncoder::new(writer1, Compression::default());
            let read_name_prefix = format!("RNEAT_generated_{}", current_block.contig);
            if config.paired_ended {
                let mut file_to_write_2 = PathBuf::from(working_dir.path());
                file_to_write_2.push(format!(
                    "temp_{}_{:010}_{:010}_r2_tmp.fastq.gz",
                    contig_name,
                    current_block.get_start()?,
                    current_block.get_end()?,
                ));
                let file2 = append_to_file(&file_to_write_2)?;
                let writer2 = BufWriter::new(&file2);
                let mut buffer2 = GzEncoder::new(writer2, Compression::default());
                debug!("Writing paired-ended contig fastq files");
                write_block_fastq(
                    block_fragments,
                    &mutated_map,
                    &current_block,
                    true,
                    &mut buffer1,
                    &mut buffer2,
                    config.read_len,
                    &read_name_prefix,
                    &quality_score_model,
                    &seq_error_model,
                    rng,
                    bam_writer.as_mut(),
                )?;
                contig_files_r1.push(file_to_write_1);
                contig_files_r2.push(file_to_write_2);
            } else {
                debug!("Writing single-ended contig fastq file");
                let dummy_data: VectorBuffer = VectorBuffer::new();
                let mut buffer2 = GzEncoder::new(dummy_data, Compression::default());
                write_block_fastq(
                    block_fragments,
                    &mutated_map,
                    &current_block,
                    false,
                    &mut buffer1,
                    &mut buffer2,
                    config.read_len,
                    &read_name_prefix,
                    &quality_score_model,
                    &seq_error_model,
                    rng,
                    bam_writer.as_mut(),
                )?;
                contig_files_r1.push(file_to_write_1);
            }
        }

        if let Some(bam) = bam_writer.as_mut() {
            bam.flush_all()?;
        }
        mutated_maps.insert(contig_name.clone(), vec![mutated_map]);
        all_fastq_files.insert(contig_name, (contig_files_r1, contig_files_r2));
    }

    info!("Read generation complete, producing output files");

    if config.produce_fastq {
        info!("Producing final fastq(s) file(s)");

        // Warn when the genome is large enough that loading all reads into memory
        // for the global shuffle may be impractical. For human-scale genomes (>500 Mbp)
        // consider post-processing with `seqkit shuffle` instead.
        let total_genome_bp: usize = fasta_lengths.values().sum();
        if total_genome_bp > 500_000_000 {
            warn!(
                "Genome is {:.1} Gbp. The global FASTQ shuffle loads all reads into memory. \
                 For large genomes, consider running `seqkit shuffle` on the output instead.",
                total_genome_bp as f64 / 1e9
            );
        }

        // Collect all per-block temp files across all contigs so we can shuffle
        // reads globally rather than per-contig.
        let mut all_r1: Vec<PathBuf> = Vec::new();
        let mut all_r2: Vec<PathBuf> = Vec::new();
        for (r1_files, r2_files) in all_fastq_files.into_values() {
            all_r1.extend(r1_files);
            all_r2.extend(r2_files);
        }

        match &config.output_fastq_1 {
            Some(filename1) => {
                create_output_file(filename1, config.overwrite_output)?;
                if config.paired_ended {
                    match &config.output_fastq_2 {
                        Some(filename2) => {
                            create_output_file(filename2, config.overwrite_output)?;
                            combine_temp_fastqs(
                                all_r1,
                                all_r2,
                                filename1,
                                Some(filename2),
                                config.shuffle_fastq,
                                rng,
                            )?;
                        },
                        None => {
                            error!("Produce fastq true and paired-ended true, but output_fastq_2 was missing.");
                            return Err(GenerateReadsError::ConfigError)
                        }
                    }
                } else {
                    combine_temp_fastqs(
                        all_r1,
                        vec![],
                        filename1,
                        None,
                        config.shuffle_fastq,
                        rng,
                    )?;
                }
            },
            None => {
                error!("Produce fastq true but output_fastq_1 was missing.");
                return Err(GenerateReadsError::ConfigError)
            },
        }
    }

    let mut files_written = Vec::new();
    if config.paired_ended {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", &filename1);
            files_written.push(filename1.clone());
            if let Some(filename2) = &config.output_fastq_2 {
                info!("Successfully wrote fastq file: {:?}", &filename2);
                files_written.push(filename2.clone());
            }
        }
    } else {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", &filename1);
            files_written.push(filename1.clone());
        }
    }

    if let Some(bam_path) = &config.output_bam {
        if config.produce_bam {
            info!("Successfully wrote BAM file: {:?}", bam_path);
            files_written.push(bam_path.clone());
        }
    }

    if let Some(filename) = &config.output_vcf {
        info!("Writing output vcf file");
        let result = write_vcf(
            &mutated_maps,
            &contig_order,
            &fasta_lengths,
            &config.reference,
            config.overwrite_output,
            &filename,
        );
        match result {
            Ok(()) => {
                info!("Successfully wrote vcf file: {:?}", filename);
                files_written.push(filename.clone());
            },
            Err(error) => {
                error!("Error writing vcf file!");
                return Err(GenerateReadsError::IoError(error))
            },
        }
    }
    Ok(files_written.clone())
}

/// Removes Complex variants from the input map, emitting a warning for each one.
/// SNPs, insertions, and deletions are kept as-is.
fn filter_input_vcf(
    raw: HashMap<String, Vec<Variant>>,
) -> HashMap<String, Vec<Variant>> {
    let mut out: HashMap<String, Vec<Variant>> = HashMap::new();
    for (contig, variants) in raw {
        let mut kept = Vec::new();
        for v in variants {
            if v.variant_type == VariantType::Complex {
                warn!(
                    "Skipping complex variant at {}:{} (multi-base REF and ALT that is not \
                     a simple indel) — not yet supported",
                    contig, v.location
                );
            } else {
                kept.push(v);
            }
        }
        if !kept.is_empty() {
            out.insert(contig, kept);
        }
    }
    out
}

/// Intersects a set of non-N regions (block-local coordinates) with BED records
/// (global contig coordinates) and returns only the overlapping sub-intervals,
/// still expressed in block-local coordinates.
///
/// `block_offset` is `SequenceBlock::ref_start` — the global contig position at
/// which this block begins.
fn intersect_with_bed(
    regions: &[&SequenceMap],
    bed_records: &[BedRecord],
    block_offset: usize,
) -> Vec<SequenceMap> {
    let mut out = Vec::new();
    for region in regions {
        let global_start = region.start + block_offset;
        let global_end = region.end + block_offset;
        for bed in bed_records {
            let isect_start = global_start.max(bed.start);
            let isect_end = global_end.min(bed.end);
            if isect_start < isect_end {
                out.push(SequenceMap::from(
                    RegionType::NonNRegion,
                    isect_start - block_offset,
                    isect_end - block_offset,
                ));
            }
        }
    }
    out
}