use rayon::prelude::*;
use common::file_tools::file_io::create_output_file;
use common::structs::variants::VariantType;
use tempfile;
use log::{info, debug, warn, error};
use common::rng::NeatRng;
use std::collections::HashMap;
use std::path::PathBuf;
use std::io::BufWriter;
use std::sync::Arc;

use flate2::{
    Compression,
    write::GzEncoder
};
use common::models::gc_bias_model::GcBiasModel;
use crate::{
    common::{
        file_tools::{
            bam_writer::{BamBodyWriter, BamContext, BamRecordStager, concat_temp_bams},
            bed_reader::read_bed,
            fasta_stream::{FastaStream, map_buffer, apply_n_substitution, scan_fasta_lengths},
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
            sequence_block::{SequenceBlock, SequenceMap, RegionType},
            mutated_map::MutatedMap,
            nucleotides::{NucleotideSelector, Nucleotide},
            variants::Variant
        }
    },
    gen_reads::{
        errors::GenerateReadsError,
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

struct ContigContext<'a> {
    config: &'a RunConfiguration,
    target_bed: &'a Option<HashMap<String, Vec<BedRecord>>>,
    input_variants: &'a Option<HashMap<String, Vec<Variant>>>,
    nuc_sub_model: &'a NucleotideSelector,
    mutation_model: &'a MutationModel,
    run_mutation_rate: f64,
    fragment_length_model: &'a FragmentLengthModel,
    gc_bias_model: &'a GcBiasModel,
    quality_score_model: &'a QualityScoreModel,
    seq_error_model: &'a SequencingErrorModel,
    working_dir: &'a std::path::Path,
    base_rng: NeatRng,
    bam_context: Option<Arc<BamContext>>,
}

struct ContigResult {
    idx: usize,
    name: String,
    len: usize,
    data: Option<ProcessedContigData>,
}

struct ProcessedContigData {
    mutated_map: MutatedMap,
    r1_files: Vec<PathBuf>,
    r2_files: Vec<PathBuf>,
    bam_body_file: Option<PathBuf>,
}

pub fn run_neat(config: &Box<RunConfiguration>, rng: &mut NeatRng) -> Result<Vec<PathBuf>, GenerateReadsError> {
    let working_dir = tempfile::tempdir().unwrap();
    info!("Created temp dir at {:?}", working_dir);

    info!("Generate mutation model");
    let mutation_model = {
        match &config.mutation_model {
            Some(filename) => MutationModel::from_file(&filename)?,
            None => MutationModel::default()?,
        }
    };
    let run_mutation_rate = match config.mutation_rate {
        Some(rate) => rate,
        None => mutation_model.mutation_rate
    };

    info!("Generate fragment length model");
    let fragment_length_model: FragmentLengthModel = {
        match &config.fragment_model {
            Some(filename) => FragmentLengthModel::discrete_from_file(&filename)?.into(),
            None => {
                match config.fragment_mean {
                    Some(mean) => {
                        FragmentLengthModel::new_normal(
                            mean,
                            config.fragment_st_dev.unwrap()
                        )?.into()
                    },
                    None => FragmentLengthModel::default()?.into(),
                }
            }
        }
    };

    info!("Generate sequencing error model");
    let seq_error_model: SequencingErrorModel = {
        match &config.sequence_error_model {
            Some(filename) => SequencingErrorModel::from_file(&filename)?,
            None => SequencingErrorModel::default()?,
        }
    };

    info!("Generate quality score model");
    let quality_score_model: QualityScoreModel = {
        match &config.quality_score_model {
            Some(filename) => QualityScoreModel::from_file(filename)?,
            None => QualityScoreModel::default()?,
        }
    };

    let gc_bias_model = match &config.gc_bias_model {
        Some(path) => {
            info!("Loading GC Bias model: {}", path.display());
            GcBiasModel::from_file(path)?
        },
        None => GcBiasModel::default(),
    };

    info!("Initialize Nucleotide selector");
    let nuc_sub_model: NucleotideSelector = NucleotideSelector::new();

    let target_bed = match &config.target_bed {
        Some(path) => {
            info!("Loading target BED: {}", path.display());
            Some(read_bed(path)?)
        },
        None => None,
    };

    let input_variants: Option<HashMap<String, Vec<Variant>>> = match &config.input_vcf {
        Some(path) => {
            info!("Loading input VCF: {}", path.display());
            let raw = read_vcf(path.to_path_buf())?;
            Some(filter_input_vcf(raw))
        },
        None => None,
    };

    info!("Reading fasta file: {}", &config.reference.display());
    let bam_context: Option<Arc<BamContext>> = if config.produce_bam {
        let contigs = scan_fasta_lengths(&config.reference)?;
        Some(Arc::new(BamContext::new(&contigs)))
    } else {
        None
    };

    let ctx = ContigContext {
        config,
        target_bed: &target_bed,
        input_variants: &input_variants,
        nuc_sub_model: &nuc_sub_model,
        mutation_model: &mutation_model,
        run_mutation_rate,
        fragment_length_model: &fragment_length_model,
        gc_bias_model: &gc_bias_model,
        quality_score_model: &quality_score_model,
        seq_error_model: &seq_error_model,
        working_dir: working_dir.path(),
        base_rng: *rng,
        bam_context,
    };

    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    let mut all_fastq_files: HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)> = HashMap::new();
    let mut contig_order: Vec<String> = Vec::new();
    let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
    let mut bam_body_files: HashMap<String, PathBuf> = HashMap::new();

    info!("Generating simulated dataset");

    // All contigs processed in parallel; BAM workers each write a temp body file.
    let fasta = FastaStream::open(&config.reference)?;
    let parallel_iter = fasta
        .enumerate()
        .par_bridge()
        .map(|(idx, result)| -> Result<ContigResult, GenerateReadsError> {
            let (name, seq) = result?;
            let child_rng = ctx.base_rng.derive_child(idx as u64);
            process_contig(idx, name, seq, &ctx, child_rng)
        });
    let collected: Result<Vec<ContigResult>, _> = match config.num_threads {
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|e| GenerateReadsError::CliError(format!("Failed to build thread pool: {}", e)))?
            .install(|| parallel_iter.collect()),
        None => parallel_iter.collect(),
    };
    // par_bridge does not preserve order — restore contig order by idx.
    let mut results = collected?;
    results.sort_unstable_by_key(|r| r.idx);
    for cr in results {
        collect_contig_result(cr, &mut contig_order, &mut fasta_lengths, &mut mutated_maps, &mut all_fastq_files, &mut bam_body_files);
    }

    info!("Read generation complete, producing output files");

    if config.produce_fastq {
        info!("Producing final fastq(s) file(s)");

        let total_genome_bp: usize = fasta_lengths.values().sum();
        if total_genome_bp > 500_000_000 {
            warn!(
                "Genome is {:.1} Gbp. The global FASTQ shuffle loads all reads into memory. \
                 For large genomes, consider running `seqkit shuffle` on the output instead.",
                total_genome_bp as f64 / 1e9
            );
        }

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

    if config.produce_bam {
        if let (Some(bam_ctx), Some(bam_path)) = (ctx.bam_context.as_ref(), &config.output_bam) {
            info!("Assembling BAM from {} temp body file(s)", bam_body_files.len());
            let ordered_bodies: Vec<PathBuf> = contig_order.iter()
                .filter_map(|name| bam_body_files.remove(name))
                .collect();
            concat_temp_bams(bam_ctx, &ordered_bodies, bam_path)?;
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

fn process_contig(
    idx: usize,
    contig_name: String,
    mut sequence: Vec<Nucleotide>,
    ctx: &ContigContext,
    mut rng: NeatRng,
) -> Result<ContigResult, GenerateReadsError> {
    let contig_len = sequence.len();
    debug!("Processing {}", contig_name);

    if let Some(bed) = ctx.target_bed {
        if !bed.contains_key(&contig_name) {
            debug!("Skipping {} — not in target BED", contig_name);
            return Ok(ContigResult { idx, name: contig_name, len: contig_len, data: None });
        }
    }

    if sequence.is_empty() {
        warn!("Contig {} has empty sequence, skipping", contig_name);
        return Ok(ContigResult { idx, name: contig_name, len: contig_len, data: None });
    }

    apply_n_substitution(&mut sequence, ctx.nuc_sub_model, &mut rng)?;
    let sequence_map = map_buffer(&sequence);
    let current_block = SequenceBlock {
        contig: contig_name.clone(),
        ref_start: 0,
        ref_end: contig_len,
        sequence,
        sequence_map,
    };

    debug!("    > Generating bias map.");
    let raw_regions = current_block.get_non_n_regions();
    let regions_of_interest: Vec<SequenceMap> = if let Some(bed) = ctx.target_bed {
        let contig_beds = bed.get(&contig_name).map(|v| v.as_slice()).unwrap_or(&[]);
        intersect_with_bed(&raw_regions, contig_beds, 0)
    } else {
        raw_regions.into_iter().cloned().collect()
    };
    if regions_of_interest.is_empty() {
        return Ok(ContigResult { idx, name: contig_name, len: contig_len, data: None });
    }

    let mut bias_map: Vec<f64> = vec![0.0; current_block.get_len()];
    let mut num_mutable_area = 0.0_f64;
    for region in &regions_of_interest {
        for j in region.start..region.end {
            bias_map[j] = 1.0;
            num_mutable_area += 1.0;
        }
    }

    let block_end = contig_len;
    let mut block_variants: Vec<Variant> = if let Some(iv) = ctx.input_variants {
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

    let num_mutations = (num_mutable_area * ctx.run_mutation_rate).trunc() as usize;
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
            ctx.mutation_model,
            num_mutations,
            ctx.config.ploidy,
            &mut rng,
        )?;
        if let Some(vec) = result {
            for variant in vec {
                if variant.variant_type == VariantType::Deletion {
                    if variant.reference.len() - 1 > max_del_len {
                        max_del_len = variant.reference.len() - 1;
                    }
                }
                block_variants.push(variant);
            }
        }
    }

    let block_fragments: Vec<(usize, usize)> = {
        let mut block_frags = Vec::new();
        for (frag_start, frag_end) in regions_of_interest.into_iter().map(|r| (r.start, r.end)) {
            let frags = if ctx.gc_bias_model.is_uniform() {
                generate_fragments(
                    frag_end - frag_start,
                    ctx.config.read_len,
                    max_del_len,
                    frag_start,
                    ctx.config.coverage,
                    ctx.config.paired_ended,
                    ctx.config.long_reads,
                    ctx.fragment_length_model,
                    &mut rng,
                )?
            } else {
                generate_weighted_fragments(
                    &current_block,
                    frag_start,
                    frag_end,
                    ctx.config.read_len,
                    max_del_len,
                    ctx.config.coverage,
                    ctx.gc_bias_model,
                    ctx.fragment_length_model,
                    ctx.config.gc_bias_normalize_coverage,
                    ctx.config.long_reads,
                    &mut rng,
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

    // Create a per-contig BAM body writer if BAM output is requested.
    let mut bam_body_writer: Option<BamBodyWriter> = if let Some(bam_ctx) = &ctx.bam_context {
        let bam_temp_path = PathBuf::from(ctx.working_dir).join(format!(
            "temp_bam_{:06}_{}.bam", idx, contig_name
        ));
        Some(BamBodyWriter::new(bam_temp_path, Arc::clone(bam_ctx))?)
    } else {
        None
    };

    let read_name_prefix = format!("RNEAT_generated_{}", current_block.contig);

    if ctx.config.produce_fastq {
        let mut file_to_write_1 = PathBuf::from(ctx.working_dir);
        file_to_write_1.push(format!(
            "temp_{}_{:010}_{:010}_r1_tmp.fastq.gz",
            contig_name,
            current_block.ref_start,
            current_block.ref_end,
        ));
        let file1 = append_to_file(&file_to_write_1)?;
        let writer1 = BufWriter::new(&file1);
        let mut buffer1 = GzEncoder::new(writer1, Compression::default());
        let bam_stager: Option<&mut dyn BamRecordStager> =
            bam_body_writer.as_mut().map(|w| w as &mut dyn BamRecordStager);
        if ctx.config.paired_ended {
            let mut file_to_write_2 = PathBuf::from(ctx.working_dir);
            file_to_write_2.push(format!(
                "temp_{}_{:010}_{:010}_r2_tmp.fastq.gz",
                contig_name,
                current_block.ref_start,
                current_block.ref_end,
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
                ctx.config.read_len,
                ctx.config.long_reads,
                &read_name_prefix,
                ctx.quality_score_model,
                ctx.seq_error_model,
                &mut rng,
                bam_stager,
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
                ctx.config.read_len,
                ctx.config.long_reads,
                &read_name_prefix,
                ctx.quality_score_model,
                ctx.seq_error_model,
                &mut rng,
                bam_stager,
            )?;
            contig_files_r1.push(file_to_write_1);
        }
    } else if ctx.config.produce_bam {
        // BAM-only: generate reads and stage them into the BAM body writer.
        // The FASTQ buffers drain into null sinks and are discarded.
        let bam_stager: Option<&mut dyn BamRecordStager> =
            bam_body_writer.as_mut().map(|w| w as &mut dyn BamRecordStager);
        let null1: VectorBuffer = VectorBuffer::new();
        let null2: VectorBuffer = VectorBuffer::new();
        let mut buf1 = GzEncoder::new(null1, Compression::default());
        let mut buf2 = GzEncoder::new(null2, Compression::default());
        debug!("BAM-only: generating reads for {}", contig_name);
        write_block_fastq(
            block_fragments,
            &mutated_map,
            &current_block,
            ctx.config.paired_ended,
            &mut buf1,
            &mut buf2,
            ctx.config.read_len,
            ctx.config.long_reads,
            &read_name_prefix,
            ctx.quality_score_model,
            ctx.seq_error_model,
            &mut rng,
            bam_stager,
        )?;
    }

    // Flush and finalize the BAM body file; the bgzf EOF is written on drop.
    let bam_body_file = if let Some(mut bw) = bam_body_writer {
        bw.flush_all()?;
        Some(bw.path.clone())
    } else {
        None
    };

    Ok(ContigResult {
        idx,
        name: contig_name,
        len: contig_len,
        data: Some(ProcessedContigData {
            mutated_map,
            r1_files: contig_files_r1,
            r2_files: contig_files_r2,
            bam_body_file,
        }),
    })
}

fn collect_contig_result(
    cr: ContigResult,
    contig_order: &mut Vec<String>,
    fasta_lengths: &mut HashMap<String, usize>,
    mutated_maps: &mut HashMap<String, Vec<MutatedMap>>,
    all_fastq_files: &mut HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)>,
    bam_body_files: &mut HashMap<String, PathBuf>,
) {
    contig_order.push(cr.name.clone());
    fasta_lengths.insert(cr.name.clone(), cr.len);
    if let Some(data) = cr.data {
        mutated_maps.insert(cr.name.clone(), vec![data.mutated_map]);
        all_fastq_files.insert(cr.name.clone(), (data.r1_files, data.r2_files));
        if let Some(bam_path) = data.bam_body_file {
            bam_body_files.insert(cr.name, bam_path);
        }
    }
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