use rayon::prelude::*;
use common::file_tools::file_io::create_output_file;
use common::structs::variants::VariantType;
use tempfile;
use log::{info, debug, warn, error};
use common::rng::NeatRng;
use std::collections::{HashMap, HashSet};
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
            fasta_stream::{FastaStream, map_buffer, apply_n_substitution, scan_fasta_lengths, resolve_iupac_bases},
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
    mutation_regions: &'a Option<HashMap<String, Vec<BedRecord>>>,
    input_variants: &'a Option<HashMap<String, Vec<Variant>>>,
    nuc_sub_model: &'a NucleotideSelector,
    mutation_model: &'a MutationModel,
    default_run_mutation_rate: f64,
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
    let mutation_regions = match &config.mutation_regions {
        Some(path) => {
            info!("Loading mutation regions BED: {:?}", path);
            Some(read_bed(path, true)?)
        },
        None => None,
    };
    let default_run_mutation_rate = match config.mutation_rate {
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
            // Fall through to the quality model embedded in the sequencing-error model
            // when the user hasn't supplied an explicit override. This matches the
            // documented contract of `sequence_error_model:` and ensures binned-quality
            // training (gen-seq-error-model with binned_quality_bins) actually drives
            // gen-reads sampling — otherwise the binned QSM sits unused inside the
            // SeqErrorModel while the default continuous model is used here.
            None => seq_error_model.quality_score_model().clone(),
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
            info!("Loading target BED: {:?}", path);
            Some(read_bed(path, false)?)
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
        mutation_regions: &mutation_regions,
        input_variants: &input_variants,
        nuc_sub_model: &nuc_sub_model,
        mutation_model: &mutation_model,
        default_run_mutation_rate,
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
            let (name, raw) = result?;
            let mut child_rng = ctx.base_rng.derive_child(idx as u64);
            let (seq, iupac_count) = resolve_iupac_bases(&raw, &mut child_rng)?;
            if iupac_count > 0 {
                warn!("Contig {}: resolved {} IUPAC ambiguity base(s) to ACGT", name, iupac_count);
            }
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
    let regions_of_interest: Vec<SequenceMap> =
        if let Some(bed) = ctx.target_bed {
            let contig_beds = bed
                .get(&contig_name)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            intersect_with_bed(&raw_regions, contig_beds, 0)
        } else {
            raw_regions.into_iter().cloned().collect()
        };
    if regions_of_interest.is_empty() {
        return Ok(ContigResult { idx, name: contig_name, len: contig_len, data: None });
    }

    // Build a compact segment list instead of a per-position Vec<f64>.
    // Each segment is (start, end, rate); N-regions and gaps are simply absent.
    // This replaces an O(chromosome_length) allocation with O(regions + BED_records).
    let mut rate_segments: Vec<(usize, usize, f64)> = regions_of_interest
        .iter()
        .map(|r| (r.start, r.end, ctx.default_run_mutation_rate))
        .collect();

    if let Some(mut_beds) = ctx.mutation_regions {
        if let Some(records) = mut_beds.get(&contig_name) {
            for rec in records {
                if let Some(custom_rate) = rec.mut_rate {
                    rate_segments = apply_rate_override(
                        rate_segments, rec.start, rec.end, custom_rate,
                    );
                }
            }
        }
    }

    let mut num_mutations_sum: f64 = rate_segments
        .iter()
        .map(|&(s, e, r)| (e - s) as f64 * r)
        .sum();

    let block_end = contig_len;
    let mut block_variants: Vec<Variant> = Vec::new();
    if let Some(iv) = ctx.input_variants {
        if let Some(vs) = iv.get(&contig_name) {
            let mut excluded: Vec<usize> = Vec::new();
            let mut seen: HashSet<usize> = HashSet::new();
            for v in vs {
                let pos0 = v.location.saturating_sub(1);
                if pos0 >= block_end {
                    continue;
                }
                let local_pos = v.location - 1;
                if seen.insert(local_pos) {
                    let rate = rate_at(&rate_segments, local_pos);
                    if rate > 0.0 {
                        num_mutations_sum -= rate;
                        excluded.push(local_pos);
                    }
                }
                let mut v2 = v.clone();
                v2.location = local_pos;
                block_variants.push(v2);
            }
            if !excluded.is_empty() {
                excluded.sort_unstable();
                rate_segments = exclude_positions(rate_segments, &excluded);
            }
        }
    }
    debug!("Seeded {} user variant(s) into contig {}", block_variants.len(), contig_name);
    let num_mutations = num_mutations_sum.trunc() as usize;
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
            &rate_segments,
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
                    ctx.config.paired_ended,
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

/// Overrides the mutation rate in [ovr_start, ovr_end) within existing segments,
/// splitting segment boundaries where needed. Positions outside existing segments
/// (N-regions, gaps) are not affected.
fn apply_rate_override(
    segments: Vec<(usize, usize, f64)>,
    ovr_start: usize,
    ovr_end: usize,
    ovr_rate: f64,
) -> Vec<(usize, usize, f64)> {
    let mut result = Vec::with_capacity(segments.len() + 2);
    for (s, e, rate) in segments {
        if ovr_end <= s || ovr_start >= e {
            result.push((s, e, rate));
            continue;
        }
        let isect_s = s.max(ovr_start);
        let isect_e = e.min(ovr_end);
        if s < isect_s { result.push((s, isect_s, rate)); }
        result.push((isect_s, isect_e, ovr_rate));
        if isect_e < e { result.push((isect_e, e, rate)); }
    }
    result
}

/// Returns the mutation rate at `pos`, or 0.0 if the position falls in an N-region
/// or gap. Segments must be sorted by start and non-overlapping.
fn rate_at(segments: &[(usize, usize, f64)], pos: usize) -> f64 {
    let idx = segments.partition_point(|&(s, _, _)| s <= pos);
    if idx == 0 { return 0.0; }
    let (_, e, rate) = segments[idx - 1];
    if pos < e { rate } else { 0.0 }
}

/// Splits segments to remove individual excluded positions (e.g. positions already
/// occupied by input variants). `excluded` must be sorted.
fn exclude_positions(
    segments: Vec<(usize, usize, f64)>,
    excluded: &[usize],
) -> Vec<(usize, usize, f64)> {
    if excluded.is_empty() { return segments; }
    let mut result = Vec::with_capacity(segments.len() + excluded.len());
    let mut ei = 0;
    for (s, e, rate) in segments {
        let mut cur = s;
        while ei < excluded.len() && excluded[ei] < e {
            let pos = excluded[ei];
            ei += 1;
            if pos < cur { continue; }
            if cur < pos { result.push((cur, pos, rate)); }
            cur = pos + 1;
        }
        if cur < e { result.push((cur, e, rate)); }
    }
    result
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

#[cfg(test)]
mod tests {
    use super::*;
    use common::structs::sequence_block::{RegionType, SequenceMap};
    use common::structs::bed_record::BedRecord;

    #[test]
    fn test_intersect_with_bed() {
        let r1 = SequenceMap::from(RegionType::NonNRegion, 100, 200);
        let r2 = SequenceMap::from(RegionType::NonNRegion, 300, 400);
        let regions = vec![&r1, &r2];

        let b1 = BedRecord::new_bed_record("chr1".to_string(), 150, 350).unwrap();
        let bed_records = vec![b1];

        // block_offset = 0
        let result = intersect_with_bed(&regions, &bed_records, 0);
        assert_eq!(result.len(), 2);
        // Intersection with r1 [100, 200] and b1 [150, 350] -> [150, 200]
        assert_eq!(result[0].start, 150);
        assert_eq!(result[0].end, 200);
        // Intersection with r2 [300, 400] and b1 [150, 350] -> [300, 350]
        assert_eq!(result[1].start, 300);
        assert_eq!(result[1].end, 350);

        // block_offset = 1000
        // r1 global [1100, 1200], r2 global [1300, 1400]
        let b2 = BedRecord::new_bed_record("chr1".to_string(), 1150, 1350).unwrap();
        let result2 = intersect_with_bed(&regions, &vec![b2], 1000);
        assert_eq!(result2.len(), 2);
        assert_eq!(result2[0].start, 150); // global 1150 - 1000
        assert_eq!(result2[0].end, 200);   // global 1200 - 1000
        assert_eq!(result2[1].start, 300); // global 1300 - 1000
        assert_eq!(result2[1].end, 350);   // global 1350 - 1000
    }

    #[test]
    fn test_filter_input_vcf() {
        use crate::common::structs::variants::{VariantType, Genotype};
        use common::structs::nucleotides::Nucleotide;
        let mut raw = HashMap::new();
        let v1 = Variant {
            location: 100,
            reference: vec![Nucleotide::A],
            alternate: vec![Nucleotide::T],
            variant_type: VariantType::SNP,
            genotype: Genotype::Homozygous,
            genotype_str: "1/1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
        };
        let v2 = Variant {
            location: 200,
            reference: vec![Nucleotide::A, Nucleotide::T],
            alternate: vec![Nucleotide::C, Nucleotide::G],
            variant_type: VariantType::Complex,
            genotype: Genotype::Homozygous,
            genotype_str: "1/1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
        };
        raw.insert("chr1".to_string(), vec![v1.clone(), v2]);

        let filtered = filter_input_vcf(raw);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered["chr1"].len(), 1);
        assert_eq!(filtered["chr1"][0].location, 100);
    }

    #[test]
    fn test_apply_rate_override() {
        let segs = vec![(0usize, 100usize, 0.001f64), (200, 400, 0.001)];

        // No overlap: override entirely before first segment
        let result = apply_rate_override(segs.clone(), 0, 0, 0.01);
        assert_eq!(result, segs);

        // No overlap: override entirely after last segment
        let result = apply_rate_override(segs.clone(), 500, 600, 0.01);
        assert_eq!(result, segs);

        // Partial overlap at start of first segment: [0,50) gets new rate, [50,100) keeps old
        let result = apply_rate_override(segs.clone(), 0, 50, 0.01);
        assert_eq!(result, vec![(0, 50, 0.01), (50, 100, 0.001), (200, 400, 0.001)]);

        // Partial overlap at end of first segment: [0,80) keeps old, [80,100) gets new rate
        let result = apply_rate_override(segs.clone(), 80, 150, 0.01);
        assert_eq!(result, vec![(0, 80, 0.001), (80, 100, 0.01), (200, 400, 0.001)]);

        // Full containment of first segment: entire [0,100) replaced
        let result = apply_rate_override(segs.clone(), 0, 100, 0.02);
        assert_eq!(result, vec![(0, 100, 0.02), (200, 400, 0.001)]);

        // Override spanning both segments (gap between them is unaffected)
        let result = apply_rate_override(segs.clone(), 50, 300, 0.05);
        assert_eq!(result, vec![(0, 50, 0.001), (50, 100, 0.05), (200, 300, 0.05), (300, 400, 0.001)]);
    }

    #[test]
    fn test_rate_at() {
        let segs = vec![(10usize, 50usize, 0.001f64), (100, 200, 0.005)];

        // Inside first segment
        assert_eq!(rate_at(&segs, 25), 0.001);
        // At start of first segment (inclusive)
        assert_eq!(rate_at(&segs, 10), 0.001);
        // At end of first segment (exclusive — gap)
        assert_eq!(rate_at(&segs, 50), 0.0);
        // In gap between segments
        assert_eq!(rate_at(&segs, 75), 0.0);
        // Before all segments
        assert_eq!(rate_at(&segs, 0), 0.0);
        // Inside second segment
        assert_eq!(rate_at(&segs, 150), 0.005);
        // At end of second segment (exclusive)
        assert_eq!(rate_at(&segs, 200), 0.0);
    }

    #[test]
    fn test_exclude_positions() {
        let segs = vec![(0usize, 100usize, 0.001f64), (200, 400, 0.001)];

        // Empty excluded list — segments unchanged
        assert_eq!(exclude_positions(segs.clone(), &[]), segs);

        // Exclude middle of first segment — splits into two
        let result = exclude_positions(segs.clone(), &[50]);
        assert_eq!(result, vec![(0, 50, 0.001), (51, 100, 0.001), (200, 400, 0.001)]);

        // Exclude start of segment — trims left boundary
        let result = exclude_positions(segs.clone(), &[0]);
        assert_eq!(result, vec![(1, 100, 0.001), (200, 400, 0.001)]);

        // Exclude last position of segment — trims right boundary
        let result = exclude_positions(segs.clone(), &[99]);
        assert_eq!(result, vec![(0, 99, 0.001), (200, 400, 0.001)]);

        // Exclude position in gap — no change to segments
        let result = exclude_positions(segs.clone(), &[150]);
        assert_eq!(result, segs);

        // Exclude multiple positions across both segments
        let result = exclude_positions(segs.clone(), &[50, 250]);
        assert_eq!(result, vec![
            (0, 50, 0.001), (51, 100, 0.001),
            (200, 250, 0.001), (251, 400, 0.001),
        ]);
    }
}
