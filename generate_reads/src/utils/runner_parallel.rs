use crate::errors;
use crate::errors::GenerateReadsErrors;
use crate::utils;
use common;
use common::models::fragment_length::DiscreteFragmentLengthModel;
use common::models::fragment_length::FragmentLengthModel;
use common::models::fragment_length::NormalFragmentLengthModel;
use common::models::quality_scores::QualityScoreModel;
use tempfile;

use std::thread;
use std::time;

use log::info;
use std::collections::{HashMap, VecDeque};
use std::sync::{Arc, Mutex};

use utils::config::RunConfiguration;
use utils::generate_reads::generate_reads;
use utils::generate_variants::generate_variants;
use common::file_tools::fasta_tools::read_fasta;
use common::structs::variants::Variant;
use common::structs::fasta_map::FastaMap;
use common::models::mutation_model::MutationModel;
use simple_rng::NeatRng;
use utils::mutate_fasta::apply_mutations;

pub fn run_neat(config: &Box<RunConfiguration>, rng: NeatRng) -> Result<(), GenerateReadsErrors> {
    let start = time::Instant::now();
    info!("////////////// Welcome to rusty-neat read generator!");
    info!("Processing started!");
    // Create the prefix of the files to write
    let output_file_prefix = format!("{}/{}", config.output_dir.display(), config.output_prefix);
    // We'll need a temp dir to store file fragments
    let working_dir = tempfile::tempdir().unwrap();

    // Load models that will be used for the runs.
    //
    // Quality score model
    let quality_score_model = {
        match config.quality_score_data {
            Some(filename) => {
                 QualityScoreModel::from_file(&filename)?
            },
            None => {
                QualityScoreModel::default()?
            }
        }
    };

    // load mutation model
    let mutation_model = {
        match config.mutation_model_data {
            Some(filename) => {
                MutationModel::from_file(&filename)?
            },
            None => {
                MutationModel::default()?
            }
        }
    };

    // Fragment Length Model
    let fragment_length_model: FragmentLengthModel = {
        match config.fragment_model_data {
            Some(filename) => {
                DiscreteFragmentLengthModel::discrete_from_file(&filename)?.into()
            },
            None => {
                match config.fragment_mean {
                    Some(mean) => {
                        NormalFragmentLengthModel::new_from_mean(
                            mean,
                            // Config should already have caught issues with missing data
                            config.fragment_st_dev.unwrap()
                        )?.into()
                    },
                    None => {
                        NormalFragmentLengthModel::default()?.into()
                    },
                }
            }
        }
    };

    // Load Indel model
    let indel_model = None;

    // The FastaMap struct consists of a vector of Contig objects, each describing a
    // block of Sequence, broken up into chunks in the SequenceBlock objects. It also
    // holds a name_map hashmap that links tho contig name from the Fasta with the shortname
    let fasta_map = read_fasta(
        &config.reference,
        config.read_len,
        &working_dir,
    ).unwrap();

    // Mutating the reference and recording the variant locations.
    info!("Generating simulated dataset");

    // all variants will be a hashmap of the contig name indexing a hashmap of variants keyed by
    // location. Basically, the Box moves the map to the heap, the mutex is a way to lock the file
    // and arc handles the communication between threads. From the rust book.
    let all_variants_mutex: Arc<Mutex<Box<HashMap<String, HashMap<usize, Variant>>>>> =
        Arc::new(Mutex::new(Box::new(HashMap::new())));
    let all_reads_mutex: Arc<Mutex<Box<Vec<(String, usize, usize)>>>> =
        Arc::new(Mutex::new(Box::new(Vec::new())));
    let fasta_order_mutex: Arc<Mutex<VecDeque<String>>> = Arc::new(Mutex::new(fasta_order));
    let fasta_map_mutex: Arc<Mutex<Box<FastaMap>>> =
        Arc::new(Mutex::new(Box::new(fasta_map)));
    let mutation_model_mutex: Arc<Mutex<MutationModel>> = Arc::new(Mutex::new(mutation_model));
    let confix_mutex: Arc<Mutex<Box<RunConfiguration>>> = Arc::new(Mutex::new(config));
    let local_rng_mutex: Arc<Mutex<ChaCha20Rng>> = Arc::new(Mutex::new(rng.clone()));

    let mut threads = Vec::with_capacity(n);
    (0..n).for_each(|_| {
        let all_variants_mutex_clone = Arc::clone(&all_variants_mutex);
        let all_reads_mutex_clone = Arc::clone(&all_reads_mutex);
        let fasta_order_mutex_clone = Arc::clone(&fasta_order_mutex);
        let fasta_map_mutex_clone = Arc::clone(&fasta_map_mutex);
        let mutation_model_mutex_clone = Arc::clone(&mutation_model_mutex);
        let config_mutex_clone = Arc::clone(&confix_mutex);
        let local_rng_mutex_clone = Arc::clone(&local_rng_mutex);

        threads.push(thread::spawn(move || {
            let result: Result<(), RunNeatError> = {
                let contig = {
                    let mut fasta_order = fasta_order_mutex_clone
                        .lock()
                        .unwrap();
                    let contig = fasta_order.pop_front().unwrap().clone();
                    fasta_order.push_back(contig.clone());
                    contig
                }.clone();
                info!("Generating variants for {}", contig);
                let contig_variants = {
                    let config = config_mutex_clone.lock().unwrap();
                    let local_rng =
                        local_rng_mutex_clone
                            .lock()
                            .unwrap();
                    let fasta_map = fasta_map_mutex_clone
                        .lock()
                        .unwrap();
                    // todo update to fasta map syntax
                    let contig_map = fasta_map.map.get(&contig).unwrap();
                    let mutation_model =
                        mutation_model_mutex_clone.lock().unwrap();
                    generate_variants(
                        contig_map.clone(),
                        &fasta_sequences[&contig],
                        &mut mutation_model.clone(),
                        config.ploidy,
                        local_rng.clone()
                    )
                };
                let _ = {
                    let mut all_variants =
                        all_variants_mutex_clone.lock().unwrap();
                    all_variants.insert(contig.to_owned(), contig_variants)
                };

                // This all means if this run is fasta only, we can skip generating reads
                let config = config_mutex_clone.lock().unwrap();
                if config.produce_fastq ||
                    config.produce_vcf ||
                    config.produce_bam {
                    info!("Generating reads for the output for {}", &contig);
                    let contig_sequence_len = {
                        let fasta_map =
                            fasta_map_mutex_clone
                                .lock()
                                .unwrap();
                        fasta_map[&contig].len()
                    }.clone();
                    let local_rng =
                        local_rng_mutex_clone
                            .lock()
                            .unwrap();
                    let contig_reads = generate_reads(
                        contig_sequence_len,
                        config.read_len,
                        config.coverage,
                        config.paired_ended,
                        config.fragment_mean,
                        config.fragment_st_dev,
                        local_rng.clone(),
                    ).expect("Error generating reads");
                    // append contig name to all reads
                    let mut all_reads =
                        all_reads_mutex_clone.lock().unwrap();
                    for read in contig_reads{
                        all_reads.push((contig.clone(), read.0, read.1))
                    }
                };

                // For vcf only, we don't need to waste time mutating the contig.
                if config.produce_vcf &&
                    !(config.produce_fasta || config.produce_bam || config.produce_fastq) {
                    info!("Mutating the reference sequence of {}", &contig);
                    let mut sequence_to_mutate = {
                        let fasta_map =
                            fasta_map_mutex_clone
                                .lock()
                                .unwrap();
                        fasta_map[&contig].clone()
                    };

                    let contig_variants = {
                        let mut all_variants =
                            all_variants_mutex_clone.lock().unwrap();
                        all_variants.get(&contig).unwrap().clone()
                    };
                    sequence_to_mutate = apply_mutations(
                        &sequence_to_mutate,
                        &contig_variants,
                    ).expect("Error applying mutations");
                    // update our fasta map with the mutated sequence.
                    let _ = {
                        let mut fasta_map =
                            fasta_map_mutex_clone
                                .lock()
                                .unwrap();
                        // I believe doing this will overwrite the reference sequence, which is what
                        // we're trying to accomplish here.
                        fasta_map.insert(contig, sequence_to_mutate);
                    };

                }
                Ok(())
            };
            result.unwrap();
        }));
    });

    for handle in threads {
        handle.join().unwrap();
    }

    // I'm not sure if this is necessary, but going to grab the lock on these files for
    // the remainder.
    let all_reads_clone = Arc::clone(&all_reads_mutex);
    let mut all_reads = all_reads_clone.lock().unwrap();
    let fast_order_mutex_clone = Arc::clone(&fasta_order_mutex);
    let fasta_order = fast_order_mutex_clone.lock().unwrap();
    let fasta_map_mutex_clone = Arc::clone(&fasta_map_mutex);
    // at this point the sequences are mutated, so we rename the variable.
    let mutated_fasta_map = fasta_map_mutex_clone.lock().unwrap();
    let config_mutex_clone = Arc::clone(&confix_mutex);
    let config = config_mutex_clone.lock().unwrap();

    if config.produce_fasta {
        info!("Producing the fasta file.");
        write_fasta(
            &mutated_fasta_map,
            &fasta_order,
            config.overwrite_output,
            &output_file_prefix,
        ).expect("Error writing fasta file!")
    }

    if config.produce_fastq {
        write_fastq(
            &mutated_fasta_map,
            &mut all_reads,
            config.read_len,
            config.overwrite_output,
            &output_file_prefix,
            config.paired_ended,
            quality_score_model,
            rng.clone()
        ).expect("Error writing fastq file(s)!")
    }

    // if config.produce_vcf {
    //     write_vcf(
    //         &all_variants,
    //         &fasta_order,
    //         &fasta_lengths,
    //         &config.reference,
    //         config.overwrite_output,
    //         &output_file,
    //     ).expect("Error writing vcf file!")
    // }

    // let (mutated_map, variant_locations) =
    //     mutate_fasta(&fasta_map, config.minimum_mutations, &mut rng);
    //
    // if config.produce_fasta {
    //     info!("Outputting fasta file");
    //     write_fasta(
    //         &mutated_map,
    //         &fasta_order,
    //         config.overwrite_output,
    //         &output_file,
    //     )
    //     .unwrap();
    // }
    //
    // if config.produce_vcf {
    //     info!("Writing vcf file");
    //     write_vcf(
    //         &variant_locations,
    //         &fasta_order,
    //         config.ploidy,
    //         &config.reference,
    //         config.overwrite_output,
    //         &output_file,
    //         &mut rng,
    //     )
    //     .unwrap();
    // }
    //
    // let mut read_sets: HashSet<Vec<Nuc>> = HashSet::new();
    // info!("Generating reads");
    // for (_name, sequence) in mutated_map.iter() {
    //     // defined as a set of read sequences that should cover
    //     // the mutated sequence `coverage` number of times
    //     let data_set = generate_reads(
    //         sequence,
    //         &config.read_len,
    //         &config.coverage,
    //         config.paired_ended,
    //         config.fragment_mean,
    //         config.fragment_st_dev,
    //         &mut rng,
    //     )
    //     .unwrap();
    //
    //     read_sets.extend(*data_set);
    // }
    //
    // if config.produce_fastq {
    //     info!("Shuffling output fastq data");
    //     let mut outsets: Box<Vec<&Vec<Nuc>>> = Box::new(read_sets.iter().collect());
    //     outsets.shuffle(&mut rng);
    //
    //     info!("Writing fastq");
    //     write_fastq(
    //         &output_file,
    //         config.overwrite_output,
    //         config.paired_ended,
    //         *outsets,
    //         quality_score_model,
    //         rng,
    //     )
    //     .unwrap();
    //     info!("Processing complete")
    // }
    let elapsed_time = time::Instant::now() - start;
    info!("Processing finished in {} milliseconds", elapsed_time.as_millis());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::SeedableRng;
    use std::fs;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use utils::config::ConfigBuilder;
    use simplelog;
    use simplelog::*;

    #[test]
    fn test_runner() {
        let mut config = RunConfiguration::build();
        config.reference = Some("test_data/references/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test_runner");
        fs::create_dir("test_runner").unwrap();
        let config = config.build();
        let _ = run_neat(Box::new(config), ChaCha20Rng::seed_from_u64(0)).unwrap();
        fs::remove_dir_all("test_runner").unwrap();
    }

    #[test]
    fn test_runner_files_message() {
        let mut config = ConfigBuilder::new();
        config.reference = Some("test_data/references/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.produce_vcf = true;
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test_run_output");

        TermLogger::init(
            LevelFilter::Trace,
            Config::default(),
            TerminalMode::Stdout,
            ColorChoice::Auto,
        )
            .unwrap();

        fs::create_dir("test_run_output").unwrap();
        let config = Box::new(config.build());
        let rng = ChaCha20Rng::seed_from_u64(0);
        run_neat(config, rng.clone()).unwrap();
        let file_path = "test_run_output/neat_out.fasta";
        let input = File::open(file_path).unwrap();
        let buffered = BufReader::new(input);
        let line_count = buffered.lines().count();
        assert!(line_count > 0);
        fs::remove_dir_all("test_run_output").unwrap();
    }
}
