use crate::{
    aligner::{AlignmentScope, MemoryModel, WFAligner, WFAlignerGapAffine2Pieces},
    locus::{self, Locus},
    model::{AlnScoring, Params, QuickMode},
    readers::{open_catalog_reader, open_genome_reader},
};
use anyhow::{anyhow, Result};
use crossbeam_channel::{self, bounded, unbounded, Receiver, Sender};
use csv::{Writer, WriterBuilder};
use log;
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};
use serde::Serialize;
use std::{io::Write, path::Path, sync::Arc, sync::LazyLock, thread};

pub trait Args {
    fn no_clip_aln(&self) -> bool;
    fn flank_len(&self) -> usize;
    fn p_quantile(&self) -> f64;
    fn partition_by_alignment(&self) -> bool;
    fn quick(&self) -> &Option<QuickMode>;
    fn bed_filename(&self) -> &Path;
    fn reference_filename(&self) -> &Path;
    fn trid(&self) -> &Option<String>;
    fn output_path(&self) -> &str;
    fn num_threads(&self) -> usize;
    fn readids_prefix(&self) -> &str;
    fn mode_name(&self) -> &str;
    fn preflight_check(&self) -> Result<()>;
}

pub trait AlleleResultExt: Serialize + Send + 'static + Clone {
    fn read_ids(&self) -> &Option<Vec<String>>;
    fn trid(&self) -> &str;
}

static ALN_SCORING: LazyLock<AlnScoring> = LazyLock::new(AlnScoring::default);

pub fn create_aligner_with_scoring() -> WFAligner {
    WFAlignerGapAffine2Pieces::create_aligner(
        ALN_SCORING.mismatch,
        ALN_SCORING.gap_opening1,
        ALN_SCORING.gap_extension1,
        ALN_SCORING.gap_opening2,
        ALN_SCORING.gap_extension2,
        AlignmentScope::Alignment,
        MemoryModel::MemoryLow, // TODO: change to MemoryHigh?
    )
}

pub fn run<A, R, F>(args: A, process_locus_fn: F) -> Result<()>
where
    A: Args + Sync,
    R: AlleleResultExt,
    F: Fn(&Locus, &A, &Arc<Params>, &Sender<Vec<R>>) + Sync + Send,
{
    let clip_len = if args.no_clip_aln() {
        0
    } else {
        args.flank_len()
    };

    let params_arc = Arc::new(Params {
        clip_len,
        parent_quantile: args.p_quantile(),
        partition_by_alignment: args.partition_by_alignment(),
        quick_mode: *args.quick(),
    });

    let mut catalog_reader = open_catalog_reader(args.bed_filename())?;
    let genome_reader = open_genome_reader(args.reference_filename())?;

    // Check if BAM/VCF files can be opened (index validity), better to do it here before spawning threads
    args.preflight_check()?;

    let name = std::path::Path::new(args.readids_prefix())
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_else(|| args.readids_prefix());
    let readids_filename = format!("{}_{}_denovo_reads.txt", name, args.mode_name());
    let readids_writer = std::fs::File::create(&readids_filename)
        .map_err(|e| anyhow!("Failed to create read IDs file {}: {}", readids_filename, e))?;

    match args.trid() {
        Some(ref trid) => {
            let locus =
                locus::get_locus(&genome_reader, &mut catalog_reader, trid, args.flank_len())?;

            let tsv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());

            let (sender_result, receiver_result) = unbounded();
            let writer_thread = process_writer_thread(tsv_writer, readids_writer, receiver_result);

            process_locus_fn(&locus, &args, &params_arc, &sender_result);

            drop(sender_result);
            writer_thread.join().unwrap();
        }
        None => {
            let bed_filename = args.bed_filename().to_path_buf();
            let reference_filename = args.reference_filename().to_path_buf();
            let flank_len = args.flank_len();
            let (sender_locus, receiver_locus) = bounded(2048);
            let locus_stream_thread = thread::spawn(move || {
                locus::stream_loci_into_channel(
                    bed_filename,
                    reference_filename,
                    flank_len,
                    sender_locus,
                )
            });

            let tsv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(args.output_path())?;
            let (sender_result, receiver_result) = unbounded();
            let writer_thread = process_writer_thread(tsv_writer, readids_writer, receiver_result);

            if args.num_threads() == 1 {
                log::debug!("Single-threaded mode");
                for locus in receiver_locus {
                    match locus {
                        Ok(locus) => process_locus_fn(&locus, &args, &params_arc, &sender_result),
                        Err(err) => log::error!("Locus Processing: {:#}", err),
                    }
                }
            } else {
                log::debug!(
                    "Multi-threaded mode: estimated available cores: {}",
                    thread::available_parallelism().unwrap().get()
                );
                let pool = initialize_thread_pool(args.num_threads())?;
                pool.install(|| {
                    receiver_locus.into_iter().par_bridge().for_each_with(
                        &sender_result,
                        |s, result| match result {
                            Ok(locus) => process_locus_fn(&locus, &args, &params_arc, s),
                            Err(err) => log::error!("Locus Processing: {:#}", err),
                        },
                    );
                });
            }
            drop(sender_result);
            writer_thread.join().unwrap();
            match locus_stream_thread
                .join()
                .expect("Locus stream thread panicked")
            {
                Ok(_) => log::trace!("Locus stream thread finished"),
                Err(e) => log::error!("Locus streaming failed: {}", e),
            }
        }
    }
    Ok(())
}

fn process_writer_thread<T: Write + Send + 'static, R: AlleleResultExt>(
    mut tsv_writer: Writer<T>,
    mut readids_writer: std::fs::File,
    receiver: Receiver<Vec<R>>,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        for results in &receiver {
            for (i, row) in results.iter().enumerate() {
                if let Some(read_ids) = row.read_ids() {
                    if !read_ids.is_empty() {
                        writeln!(
                            readids_writer,
                            ">TRID={}\tALLELE={}\tN={}\n{}",
                            row.trid(),
                            i,
                            read_ids.len(),
                            read_ids.join("\n")
                        )
                        .unwrap_or_else(|e| {
                            log::error!("Failed to write read IDs: {}", e);
                        });
                    }
                }

                if let Err(err) = tsv_writer.serialize(row) {
                    log::error!("Failed to write record: {}", err);
                }
                tsv_writer.flush().unwrap();
            }
        }
        if receiver.recv().is_err() {
            log::debug!("All data processed, exiting writer thread.");
        }
    })
}

fn initialize_thread_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    log::info!("Starting job pool with {} thread(s)...", num_threads);
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| anyhow!("Failed to initialize thread pool: {}", e))
}
