use crate::{
    aligner::{AlignmentScope, MemoryModel, WFAligner, WFAlignerGapAffine2Pieces},
    cli::DuoArgs,
    duo::allele::{process_alleles, AlleleResult},
    handles::DuoLocalData,
    locus::{self, Locus},
    util::{AlnScoring, Params},
};
use anyhow::{anyhow, Result};
use crossbeam_channel::{self, bounded, unbounded, Receiver, Sender};
use csv::{Writer, WriterBuilder};
use log;
use noodles::{bed, fasta};
use once_cell::sync::OnceCell;
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};
use std::{
    cell::RefCell,
    fs,
    io::{BufReader, Write},
    sync::Arc,
    thread,
};

static ALN_SCORING: OnceCell<AlnScoring> = OnceCell::new();

fn create_aligner_with_scoring() -> WFAligner {
    let aln_scoring = ALN_SCORING.get().expect("AlnScoring not initialized");
    WFAlignerGapAffine2Pieces::create_aligner(
        aln_scoring.mismatch,
        aln_scoring.gap_opening1,
        aln_scoring.gap_extension1,
        aln_scoring.gap_opening2,
        aln_scoring.gap_extension2,
        AlignmentScope::Alignment,
        MemoryModel::MemoryLow, // TODO: change to MemoryHigh?
    )
}

thread_local! {
    static ALIGNER: RefCell<WFAligner> = RefCell::new(create_aligner_with_scoring());
    static LOCAL_DUO_DATA: RefCell<Option<DuoLocalData>> = const { RefCell::new(None) };
}

pub fn duo(args: DuoArgs) -> Result<()> {
    ALN_SCORING
        .set(args.aln_scoring)
        .map_err(|_| anyhow!("AlnScoring was already set"))?;

    let clip_len = if args.no_clip_aln { 0 } else { args.flank_len };

    let params_arc = Arc::new(Params {
        clip_len,
        parent_quantile: args.p_quantile,
        partition_by_alignment: args.partition_by_alignment,
        quick_mode: args.quick,
    });

    let mut catalog_reader = fs::File::open(args.bed_filename.clone())
        .map(BufReader::new)
        .map(bed::Reader::new)?;

    let mut genome_reader =
        fasta::indexed_reader::Builder::default().build_from_path(&args.reference_filename)?;

    // Check if BAM/VCF files can be opened (index validity), better to do it here before spawning threads
    DuoLocalData::new(&args.a_prefix, &args.b_prefix)?;

    let a_name = std::path::Path::new(&args.a_prefix)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&args.a_prefix);
    let readids_filename = format!("{}_duo_denovo_reads.txt", a_name);
    let readids_writer = std::fs::File::create(&readids_filename)
        .map_err(|e| anyhow!("Failed to create read IDs file {}: {}", readids_filename, e))?;

    match args.trid {
        Some(ref trid) => {
            let locus = locus::get_locus(
                &mut genome_reader,
                &mut catalog_reader,
                trid,
                args.flank_len,
            )?;

            let tsv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());

            let (sender_result, receiver_result) = unbounded();
            let writer_thread = process_writer_thread(tsv_writer, readids_writer, receiver_result);

            process_locus(&locus, &args, &params_arc, &sender_result);

            drop(sender_result);
            writer_thread.join().unwrap();
        }

        None => {
            let bed_filename = args.bed_filename.clone();
            let reference_filename = args.reference_filename.clone();
            let flank_len = args.flank_len;
            let (sender_locus, receiver_locus) = bounded(2048);
            let locus_stream_thread = thread::spawn(move || {
                locus::stream_loci_into_channel(
                    bed_filename,
                    reference_filename,
                    flank_len,
                    sender_locus,
                );
            });

            let tsv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&args.output_path)?;
            let (sender_result, receiver_result) = unbounded();
            let writer_thread = process_writer_thread(tsv_writer, readids_writer, receiver_result);

            if args.num_threads == 1 {
                log::debug!("Single-threaded mode");
                for locus in receiver_locus {
                    match locus {
                        Ok(locus) => process_locus(&locus, &args, &params_arc, &sender_result),
                        Err(err) => log::error!("Locus Processing: {:#}", err),
                    }
                }
            } else {
                log::debug!(
                    "Multi-threaded mode: estimated available cores: {}",
                    thread::available_parallelism().unwrap().get()
                );
                let pool = initialize_thread_pool(args.num_threads)?;
                pool.install(|| {
                    receiver_locus.into_iter().par_bridge().for_each_with(
                        &sender_result,
                        |s, result| match result {
                            Ok(locus) => process_locus(&locus, &args, &params_arc, s),
                            Err(err) => log::error!("Locus Processing: {:#}", err),
                        },
                    );
                });
            }
            drop(sender_result);
            writer_thread.join().unwrap();
            locus_stream_thread.join().unwrap();
        }
    }
    Ok(())
}

fn process_writer_thread<T: Write + Send + 'static>(
    mut tsv_writer: Writer<T>,
    mut readids_writer: std::fs::File,
    receiver: Receiver<Vec<AlleleResult>>,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        for results in &receiver {
            for (i, row) in results.iter().enumerate() {
                if let Some(read_ids) = &row.read_ids {
                    if !read_ids.is_empty() {
                        writeln!(
                            readids_writer,
                            ">TRID={}\tALLELE={}\tN={}\n{}",
                            row.trid,
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

fn process_locus(
    locus: &Locus,
    args: &DuoArgs,
    params_arc: &Arc<Params>,
    sender_result: &Sender<Vec<AlleleResult>>,
) {
    ALIGNER.with(|aligner| {
        LOCAL_DUO_DATA.with(|local_family_data| {
            let mut duo_data = local_family_data.borrow_mut();
            if duo_data.is_none() {
                *duo_data = Some(
                    DuoLocalData::new(&args.a_prefix, &args.b_prefix)
                        .expect("Failed to initialize DuoLocalData"),
                );
            }
            if let Ok(result) = process_alleles(
                locus,
                duo_data.as_mut().unwrap(),
                params_arc,
                &mut aligner.borrow_mut(),
            ) {
                sender_result.send(result).unwrap();
            }
        });
    });
}

fn initialize_thread_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    log::info!("Starting job pool with {} thread(s)...", num_threads);
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| anyhow!("Failed to initialize thread pool: {}", e))
}
