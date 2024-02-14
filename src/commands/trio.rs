use crate::{
    aligner::{AlignmentScope, MemoryModel, WFAligner, WFAlignerGapAffine2Pieces},
    allele,
    cli::TrioArgs,
    handles, locus,
    util::{self, AlnScoring, Params},
};
use anyhow::{anyhow, Result};
use csv::WriterBuilder;
use log;
use noodles::{
    bed,
    fasta::{self, io::BufReadSeek},
};
use once_cell::sync::OnceCell;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::{
    cell::RefCell,
    fs,
    io::BufReader,
    sync::{mpsc::channel, Arc},
    thread,
    time::Instant,
};

fn create_aligner_with_scoring() -> WFAligner {
    let aln_scoring = ALN_SCORING.get().expect("AlnScoring not initialized");
    WFAlignerGapAffine2Pieces::create_aligner(
        aln_scoring.mismatch,
        aln_scoring.gap_opening1,
        aln_scoring.gap_extension1,
        aln_scoring.gap_opening2,
        aln_scoring.gap_extension2,
        AlignmentScope::Alignment,
        MemoryModel::MemoryLow,
    )
}

static ALN_SCORING: OnceCell<AlnScoring> = OnceCell::new();

thread_local! {
    static ALIGNER: RefCell<WFAligner> = RefCell::new(create_aligner_with_scoring());
}

pub fn trio(args: TrioArgs) -> Result<()> {
    log::info!(
        "{}-{} trio start",
        env!("CARGO_PKG_NAME"),
        *crate::cli::FULL_VERSION
    );
    let start_timer = Instant::now();

    ALN_SCORING
        .set(args.aln_scoring)
        .map_err(|_| anyhow!("AlnScoring was already set"))?;

    let clip_len = if args.no_clip_aln { 0 } else { args.flank_len };

    let params = Params {
        clip_len,
        parent_quantile: args.parent_quantile,
        partition_by_alignment: args.partition_by_alignment,
    };
    let params_arc = Arc::new(params);

    let handles =
        handles::Handles::new(&args.mother_prefix, &args.father_prefix, &args.child_prefix)
            .unwrap_or_else(|err| util::handle_error_and_exit(err));
    let handles_arc = Arc::new(handles);

    let mut catalog_reader = fs::File::open(args.bed_filename)
        .map(BufReader::new)
        .map(bed::Reader::new)?;

    let mut genome_reader: fasta::IndexedReader<Box<dyn BufReadSeek>> =
        fasta::indexed_reader::Builder::default()
            .build_from_path(&args.reference_filename)
            .unwrap_or_else(|err| util::handle_error_and_exit(err.into()));

    match args.trid {
        Some(trid) => {
            let locus = locus::get_locus(
                &mut genome_reader,
                &mut catalog_reader,
                &trid,
                args.flank_len,
            )
            .unwrap_or_else(|err| util::handle_error_and_exit(err));

            let mut csv_wtr = WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());

            ALIGNER.with(|aligner| {
                let mut aligner = aligner.borrow_mut();
                if let Ok(result) =
                    allele::process_alleles(&locus, handles_arc, &params_arc, &mut aligner)
                {
                    for row in result {
                        if let Err(err) = csv_wtr.serialize(row) {
                            log::error!("Failed to write record: {}", err);
                        }
                        csv_wtr.flush().unwrap();
                    }
                }
            });
        }
        None => {
            let all_loci = locus::get_loci(&mut genome_reader, &mut catalog_reader, args.flank_len)
                .collect::<Result<Vec<_>>>()
                .unwrap_or_else(|err| util::handle_error_and_exit(err));

            let mut csv_wtr = WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&args.output_path)?;

            let (sender, receiver) = channel();
            let writer_thread = thread::spawn(move || {
                for results in &receiver {
                    for row in results {
                        if let Err(err) = csv_wtr.serialize(row) {
                            log::error!("Failed to write record: {}", err);
                        }
                        csv_wtr.flush().unwrap();
                    }
                }
            });

            log::info!("Starting job pool with {} thread(s)...", args.num_threads);
            let pool = ThreadPoolBuilder::new()
                .num_threads(args.num_threads)
                .build()
                .unwrap();

            pool.install(|| {
                all_loci
                    .into_iter()
                    .par_bridge()
                    .for_each_with(sender, |s, locus| {
                        ALIGNER.with(|aligner| {
                            let mut aligner = aligner.borrow_mut();
                            if let Ok(result) = allele::process_alleles(
                                &locus,
                                handles_arc.clone(),
                                &params_arc.clone(),
                                &mut aligner,
                            ) {
                                s.send(result).unwrap();
                            }
                        });
                    });
            });
            writer_thread.join().unwrap();
        }
    }
    log::info!("Total execution time: {:?}", start_timer.elapsed());
    log::info!("{} trio end", env!("CARGO_PKG_NAME"));

    Ok(())
}
