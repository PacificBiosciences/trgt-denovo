use crate::{
    alignment_scores::ScoreOutput,
    locus::{self, Locus},
    model::{AlnScoring, Params, QuickMode},
    readers::{open_catalog_reader, open_genome_reader},
};
use anyhow::{Result, anyhow};
use crossbeam_channel::{self, Receiver, Sender, bounded, unbounded};
use csv::{Writer, WriterBuilder};
use log;
use rayon::{
    ThreadPoolBuilder,
    iter::{ParallelBridge, ParallelIterator},
};
use rust_wfa2::aligner::{AlignmentScope, Heuristics, MemoryModel, WFAligner};
use serde::Serialize;
use std::{io::Write, path::Path, sync::Arc, thread};

pub trait Args {
    fn no_clip_aln(&self) -> bool;
    fn flank_len(&self) -> usize;
    fn p_quantile(&self) -> f64;
    fn partition_by_alignment(&self) -> bool;
    fn skip_tr_check(&self) -> bool;
    fn quick(&self) -> &Option<QuickMode>;
    fn bed_filename(&self) -> &Path;
    fn reference_filename(&self) -> &Path;
    fn trid(&self) -> &Option<String>;
    fn output_path(&self) -> Option<&str>;
    fn num_threads(&self) -> usize;
    fn readids_prefix(&self) -> &str;
    fn mode_name(&self) -> &str;
    fn preflight_check(&self) -> Result<()>;
}

pub trait AlleleResultExt: Serialize + Send + 'static + Clone {
    fn read_ids(&self) -> &Option<Vec<String>>;
    fn trid(&self) -> &str;
}

pub fn create_aligner_with_scoring(scoring: AlnScoring) -> Result<WFAligner> {
    WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine2p(
            scoring.mismatch,
            scoring.gap_opening1,
            scoring.gap_extension1,
            scoring.gap_opening2,
            scoring.gap_extension2,
        )
        // Preserve the pruning used by the previous wrapper's C defaults.
        .with_heuristics(Heuristics::wfa2_default())
        .build()
        .map_err(|e| anyhow!("Invalid two-piece affine alignment configuration: {e}"))
}

pub(crate) fn configured_aligner(
    cache: &mut Option<(AlnScoring, WFAligner)>,
    scoring: AlnScoring,
) -> Result<&mut WFAligner> {
    if cache
        .as_ref()
        .is_none_or(|(current, _)| *current != scoring)
    {
        *cache = Some((scoring, create_aligner_with_scoring(scoring)?));
    }
    Ok(&mut cache.as_mut().unwrap().1)
}

pub fn run<A, R, F, O>(args: A, process_locus_fn: F, score_output: O) -> Result<()>
where
    A: Args + Sync,
    R: AlleleResultExt,
    O: ScoreOutput,
    F: for<'a> Fn(&'a Locus, &A, &Arc<Params>, &Sender<Vec<R>>, &mut O::Recorder<'a>) -> Result<()>
        + Sync
        + Send,
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
        skip_tr_check: args.skip_tr_check(),
        quick_mode: *args.quick(),
    });

    let mut catalog_reader = open_catalog_reader(args.bed_filename())?;
    let genome_reader = open_genome_reader(args.reference_filename())?;

    // Check if BAM/VCF files can be opened (index validity), better to do it here before spawning threads
    args.preflight_check()?;

    let output_file = args.output_path().map(std::fs::File::create).transpose()?;

    let name = std::path::Path::new(args.readids_prefix())
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_else(|| args.readids_prefix());
    let readids_filename = format!("{}_{}_denovo_reads.txt", name, args.mode_name());
    let readids_writer = std::fs::File::create(&readids_filename)
        .map_err(|e| anyhow!("Failed to create read IDs file {}: {}", readids_filename, e))?;

    let (sender_result, receiver_result) = unbounded();
    let writer_thread = match output_file {
        Some(file) => process_writer_thread(
            WriterBuilder::new().delimiter(b'\t').from_writer(file),
            readids_writer,
            receiver_result,
        ),
        None => process_writer_thread(
            WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout()),
            readids_writer,
            receiver_result,
        ),
    };

    let process_locus = |locus: &Locus, sender: &Sender<Vec<R>>| -> Result<()> {
        let mut recorder = score_output.recorder(locus);
        process_locus_fn(locus, &args, &params_arc, sender, &mut recorder)?;
        score_output.write_locus(recorder)
    };

    // Capture errors until all producer/worker threads and the result writer have joined.
    let processing_result = (|| -> Result<()> {
        match args.trid() {
            Some(trid) => {
                let locus =
                    locus::get_locus(&genome_reader, &mut catalog_reader, trid, args.flank_len())?;
                process_locus(&locus, &sender_result)?;
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

                // Owning the receiver in this scope unblocks the producer on any error.
                let locus_result = (|| -> Result<()> {
                    if args.num_threads() == 1 {
                        log::debug!("Single-threaded mode");
                        for locus in receiver_locus {
                            match locus {
                                Ok(locus) => process_locus(&locus, &sender_result)?,
                                Err(err) => log::error!("Locus Processing: {err:#}"),
                            }
                        }
                    } else {
                        log::debug!(
                            "Multi-threaded mode: estimated available cores: {}",
                            thread::available_parallelism().unwrap().get()
                        );
                        let pool = initialize_thread_pool(args.num_threads())?;
                        pool.install(|| {
                            receiver_locus.into_iter().par_bridge().try_for_each_with(
                                &sender_result,
                                |sender, result| match result {
                                    Ok(locus) => process_locus(&locus, sender),
                                    Err(err) => {
                                        log::error!("Locus Processing: {err:#}");
                                        Ok(())
                                    }
                                },
                            )
                        })?;
                    }
                    Ok(())
                })();
                let stream_result = locus_stream_thread
                    .join()
                    .map_err(|_| anyhow!("Locus stream thread panicked"));
                locus_result?;
                match stream_result? {
                    Ok(_) => log::trace!("Locus stream thread finished"),
                    Err(e) => log::error!("Locus streaming failed: {e}"),
                }
            }
        }
        Ok(())
    })();
    drop(sender_result);
    let writer_result = writer_thread
        .join()
        .map_err(|_| anyhow!("Result writer thread panicked"));
    processing_result?;
    writer_result?;
    score_output.finish()
}

fn process_writer_thread<T: Write + Send + 'static, R: AlleleResultExt>(
    mut tsv_writer: Writer<T>,
    mut readids_writer: std::fs::File,
    receiver: Receiver<Vec<R>>,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        for results in &receiver {
            for (i, row) in results.iter().enumerate() {
                if let Some(read_ids) = row.read_ids()
                    && !read_ids.is_empty()
                {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;
    use tempfile::NamedTempFile;

    #[test]
    fn test_cli_scoring_changes_alignment_penalties() {
        let reference = NamedTempFile::new().unwrap();
        let bed = NamedTempFile::new().unwrap();
        let cli = Cli::try_parse_from([
            "trgt-denovo",
            "duo",
            "-r",
            reference.path().to_str().unwrap(),
            "-b",
            bed.path().to_str().unwrap(),
            "-o",
            "out.tsv",
            "-1",
            "a",
            "-2",
            "b",
            "--aln-scoring",
            "12,4,2,24,1",
        ])
        .unwrap();
        let Command::Duo(args) = cli.command else {
            panic!("Expected duo command");
        };
        let mut aligner = create_aligner_with_scoring(args.aln_scoring).unwrap();
        aligner.align_end_to_end(b"ACGT", b"AGGT");
        assert_eq!(aligner.cigar_score_clipped(0), -12);
    }

    #[test]
    fn test_trid_writes_selected_locus_to_output_file() {
        use rust_htslib::{bam, bcf};

        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        let bed = directory.path().join("repeats.bed");
        let output = directory.path().join("out.tsv");
        std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(400))).unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t400\t6\t400\t401\n",
        )
        .unwrap();
        std::fs::write(
            &bed,
            b"chr1\t100\t120\tID=selected;MOTIFS=A;STRUC=(A)n\n\
              chr1\t200\t220\tID=other;MOTIFS=A;STRUC=(A)n\n",
        )
        .unwrap();

        // Own the command's uniquely named sidecar so it is removed even on failure.
        let readids = tempfile::Builder::new()
            .suffix("_duo_denovo_reads.txt")
            .tempfile_in(".")
            .unwrap();
        let sample_name = readids
            .path()
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .strip_suffix("_duo_denovo_reads.txt")
            .unwrap();
        let bam_path = directory.path().join(sample_name);
        let mut bam_header = bam::Header::new();
        bam_header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 400),
        );
        drop(bam::Writer::from_path(&bam_path, &bam_header, bam::Format::Bam).unwrap());
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let vcf_path = directory.path().join("sample.vcf.gz");
        let mut vcf_header = bcf::Header::new();
        vcf_header.push_record(b"##contig=<ID=chr1,length=400>");
        drop(bcf::Writer::from_path(&vcf_path, &vcf_header, false, bcf::Format::Vcf).unwrap());
        bcf::index::build(&vcf_path, None, 1, bcf::index::Type::Csi(14)).unwrap();

        let cli = Cli::try_parse_from([
            "trgt-denovo",
            "duo",
            "-r",
            reference.to_str().unwrap(),
            "-b",
            bed.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--sample-a-vcf",
            vcf_path.to_str().unwrap(),
            "--sample-a-bam",
            bam_path.to_str().unwrap(),
            "--sample-b-vcf",
            vcf_path.to_str().unwrap(),
            "--sample-b-bam",
            bam_path.to_str().unwrap(),
            "--trid",
            "selected",
        ])
        .unwrap();
        let Command::Duo(args) = cli.command else {
            panic!("Expected duo command");
        };
        crate::commands::duo(args).unwrap();

        #[derive(Debug, serde::Deserialize, PartialEq)]
        struct OutputLocus {
            chrom: String,
            start: u32,
            end: u32,
            trid: String,
        }
        let loci: Vec<OutputLocus> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output)
            .unwrap()
            .deserialize()
            .collect::<std::result::Result<_, _>>()
            .unwrap();
        assert_eq!(
            loci,
            [OutputLocus {
                chrom: "chr1".to_owned(),
                start: 100,
                end: 120,
                trid: "selected".to_owned(),
            }]
        );
    }

    fn alignment_scores_fixture(mode: &str) -> (tempfile::TempDir, tempfile::NamedTempFile) {
        use rust_htslib::{
            bam::{
                self,
                record::{Aux, Cigar, CigarString},
            },
            bcf::{self, record::GenotypeAllele},
        };

        let directory = tempfile::tempdir().unwrap();
        let root = directory.path();
        std::fs::write(
            root.join("reference.fa"),
            format!(">chr1\n{}\n", "A".repeat(400)),
        )
        .unwrap();
        std::fs::write(root.join("reference.fa.fai"), b"chr1\t400\t6\t400\t401\n").unwrap();
        std::fs::write(
            root.join("repeats.bed"),
            b"chr1\t100\t104\tID=first;MOTIFS=ACGT;STRUC=(ACGT)n\n\
              chr1\t200\t204\tID=second;MOTIFS=ACGT;STRUC=(ACGT)n\n",
        )
        .unwrap();

        let suffix = format!("_{mode}_denovo_reads.txt");
        let readids = tempfile::Builder::new()
            .suffix(&suffix)
            .tempfile_in(".")
            .unwrap();
        let sample = readids
            .path()
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .strip_suffix(&suffix)
            .unwrap();
        let bam_path = root.join(sample);
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 400),
        );
        let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for (start, trid) in [(100, "first"), (200, "second")] {
            for (index, bases) in [b"AAAAACGTAAAA", b"AAAAATGTAAAA"].into_iter().enumerate() {
                let mut record = bam::Record::new();
                record.set(
                    format!("r{index}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(12)])),
                    bases,
                    &[30; 12],
                );
                record.set_tid(0);
                record.set_pos(start - 4);
                record.set_flags(0);
                record.push_aux(b"TR", Aux::String(trid)).unwrap();
                record.push_aux(b"AL", Aux::U8(index as u8)).unwrap();
                writer.write(&record).unwrap();
            }
        }
        drop(writer);
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
        // The explicit input path remains stable without depending on the sidecar's random name.
        std::fs::write(root.join("bam-path"), bam_path.to_str().unwrap()).unwrap();

        let mut header = bcf::Header::new();
        for line in [
            "##contig=<ID=chr1,length=400>",
            "##INFO=<ID=TRID,Number=1,Type=String,Description=\"Repeat ID\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=MC,Number=1,Type=String,Description=\"Motif counts\">",
            "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Allele lengths\">",
        ] {
            header.push_record(line.as_bytes());
        }
        header.push_sample(b"sample");
        let vcf_path = root.join("sample.vcf.gz");
        let mut writer =
            bcf::Writer::from_path(&vcf_path, &header, false, bcf::Format::Vcf).unwrap();
        for (start, trid) in [(100, b"first".as_slice()), (200, b"second".as_slice())] {
            let mut record = writer.empty_record();
            record.set_rid(Some(0));
            record.set_pos(start);
            record.set_alleles(&[b"ACGT", b"ATGT"]).unwrap();
            record.push_info_string(b"TRID", &[trid]).unwrap();
            record
                .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                .unwrap();
            record
                .push_format_string(b"MC", &[b"1,1".as_slice()])
                .unwrap();
            record.push_format_integer(b"AL", &[4, 4]).unwrap();
            writer.write(&record).unwrap();
        }
        drop(writer);
        bcf::index::build(&vcf_path, None, 1, bcf::index::Type::Csi(14)).unwrap();
        (directory, readids)
    }

    fn run_alignment_scores_fixture(
        root: &Path,
        mode: &str,
        extra: &[&str],
        output: &Path,
    ) -> Result<()> {
        let reference = root.join("reference.fa");
        let bed = root.join("repeats.bed");
        let vcf = root.join("sample.vcf.gz");
        let bam = std::fs::read_to_string(root.join("bam-path")).unwrap();
        let mut arguments = vec![
            "trgt-denovo",
            "--quiet",
            mode,
            "--reference",
            reference.to_str().unwrap(),
            "--bed",
            bed.to_str().unwrap(),
            "--out",
            output.to_str().unwrap(),
            "--flank-len",
            "4",
        ];
        let roles = if mode == "trio" {
            vec![
                ("--mother-vcf", "--mother-bam"),
                ("--father-vcf", "--father-bam"),
                ("--child-vcf", "--child-bam"),
            ]
        } else {
            vec![
                ("--sample-a-vcf", "--sample-a-bam"),
                ("--sample-b-vcf", "--sample-b-bam"),
            ]
        };
        for (vcf_flag, bam_flag) in roles {
            arguments.extend([vcf_flag, vcf.to_str().unwrap(), bam_flag, &bam]);
        }
        arguments.extend(extra);
        match Cli::try_parse_from(arguments)?.command {
            Command::Trio(args) => crate::commands::trio(args),
            Command::Duo(args) => crate::commands::duo(args),
        }
    }

    #[test]
    fn test_alignment_scores_exports_comparisons_without_changing_calls() {
        for mode in ["trio", "duo"] {
            let (directory, _sidecar) = alignment_scores_fixture(mode);
            let root = directory.path();
            let baseline = root.join(format!("{mode}-baseline.tsv"));
            let output = root.join(format!("{mode}-calls.tsv"));
            let scores = root.join(format!("{mode}-scores.tsv"));
            run_alignment_scores_fixture(root, mode, &[], &baseline).unwrap();
            run_alignment_scores_fixture(
                root,
                mode,
                &["--alignment-scores", scores.to_str().unwrap()],
                &output,
            )
            .unwrap();
            assert_eq!(
                std::fs::read(&baseline).unwrap(),
                std::fs::read(&output).unwrap()
            );
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(&scores)
                .unwrap();
            assert_eq!(
                reader.headers().unwrap(),
                &csv::StringRecord::from(vec![
                    "chrom",
                    "start",
                    "end",
                    "trid",
                    "candidate_genotype",
                    "candidate_index",
                    "sample_role",
                    "source_genotype",
                    "source_index",
                    "read_id",
                    "score",
                ])
            );
            let mut actual: Vec<_> = reader.records().map(std::result::Result::unwrap).collect();
            let mut expected = Vec::new();
            for (start, end, trid) in [("100", "104", "first"), ("200", "204", "second")] {
                for candidate in ["0", "1"] {
                    let roles = if mode == "trio" {
                        vec!["mother", "father", "child"]
                    } else {
                        vec!["b", "a"]
                    };
                    for role in roles {
                        for (source, read_id) in [("0", "r0"), ("1", "r1")] {
                            if (role == "child" || role == "a") && source != candidate {
                                continue;
                            }
                            expected.push(csv::StringRecord::from(vec![
                                "chr1",
                                start,
                                end,
                                trid,
                                candidate,
                                candidate,
                                role,
                                source,
                                source,
                                read_id,
                                if source == candidate { "0" } else { "-8" },
                            ]));
                        }
                    }
                }
            }
            actual.sort_by(|a, b| a.iter().cmp(b.iter()));
            expected.sort_by(|a, b| a.iter().cmp(b.iter()));
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_alignment_scores_parallel_batches_match_serial_comparisons() {
        for mode in ["trio", "duo"] {
            let (directory, _sidecar) = alignment_scores_fixture(mode);
            let root = directory.path();
            let serial_scores = root.join("serial-scores.tsv");
            let parallel_scores = root.join("parallel-scores.tsv");
            for (threads, scores, calls) in [
                ("1", &serial_scores, root.join("serial-calls.tsv")),
                ("4", &parallel_scores, root.join("parallel-calls.tsv")),
            ] {
                run_alignment_scores_fixture(
                    root,
                    mode,
                    &[
                        "-@",
                        threads,
                        "--alignment-scores",
                        scores.to_str().unwrap(),
                    ],
                    &calls,
                )
                .unwrap();
            }
            let read_rows = |path: &Path| {
                csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .from_path(path)
                    .unwrap()
                    .records()
                    .map(std::result::Result::unwrap)
                    .collect::<Vec<_>>()
            };
            let mut serial = read_rows(&serial_scores);
            let mut parallel = read_rows(&parallel_scores);
            serial.sort_by(|a, b| a.iter().cmp(b.iter()));
            parallel.sort_by(|a, b| a.iter().cmp(b.iter()));
            assert_eq!(parallel, serial);
        }
    }

    #[test]
    fn test_alignment_scores_reject_output_aliases_before_truncation() {
        for mode in ["trio", "duo"] {
            let (fixture, _sidecar) = alignment_scores_fixture(mode);
            let output_directory = tempfile::tempdir_in(".").unwrap();
            let relative_directory =
                std::path::PathBuf::from(output_directory.path().file_name().unwrap());
            let relative = relative_directory.join("calls.tsv");
            let absolute = std::env::current_dir().unwrap().join(&relative);
            let mut aliases = vec![absolute];
            #[cfg(unix)]
            {
                let directory_alias = relative_directory.join("directory-alias");
                std::os::unix::fs::symlink(".", &directory_alias).unwrap();
                aliases.push(directory_alias.join("calls.tsv"));
                let file_alias = relative_directory.join("file-alias.tsv");
                std::os::unix::fs::symlink("calls.tsv", &file_alias).unwrap();
                aliases.push(file_alias);
            }
            for existing in [false, true] {
                if existing {
                    std::fs::write(&relative, b"original results").unwrap();
                }
                for alias in &aliases {
                    let result = run_alignment_scores_fixture(
                        fixture.path(),
                        mode,
                        &["--alignment-scores", alias.to_str().unwrap()],
                        &relative,
                    );
                    assert_eq!(
                        (result.is_err(), std::fs::read(&relative).ok()),
                        (true, existing.then(|| b"original results".to_vec())),
                        "{mode}: aliased outputs must be rejected before either file is created",
                    );
                }
            }
        }
    }

    #[test]
    fn test_alignment_scores_quick_skips_leave_header_only_export() {
        const HEADER: &str = "chrom\tstart\tend\ttrid\tcandidate_genotype\tcandidate_index\tsample_role\tsource_genotype\tsource_index\tread_id\tscore\n";
        for mode in ["trio", "duo"] {
            let (directory, _sidecar) = alignment_scores_fixture(mode);
            let root = directory.path();
            let scores = root.join("quick-scores.tsv");
            let baseline = root.join("quick-baseline.tsv");
            let calls = root.join("quick-calls.tsv");
            run_alignment_scores_fixture(root, mode, &["--quick", "AL"], &baseline).unwrap();
            run_alignment_scores_fixture(
                root,
                mode,
                &[
                    "--quick",
                    "AL",
                    "--alignment-scores",
                    scores.to_str().unwrap(),
                ],
                &calls,
            )
            .unwrap();
            assert_eq!(
                std::fs::read(&calls).unwrap(),
                std::fs::read(&baseline).unwrap()
            );
            assert_eq!(std::fs::read_to_string(&scores).unwrap(), HEADER);
        }
    }
}
