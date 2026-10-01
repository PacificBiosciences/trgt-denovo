use crate::{allele::Allele, locus::Locus, read::TrgtRead};
use anyhow::{Context, Result, anyhow, ensure};
use csv::{Writer, WriterBuilder};
use serde::Serialize;
use std::{
    fs::File,
    io::{BufWriter, ErrorKind, Write},
    path::{Path, PathBuf},
    sync::Mutex,
};

pub trait ScoreRecorder {
    fn record(
        &mut self,
        sample_role: &'static str,
        candidate: &Allele,
        source: &Allele,
        read: &TrgtRead,
        score: i32,
    );
}

#[derive(Default)]
pub struct NoScores;

impl ScoreRecorder for NoScores {
    #[inline(always)]
    fn record(
        &mut self,
        _sample_role: &'static str,
        _candidate: &Allele,
        _source: &Allele,
        _read: &TrgtRead,
        _score: i32,
    ) {
    }
}

pub trait ScoreOutput: Sync {
    type Recorder<'a>: ScoreRecorder
    where
        Self: 'a;

    fn recorder<'a>(&'a self, locus: &'a Locus) -> Self::Recorder<'a>;
    fn write_locus(&self, recorder: Self::Recorder<'_>) -> Result<()>;
    fn finish(&self) -> Result<()>;
}

impl ScoreOutput for NoScores {
    type Recorder<'a> = NoScores;

    #[inline(always)]
    fn recorder<'a>(&'a self, _locus: &'a Locus) -> NoScores {
        NoScores
    }

    #[inline(always)]
    fn write_locus(&self, _recorder: NoScores) -> Result<()> {
        Ok(())
    }

    #[inline(always)]
    fn finish(&self) -> Result<()> {
        Ok(())
    }
}

pub struct AlignmentScoreFile {
    writer: Mutex<BufWriter<File>>,
}

impl AlignmentScoreFile {
    pub fn create(path: impl AsRef<Path>, main_output: Option<&Path>) -> Result<Self> {
        let path = path.as_ref();
        if let Some(main_output) = main_output {
            ensure!(
                output_destination(path)? != output_destination(main_output)?,
                "--alignment-scores must use a different path from --out"
            );
        }
        let file = File::create(path).with_context(|| {
            format!("Failed to create alignment scores file {}", path.display())
        })?;
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(file);
        writer
            .write_record([
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
            .context("Failed to write alignment scores header")?;
        let file = writer
            .into_inner()
            .map_err(|error| error.into_error())
            .context("Failed to flush alignment scores header")?;
        Ok(Self {
            writer: Mutex::new(BufWriter::new(file)),
        })
    }
}

fn output_destination(path: &Path) -> Result<PathBuf> {
    let destination = match path.symlink_metadata() {
        Ok(_) => path.canonicalize(),
        Err(error) if error.kind() == ErrorKind::NotFound => {
            let parent = path
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            let filename = path
                .file_name()
                .ok_or_else(|| anyhow!("Output path has no filename: {}", path.display()))?;
            parent.canonicalize().map(|parent| parent.join(filename))
        }
        Err(error) => Err(error),
    };
    destination.with_context(|| format!("Failed to resolve output path {}", path.display()))
}

/// Serializes borrowed identities immediately, retaining only one locus's TSV bytes.
pub struct LocusScores<'a> {
    locus: &'a Locus,
    writer: Writer<Vec<u8>>,
    error: Option<csv::Error>,
}

#[derive(Serialize)]
struct ScoreRow<'a> {
    chrom: &'a str,
    start: u32,
    end: u32,
    trid: &'a str,
    candidate_genotype: usize,
    candidate_index: usize,
    sample_role: &'static str,
    source_genotype: usize,
    source_index: usize,
    read_id: &'a str,
    score: i32,
}

impl ScoreRecorder for LocusScores<'_> {
    fn record(
        &mut self,
        sample_role: &'static str,
        candidate: &Allele,
        source: &Allele,
        read: &TrgtRead,
        score: i32,
    ) {
        if self.error.is_some() {
            return;
        }
        if let Err(error) = self.writer.serialize(ScoreRow {
            chrom: &self.locus.region.contig,
            start: self.locus.region.start,
            end: self.locus.region.end,
            trid: &self.locus.id,
            candidate_genotype: candidate.genotype,
            candidate_index: candidate.index,
            sample_role,
            source_genotype: source.genotype,
            source_index: source.index,
            read_id: &read.name,
            score,
        }) {
            self.error = Some(error);
        }
    }
}

impl ScoreOutput for AlignmentScoreFile {
    type Recorder<'a> = LocusScores<'a>;

    fn recorder<'a>(&'a self, locus: &'a Locus) -> LocusScores<'a> {
        LocusScores {
            locus,
            writer: WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_writer(Vec::new()),
            error: None,
        }
    }

    fn write_locus(&self, recorder: LocusScores<'_>) -> Result<()> {
        if let Some(error) = recorder.error {
            return Err(error).context("Failed to serialize alignment scores");
        }
        let bytes = recorder
            .writer
            .into_inner()
            .map_err(|error| error.into_error())
            .context("Failed to finalize alignment scores batch")?;
        if bytes.is_empty() {
            return Ok(());
        }
        self.writer
            .lock()
            .map_err(|_| anyhow!("Alignment scores writer lock poisoned"))?
            .write_all(&bytes)
            .context("Failed to write alignment scores batch")
    }

    fn finish(&self) -> Result<()> {
        self.writer
            .lock()
            .map_err(|_| anyhow!("Alignment scores writer lock poisoned"))?
            .flush()
            .context("Failed to flush alignment scores file")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::GenomicRegion;
    use tempfile::tempdir;

    const HEADER: &str = "chrom\tstart\tend\ttrid\tcandidate_genotype\tcandidate_index\tsample_role\tsource_genotype\tsource_index\tread_id\tscore\n";

    #[test]
    fn empty_export_contains_header() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("scores.tsv");
        let output = AlignmentScoreFile::create(&path, None).unwrap();
        output.finish().unwrap();
        assert_eq!(std::fs::read_to_string(path).unwrap(), HEADER);
    }

    #[test]
    fn export_escapes_metadata_and_preserves_score_identity() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("scores.tsv");
        let output = AlignmentScoreFile::create(&path, None).unwrap();
        let locus = Locus {
            id: "repeat\"\nidentifier".into(),
            motifs: vec![],
            left_flank: vec![],
            right_flank: vec![],
            region: GenomicRegion::new("chr\t1", 100, 120).unwrap(),
        };
        let candidate = Allele {
            genotype: 3,
            index: 1,
            ..Allele::dummy(vec![])
        };
        let source = Allele {
            genotype: 7,
            index: 0,
            ..Allele::dummy(vec![])
        };
        let read = TrgtRead {
            name: "read\t\"quoted\"\nnext".into(),
            bases: Box::new([]),
            classification: None,
            haplotype: None,
            start_offset: None,
            end_offset: None,
            mismatch_offsets: None,
            pos: 0,
        };
        let mut recorder = output.recorder(&locus);
        recorder.record("parent", &candidate, &source, &read, -42);
        output.write_locus(recorder).unwrap();
        output.finish().unwrap();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .unwrap();
        assert_eq!(
            reader.headers().unwrap().iter().collect::<Vec<_>>(),
            HEADER.trim_end().split('\t').collect::<Vec<_>>()
        );
        let rows = reader
            .records()
            .map(|row| row.unwrap().iter().map(str::to_owned).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(
            rows,
            vec![vec![
                "chr\t1",
                "100",
                "120",
                "repeat\"\nidentifier",
                "3",
                "1",
                "parent",
                "7",
                "0",
                "read\t\"quoted\"\nnext",
                "-42"
            ]]
        );
    }

    #[test]
    fn export_creation_errors_are_returned() {
        let dir = tempdir().unwrap();
        assert!(AlignmentScoreFile::create(dir.path(), None).is_err());
        assert!(
            AlignmentScoreFile::create(dir.path().join("missing").join("scores.tsv"), None)
                .is_err()
        );
    }

    #[test]
    fn export_batch_write_and_final_flush_errors_are_returned() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let locus = Locus {
            id: "repeat".into(),
            motifs: vec![],
            left_flank: vec![],
            right_flank: vec![],
            region: GenomicRegion::new("chr1", 100, 120).unwrap(),
        };
        let allele = Allele::dummy(vec![]);
        let read = crate::read::TrgtReadBuilder::default().build();
        for capacity in [0, 8192] {
            let output = AlignmentScoreFile {
                writer: Mutex::new(BufWriter::with_capacity(
                    capacity,
                    File::open(file.path()).unwrap(),
                )),
            };
            let mut recorder = output.recorder(&locus);
            recorder.record("parent", &allele, &allele, &read, -42);
            if capacity == 0 {
                assert!(output.write_locus(recorder).is_err());
            } else {
                output.write_locus(recorder).unwrap();
                assert!(output.finish().is_err());
            }
        }
    }
}
