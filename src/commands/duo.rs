use crate::{
    alignment_scores::{AlignmentScoreFile, NoScores, ScoreRecorder},
    cli::DuoArgs,
    commands::shared::{self, AlleleResultExt, Args},
    duo::allele::{self, AlleleResult},
    handles::{DuoLocalData, SampleInput},
    locus::Locus,
    model::{Params, QuickMode},
};
use anyhow::Result;
use crossbeam_channel::Sender;
use rust_wfa2::aligner::WFAligner;
use std::{cell::RefCell, path::Path, sync::Arc};

impl DuoArgs {
    /// Returns the `SampleInput` for sample A
    pub fn sample_a_input(&self) -> Result<SampleInput<'_>> {
        SampleInput::from_paths(
            self.a_prefix.as_deref(),
            self.a_vcf.as_deref(),
            self.a_bam.as_deref(),
        )
    }

    /// Returns the `SampleInput` for sample B
    pub fn sample_b_input(&self) -> Result<SampleInput<'_>> {
        SampleInput::from_paths(
            self.b_prefix.as_deref(),
            self.b_vcf.as_deref(),
            self.b_bam.as_deref(),
        )
    }

    /// Returns a string identifier for sample A (used for output naming)
    fn sample_a_identifier(&self) -> &str {
        self.a_prefix
            .as_deref()
            .or_else(|| self.a_bam.as_deref().and_then(Path::to_str))
            .unwrap_or("sample_a")
    }
}

impl Args for DuoArgs {
    fn no_clip_aln(&self) -> bool {
        self.no_clip_aln
    }
    fn flank_len(&self) -> usize {
        self.flank_len
    }
    fn p_quantile(&self) -> f64 {
        self.p_quantile
    }
    fn partition_by_alignment(&self) -> bool {
        self.partition_by_alignment
    }
    fn skip_tr_check(&self) -> bool {
        self.skip_tr_check
    }
    fn quick(&self) -> &Option<QuickMode> {
        &self.quick
    }
    fn bed_filename(&self) -> &Path {
        &self.bed_filename
    }
    fn reference_filename(&self) -> &Path {
        &self.reference_filename
    }
    fn trid(&self) -> &Option<String> {
        &self.trid
    }
    fn output_path(&self) -> Option<&str> {
        self.output_path.as_deref()
    }
    fn num_threads(&self) -> usize {
        self.num_threads
    }
    fn readids_prefix(&self) -> &str {
        self.sample_a_identifier()
    }
    fn mode_name(&self) -> &str {
        "duo"
    }
    fn preflight_check(&self) -> Result<()> {
        DuoLocalData::new(self.sample_a_input()?, self.sample_b_input()?)?;
        Ok(())
    }
}

impl AlleleResultExt for AlleleResult {
    fn read_ids(&self) -> &Option<Vec<String>> {
        &self.read_ids
    }
    fn trid(&self) -> &str {
        &self.trid
    }
}

thread_local! {
    static ALIGNER: RefCell<Option<(crate::model::AlnScoring, WFAligner)>> = const { RefCell::new(None) };
    static LOCAL_DUO_DATA: RefCell<Option<DuoLocalData>> = const { RefCell::new(None) };
}

pub fn duo(args: DuoArgs) -> Result<()> {
    if let Some(path) = args.alignment_scores.as_deref() {
        let file = AlignmentScoreFile::create(path, args.output_path.as_deref().map(Path::new))?;
        shared::run(
            args,
            |locus, args, params, sender, scores| {
                process_locus(locus, args, params, sender, scores)
            },
            file,
        )
    } else {
        shared::run(args, process_locus::<NoScores>, NoScores)
    }
}

fn process_locus<S: ScoreRecorder>(
    locus: &Locus,
    args: &DuoArgs,
    params_arc: &Arc<Params>,
    sender_result: &Sender<Vec<AlleleResult>>,
    scores: &mut S,
) -> Result<()> {
    ALIGNER.with(|aligner| {
        let mut aligner = aligner.borrow_mut();
        let aligner = shared::configured_aligner(&mut aligner, args.aln_scoring)?;
        LOCAL_DUO_DATA.with(|local_family_data| {
            let mut duo_data = local_family_data.borrow_mut();
            if duo_data.is_none() {
                *duo_data = Some(DuoLocalData::new(
                    args.sample_a_input()?,
                    args.sample_b_input()?,
                )?);
            }
            if let Ok(result) = allele::process_alleles(
                locus,
                duo_data.as_mut().unwrap(),
                params_arc,
                aligner,
                scores,
            ) {
                sender_result.send(result).unwrap();
            }
            Ok(())
        })
    })
}
