use crate::{
    alignment_scores::{AlignmentScoreFile, NoScores, ScoreRecorder},
    cli::TrioArgs,
    commands::shared::{self, AlleleResultExt, Args},
    handles::{SampleInput, TrioLocalData},
    locus::Locus,
    model::{Params, QuickMode},
    trio::allele::{self, AlleleResult},
};
use anyhow::Result;
use crossbeam_channel::Sender;
use rust_wfa2::aligner::WFAligner;
use std::{cell::RefCell, path::Path, sync::Arc};

impl TrioArgs {
    /// Returns the `SampleInput` for the mother sample
    pub fn mother_input(&self) -> Result<SampleInput<'_>> {
        SampleInput::from_paths(
            self.mother_prefix.as_deref(),
            self.mother_vcf.as_deref(),
            self.mother_bam.as_deref(),
        )
    }

    /// Returns the `SampleInput` for the father sample
    pub fn father_input(&self) -> Result<SampleInput<'_>> {
        SampleInput::from_paths(
            self.father_prefix.as_deref(),
            self.father_vcf.as_deref(),
            self.father_bam.as_deref(),
        )
    }

    /// Returns the `SampleInput` for the child sample
    pub fn child_input(&self) -> Result<SampleInput<'_>> {
        SampleInput::from_paths(
            self.child_prefix.as_deref(),
            self.child_vcf.as_deref(),
            self.child_bam.as_deref(),
        )
    }

    /// Returns a string identifier for the child sample (used for output naming)
    fn child_identifier(&self) -> &str {
        self.child_prefix
            .as_deref()
            .or_else(|| self.child_bam.as_deref().and_then(Path::to_str))
            .unwrap_or("child")
    }
}

impl Args for TrioArgs {
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
        self.child_identifier()
    }
    fn mode_name(&self) -> &str {
        "trio"
    }
    fn preflight_check(&self) -> Result<()> {
        TrioLocalData::new(
            self.mother_input()?,
            self.father_input()?,
            self.child_input()?,
        )?;
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
    static LOCAL_FAMILY_DATA: RefCell<Option<TrioLocalData>> = const { RefCell::new(None) };
}

pub fn trio(args: TrioArgs) -> Result<()> {
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
    args: &TrioArgs,
    params_arc: &Arc<Params>,
    sender_result: &Sender<Vec<AlleleResult>>,
    scores: &mut S,
) -> Result<()> {
    ALIGNER.with(|aligner| {
        let mut aligner = aligner.borrow_mut();
        let aligner = shared::configured_aligner(&mut aligner, args.aln_scoring)?;
        LOCAL_FAMILY_DATA.with(|local_family_data| {
            let mut family_data = local_family_data.borrow_mut();
            if family_data.is_none() {
                *family_data = Some(TrioLocalData::new(
                    args.mother_input()?,
                    args.father_input()?,
                    args.child_input()?,
                )?);
            }
            if let Ok(result) = allele::process_alleles(
                locus,
                family_data.as_mut().unwrap(),
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
