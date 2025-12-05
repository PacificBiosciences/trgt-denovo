use crate::{
    aligner::WFAligner,
    cli::DuoArgs,
    commands::shared::{self, AlleleResultExt, Args},
    duo::allele::{self, AlleleResult},
    handles::{DuoLocalData, SampleInput},
    locus::Locus,
    model::{Params, QuickMode},
};
use anyhow::Result;
use crossbeam_channel::Sender;
use std::{cell::RefCell, path::Path, sync::Arc};

impl DuoArgs {
    /// Returns the `SampleInput` for sample A
    pub fn sample_a_input(&self) -> SampleInput {
        if let Some(ref prefix) = self.a_prefix {
            SampleInput::Prefix(prefix)
        } else {
            SampleInput::Explicit {
                vcf: self.a_vcf.as_ref().unwrap(),
                bam: self.a_bam.as_ref().unwrap(),
            }
        }
    }

    /// Returns the `SampleInput` for sample B
    pub fn sample_b_input(&self) -> SampleInput {
        if let Some(ref prefix) = self.b_prefix {
            SampleInput::Prefix(prefix)
        } else {
            SampleInput::Explicit {
                vcf: self.b_vcf.as_ref().unwrap(),
                bam: self.b_bam.as_ref().unwrap(),
            }
        }
    }

    /// Returns a string identifier for sample A (used for output naming)
    fn sample_a_identifier(&self) -> &str {
        if let Some(ref prefix) = self.a_prefix {
            prefix
        } else {
            self.a_bam.as_ref().unwrap().to_str().unwrap_or("sample_a")
        }
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
    fn output_path(&self) -> &str {
        &self.output_path
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
        DuoLocalData::new(self.sample_a_input(), self.sample_b_input())?;
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
    static ALIGNER: RefCell<WFAligner> = RefCell::new(shared::create_aligner_with_scoring());
    static LOCAL_DUO_DATA: RefCell<Option<DuoLocalData>> = const { RefCell::new(None) };
}

pub fn duo(args: DuoArgs) -> Result<()> {
    shared::run(args, process_locus)
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
                    DuoLocalData::new(args.sample_a_input(), args.sample_b_input())
                        .expect("Failed to initialize DuoLocalData"),
                );
            }
            if let Ok(result) = allele::process_alleles(
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
