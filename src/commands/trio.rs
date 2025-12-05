use crate::{
    aligner::WFAligner,
    cli::TrioArgs,
    commands::shared::{self, AlleleResultExt, Args},
    handles::{SampleInput, TrioLocalData},
    locus::Locus,
    model::{Params, QuickMode},
    trio::allele::{self, AlleleResult},
};
use anyhow::Result;
use crossbeam_channel::Sender;
use std::{cell::RefCell, path::Path, sync::Arc};

impl TrioArgs {
    /// Returns the `SampleInput` for the mother sample
    pub fn mother_input(&self) -> SampleInput {
        if let Some(ref prefix) = self.mother_prefix {
            SampleInput::Prefix(prefix)
        } else {
            SampleInput::Explicit {
                vcf: self.mother_vcf.as_ref().unwrap(),
                bam: self.mother_bam.as_ref().unwrap(),
            }
        }
    }

    /// Returns the `SampleInput` for the father sample
    pub fn father_input(&self) -> SampleInput {
        if let Some(ref prefix) = self.father_prefix {
            SampleInput::Prefix(prefix)
        } else {
            SampleInput::Explicit {
                vcf: self.father_vcf.as_ref().unwrap(),
                bam: self.father_bam.as_ref().unwrap(),
            }
        }
    }

    /// Returns the `SampleInput` for the child sample
    pub fn child_input(&self) -> SampleInput {
        if let Some(ref prefix) = self.child_prefix {
            SampleInput::Prefix(prefix)
        } else {
            SampleInput::Explicit {
                vcf: self.child_vcf.as_ref().unwrap(),
                bam: self.child_bam.as_ref().unwrap(),
            }
        }
    }

    /// Returns a string identifier for the child sample (used for output naming)
    fn child_identifier(&self) -> &str {
        if let Some(ref prefix) = self.child_prefix {
            prefix
        } else {
            self.child_bam.as_ref().unwrap().to_str().unwrap_or("child")
        }
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
        self.child_identifier()
    }
    fn mode_name(&self) -> &str {
        "trio"
    }
    fn preflight_check(&self) -> Result<()> {
        TrioLocalData::new(self.mother_input(), self.father_input(), self.child_input())?;
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
    static LOCAL_FAMILY_DATA: RefCell<Option<TrioLocalData>> = const { RefCell::new(None) };
}

pub fn trio(args: TrioArgs) -> Result<()> {
    shared::run(args, process_locus)
}

fn process_locus(
    locus: &Locus,
    args: &TrioArgs,
    params_arc: &Arc<Params>,
    sender_result: &Sender<Vec<AlleleResult>>,
) {
    ALIGNER.with(|aligner| {
        LOCAL_FAMILY_DATA.with(|local_family_data| {
            let mut family_data = local_family_data.borrow_mut();
            if family_data.is_none() {
                *family_data = Some(
                    TrioLocalData::new(
                        args.mother_input(),
                        args.father_input(),
                        args.child_input(),
                    )
                    .expect("Failed to initialize FamilyLocalData"),
                );
            }
            if let Ok(result) = allele::process_alleles(
                locus,
                family_data.as_mut().unwrap(),
                params_arc,
                &mut aligner.borrow_mut(),
            ) {
                sender_result.send(result).unwrap();
            }
        });
    });
}
