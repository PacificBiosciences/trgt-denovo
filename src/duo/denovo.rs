use crate::aligner::WFAligner;
use crate::allele::AlleleSet;
use crate::denovo::{
    align_allele, align_alleleset, get_overlap_coverage, get_score_count_diff, get_top_other_score,
};
use crate::math;
use crate::model::{DenovoStatus, DenovoType, Params};

/// Represents a de novo allele event with associated scoring and classification information.
#[derive(Debug)]
pub struct DenovoAllele {
    /// TRGT genotype index of the allele.
    pub genotype: usize,
    /// The number of A reads with de novo signal relative to the B alleles
    pub denovo_coverage: usize,
    /// The number of reads covering this allele in A.
    pub a_coverage: usize,
    /// The number of reads covering this locus.
    pub allele_coverage: usize,
    /// Average difference in alignment scores compared to B alleles.
    pub mean_diff_b: f32,
    /// Status indicating whether the allele is de novo and its type.
    pub denovo_status: DenovoStatus,
    /// Index used to identify the allele.
    pub index: usize,
    /// The number of B reads per allele that overlap with the A allele
    pub b_overlap_coverage: Vec<i32>,
    /// The read IDs that contribute to de novo coverage
    pub read_ids: Option<Vec<String>>,
}

/// Assesses de novo alleles by comparing between two samples.
///
/// Calculates de novo coverage and type of potential de novo events.
///
/// # Arguments
///
/// * `b_gts` - A slice of alleles from the mother.
/// * `a_gts` - A slice of alleles from the father.
/// * `params` - Parameters used by downstream functions.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// An iterator over `DenovoAllele` instances with updated de novo information.
pub fn assess_denovo<'a>(
    a_gts: &'a AlleleSet,
    b_gts: &'a AlleleSet,
    params: &Params,
    aligner: &mut WFAligner,
) -> impl Iterator<Item = DenovoAllele> + 'a {
    let mut dnrs = Vec::with_capacity(a_gts.len());

    for denovo_allele in a_gts.iter() {
        let b_align_scores = align_alleleset(b_gts, &denovo_allele.seq, params.clip_len, aligner);

        let a_align_scores =
            align_allele(denovo_allele, &denovo_allele.seq, params.clip_len, aligner);

        let child_read_scores: Vec<(String, i32)> = a_align_scores
            .iter()
            .zip(&denovo_allele.read_aligns)
            .map(|(&score, (read_info, _))| (read_info.name.clone(), score))
            .collect();

        let (denovo_coverage, mean_diff_b, read_ids) = get_denovo_coverage(
            &b_align_scores,
            &a_align_scores,
            &child_read_scores,
            params.parent_quantile,
        );

        let child_score_threshold = math::median(&a_align_scores).unwrap_or(f64::MAX);
        let b_overlap_coverage = get_overlap_coverage(child_score_threshold, &b_align_scores);

        let denovo_status = match denovo_coverage {
            0 => DenovoStatus::NotDenovo,
            _ => DenovoStatus::Denovo(DenovoType::Unclear),
        };

        dnrs.push(DenovoAllele {
            genotype: denovo_allele.genotype,
            denovo_coverage,
            a_coverage: a_gts.iter().map(|vec| vec.read_aligns.len()).sum(),
            allele_coverage: denovo_allele.read_aligns.len(),
            denovo_status,
            mean_diff_b,
            index: denovo_allele.index,
            b_overlap_coverage,
            read_ids: if read_ids.is_empty() {
                None
            } else {
                Some(read_ids)
            },
        });
    }

    dnrs.into_iter().filter_map(Some)
}

fn get_denovo_coverage(
    b_align_scores: &[Vec<i32>],
    a_align_scores: &[i32],
    child_reads: &[(String, i32)],
    p_quantile: f64,
) -> (usize, f32, Vec<String>) {
    let top_b_score = get_top_other_score(b_align_scores, p_quantile).unwrap();

    let mut denovo_read_ids = Vec::new();
    for (read_name, score) in child_reads {
        if (*score as f64) > top_b_score {
            denovo_read_ids.push(read_name.clone());
        }
    }

    let (score, mean_diff_b) = get_score_count_diff(top_b_score, a_align_scores);

    (score, mean_diff_b, denovo_read_ids)
}
