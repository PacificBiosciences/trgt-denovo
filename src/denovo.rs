use crate::{
    aligner::{AlignmentStatus, WFAligner},
    allele::{Allele, AlleleSet},
    math,
};

/// Aligns reads from alleles to a target sequence and calculates alignment scores.
///
/// This function performs end-to-end alignments of reads from alleles to a given target sequence
/// and returns the alignment scores.
///
/// # Arguments
///
/// * `gts` - A slice of alleles containing the reads to align.
/// * `target` - A slice of the target sequence to align to.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// A vector of vectors containing alignment scores for each allele.
pub fn align_alleleset(
    gts: &AlleleSet,
    target: &[u8],
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Vec<Vec<i32>> {
    let mut align_scores = vec![vec![]; gts.len()];
    for (i, allele) in gts.iter().enumerate() {
        for (read, _align) in &allele.read_aligns {
            if let AlignmentStatus::StatusAlgCompleted =
                aligner.align_end_to_end(&read.bases, target)
            {
                align_scores[i].push(aligner.cigar_score_clipped(clip_len));
            }
        }
    }
    align_scores
}

pub fn align_allele(
    allele: &Allele,
    target: &[u8],
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Vec<i32> {
    let mut align_scores = vec![];
    for (read, _align) in &allele.read_aligns {
        if let AlignmentStatus::StatusAlgCompleted = aligner.align_end_to_end(&read.bases, target) {
            align_scores.push(aligner.cigar_score_clipped(clip_len));
        }
    }
    align_scores
}

/// Count the number of reads for each allele that exceed a given score threshold obtained in a target de novo allele.
///
/// # Arguments
///
/// * `score_threshold` - The score threshold above which alignments are counted.
/// * `aligns` - A slice of vectors, each containing alignment scores of a sample against
///   the putative de novo allele.
///
/// # Returns
///
/// A vector of integers, where each integer represents the count per allele of parental alignments exceeding
/// the score threshold.
pub fn get_overlap_coverage(score_threshold: f64, aligns: &[Vec<i32>]) -> Vec<i32> {
    aligns
        .iter()
        .map(|scores| {
            scores
                .iter()
                .filter(|&&score| score as f64 >= score_threshold)
                .count() as i32
        })
        .collect()
}

/// Determines the top alignment score among other alleles at a given quantile.
///
/// Finds the highest alignment score at a specified quantile across all other
/// alleles. It is used to establish a threshold for comparing to de novo allele alignment scores.
///
/// # Arguments
///
/// * `align_scores`: A slice of vectors containing alignment scores for other alleles.
/// * `quantile`: The quantile used to determine the top score.
///
/// # Returns
///
/// Returns an `Option<f64>` representing the top alignment score at the given quantile, if scores are available.
pub fn get_top_other_score(align_scores: &[Vec<i32>], quantile: f64) -> Option<f64> {
    align_scores
        .iter()
        .filter_map(|scores| {
            let mut scores_f64: Vec<f64> = scores.iter().map(|&x| x as f64).collect();
            math::quantile(&mut scores_f64, quantile)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

/// Calculates the count and mean difference of alignments that exceed the top score.
///
/// Determines the number of alignments with scores higher than the top score
/// and calculates the mean difference of these scores with respect to the top score.
///
/// # Arguments
///
/// * `top_score`: The highest alignment score.
/// * `aligns`: A slice of vectors containing alignment scores for some alleles.
///
/// # Returns
///
/// Returns a tuple containing the count of alignments exceeding the top score and the mean
/// difference of these alignments from the top score.
pub fn get_score_count_diff(top_score: f64, aligns: &[i32]) -> (usize, f32) {
    let (count, sum) = aligns
        .iter()
        .map(|&a| a as f64)
        .filter(|a| a > &top_score)
        .fold((0, 0.0), |(count, sum), a| (count + 1, sum + a));
    let mean_diff = if count > 0 {
        (sum / count as f64 - top_score).abs() as f32
    } else {
        0.0
    };
    (count, mean_diff)
}
