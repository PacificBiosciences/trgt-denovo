use crate::{
    allele::{Allele, AlleleSet},
    math,
    read::TrgtRead,
};
use rust_wfa2::aligner::{AlignmentStatus, WFAligner};

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
/// * `on_score` - Receives each successful alignment's source allele, read, and exact score.
///
/// # Returns
///
/// A vector of vectors containing alignment scores for each allele.
pub fn align_alleleset(
    gts: &AlleleSet,
    target: &[u8],
    clip_len: usize,
    aligner: &mut WFAligner,
    mut on_score: impl FnMut(&Allele, &TrgtRead, i32),
) -> Vec<Vec<i32>> {
    let mut align_scores = vec![vec![]; gts.len()];
    for (i, allele) in gts.iter().enumerate() {
        for (read, _align) in &allele.read_aligns {
            if let AlignmentStatus::StatusAlgCompleted =
                aligner.align_end_to_end(&read.bases, target).status
            {
                let score = aligner.cigar_score_clipped(clip_len);
                align_scores[i].push(score);
                on_score(allele, read, score);
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
    mut on_score: impl FnMut(&Allele, &TrgtRead, i32),
) -> Vec<i32> {
    let mut align_scores = vec![];
    for (read, _align) in &allele.read_aligns {
        if let AlignmentStatus::StatusAlgCompleted =
            aligner.align_end_to_end(&read.bases, target).status
        {
            let score = aligner.cigar_score_clipped(clip_len);
            align_scores.push(score);
            on_score(allele, read, score);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{commands::shared::create_aligner_with_scoring, read::TrgtReadBuilder};

    #[test]
    fn test_upstream_aligner_preserves_clipped_read_scores() {
        let left = b"AAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let right = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT";
        let read_left = b"AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let read_right = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT";
        let target = [left.as_slice(), &b"ATTT".repeat(10), right].concat();
        let read = [read_left.as_slice(), &b"ATTT".repeat(8), read_right].concat();
        let allele = Allele::dummy(vec![
            TrgtReadBuilder::default().with_bases(read).build(),
            TrgtReadBuilder::default()
                .with_bases(target.clone())
                .build(),
        ]);
        let mut aligner = create_aligner_with_scoring(crate::model::AlnScoring::default()).unwrap();

        // Clipping removes the two flank mismatches, but retains the 8-base gap.
        assert_eq!(
            align_allele(&allele, &target, 0, &mut aligner, |_, _, _| {}),
            [-36, 0]
        );
        assert_eq!(
            align_allele(&allele, &target, 50, &mut aligner, |_, _, _| {}),
            [-20, 0]
        );
        let alleles = AlleleSet {
            alleles: vec![allele, Allele::dummy(vec![])],
            hp_counts: [0; 3],
        };
        assert_eq!(
            align_alleleset(&alleles, &target, 50, &mut aligner, |_, _, _| {}),
            [vec![-20, 0], vec![]]
        );
    }

    #[test]
    fn test_failed_alignments_do_not_contribute_scores() {
        use rust_wfa2::aligner::{AlignmentScope, Heuristics, MemoryModel};

        let target = b"GGGGACGTCCCC";
        let allele = Allele::dummy(
            [target.as_slice(), b"GGGGATGTCCCC", target]
                .into_iter()
                .map(|bases| TrgtReadBuilder::default().with_bases(bases).build())
                .collect(),
        );
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p(8, 4, 2, 24, 1)
            .with_heuristics(Heuristics::wfa2_default())
            .with_max_alignment_steps(1)
            .build()
            .unwrap();

        assert_eq!(
            align_allele(&allele, target, 4, &mut aligner, |_, _, _| {}),
            [0, 0]
        );
        let alleles = AlleleSet {
            alleles: vec![allele],
            hp_counts: [0; 3],
        };
        assert_eq!(
            align_alleleset(&alleles, target, 4, &mut aligner, |_, _, _| {}),
            [vec![0, 0]]
        );
    }

    #[test]
    fn test_alignment_scores_preserve_read_identity_after_failed_alignment() {
        use rust_wfa2::aligner::{AlignmentScope, Heuristics, MemoryModel};

        let target = b"GGGGACGTCCCC";
        let allele = Allele::dummy(
            [
                ("first", target.as_slice()),
                ("failed", b"GGGGATGTCCCC"),
                ("last", target),
            ]
            .into_iter()
            .map(|(name, bases)| {
                let mut read = TrgtReadBuilder::default().with_bases(bases).build();
                read.name = name.to_owned();
                read
            })
            .collect(),
        );
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p(8, 4, 2, 24, 1)
            .with_heuristics(Heuristics::wfa2_default())
            .with_max_alignment_steps(1)
            .build()
            .unwrap();
        let mut records = Vec::new();
        let scores = align_allele(&allele, target, 4, &mut aligner, |source, read, score| {
            records.push((source.index, read.name.clone(), score));
        });
        assert_eq!(
            (scores, records),
            (
                vec![0, 0],
                vec![(0, "first".to_owned(), 0), (0, "last".to_owned(), 0)],
            ),
        );
    }
}
