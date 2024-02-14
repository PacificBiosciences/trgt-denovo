//! Provides functionality for assessing de novo allele events.
//!
//! This module contains structures and functions for identifying and classifying de novo
//! mutations in genomic data. It includes methods for comparing allele sequences, determining
//! the origin of alleles, and calculating de novo scores.

use crate::{
    aligner::{AlignmentStatus, WFAligner},
    allele::{Allele, AlleleSet},
    math,
    util::Params,
};
use ndarray::{Array2, ArrayBase, Dim, OwnedRepr};
use serde::Serialize;
use std::{cmp::Ordering, fmt};

/// Represents a de novo allele event with associated scoring and classification information.
#[derive(Debug)]
pub struct DenovoAllele {
    /// TRGT genotype index of the allele.
    pub genotype: usize,
    /// The number of child reads with de novo signal relative to the parental alleles
    pub denovo_coverage: usize,
    /// The number of reads covering this allele in the child.
    pub child_coverage: usize,
    /// The number of reads covering this allele.
    pub allele_coverage: usize,
    /// Average difference in alignment scores compared to the father's alleles.
    pub mean_diff_father: f32,
    /// Average difference in alignment scores compared to the mother's alleles.
    pub mean_diff_mother: f32,
    /// Inferred origin of the allele (father or mother).
    pub allele_origin: AlleleOrigin,
    /// Status indicating whether the allele is de novo and its type.
    pub denovo_status: DenovoStatus,
    /// TRGT count of the motif present in the allele.
    pub motif_count: String,
    /// Index used to identify the allele.
    pub index: usize,
}

/// Enumerates the types of de novo events that can occur.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoType {
    Expansion,
    Contraction,
    Substitution,
    Unclear,
}

/// Represents the status of a de novo event, indicating whether it is de novo and its type.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoStatus {
    /// Indicates a de novo event with a specified type.
    Denovo(DenovoType),
    NotDenovo,
}

impl std::fmt::Display for DenovoType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DenovoType::Expansion => write!(f, "+"),
            DenovoType::Contraction => write!(f, "-"),
            DenovoType::Substitution => write!(f, "="),
            DenovoType::Unclear => write!(f, "?"),
        }
    }
}

impl std::fmt::Display for DenovoStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DenovoStatus::Denovo(denovo_type) => write!(f, "Y:{}", denovo_type),
            DenovoStatus::NotDenovo => write!(f, "X"),
        }
    }
}

/// Represents the number of alleles inherited from a parent.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum AlleleNum {
    One,
    Two,
    Unclear,
}

/// Represents the origin of an allele, specifying which parent it was inherited from.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum AlleleOrigin {
    Father { allele: AlleleNum },
    Mother { allele: AlleleNum },
    Unclear,
}

impl AlleleOrigin {
    pub fn new(value: usize) -> Option<Self> {
        match value {
            0 => Some(AlleleOrigin::Mother {
                allele: AlleleNum::One,
            }),
            1 => Some(AlleleOrigin::Mother {
                allele: AlleleNum::Two,
            }),
            2 => Some(AlleleOrigin::Father {
                allele: AlleleNum::One,
            }),
            3 => Some(AlleleOrigin::Father {
                allele: AlleleNum::Two,
            }),
            _ => None,
        }
    }
}

impl fmt::Display for AlleleNum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlleleNum::One => write!(f, "1"),
            AlleleNum::Two => write!(f, "2"),
            AlleleNum::Unclear => write!(f, "?"),
        }
    }
}

impl fmt::Display for AlleleOrigin {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlleleOrigin::Father { allele } => write!(f, "F:{}", allele),
            AlleleOrigin::Mother { allele } => write!(f, "M:{}", allele),
            AlleleOrigin::Unclear => write!(f, "?"),
        }
    }
}

// Valid allele combinations
static COMBS_1D: [[(usize, usize); 2]; 4] = [
    [(0, 0), (0, 0)],
    [(1, 0), (1, 0)],
    [(2, 0), (2, 0)],
    [(3, 0), (3, 0)],
];
static COMBS_2D: [[(usize, usize); 2]; 8] = [
    [(0, 0), (2, 1)],
    [(0, 0), (3, 1)],
    [(1, 0), (2, 1)],
    [(1, 0), (3, 1)],
    [(2, 0), (0, 1)],
    [(3, 0), (0, 1)],
    [(2, 0), (1, 1)],
    [(3, 0), (1, 1)],
];

/// Assesses de novo alleles by comparing child alleles to parental alleles.
///
/// Calculates de novo coverage and determines the origin and type of potential
/// de novo events.
///
/// # Arguments
///
/// * `mother_gts` - A slice of alleles from the mother.
/// * `father_gts` - A slice of alleles from the father.
/// * `child_gts` - A slice of alleles from the child.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `parent_quantile` - The quantile used for parental allele frequency calculations.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// An iterator over `DenovoAllele` instances with updated de novo information.
pub fn assess_denovo<'a>(
    mother_gts: &'a AlleleSet,
    father_gts: &'a AlleleSet,
    child_gts: &'a AlleleSet,
    params: &Params,
    aligner: &mut WFAligner,
) -> impl Iterator<Item = DenovoAllele> + 'a {
    let mut matrix = Array2::from_elem((4, 2), f64::MIN);
    let mut dnrs = Vec::with_capacity(child_gts.len());

    for (index, denovo_allele) in child_gts.iter().enumerate() {
        let mother_align_scores =
            align_alleleset(mother_gts, &denovo_allele.seq, params.clip_len, aligner);
        let father_align_scores =
            align_alleleset(father_gts, &denovo_allele.seq, params.clip_len, aligner);

        let child_align_scores =
            align_allele(denovo_allele, &denovo_allele.seq, params.clip_len, aligner);

        let (denovo_coverage, mean_diff_father, mean_diff_mother) = get_denovo_coverage(
            &mother_align_scores,
            &father_align_scores,
            &child_align_scores,
            params.parent_quantile,
        );

        update_mean_matrix(
            &mut matrix,
            index,
            &mother_align_scores,
            &father_align_scores,
        );

        dnrs.push(DenovoAllele {
            genotype: denovo_allele.genotype,
            denovo_coverage,
            child_coverage: child_gts.iter().map(|vec| vec.read_aligns.len()).sum(),
            allele_coverage: denovo_allele.read_aligns.len(),
            mean_diff_mother,
            mean_diff_father,
            allele_origin: AlleleOrigin::Unclear,
            denovo_status: DenovoStatus::NotDenovo,
            motif_count: denovo_allele.motif_count.clone(),
            index: denovo_allele.index,
        });
    }
    log::trace!("Mean matrix={:?}", matrix);

    // For each valid allele inheritance pattern calculate the scores
    let mut combs_score: Vec<(f64, [(usize, usize); 2])> = Vec::new();
    if child_gts.len() > 1 {
        for c in COMBS_2D.iter() {
            combs_score.push((c.iter().map(|&(i, j)| matrix[[i, j]]).sum::<f64>(), *c));
        }
    } else {
        for c in COMBS_1D.iter() {
            combs_score.push((matrix[[c[0].0, c[0].1]], *c));
        }
    }
    // Sort the scores (higher more similar)
    combs_score.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    // Absolute difference between top two assignments
    let comb_diff = (combs_score[0].0 - combs_score[1].0).abs();
    log::trace!("comb_diff: {}", comb_diff);
    // Update allele origin and de novo type
    let comb_0 = combs_score[0].1;
    for (index, denovo_allele) in child_gts.iter().enumerate() {
        let dna = &mut dnrs[index];
        if comb_diff < 1.0 {
            dna.allele_origin = AlleleOrigin::Unclear;
        } else {
            dna.allele_origin = AlleleOrigin::new(comb_0[index].0).unwrap();
        }

        dna.denovo_status = match dna.denovo_coverage {
            0 => DenovoStatus::NotDenovo,
            _ => DenovoStatus::Denovo(get_denovo_type(
                mother_gts,
                father_gts,
                &dna.allele_origin,
                &denovo_allele.seq,
            )),
        };
    }
    dnrs.into_iter().filter_map(Some)
}

/// Updates the mean matrix with alignment scores for parental alleles.
///
/// Calculates the mean alignment scores for each parent and updates the matrix
/// used to determine the most likely allele origin.
///
/// # Arguments
///
/// * `matrix` - A mutable reference to the mean matrix.
/// * `index` - The index of the current allele being assessed.
/// * `mother_align_scores` - A slice of alignment scores for the mother's alleles.
/// * `father_align_scores` - A slice of alignment scores for the father's alleles.
pub fn update_mean_matrix(
    matrix: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
    index: usize,
    mother_align_scores: &[Vec<i32>],
    father_align_scores: &[Vec<i32>],
) {
    // Calculates the average of each non-empty vector in a slice of integer vectors,
    // and pairs it with the vector's index. For empty vectors, it pairs negative infinity with their index.
    // Returns a vector of these pairs (average or negative infinity, index).
    let get_mean = |aligns: &[Vec<i32>]| -> Vec<(f64, usize)> {
        aligns
            .iter()
            .enumerate()
            .map(|(i, v)| {
                if v.is_empty() {
                    (f64::NEG_INFINITY, i)
                } else {
                    (v.iter().sum::<i32>() as f64 / v.len() as f64, i)
                }
            })
            .collect()
    };

    let mother_mean: Vec<(f64, usize)> = get_mean(mother_align_scores);
    let father_mean: Vec<(f64, usize)> = get_mean(father_align_scores);

    matrix[[0, index]] = mother_mean[0].0;
    matrix[[1, index]] = mother_mean.get(1).map_or(f64::MIN, |v| v.0);
    matrix[[2, index]] = father_mean[0].0;
    matrix[[3, index]] = father_mean.get(1).map_or(f64::MIN, |v| v.0);
}

/// Compares the lengths of allele sequences to determine the type of de novo event.
///
/// Compares the lengths of the sequences of the parental alleles to the child's allele
/// to classify the de novo event as an expansion, contraction, or substitution.
///
/// # Arguments
///
/// * `gts`: A slice of alleles from a parent.
/// * `allele`: The number of alleles inherited from the parent.
/// * `child_seq`: The sequence of the child's allele.
///
/// # Returns
///
/// Returns an `Option<DenovoType>` indicating the type of de novo event, if determinable.
fn compare_seq_lengths(
    gts: &AlleleSet,
    allele: &AlleleNum,
    child_seq: &[u8],
) -> Option<DenovoType> {
    match allele {
        AlleleNum::One => compare_sequences(&gts[0].seq, child_seq),
        AlleleNum::Two => compare_sequences(&gts[1].seq, child_seq),
        _ => None,
    }
}

/// Compares two sequences to determine the type of de novo event based on length differences.
///
/// Compares the length of a parental allele sequence to the child's allele sequence
/// to classify the de novo event as an expansion, contraction, or substitution.
///
/// # Arguments
///
/// * `seq`: The sequence of the parental allele.
/// * `child_seq`: The sequence of the child's allele.
///
/// # Returns
///
/// Returns an `Option<DenovoType>` indicating the type of de novo event, if determinable.
fn compare_sequences(seq: &[u8], child_seq: &[u8]) -> Option<DenovoType> {
    match seq.len().cmp(&child_seq.len()) {
        Ordering::Less => Some(DenovoType::Expansion),
        Ordering::Greater => Some(DenovoType::Contraction),
        Ordering::Equal => Some(DenovoType::Substitution),
    }
}

/// Determines the type of de novo event based on allele sequences.
///
/// Compare the lengths of the child's allele sequence to the parental allele
/// sequences to classify the de novo event as an expansion, contraction, substitution, or unclear.
///
/// # Arguments
///
/// * `mother_gts` - A slice of alleles from the mother.
/// * `father_gts` - A slice of alleles from the father.
/// * `allele_origin` - A reference to the `AlleleOrigin` indicating the parental origin.
/// * `child_seq` - A slice of the child's allele sequence.
///
/// # Returns
///
/// The `DenovoType` representing the type of de novo event.
fn get_denovo_type(
    mother_gts: &AlleleSet,
    father_gts: &AlleleSet,
    allele_origin: &AlleleOrigin,
    child_seq: &[u8],
) -> DenovoType {
    match allele_origin {
        AlleleOrigin::Father { allele } => {
            compare_seq_lengths(father_gts, allele, child_seq).unwrap_or(DenovoType::Unclear)
        }
        AlleleOrigin::Mother { allele } => {
            compare_seq_lengths(mother_gts, allele, child_seq).unwrap_or(DenovoType::Unclear)
        }
        _ => DenovoType::Unclear,
    }
}

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
fn align_alleleset(
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

fn align_allele(
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

/// Calculates the count and mean difference of child alignments that exceed the top parental score.
///
/// Determines the number of child alignments with scores higher than the top score
/// of a parent and calculates the mean difference of these scores from the top parental score.
///
/// # Arguments
///
/// * `top_parent_score`: The highest alignment score of a parent.
/// * `child_aligns`: A slice of vectors containing alignment scores for the child's alleles.
///
/// # Returns
///
/// Returns a tuple containing the count of alignments exceeding the top parental score and the mean
/// difference of these alignments from the top score.
fn get_parental_count_diff(top_parent_score: f64, child_aligns: &[i32]) -> (usize, f32) {
    let (count, sum) = child_aligns
        .iter()
        .map(|&a| a as f64)
        .filter(|a| a > &top_parent_score)
        .fold((0, 0.0), |(count, sum), a| (count + 1, sum + a));
    let mean_diff = if count > 0 {
        (sum / count as f64 - top_parent_score).abs() as f32
    } else {
        0.0
    };
    (count, mean_diff)
}

/// Determines the top alignment score among parental alleles at a given quantile.
///
/// Finds the highest alignment score at a specified quantile across all parental
/// alleles. It is used to establish a threshold for comparing child alignment scores.
///
/// # Arguments
///
/// * `align_scores`: A slice of vectors containing alignment scores for parental alleles.
/// * `quantile`: The quantile used to determine the top score.
///
/// # Returns
///
/// Returns an `Option<f64>` representing the top alignment score at the given quantile, if scores are available.
fn get_top_parent_score(align_scores: &[Vec<i32>], quantile: f64) -> Option<f64> {
    align_scores
        .iter()
        .filter_map(|scores| {
            let mut scores_f64: Vec<f64> = scores.iter().map(|&x| x as f64).collect();
            let quantile_score = math::quantile(&mut scores_f64, quantile);
            log::trace!(
                "Quantile score for {:?} at quantile {}: {}",
                scores_f64,
                quantile,
                quantile_score.unwrap_or_default()
            );
            quantile_score
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

/// Calculates the de novo coverage and mean difference from parental alignment scores.
///
/// Determines the number of child alignments that exceed the top parental alignment
/// scores and calculates the mean difference in scores.
///
/// # Arguments
///
/// * `mother_align_scores` - A slice of alignment scores for the mother's alleles.
/// * `father_align_scores` - A slice of alignment scores for the father's alleles.
/// * `child_align_scores` - A slice of alignment scores for the child's alleles.
/// * `parent_quantile` - The quantile used for parental allele frequency calculations.
///
/// # Returns
///
/// A tuple containing the de novo coverage, mean difference from the father's scores, and mean
/// difference from the mother's scores.
fn get_denovo_coverage(
    mother_align_scores: &[Vec<i32>],
    father_align_scores: &[Vec<i32>],
    child_align_scores: &[i32],
    parent_quantile: f64,
) -> (usize, f32, f32) {
    let top_mother_score = get_top_parent_score(mother_align_scores, parent_quantile).unwrap();
    let top_father_score = get_top_parent_score(father_align_scores, parent_quantile).unwrap();

    let (mother_count, mean_diff_mother) =
        get_parental_count_diff(top_mother_score, child_align_scores);
    let (father_count, mean_diff_father) =
        get_parental_count_diff(top_father_score, child_align_scores);

    let top_score = if top_mother_score >= top_father_score {
        mother_count
    } else {
        father_count
    };

    (top_score, mean_diff_father, mean_diff_mother)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::read::ReadInfoBuilder;

    #[test]
    fn test_parent_origin_matrix_1() {
        let mut matrix = Array2::from_elem((4, 2), f64::MIN);
        matrix[[0, 0]] = -25.0;
        matrix[[0, 1]] = -50.0;
        matrix[[1, 0]] = -40.0;
        matrix[[1, 1]] = -30.0;
        matrix[[2, 0]] = 0.0;
        matrix[[2, 1]] = -50.0;
        matrix[[3, 0]] = -10.0;
        matrix[[3, 1]] = -20.0;

        let mut combs_score: Vec<(f64, [(usize, usize); 2])> = Vec::new();
        for c in COMBS_2D.iter() {
            combs_score.push((c.iter().map(|&(i, j)| matrix[[i, j]]).sum::<f64>(), *c));
        }
        combs_score.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

        let expected_result = (-30.0, [(2, 0), (1, 1)]);
        assert_eq!(combs_score[0], expected_result);

        assert_eq!(
            AlleleOrigin::new(combs_score[0].1[0].0).unwrap(),
            AlleleOrigin::Father {
                allele: AlleleNum::One,
            },
        );

        assert_eq!(
            AlleleOrigin::new(combs_score[0].1[1].0).unwrap(),
            AlleleOrigin::Mother {
                allele: AlleleNum::Two,
            },
        );
    }

    #[test]
    fn test_get_parental_count_diff() {
        let top_parent_score = -8.0;
        let child_aligns = vec![vec![-30, -20, -30], vec![-6, 3, 0]];
        assert_eq!(
            get_parental_count_diff(top_parent_score, &child_aligns[1]),
            (3, 7.0)
        );
    }

    #[test]
    fn test_get_denovo_coverage() {
        let mother_aligns = vec![vec![-24, -24, -24], vec![-24, -12, -24]];
        let father_aligns = vec![vec![-8, -8, -12], vec![-8, -8, -20]];
        let child_aligns = vec![vec![-30, -20, -30], vec![-6, 0, 0]];
        assert_eq!(
            get_denovo_coverage(&mother_aligns, &father_aligns, &child_aligns[1], 1.0),
            (3, 6.0, 10.0)
        );
    }

    #[test]
    fn test_get_denovo_coverage_single_allele() {
        let mother_aligns = vec![vec![-24, -24, -24], vec![-24, -12, -24]];
        let father_aligns = vec![vec![-8, -8, -12], vec![-8, -8, -20]];
        let child_aligns = vec![vec![-30, 0, 0], vec![-6, 0, 0]];
        assert_eq!(
            get_denovo_coverage(&mother_aligns, &father_aligns, &child_aligns[1], 1.0),
            (3, 6.0, 10.0)
        );
    }

    #[test]
    fn test_compare_seq_lengths_substitution() {
        let gts = vec![
            Allele {
                seq: b"ATATATAT".to_vec(),
                genotype: 0,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 0,
            },
            Allele {
                seq: b"ATATATAT".to_vec(),
                genotype: 1,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 1,
            },
        ];
        let gts = AlleleSet {
            alleles: gts,
            hp_counts: [0; 3],
        };
        let allele = AlleleNum::One;
        let child_seq = b"ATCGATAT".to_vec();
        assert_eq!(
            compare_seq_lengths(&gts, &allele, &child_seq).unwrap(),
            DenovoType::Substitution
        );
    }

    #[test]
    fn test_compare_seq_lengths_contraction() {
        let gts = vec![
            Allele {
                seq: b"ATATATAT".to_vec(),
                genotype: 0,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 0,
            },
            Allele {
                seq: b"ATA".to_vec(),
                genotype: 1,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 1,
            },
        ];
        let gts = AlleleSet {
            alleles: gts,
            hp_counts: [0; 3],
        };
        let allele = AlleleNum::One;
        let child_seq = b"ATATAT".to_vec();
        assert_eq!(
            compare_seq_lengths(&gts, &allele, &child_seq).unwrap(),
            DenovoType::Contraction
        );
    }

    #[test]
    fn test_compare_seq_lengths_expansion() {
        let gts = vec![
            Allele {
                seq: b"ATATA".to_vec(),
                genotype: 0,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 0,
            },
            Allele {
                seq: b"ATATATAT".to_vec(),
                genotype: 1,
                read_aligns: vec![
                    (ReadInfoBuilder::default().build().unwrap(), 0),
                    (ReadInfoBuilder::default().build().unwrap(), 1),
                ],
                motif_count: "".to_string(),
                index: 1,
            },
        ];
        let gts = AlleleSet {
            alleles: gts,
            hp_counts: [0; 3],
        };
        let allele = AlleleNum::Two;
        let child_seq = b"ATATATATAT".to_vec();
        assert_eq!(
            compare_seq_lengths(&gts, &allele, &child_seq).unwrap(),
            DenovoType::Expansion
        );
    }
}
