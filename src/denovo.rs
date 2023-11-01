use crate::aligner::WFAligner;
use crate::allele::Allele;
use ndarray::{Array2, ArrayBase, Dim, OwnedRepr};
use serde::Serialize;
use std::{cmp::Ordering, fmt};

#[derive(Debug)]
pub struct DenovoAllele {
    pub genotype: usize,
    pub denovo_score: usize,
    pub child_coverage: usize,
    pub allele_coverage: usize,
    pub mean_diff_father: f32,
    pub mean_diff_mother: f32,
    pub allele_origin: AlleleOrigin,
    pub denovo_status: DenovoStatus,
    pub motif_count: String,
    pub index: usize,
}

#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoType {
    Expansion,
    Contraction,
    Substitution,
    Unclear,
}

#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoStatus {
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

#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum AlleleNum {
    One,
    Two,
    Unclear,
}

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

pub fn assess_denovo<'a>(
    mother_gts: &'a [Allele],
    father_gts: &'a [Allele],
    child_gts: &'a [Allele],
    clip_len: usize,
    parent_quantile: f64,
    aligner: &mut WFAligner,
) -> impl Iterator<Item = DenovoAllele> + 'a {
    let mut matrix = Array2::from_elem((4, 2), f64::MIN);
    let mut dnrs = Vec::with_capacity(child_gts.len());

    for (index, denovo_allele) in child_gts.iter().enumerate() {
        let mother_align_scores = align(mother_gts, &denovo_allele.seq, clip_len, aligner);
        let father_align_scores = align(father_gts, &denovo_allele.seq, clip_len, aligner);
        let child_align_scores = align(child_gts, &denovo_allele.seq, clip_len, aligner);

        let (denovo_coverage, mean_diff_father, mean_diff_mother) = get_denovo_coverage(
            &mother_align_scores,
            &father_align_scores,
            &child_align_scores,
            parent_quantile,
        );

        update_mean_matrix(
            &mut matrix,
            index,
            &mother_align_scores,
            &father_align_scores,
        );

        dnrs.push(DenovoAllele {
            genotype: denovo_allele.genotype,
            denovo_score: denovo_coverage,
            child_coverage: child_align_scores.iter().map(|vec| vec.len()).sum(),
            allele_coverage: denovo_allele.read_aligns.len(),
            mean_diff_mother,
            mean_diff_father,
            allele_origin: AlleleOrigin::Unclear,
            denovo_status: DenovoStatus::NotDenovo,
            motif_count: denovo_allele.motif_count.clone(),
            index: denovo_allele.index,
        });
    }

    // Get best allele combinations
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
    combs_score.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    // Update allele origin and de novo type
    let comb_diff = (combs_score[0].0 - combs_score[1].0).abs();
    let comb_0 = combs_score[0].1;
    for (index, denovo_allele) in child_gts.iter().enumerate() {
        let dna = &mut dnrs[index];
        // TODO: Alignment based: while we might not know the allele origin, may still figure out parental origin
        if comb_diff < 1.0 {
            dna.allele_origin = AlleleOrigin::Unclear;
        } else {
            dna.allele_origin = AlleleOrigin::new(comb_0[index].0).unwrap();
        }

        dna.denovo_status = match dna.denovo_score {
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

pub fn update_mean_matrix(
    matrix: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
    index: usize,
    mother_align_scores: &[Vec<i32>],
    father_align_scores: &[Vec<i32>],
) {
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

fn compare_seq_lengths(gts: &[Allele], allele: &AlleleNum, child_seq: &[u8]) -> Option<DenovoType> {
    match allele {
        AlleleNum::One => compare_sequences(&gts[0].seq, child_seq),
        AlleleNum::Two => compare_sequences(&gts[1].seq, child_seq),
        _ => None,
    }
}

fn compare_sequences(seq: &[u8], child_seq: &[u8]) -> Option<DenovoType> {
    match seq.len().cmp(&child_seq.len()) {
        Ordering::Less => Some(DenovoType::Expansion),
        Ordering::Greater => Some(DenovoType::Contraction),
        Ordering::Equal => Some(DenovoType::Substitution),
    }
}

fn get_denovo_type(
    mother_gts: &[Allele],
    father_gts: &[Allele],
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

fn align(gts: &[Allele], target: &[u8], clip_len: usize, aligner: &mut WFAligner) -> Vec<Vec<i32>> {
    let mut align_scores = vec![vec![]; gts.len()];
    for (i, allele) in gts.iter().enumerate() {
        for (read, _align) in &allele.read_aligns {
            let _status = aligner.align_end_to_end(&read.bases, target);
            align_scores[i].push(aligner.cigar_score_clipped(clip_len));
        }
    }
    align_scores
}

fn get_parental_count_diff(top_parent_score: f64, child_aligns: &[Vec<i32>]) -> (usize, f32) {
    let (count, sum) = child_aligns
        .iter()
        .flatten()
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

fn get_score_quantile(xs: &mut [f64], q: f64) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let index = q * (xs.len() as f64 - 1.0);
    let lower_index = index.floor() as usize;
    let upper_index = lower_index + 1;
    let fraction = index - lower_index as f64;

    if upper_index >= xs.len() {
        return Some(xs[xs.len() - 1]);
    }

    Some(xs[lower_index] + (xs[upper_index] - xs[lower_index]) * fraction)
}

fn get_top_parent_score(align_scores: &[Vec<i32>], quantile: f64) -> Option<f64> {
    align_scores
        .iter()
        .filter_map(|scores| {
            let mut scores_f64: Vec<f64> = scores.iter().map(|&x| x as f64).collect();
            get_score_quantile(&mut scores_f64, quantile)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

fn get_denovo_coverage(
    mother_align_scores: &[Vec<i32>],
    father_align_scores: &[Vec<i32>],
    child_align_scores: &[Vec<i32>],
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
    use crate::read::ReadInfoBuilder;

    use super::*;

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
            get_parental_count_diff(top_parent_score, &child_aligns),
            (3, 7.0)
        );
    }

    #[test]
    fn test_get_denovo_coverage() {
        let mother_aligns = vec![vec![-24, -24, -24], vec![-24, -12, -24]];
        let father_aligns = vec![vec![-8, -8, -12], vec![-8, -8, -20]];
        let child_aligns = vec![vec![-30, -20, -30], vec![-6, 0, 0]];
        assert_eq!(
            get_denovo_coverage(&mother_aligns, &father_aligns, &child_aligns, 1.0),
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
        let allele = AlleleNum::Two;
        let child_seq = b"ATATATATAT".to_vec();
        assert_eq!(
            compare_seq_lengths(&gts, &allele, &child_seq).unwrap(),
            DenovoType::Expansion
        );
    }
}
