use crate::{allele::Allele, read::ReadInfo};
use ndarray::{s, Array2, ArrayView1, ArrayView2};
use std::collections::{HashMap, HashSet};
use std::iter;

#[cfg(not(test))]
mod constants {
    pub const MAX_SNP_DIFF: i32 = 6000;
    pub const MIN_SNP_FREQ: f64 = 0.2;
}

#[cfg(test)]
mod constants {
    pub const MAX_SNP_DIFF: i32 = 4000;
    pub const MIN_SNP_FREQ: f64 = 0.2;
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FamilyMember {
    Father,
    Mother,
    Child,
}

impl FamilyMember {
    fn index(&self) -> usize {
        match self {
            FamilyMember::Father => 0,
            FamilyMember::Mother => 1,
            FamilyMember::Child => 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct FamilyOffsets {
    pub start_offset: usize,
    pub end_offset: usize,
    pub allele_offsets: Vec<usize>,
}

#[derive(Debug, PartialEq)]
pub struct TrinaryMatrix {
    pub matrix: Array2<u8>,
    pub offsets: [FamilyOffsets; 3],
    pub mismatch_offsets: Vec<i32>,
}

impl TrinaryMatrix {
    pub fn new(
        child: &Vec<Allele>,
        father: &Vec<Allele>,
        mother: &Vec<Allele>,
    ) -> Option<TrinaryMatrix> {
        let mut n_reads = 0;
        let mut all_pois: HashSet<i32> = HashSet::new();
        let family_members = [father, mother, child];

        for family_member in &family_members {
            for allele in *family_member {
                n_reads += allele.read_aligns.len();
                for (read, _) in &allele.read_aligns {
                    all_pois.extend(read.mismatch_offsets.as_ref().unwrap());
                }
            }
        }

        if all_pois.is_empty() {
            return None;
        }

        let mut all_pois: Vec<i32> = all_pois.into_iter().collect();
        all_pois.sort_unstable();

        let mut matrix = Array2::from_elem((n_reads, all_pois.len()), 2);

        // TODO: ideally want: offsets = [FamilyOffsets::default(); 3], can't because of vec, maybe box Vec?
        let mut offsets = [
            (FamilyOffsets {
                start_offset: 0,
                end_offset: 0,
                allele_offsets: Vec::new(),
            }),
            (FamilyOffsets {
                start_offset: 0,
                end_offset: 0,
                allele_offsets: Vec::new(),
            }),
            (FamilyOffsets {
                start_offset: 0,
                end_offset: 0,
                allele_offsets: Vec::new(),
            }),
        ];

        let mut row_idx = 0;
        let mut allele_offsets = Vec::new();
        for (member_idx, family_member) in family_members.iter().enumerate() {
            let start_offset = row_idx;
            for (allele_idx, allele) in family_member.iter().enumerate() {
                if allele_idx > 0 {
                    allele_offsets.push(row_idx);
                }
                for (read, _) in &allele.read_aligns {
                    let mismatch_offsets = read.mismatch_offsets.as_ref().unwrap();
                    let start_col_idx = all_pois
                        .binary_search(&read.start_offset.unwrap())
                        .unwrap_or_else(|x| x);
                    let end_col_idx = all_pois
                        .binary_search(&read.end_offset.unwrap())
                        .unwrap_or_else(|x| x);
                    for col_idx in start_col_idx..end_col_idx {
                        let offset = all_pois[col_idx];
                        matrix[[row_idx, col_idx]] =
                            mismatch_offsets.binary_search(&offset).is_ok() as u8;
                    }
                    row_idx += 1;
                }
            }
            offsets[member_idx] = FamilyOffsets {
                start_offset,
                end_offset: row_idx - 1,
                allele_offsets: allele_offsets.clone(),
            };
            allele_offsets.clear();
        }

        Some(TrinaryMatrix {
            matrix,
            offsets,
            mismatch_offsets: all_pois,
        })
    }

    pub fn family_submatrix(&self, member: FamilyMember) -> ArrayView2<u8> {
        let offset = &self.offsets[member.index()];
        self.matrix
            .slice(s![offset.start_offset..=offset.end_offset, ..])
    }

    pub fn allele_submatrix(
        &self,
        member: FamilyMember,
        allele_idx: usize,
    ) -> Option<ArrayView2<u8>> {
        let offset = &self.offsets[member.index()];

        if allele_idx > offset.allele_offsets.len() {
            return None;
        }

        let start_offset = if allele_idx == 0 {
            offset.start_offset
        } else {
            offset.allele_offsets[allele_idx - 1]
        };

        let end_offset = if allele_idx < offset.allele_offsets.len() {
            offset.allele_offsets[allele_idx] - 1
        } else {
            offset.end_offset
        };

        Some(self.matrix.slice(s![start_offset..=end_offset, ..]))
    }

    pub fn iter_alleles<'a>(
        &'a self,
        member: FamilyMember,
    ) -> impl Iterator<Item = ArrayView2<u8>> + 'a {
        let offset = &self.offsets[member.index()];
        let start_offsets =
            iter::once(offset.start_offset).chain(offset.allele_offsets.iter().cloned());
        let end_offsets = offset
            .allele_offsets
            .iter()
            .cloned()
            .chain(iter::once(offset.end_offset + 1));

        start_offsets
            .zip(end_offsets)
            .map(move |(start, end)| self.matrix.slice(s![start..end, ..]))
    }

    pub fn count_alleles(&self, member: FamilyMember) -> usize {
        let offset = &self.offsets[member.index()];
        offset.allele_offsets.len() + 1
    }
}

pub trait ReadFilter {
    fn filter(&self, reads: &mut Vec<ReadInfo>);
}

pub struct FilterByDist;
impl ReadFilter for FilterByDist {
    fn filter(&self, reads: &mut Vec<ReadInfo>) {
        for read in reads.iter_mut() {
            if let Some(offsets) = &mut read.mismatch_offsets {
                offsets.drain(
                    ..offsets
                        .binary_search(&(-constants::MAX_SNP_DIFF - 1))
                        .unwrap_or_else(|x| x),
                );
                offsets.drain(
                    offsets
                        .binary_search(&(constants::MAX_SNP_DIFF + 1))
                        .unwrap_or_else(|x| x)..,
                );
            }
        }
    }
}

fn calc_similarity(child_row: &ArrayView1<u8>, parent_row: &ArrayView1<u8>) -> f64 {
    let mut matching_positions = 0;
    let mut total_covered_positions = 0;

    for (&child_value, &parent_value) in child_row.iter().zip(parent_row.iter()) {
        if child_value == 2 || parent_value == 2 {
            continue;
        }

        total_covered_positions += 1;

        if child_value == parent_value {
            matching_positions += 1;
        }
    }

    let total_covered_positions = total_covered_positions as f64;
    let total_covered_positions = if total_covered_positions == 0.0 {
        0.00000001
    } else {
        total_covered_positions
    };

    matching_positions as f64 / total_covered_positions
}

fn get_likelihood(child_allele: ArrayView2<u8>, parent_allele: ArrayView2<u8>) -> f64 {
    let similarity_scores: Vec<f64> = child_allele
        .outer_iter()
        .flat_map(|child_row| {
            parent_allele
                .outer_iter()
                .map(move |parent_row| calc_similarity(&child_row, &parent_row))
        })
        .collect();

    let likelihood: f64 = similarity_scores.iter().sum::<f64>() / similarity_scores.len() as f64;
    likelihood
}

fn get_individual_prob(
    trinary_matrix: &TrinaryMatrix,
    child_allele_idx: usize,
) -> Option<(Vec<f64>, f64)> {
    let child_allele = trinary_matrix.allele_submatrix(FamilyMember::Child, child_allele_idx)?;
    if child_allele.is_empty() {
        return None;
    }

    let father_count = trinary_matrix.count_alleles(FamilyMember::Father) as f64;
    let mother_count = trinary_matrix.count_alleles(FamilyMember::Mother) as f64;

    let mut hypotheses = Vec::new();
    let mut priors = Vec::new();

    for father_allele in trinary_matrix.iter_alleles(FamilyMember::Father) {
        if father_allele.is_empty() {
            return None;
        }
        hypotheses.push((child_allele, father_allele));
        priors.push(0.5 / father_count);
    }

    for mother_allele in trinary_matrix.iter_alleles(FamilyMember::Mother) {
        if mother_allele.is_empty() {
            return None;
        }
        hypotheses.push((child_allele, mother_allele));
        priors.push(0.5 / mother_count);
    }

    let likelihoods: Vec<f64> = hypotheses
        .iter()
        .map(|(c, p)| get_likelihood(*c, *p))
        .collect();

    let max_likelihood = likelihoods
        .iter()
        .cloned()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    let all_probabilities: f64 = likelihoods.iter().zip(&priors).map(|(l, p)| l * p).sum();
    let posteriors: Vec<f64> = likelihoods
        .iter()
        .zip(&priors)
        .map(|(l, p)| l * p / all_probabilities)
        .collect();

    Some((posteriors, max_likelihood))
}

fn combine_probabilities(
    child_allele1_prob: &[f64],
    child_allele2_prob: &[f64],
    father_count: usize,
    mother_count: usize,
) -> HashMap<String, f64> {
    let mut combined_probabilities = HashMap::new();
    let mut total_prob = 0.0;

    for i in 0..(father_count + mother_count) {
        for j in 0..(father_count + mother_count) {
            if (i < father_count && j < father_count) || (i >= father_count && j >= father_count) {
                continue;
            }

            let probability = child_allele1_prob[i] * child_allele2_prob[j];
            combined_probabilities.insert(format!("H{}_{}", i, j), probability);
            total_prob += probability;
        }
    }

    // Normalize probabilities
    for value in combined_probabilities.values_mut() {
        *value /= total_prob;
    }

    combined_probabilities
}

pub fn inheritance_prob(
    trinary_matrix: &TrinaryMatrix,
) -> Option<(HashMap<String, f64>, (f64, f64))> {
    match trinary_matrix.count_alleles(FamilyMember::Child) {
        1 => {
            let (child_allele1_prob, max_likelihood) = get_individual_prob(trinary_matrix, 0)?;
            let mut result: HashMap<String, f64> = HashMap::new();
            for (i, prob) in child_allele1_prob.iter().enumerate() {
                result.insert(i.to_string(), *prob);
            }
            Some((result, ((max_likelihood, -1.0))))
        }
        2 => {
            let (child_allele1_prob, max_likelihood_1) = get_individual_prob(trinary_matrix, 0)?;
            let (child_allele2_prob, max_likelihood_2) = get_individual_prob(trinary_matrix, 1)?;
            Some((
                combine_probabilities(
                    &child_allele1_prob,
                    &child_allele2_prob,
                    trinary_matrix.count_alleles(FamilyMember::Father),
                    trinary_matrix.count_alleles(FamilyMember::Mother),
                ),
                (max_likelihood_1, max_likelihood_2),
            ))
        }
        _ => None,
    }
}

pub struct FilterByFreq;
impl ReadFilter for FilterByFreq {
    fn filter(&self, reads: &mut Vec<ReadInfo>) {
        let mut offset_counts: HashMap<i32, usize> = HashMap::new();
        let total_reads = reads.len();

        for read in reads.iter() {
            if let Some(offsets) = &read.mismatch_offsets {
                for &offset in offsets {
                    *offset_counts.entry(offset).or_insert(0) += 1;
                }
            }
        }
        for read in reads.iter_mut() {
            if let Some(offsets) = &mut read.mismatch_offsets {
                offsets.retain(|&offset| {
                    let count = offset_counts.get(&offset).unwrap_or(&0);
                    (*count as f64) / (total_reads as f64) >= constants::MIN_SNP_FREQ
                });
            }
        }
    }
}

pub fn apply_read_filters(reads: &mut Vec<ReadInfo>, filters: &[&(dyn ReadFilter + Sync)]) {
    for filter in filters {
        filter.filter(reads);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::ReadInfoBuilder;
    use itertools::Itertools;

    fn create_read_info(
        mismatch_offsets: Vec<i32>,
        start_offset: i32,
        end_offset: i32,
    ) -> ReadInfo {
        ReadInfoBuilder::new()
            .with_mismatch_offsets(Some(mismatch_offsets))
            .with_start_offset(Some(start_offset))
            .with_end_offset(Some(end_offset))
            .build()
            .unwrap()
    }

    #[test]
    fn test_empty_trinary_matrix() {
        let child_alleles: Vec<Allele> = vec![];
        let father_alleles: Vec<Allele> = vec![];
        let mother_alleles: Vec<Allele> = vec![];
        let trinaray_mat: Option<TrinaryMatrix> =
            TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles);
        assert_eq!(trinaray_mat, None);
    }

    // TODO: split into multiple tests
    #[test]
    fn test_trinary_matrix() {
        let father_alleles = vec![
            Allele::dummy(vec![
                create_read_info(vec![-8, -2, 3, 5, 8, 42], -12, 70),
                create_read_info(vec![-2, 3, 5, 8], -6, 30),
                create_read_info(vec![3, 5, 8, 42], 0, 50),
            ]),
            Allele::dummy(vec![
                create_read_info(vec![-8, -4, 2], -22, 30),
                create_read_info(vec![-8, -4, 2], -40, 60),
                create_read_info(vec![-8, -4, 2], -30, 90),
            ]),
        ];

        let mother_alleles = vec![
            Allele::dummy(vec![
                create_read_info(vec![-4, -1, 2, 6], -5, 7),
                create_read_info(vec![-4, -1, 2, 6, 23], -20, 70),
                create_read_info(vec![2, 6, 23], 1, 120),
            ]),
            // Allele::dummy(vec![
            //     create_read_info(vec![6, 23], -30, 40),
            //     create_read_info(vec![6, 23], -10, 90),
            //     create_read_info(vec![6, 23], -2, 120),
            // ]),
        ];

        let child_alleles = vec![
            Allele::dummy(vec![
                create_read_info(vec![-8, -2, 3, 5, 8, 42], -12, 70),
                create_read_info(vec![-8, -2, 3, 5, 8], -10, 30),
                create_read_info(vec![3, 5, 8, 42], 0, 50),
            ]),
            Allele::dummy(vec![
                create_read_info(vec![-4, -1, 2, 6, 23], -40, 90),
                create_read_info(vec![-4, -1, 2, 6], -32, 20),
                create_read_info(vec![2, 6, 23], 1, 120),
            ]),
            Allele::dummy(vec![create_read_info(vec![], -40, 90)]),
            Allele::dummy(vec![
                create_read_info(vec![-1, 2, 6], -32, 20),
                create_read_info(vec![2, 6], 1, 10),
            ]),
        ];

        let trinaray_mat =
            TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles).unwrap();

        // Test dimensionality
        assert_eq!(trinaray_mat.matrix.shape(), &[18, 11]);

        assert_eq!(
            get_likelihood(
                trinaray_mat
                    .allele_submatrix(FamilyMember::Child, 0)
                    .unwrap(),
                trinaray_mat
                    .allele_submatrix(FamilyMember::Child, 0)
                    .unwrap()
            ),
            1.0
        );

        // Counts alleles
        assert_eq!(trinaray_mat.count_alleles(FamilyMember::Child), 4);
        assert_eq!(trinaray_mat.count_alleles(FamilyMember::Mother), 1);
        assert_eq!(trinaray_mat.count_alleles(FamilyMember::Father), 2);

        // Check for correctness of mismatch offset labels
        assert_eq!(
            trinaray_mat.mismatch_offsets,
            vec![-8, -4, -2, -1, 2, 3, 5, 6, 8, 23, 42,]
        );

        // Row 0 should correspond to the first read of the father
        assert_eq!(
            trinaray_mat.matrix.row(0),
            trinaray_mat
                .matrix
                .row(trinaray_mat.offsets[FamilyMember::Father.index()].start_offset),
        );

        // Check for correctness of matrix row 0
        assert_eq!(
            trinaray_mat.matrix.row(0).to_vec(),
            vec![1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
        );

        // Check for correctness for the first read of the mother
        assert_eq!(
            trinaray_mat
                .matrix
                .row(trinaray_mat.offsets[FamilyMember::Mother.index()].start_offset)
                .to_vec(),
            vec![2, 1, 0, 1, 1, 0, 0, 1, 2, 2, 2]
        );

        // Mother has only one allele so the mother submatrix should be equivalent to allele submatrix 0
        assert_eq!(
            trinaray_mat.family_submatrix(FamilyMember::Mother),
            trinaray_mat
                .allele_submatrix(FamilyMember::Mother, 0)
                .unwrap()
        );

        // Allele submatrix 1 (and beyond) of the mother should be None
        assert_eq!(None, trinaray_mat.allele_submatrix(FamilyMember::Mother, 1));

        // Child has 4 alleles, check that allele 3 is correct
        assert_eq!(
            trinaray_mat
                .allele_submatrix(FamilyMember::Child, 2)
                .unwrap()
                .rows()
                .into_iter()
                .map(|row| row.to_vec())
                .collect_vec(),
            vec![[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
        );
    }

    #[test]
    fn test_filter_by_dist() {
        let mut reads = vec![
            create_read_info(
                vec![
                    -5000, -4000, -3000, -2000, -1000, 1000, 2000, 3000, 4000, 5000,
                ],
                0,
                0,
            ),
            create_read_info(vec![4001, 6000, 7000, 8000, 9000, 10000], 0, 0),
            create_read_info(vec![0, 1000, 2000, 3000, 3999, 4000, 4001], 0, 0),
        ];

        FilterByDist.filter(&mut reads);

        assert_eq!(
            reads[0].mismatch_offsets.as_ref().unwrap(),
            &vec![-4000, -3000, -2000, -1000, 1000, 2000, 3000, 4000]
        );
        assert_eq!(reads[1].mismatch_offsets.as_ref().unwrap(), &vec![]);
        assert_eq!(
            reads[2].mismatch_offsets.as_ref().unwrap(),
            &vec![0, 1000, 2000, 3000, 3999, 4000]
        );
    }

    #[test]
    fn test_filter_by_freq() {
        let mut reads = vec![
            create_read_info(vec![1000, 2000, 3000, 3002, 4000, 5000], 0, 0),
            create_read_info(vec![1000, 2000, 3000, 4000, 4090, 5000], 0, 0),
            create_read_info(vec![-3213, 1000, 2000, 3000, 4000, 5000], 0, 0),
            create_read_info(vec![1, 2, 3, 4, 5], 0, 0),
            create_read_info(vec![-3, -2, -1, 0, 2000], 0, 0),
            create_read_info(vec![], 0, 0),
        ];

        FilterByFreq.filter(&mut reads);

        assert_eq!(
            reads[0].mismatch_offsets.as_ref().unwrap(),
            &vec![1000, 2000, 3000, 4000, 5000]
        );
        assert_eq!(
            reads[1].mismatch_offsets.as_ref().unwrap(),
            &vec![1000, 2000, 3000, 4000, 5000]
        );
        assert_eq!(
            reads[2].mismatch_offsets.as_ref().unwrap(),
            &vec![1000, 2000, 3000, 4000, 5000]
        );
        assert_eq!(reads[3].mismatch_offsets.as_ref().unwrap(), &vec![]);
        assert_eq!(reads[4].mismatch_offsets.as_ref().unwrap(), &vec![2000]);
    }

    #[test]
    fn test_filter_by_dist_and_freq() {
        let mut reads = vec![
            create_read_info(
                vec![
                    -5000, -4000, -3000, -2000, -1000, 1000, 2000, 3000, 3002, 4000, 5000,
                ],
                0,
                0,
            ),
            create_read_info(
                vec![
                    -5000, -4000, -3000, -2000, -1000, 1000, 2000, 3000, 4000, 4090, 5000,
                ],
                0,
                0,
            ),
            create_read_info(vec![-3213, 1000, 2000, 3000, 4000, 5000], 0, 0),
            create_read_info(vec![1, 2, 3, 4, 5], 0, 0),
            create_read_info(vec![-3, -2, -1, 0, 2000], 0, 0),
            create_read_info(vec![], 0, 0),
        ];

        FilterByDist.filter(&mut reads);
        FilterByFreq.filter(&mut reads);

        assert_eq!(
            reads[0].mismatch_offsets.as_ref().unwrap(),
            &vec![-4000, -3000, -2000, -1000, 1000, 2000, 3000, 4000]
        );
        assert_eq!(
            reads[1].mismatch_offsets.as_ref().unwrap(),
            &vec![-4000, -3000, -2000, -1000, 1000, 2000, 3000, 4000]
        );
        assert_eq!(
            reads[2].mismatch_offsets.as_ref().unwrap(),
            &vec![1000, 2000, 3000, 4000]
        );
        assert_eq!(reads[3].mismatch_offsets.as_ref().unwrap(), &vec![]);
        assert_eq!(reads[4].mismatch_offsets.as_ref().unwrap(), &vec![2000]);
    }
}
