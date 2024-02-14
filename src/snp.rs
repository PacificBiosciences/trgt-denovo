//! SNP analysis module.
//!
//! This module provides functionality for analyzing SNPs within genomic data.
//! It includes structures and functions for representing family members, calculating
//! SNP similarities, and determining inheritance probabilities.

use crate::{allele::AlleleSet, read::ReadInfo};
use ndarray::{s, Array2, ArrayView1, ArrayView2};
use std::{
    collections::{HashMap, HashSet},
    iter,
};

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

/// Represents a family member in the context of SNP analysis.
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

/// Stores the offsets for SNP analysis for a family member.
///
/// This struct is used to keep track of the start and end offsets within the SNP matrix,
/// as well as the offsets for each allele.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FamilyOffsets {
    /// The start offset in the SNP matrix for this family member.
    pub start_offset: usize,
    /// The end offset in the SNP matrix for this family member.
    pub end_offset: usize,
    /// A vector of offsets for each allele of this family member.
    pub allele_offsets: Vec<usize>,
}

/// Represents a matrix used for SNP analysis, encoding the presence or absence of SNPs.
///
/// This matrix is used to compare SNP positions across family members and to calculate
/// inheritance probabilities. It uses a trinary encoding where `0` indicates absence,
/// `1` indicates presence, and `2` indicates missing information.
#[derive(Debug, PartialEq)]
pub struct TrinaryMatrix {
    /// The underlying matrix storing SNP information.
    pub matrix: Array2<u8>,
    /// Offsets for each family member within the matrix.
    pub offsets: [FamilyOffsets; 3],
    /// A vector of mismatch offsets used for SNP analysis.
    pub mismatch_offsets: Vec<i32>,
}

/// Provides methods for working with the `TrinaryMatrix`.
impl TrinaryMatrix {
    /// Creates a new `TrinaryMatrix` from the given alleles of a family trio.
    ///
    /// # Arguments
    ///
    /// * `child` - A slice of `Allele` instances for the child.
    /// * `father` - A slice of `Allele` instances for the father.
    /// * `mother` - A slice of `Allele` instances for the mother.
    ///
    /// # Returns
    ///
    /// An `Option<TrinaryMatrix>` which is `None` if no SNPs are found, or contains the constructed matrix otherwise.
    pub fn new(child: &AlleleSet, father: &AlleleSet, mother: &AlleleSet) -> Option<TrinaryMatrix> {
        let mut total_reads = 0;
        let mut pois_set: HashSet<i32> = HashSet::new();
        let family_members = [father, mother, child];

        for family_member in &family_members {
            for allele in *family_member {
                total_reads += allele.read_aligns.len();
                allele
                    .read_aligns
                    .iter()
                    .filter_map(|(read, _)| read.mismatch_offsets.as_ref())
                    .for_each(|mismatches| {
                        pois_set.extend(mismatches);
                    });
            }
        }

        if pois_set.is_empty() {
            return None;
        }

        let mut pois: Vec<i32> = pois_set.into_iter().collect();
        pois.sort_unstable();

        let mut matrix = Array2::from_elem((total_reads, pois.len()), 2);

        let mut offsets: [FamilyOffsets; 3] = std::iter::repeat_with(FamilyOffsets::default)
            .take(3)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| unreachable!());

        let mut current_row = 0;
        let mut allele_offsets = Vec::with_capacity(total_reads);
        for (member_idx, family_member) in family_members.iter().enumerate() {
            let start_offset = current_row;
            for (allele_idx, allele) in family_member.iter().enumerate() {
                if allele_idx > 0 {
                    allele_offsets.push(current_row);
                }
                for (read, _) in &allele.read_aligns {
                    let mismatch_offsets = read.mismatch_offsets.as_ref().unwrap();
                    let start_col = pois
                        .binary_search(&read.start_offset.unwrap())
                        .unwrap_or_else(|x| x);
                    let end_col = pois
                        .binary_search(&read.end_offset.unwrap())
                        .unwrap_or_else(|x| x);
                    for col in start_col..end_col {
                        let offset = pois[col];
                        matrix[[current_row, col]] =
                            mismatch_offsets.binary_search(&offset).is_ok() as u8;
                    }
                    current_row += 1;
                }
            }
            offsets[member_idx] = FamilyOffsets {
                start_offset,
                end_offset: current_row - 1,
                allele_offsets: allele_offsets.clone(),
            };
            allele_offsets.clear();
        }

        Some(TrinaryMatrix {
            matrix,
            offsets,
            mismatch_offsets: pois,
        })
    }

    /// Retrieves a submatrix view for a specific family member's SNPs.
    ///
    /// Returns a view of the SNP matrix that corresponds to the SNPs
    /// associated with the specified family member.
    ///
    /// # Arguments
    ///
    /// * `member` - The family member for which to retrieve the SNP submatrix.
    ///
    /// # Returns
    ///
    /// An `ArrayView2<u8>` representing the submatrix of SNPs for the given family member.
    pub fn family_submatrix(&self, member: FamilyMember) -> ArrayView2<u8> {
        let offset = &self.offsets[member.index()];
        self.matrix
            .slice(s![offset.start_offset..=offset.end_offset, ..])
    }

    /// Retrieves a submatrix view for a specific allele of a family member.
    ///
    /// Returns a view of the SNP matrix that corresponds to the SNPs
    /// associated with a particular allele of the specified family member.
    ///
    /// # Arguments
    ///
    /// * `member` - The family member to which the allele belongs.
    /// * `allele_idx` - The index of the allele within the family member's alleles.
    ///
    /// # Returns
    ///
    /// An `Option<ArrayView2<u8>>` representing the submatrix of SNPs for the given allele,
    /// or `None` if the allele index is out of bounds.
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

    /// Creates an iterator over submatrices for each allele of a family member.
    ///
    /// Returns an iterator that yields views of the SNP matrix for each allele
    /// of the specified family member. Each item in the iterator is a submatrix corresponding
    /// to a single allele's SNPs.
    ///
    /// # Arguments
    ///
    /// * `member` - The family member whose alleles are to be iterated over.
    ///
    /// # Returns
    ///
    /// An iterator over `ArrayView2<u8>` where each item is a submatrix for an allele of the family member.
    pub fn iter_alleles(&self, member: FamilyMember) -> impl Iterator<Item = ArrayView2<u8>> {
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

/// Calculates the similarity between a child's SNP row and a parent's SNP row.
///
/// # Arguments
///
/// * `child_row` - A view of the SNP row for the child.
/// * `parent_row` - A view of the SNP row for the parent.
///
/// # Returns
///
/// A `f64` representing the similarity score between the two rows.
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

/// Calculates the likelihood of inheritance based on the similarity of SNP rows between a child and a parent.
///
/// # Arguments
///
/// * `child_allele` - A view of the SNP matrix for the child's allele.
/// * `parent_allele` - A view of the SNP matrix for the parent's allele.
///
/// # Returns
///
/// A `f64` representing the likelihood of inheritance.
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

/// Determines the individual probability of inheritance for a given child allele.
///
/// # Arguments
///
/// * `trinary_matrix` - A reference to the `TrinaryMatrix`.
/// * `child_allele_idx` - The index of the child's allele.
///
/// # Returns
///
/// An `Option<(Vec<f64>, f64)>` containing the posterior probabilities and the maximum likelihood, or `None` if the allele is empty.
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

/// Combines the probabilities of two child alleles to determine the overall inheritance probabilities.
///
/// # Arguments
///
/// * `child_allele1_prob` - A slice of probabilities for the first child allele.
/// * `child_allele2_prob` - A slice of probabilities for the second child allele.
/// * `father_count` - The number of alleles for the father.
/// * `mother_count` - The number of alleles for the mother.
///
/// # Returns
///
/// A `HashMap<String, f64>` containing the combined probabilities for each inheritance hypothesis.
fn combine_probabilities(
    child_allele1_prob: &[f64],
    child_allele2_prob: &[f64],
    father_count: usize,
    mother_count: usize,
) -> HashMap<String, f64> {
    let mut combined_probabilities = HashMap::new();
    let mut total_prob = 0.0;
    for (i, c1_prob) in child_allele1_prob
        .iter()
        .enumerate()
        .take(father_count + mother_count)
    {
        for (j, c2_prob) in child_allele2_prob
            .iter()
            .enumerate()
            .take(father_count + mother_count)
        {
            if (i < father_count && j < father_count) || (i >= father_count && j >= father_count) {
                continue;
            }
            let probability = c1_prob * c2_prob;
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

/// Calculates the probabilities of inheritance for each possible combination of parental alleles.
///
/// # Arguments
///
/// * `trinary_matrix` - A reference to the `TrinaryMatrix`.
///
/// # Returns
///
/// An `Option<(HashMap<String, f64>, (f64, f64))>` containing the inheritance probabilities and the maximum likelihoods for each child allele, or `None` if there are more than two alleles.
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
            Some((result, (max_likelihood, -1.0)))
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

/// Defines a trait for filtering reads based on SNP criteria.
///
/// Implementors of this trait provide methods to filter reads to remove noise and improve SNP analysis accuracy.
pub trait ReadFilter {
    /// Filters a vector of `ReadInfo` instances based on specific criteria.
    ///
    /// # Arguments
    ///
    /// * `reads` - A mutable reference to a vector of `ReadInfo` instances to be filtered.
    fn filter(&self, reads: &mut Vec<ReadInfo>);
}

/// A `ReadFilter` that filters reads based on distance criteria.
///
/// This filter removes SNPs that are too far from the region of interest, as defined by a maximum distance threshold.
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

/// A `ReadFilter` that filters reads based on frequency criteria.
///
/// This filter retains SNPs that occur with a frequency above a specified threshold, ensuring that only common SNPs are considered.
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

/// Applies a set of read filters to a vector of `ReadInfo` instances.
///
/// # Arguments
///
/// * `reads` - A mutable reference to a vector of `ReadInfo` instances to be filtered.
/// * `filters` - A slice of references to objects that implement the `ReadFilter` trait.
pub fn apply_read_filters(reads: &mut Vec<ReadInfo>, filters: &[&(dyn ReadFilter + Sync)]) {
    for filter in filters {
        filter.filter(reads);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{allele::Allele, read::ReadInfoBuilder};
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
        let child_alleles = AlleleSet {
            alleles: vec![],
            hp_counts: [0; 3],
        };
        let father_alleles = AlleleSet {
            alleles: vec![],
            hp_counts: [0; 3],
        };
        let mother_alleles = AlleleSet {
            alleles: vec![],
            hp_counts: [0; 3],
        };
        let trinaray_mat: Option<TrinaryMatrix> =
            TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles);
        assert_eq!(trinaray_mat, None);
    }

    // TODO: split into multiple tests
    #[test]
    fn test_trinary_matrix() {
        let father_alleles = AlleleSet {
            alleles: vec![
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
            ],
            hp_counts: [0; 3],
        };

        let mother_alleles = AlleleSet {
            alleles: vec![
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
            ],
            hp_counts: [0; 3],
        };

        let child_alleles = AlleleSet {
            alleles: vec![
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
            ],
            hp_counts: [0; 3],
        };

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

    #[test]
    fn test_combine_probabilities() {
        let child_allele1_prob = vec![0.1, 0.2, 0.3, 0.4];
        let child_allele2_prob = vec![0.3, 0.2, 0.1, 0.4];
        let father_count = 2;
        let mother_count = 2;

        let combined_probs = combine_probabilities(
            &child_allele1_prob,
            &child_allele2_prob,
            father_count,
            mother_count,
        );

        let expected_probs = HashMap::from([
            ("H3_1".to_string(), 0.16),
            ("H2_0".to_string(), 0.18),
            ("H1_3".to_string(), 0.16),
            ("H0_3".to_string(), 0.08),
            ("H1_2".to_string(), 0.04),
            ("H3_0".to_string(), 0.24),
            ("H0_2".to_string(), 0.02),
            ("H2_1".to_string(), 0.12),
        ]);

        for (key, &expected_prob) in &expected_probs {
            let combined_prob = combined_probs.get(key).unwrap();
            assert!(
                (combined_prob - expected_prob).abs() < 1e-6,
                "Probabilities do not match for key {}: expected {}, got {}",
                key,
                expected_prob,
                combined_prob
            );
        }
    }
}
