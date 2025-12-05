//! SNP analysis module.
//!
//! This module provides functionality for analyzing SNPs within genomic data.
//! It includes structures and functions for representing members, calculating
//! SNP similarities, and determining inheritance probabilities.

use crate::{allele::AlleleSet, read::TrgtRead};
use ndarray::{s, Array2, ArrayView2};
use std::{
    collections::{HashMap, HashSet},
    fmt, iter,
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

/// Represents a member in the context of SNP analysis.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Member {
    Father,
    Mother,
    Child,
    Sample0,
    Sample1,
}

impl Member {
    fn index(&self) -> usize {
        match self {
            Member::Father => 0,
            Member::Mother => 1,
            Member::Child => 2,
            Member::Sample0 => 0,
            Member::Sample1 => 1,
        }
    }
}

/// Stores the offsets for SNP analysis for a member.
///
/// This struct is used to keep track of the start and end offsets within the SNP matrix,
/// as well as the offsets for each allele.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MemberOffsets {
    /// The start offset in the SNP matrix for this member.
    pub start_offset: usize,
    /// The end offset in the SNP matrix for this member.
    pub end_offset: usize,
    /// A vector of offsets for each allele of this member.
    pub allele_offsets: Vec<usize>,
}

/// Represents a matrix used for SNP analysis, encoding the presence or absence of SNPs.
///
/// This matrix is used to compare SNP positions across members and to calculate
/// inheritance probabilities. It uses a trinary encoding where `0` indicates absence,
/// `1` indicates presence, and `2` indicates missing information.
#[derive(Debug, PartialEq)]
pub struct TrinaryMatrix<const N: usize> {
    /// The underlying matrix storing SNP information.
    pub matrix: Array2<u8>,
    /// Offsets for each member within the matrix.
    pub offsets: [MemberOffsets; N],
    /// A vector of mismatch offsets used for SNP analysis.
    pub mismatch_offsets: Vec<i32>,
    /// The members in this matrix
    pub members: Vec<Member>,
}

impl<const N: usize> TrinaryMatrix<N> {
    /// Helper method to construct a TrinaryMatrix from a slice of AlleleSet members
    fn construct_matrix(members: &[&AlleleSet], member_types: Vec<Member>) -> Option<Self> {
        let mut total_reads = 0;
        let mut mismatch_positions_set: HashSet<i32> = HashSet::new();

        // First pass: count total reads and collect unique mismatch positions.
        for member in members {
            for allele in *member {
                total_reads += allele.read_aligns.len();
                allele
                    .read_aligns
                    .iter()
                    .filter_map(|(read, _)| read.mismatch_offsets.as_ref())
                    .for_each(|mismatches| {
                        mismatch_positions_set.extend(mismatches);
                    });
            }
        }

        if mismatch_positions_set.is_empty() {
            return None;
        }

        let mut mismatch_positions: Vec<i32> = mismatch_positions_set.into_iter().collect();
        mismatch_positions.sort_unstable();

        let mut matrix = Array2::from_elem((total_reads, mismatch_positions.len()), 2);

        let mut offsets: Vec<MemberOffsets> = std::iter::repeat_with(MemberOffsets::default)
            .take(N)
            .collect();

        let mut current_row = 0;
        let mut allele_offsets = Vec::with_capacity(total_reads);
        for (member_idx, member) in members.iter().enumerate() {
            let start_offset = current_row;
            for (allele_idx, allele) in member.iter().enumerate() {
                if allele_idx > 0 {
                    allele_offsets.push(current_row);
                }
                for (read, _) in &allele.read_aligns {
                    let read_start = read
                        .start_offset
                        .expect("Invariant: read.start_offset must be present");
                    let read_end = read
                        .end_offset
                        .expect("Invariant: read.end_offset must be present");
                    let mismatch_offsets = read
                        .mismatch_offsets
                        .as_ref()
                        .expect("Invariant: read.mismatch_offsets must be present");

                    let start_col = mismatch_positions
                        .binary_search(&read_start)
                        .unwrap_or_else(|x| x);
                    let end_col = mismatch_positions
                        .binary_search(&read_end)
                        .unwrap_or_else(|x| x);
                    for col in start_col..end_col {
                        let offset = mismatch_positions[col];
                        matrix[[current_row, col]] =
                            mismatch_offsets.binary_search(&offset).is_ok() as u8;
                    }
                    current_row += 1;
                }
            }
            offsets[member_idx] = MemberOffsets {
                start_offset,
                end_offset: current_row - 1,
                allele_offsets: allele_offsets.clone(),
            };
            allele_offsets.clear();
        }

        Some(TrinaryMatrix {
            matrix,
            offsets: offsets.try_into().expect("Wrong number of offsets"),
            mismatch_offsets: mismatch_positions,
            members: member_types,
        })
    }

    /// Retrieves a submatrix view for a specific member's SNPs.
    ///
    /// Returns a view of the SNP matrix that corresponds to the SNPs
    /// associated with the specified member.
    ///
    /// # Arguments
    ///
    /// * `member` - The member for which to retrieve the SNP submatrix.
    ///
    /// # Returns
    ///
    /// An `ArrayView2<u8>` representing the submatrix of SNPs for the given member.
    pub fn member_submatrix(&self, member: &Member) -> ArrayView2<u8> {
        let offset = &self.offsets[member.index()];
        self.matrix
            .slice(s![offset.start_offset..=offset.end_offset, ..])
    }

    /// Retrieves a submatrix view for a specific allele of a member.
    ///
    /// Returns a view of the SNP matrix that corresponds to the SNPs
    /// associated with a particular allele of the specified member.
    ///
    /// # Arguments
    ///
    /// * `member` - The member to which the allele belongs.
    /// * `allele_idx` - The index of the allele within the member's alleles.
    ///
    /// # Returns
    ///
    /// An `Option<ArrayView2<u8>>` representing the submatrix of SNPs for the given allele,
    /// or `None` if the allele index is out of bounds.
    pub fn allele_submatrix(&self, member: &Member, allele_idx: usize) -> Option<ArrayView2<u8>> {
        let offset = &self.offsets[member.index()];

        // TODO: Check
        // if allele_idx > offset.allele_offsets.len() + 1 {
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

    /// Creates an iterator over submatrices for each allele of a member.
    ///
    /// Returns an iterator that yields views of the SNP matrix for each allele
    /// of the specified member. Each item in the iterator is a submatrix corresponding
    /// to a single allele's SNPs.
    ///
    /// # Arguments
    ///
    /// * `member` - The member whose alleles are to be iterated over.
    ///
    /// # Returns
    ///
    /// An iterator over `ArrayView2<u8>` where each item is a submatrix for an allele of the member.
    pub fn iter_alleles(&self, member: Member) -> impl Iterator<Item = ArrayView2<u8>> {
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

    pub fn count_alleles(&self, member: &Member) -> usize {
        let offset = &self.offsets[member.index()];
        offset.allele_offsets.len() + 1
    }
}

impl<const N: usize> fmt::Display for TrinaryMatrix<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for member in &self.members {
            writeln!(f, "{:?}=(", member)?;
            let allele_count = self.count_alleles(member);
            for allele_idx in 0..allele_count {
                match self.allele_submatrix(member, allele_idx) {
                    Some(submatrix) => {
                        write!(f, "  [")?;
                        for row in submatrix.outer_iter() {
                            write!(f, "\n    \"")?;
                            for &val in row.iter() {
                                let ch = if val == 2 {
                                    '-'
                                } else {
                                    char::from_digit(val as u32, 10)
                                        .expect("Value should be 0 or 1")
                                };
                                write!(f, "{}", ch)?;
                            }
                            write!(f, "\",")?;
                        }
                        writeln!(f, "\n  ],")?;
                    }
                    None => {
                        writeln!(f, "  Allele {}: None", allele_idx)?;
                    }
                }
            }
            writeln!(f, ")")?;
        }
        write!(f, "")
    }
}

/// Provides methods for working with the `TrinaryMatrix`.
impl TrinaryMatrix<2> {
    /// Creates a new `TrinaryMatrix` from the given alleles of a duo.
    ///
    /// # Arguments
    ///
    /// * `sample1` - A slice of `Allele` instances for the first sample.
    /// * `sample2` - A slice of `Allele` instances for the second sample.
    ///
    /// # Returns
    ///
    /// An `Option<TrinaryMatrix>` which is `None` if no SNPs are found, or contains the constructed matrix otherwise.
    pub fn new_duo(sample1: &AlleleSet, sample2: &AlleleSet) -> Option<Self> {
        Self::construct_matrix(&[sample1, sample2], vec![Member::Sample0, Member::Sample1])
    }
}

impl TrinaryMatrix<3> {
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
    pub fn new_trio(child: &AlleleSet, father: &AlleleSet, mother: &AlleleSet) -> Option<Self> {
        Self::construct_matrix(
            &[father, mother, child],
            vec![Member::Father, Member::Mother, Member::Child],
        )
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
    fn filter(&self, reads: &mut Vec<TrgtRead>);
}

/// A `ReadFilter` that filters reads based on distance criteria.
///
/// This filter removes SNPs that are too far from the region of interest, as defined by a maximum distance threshold.
pub struct FilterByDist;
impl ReadFilter for FilterByDist {
    fn filter(&self, reads: &mut Vec<TrgtRead>) {
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
    fn filter(&self, reads: &mut Vec<TrgtRead>) {
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
pub fn apply_read_filters(reads: &mut Vec<TrgtRead>, filters: &[&(dyn ReadFilter + Sync)]) {
    for filter in filters {
        filter.filter(reads);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{allele::Allele, read::TrgtReadBuilder};
    use itertools::Itertools;

    fn create_read_info(
        mismatch_offsets: Vec<i32>,
        start_offset: i32,
        end_offset: i32,
    ) -> TrgtRead {
        TrgtReadBuilder::default()
            .with_mismatch_offsets(Some(mismatch_offsets))
            .with_start_offset(Some(start_offset))
            .with_end_offset(Some(end_offset))
            .build()
    }

    #[test]
    fn test_empty_trinary_matrix_trio() {
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
        let trinaray_mat =
            TrinaryMatrix::new_trio(&child_alleles, &father_alleles, &mother_alleles);
        assert_eq!(trinaray_mat, None);
    }

    #[test]
    fn test_empty_trinary_matrix_duo() {
        let sample0_alleles = AlleleSet {
            alleles: vec![],
            hp_counts: [0; 3],
        };
        let sample1_alleles = AlleleSet {
            alleles: vec![],
            hp_counts: [0; 3],
        };
        let trinaray_mat = TrinaryMatrix::new_duo(&sample0_alleles, &sample1_alleles);
        assert_eq!(trinaray_mat, None);
    }

    #[test]
    fn test_trinary_matrix_duo() {
        let sample0_alleles = AlleleSet {
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

        let sample1_alleles = AlleleSet {
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

        let mat = TrinaryMatrix::new_duo(&sample0_alleles, &sample1_alleles).unwrap();

        // Test dimensionality
        assert_eq!(mat.matrix.shape(), &[15, 11]);

        // Counts alleles
        assert_eq!(mat.count_alleles(&Member::Sample0), 2);
        assert_eq!(mat.count_alleles(&Member::Sample1), 4);

        // Check for correctness of mismatch offset labels
        assert_eq!(
            mat.mismatch_offsets,
            vec![-8, -4, -2, -1, 2, 3, 5, 6, 8, 23, 42,]
        );

        // Row 0 should correspond to the first read of sample0
        assert_eq!(
            mat.matrix.row(0),
            mat.matrix
                .row(mat.offsets[Member::Sample0.index()].start_offset),
        );

        // Check for correctness of matrix row 0
        assert_eq!(
            mat.matrix.row(0).to_vec(),
            vec![1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
        );

        // Check for correctness for the first read of sample1
        assert_eq!(
            mat.matrix
                .row(mat.offsets[Member::Sample1.index()].start_offset)
                .to_vec(),
            vec![1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
        );

        // Allele submatrix 4 (and beyond) of the sample1 should be None
        assert_eq!(None, mat.allele_submatrix(&Member::Sample1, 4));

        // Sample1 has 4 alleles, check that allele 3 is correct
        assert_eq!(
            mat.allele_submatrix(&Member::Sample1, 2)
                .unwrap()
                .rows()
                .into_iter()
                .map(|row| row.to_vec())
                .collect_vec(),
            vec![[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
        );
    }

    #[test]
    fn test_trinary_matrix_trio() {
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

        let mat =
            TrinaryMatrix::new_trio(&child_alleles, &father_alleles, &mother_alleles).unwrap();

        // Test dimensionality
        assert_eq!(mat.matrix.shape(), &[18, 11]);

        // Counts alleles
        assert_eq!(mat.count_alleles(&Member::Child), 4);
        assert_eq!(mat.count_alleles(&Member::Mother), 1);
        assert_eq!(mat.count_alleles(&Member::Father), 2);

        // Check for correctness of mismatch offset labels
        assert_eq!(
            mat.mismatch_offsets,
            vec![-8, -4, -2, -1, 2, 3, 5, 6, 8, 23, 42,]
        );

        // Row 0 should correspond to the first read of the father
        assert_eq!(
            mat.matrix.row(0),
            mat.matrix
                .row(mat.offsets[Member::Father.index()].start_offset),
        );

        // Check for correctness of matrix row 0
        assert_eq!(
            mat.matrix.row(0).to_vec(),
            vec![1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
        );

        // Check for correctness for the first read of the mother
        assert_eq!(
            mat.matrix
                .row(mat.offsets[Member::Mother.index()].start_offset)
                .to_vec(),
            vec![2, 1, 0, 1, 1, 0, 0, 1, 2, 2, 2]
        );

        // Mother has only one allele so the mother submatrix should be equivalent to allele submatrix 0
        assert_eq!(
            mat.member_submatrix(&Member::Mother),
            mat.allele_submatrix(&Member::Mother, 0).unwrap()
        );

        // Allele submatrix 1 (and beyond) of the mother should be None
        assert_eq!(None, mat.allele_submatrix(&Member::Mother, 1));

        // Child has 4 alleles, check that allele 3 is correct
        assert_eq!(
            mat.allele_submatrix(&Member::Child, 2)
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
