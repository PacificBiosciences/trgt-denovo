//! Statistical functions for analyzing allele data.
//!
//! This module provides functions to calculate various statistics related to alleles,
//! such as read counts, total reads, simple dropout probabilities, and allele frequencies.

use crate::allele::AlleleSet;
use std::cmp::Ordering;

/// Calculates the number of reads for each allele.
///
/// # Arguments
///
/// * `alleles` - A slice of `Allele` instances to analyze.
///
/// # Returns
///
/// A vector containing the number of reads for each allele.
pub fn get_per_allele_reads(alleles: &AlleleSet) -> Vec<usize> {
    alleles.iter().map(|a| a.read_aligns.len()).collect()
}

/// Calculates the total number of reads across all alleles.
///
/// # Arguments
///
/// * `alleles` - A slice of `Allele` instances to analyze.
///
/// # Returns
///
/// The total number of reads across all alleles.
pub fn get_total_reads(alleles: &AlleleSet) -> usize {
    get_per_allele_reads(alleles).iter().sum()
}

/// Estimates the probability of (simple) dropout for a set of alleles.
///
/// # Arguments
///
/// * `alleles` - A slice of `Allele` instances to analyze.
///
/// # Returns
///
/// The estimated probability of dropout for the alleles.
pub fn get_dropout_prob(alleles: &AlleleSet) -> f64 {
    assert!(!alleles.is_empty());
    let allele_frac: f64 = 0.5;
    let num_reads = get_total_reads(alleles) as f64;
    allele_frac.powf(num_reads)
}

/// Calculates the frequency of each allele based on read counts.
///
/// # Arguments
///
/// * `alleles` - A slice of `Allele` instances to analyze.
///
/// # Returns
///
/// A vector containing the frequency of each allele.
pub fn get_allele_freqs(alleles: &AlleleSet) -> Vec<f64> {
    let num_reads = get_total_reads(alleles) as f64;
    let allele_freqs: Vec<f64> = alleles
        .iter()
        .map(|a| a.read_aligns.len() as f64 / num_reads)
        .collect();
    allele_freqs
}

/// Calculates the quantile value of a sorted list of scores.
///
/// Computes the value at a given quantile in a list of scores. The list is sorted
/// and the quantile value is interpolated if necessary.
///
/// # Arguments
///
/// * `xs`: A mutable slice of scores.
/// * `q`: The quantile to compute (between 0 and 1).
///
/// # Returns
///
/// Returns an `Option<f64>` representing the quantile value, if the list is not empty.
pub fn quantile(xs: &mut [f64], q: f64) -> Option<f64> {
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

fn partition(data: &[i32]) -> Option<(Vec<i32>, i32, Vec<i32>)> {
    match data.len() {
        0 => None,
        _ => {
            let (pivot_slice, tail) = data.split_at(1);
            let pivot = pivot_slice[0];
            let (left, right) = tail.iter().fold((vec![], vec![]), |mut splits, next| {
                {
                    let (ref mut left, ref mut right) = &mut splits;
                    if next < &pivot {
                        left.push(*next);
                    } else {
                        right.push(*next);
                    }
                }
                splits
            });

            Some((left, pivot, right))
        }
    }
}

fn select(data: &[i32], k: usize) -> Option<i32> {
    let part = partition(data);
    match part {
        None => None,
        Some((left, pivot, right)) => {
            let pivot_idx = left.len();
            match pivot_idx.cmp(&k) {
                Ordering::Equal => Some(pivot),
                Ordering::Greater => select(&left, k),
                Ordering::Less => select(&right, k - (pivot_idx + 1)),
            }
        }
    }
}

pub fn median(data: &[i32]) -> Option<f64> {
    let size = data.len();
    if size == 0 {
        return None;
    }
    match size {
        even if even % 2 == 0 => {
            let fst_med = select(data, (even / 2) - 1);
            let snd_med = select(data, even / 2);

            match (fst_med, snd_med) {
                (Some(fst), Some(snd)) => Some((fst + snd) as f64 / 2.0),
                _ => None,
            }
        }
        odd => select(data, odd / 2).map(|x| x as f64),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_empty() {
        assert_eq!(median(&[]), None);
    }

    #[test]
    fn test_median_single() {
        assert_eq!(median(&[5]), Some(5.0));
    }

    #[test]
    fn test_median_even() {
        assert_eq!(median(&[3, 1, 4, 1, 5, 9, 2, 6, 5, 3]), Some(3.5));
    }

    #[test]
    fn test_median_odd() {
        assert_eq!(median(&[1, 3, 2, 5, 4]), Some(3.0));
    }

    #[test]
    fn test_partition_empty() {
        assert_eq!(partition(&[]), None);
    }

    #[test]
    fn test_partition_single() {
        assert_eq!(partition(&[5]), Some((vec![], 5, vec![])));
    }

    #[test]
    fn test_partition_multiple() {
        let result = partition(&[3, 1, 4, 1, 5, 9, 2, 6, 5, 3]);
        assert_eq!(result, Some((vec![1, 1, 2], 3, vec![4, 5, 9, 6, 5, 3])));
    }

    #[test]
    fn test_partition_all_equal() {
        assert_eq!(partition(&[2, 2, 2, 2]), Some((vec![], 2, vec![2, 2, 2])));
    }
}
