//! Statistical functions for analyzing allele data.
//!
//! This module provides functions to calculate various statistics related to alleles,
//! such as read counts, total reads, simple dropout probabilities, and allele frequencies.

use crate::allele::AlleleSet;

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
