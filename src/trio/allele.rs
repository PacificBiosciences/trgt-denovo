//! Provides functionality for handling alleles and related operations.
//!
//! This module includes structures and functions to represent alleles, load them from VCF and BAM
//! files, perform read alignments, and process allele information for further analysis.
//!
use super::denovo;
use crate::{
    aligner::WFAligner,
    allele::{
        join_allele_attribute, load_alleles_handle, serialize_as_display, serialize_with_precision,
        Allele, AlleleSet,
    },
    handles::TrioLocalData,
    locus::Locus,
    math,
    util::{AlleleOrigin, DenovoStatus, Params, QuickMode, Result},
};
use anyhow::anyhow;
use itertools::Itertools;
use serde::Serialize;
use std::{cmp::max, collections::HashSet};

/// Represents the result of allele processing, including various statistics and classifications.
#[derive(Debug, PartialEq, Serialize)]
#[allow(non_snake_case)]
pub struct AlleleResult {
    pub trid: String,
    pub genotype: usize,
    pub denovo_coverage: usize,
    pub allele_coverage: usize,
    #[serde(serialize_with = "serialize_with_precision")]
    pub allele_ratio: f64,
    pub child_coverage: usize,
    #[serde(serialize_with = "serialize_with_precision")]
    pub child_ratio: f64,
    #[serde(serialize_with = "serialize_with_precision")]
    pub mean_diff_father: f32,
    #[serde(serialize_with = "serialize_with_precision")]
    pub mean_diff_mother: f32,
    #[serde(serialize_with = "serialize_with_precision")]
    pub father_dropout_prob: f64,
    #[serde(serialize_with = "serialize_with_precision")]
    pub mother_dropout_prob: f64,
    #[serde(serialize_with = "serialize_as_display")]
    pub allele_origin: AlleleOrigin,
    #[serde(serialize_with = "serialize_as_display")]
    pub denovo_status: DenovoStatus,
    pub per_allele_reads_father: String,
    pub per_allele_reads_mother: String,
    pub per_allele_reads_child: String,
    pub father_dropout: String,
    pub mother_dropout: String,
    pub child_dropout: String,
    pub index: usize,
    pub father_MC: String,
    pub mother_MC: String,
    pub child_MC: String,
    pub father_AL: String,
    pub mother_AL: String,
    pub child_AL: String,
    pub father_overlap_coverage: String,
    pub mother_overlap_coverage: String,
}

pub fn check_field_equivalence<'a>(
    alleles_father: &'a AlleleSet,
    alleles_mother: &'a AlleleSet,
    alleles_child: &'a AlleleSet,
    quick_mode: &QuickMode,
) -> bool {
    let get_field_value = |allele: &'a Allele| -> &str {
        match quick_mode {
            QuickMode::AL(_) => &allele.allele_length,
            QuickMode::MC(_) => &allele.motif_count,
        }
    };

    let father_values: Vec<&str> = alleles_father.iter().map(get_field_value).collect();
    let mother_values: Vec<&str> = alleles_mother.iter().map(get_field_value).collect();
    let mut child_values: Vec<&str> = alleles_child.iter().map(get_field_value).collect();
    child_values.sort();

    if child_values.len() == 1 {
        let child_value = child_values[0];
        return father_values.contains(&child_value) || mother_values.contains(&child_value);
    }

    let mut possible_pairs = HashSet::new();
    for (f, m) in father_values.iter().cartesian_product(mother_values.iter()) {
        let mut pair = vec![*f, *m];
        pair.sort();
        possible_pairs.insert(pair);
    }
    possible_pairs.contains(&child_values)
}

fn check_field_similarity(
    alleles_father: &AlleleSet,
    alleles_mother: &AlleleSet,
    alleles_child: &AlleleSet,
    quick_mode: &QuickMode,
) -> bool {
    let get_field_value = |allele: &Allele| -> i32 {
        match quick_mode {
            QuickMode::AL(_) => allele.allele_length.parse::<i32>().unwrap(),
            QuickMode::MC(_) => allele.motif_count.parse::<i32>().unwrap(),
        }
    };

    let tolerance = match quick_mode {
        QuickMode::AL(tol) | QuickMode::MC(tol) => *tol,
    }
    .unwrap();

    let father_values: Vec<i32> = alleles_father.iter().map(get_field_value).collect();
    let mother_values: Vec<i32> = alleles_mother.iter().map(get_field_value).collect();
    let mut child_values: Vec<i32> = alleles_child.iter().map(get_field_value).collect();
    child_values.sort();

    let relative_difference = if child_values.len() == 1 {
        // Single allele case
        let child_value = child_values[0];
        let diff_father: Vec<f64> = father_values
            .iter()
            .map(|&f| (f - child_value).abs() as f64)
            .collect();
        let diff_mother: Vec<f64> = mother_values
            .iter()
            .map(|&m| (m - child_value).abs() as f64)
            .collect();

        let best_father_diff = diff_father.into_iter().fold(f64::INFINITY, f64::min);
        let best_mother_diff = diff_mother.into_iter().fold(f64::INFINITY, f64::min);

        let relative_diff_father = best_father_diff / *father_values.iter().max().unwrap() as f64;
        let relative_diff_mother = best_mother_diff / *mother_values.iter().max().unwrap() as f64;

        f64::min(relative_diff_father, relative_diff_mother)
    } else {
        // Two alleles case
        let mut best_pair: Option<(i32, i32)> = None;
        let mut best_difference_sum = f64::INFINITY;

        for (f, m) in father_values.iter().cartesian_product(mother_values.iter()) {
            let pair = if f < m { (*f, *m) } else { (*m, *f) };

            let diff_0 = (pair.0 - child_values[0]).abs() as f64;
            let diff_1 = (pair.1 - child_values[1]).abs() as f64;
            let difference_sum = diff_0 + diff_1;

            if difference_sum < best_difference_sum {
                best_difference_sum = difference_sum;
                best_pair = Some(pair);
            }
        }

        if let Some((best_f, best_m)) = best_pair {
            let diff_0 = if best_f == 0 && child_values[0] == 0 {
                0.0
            } else {
                (best_f - child_values[0]).abs() as f64 / max(best_f, child_values[0]) as f64
            };

            let diff_1 = if best_m == 0 && child_values[1] == 0 {
                0.0
            } else {
                (best_m - child_values[1]).abs() as f64 / max(best_m, child_values[1]) as f64
            };

            f64::max(diff_0, diff_1)
        } else {
            f64::INFINITY
        }
    };

    relative_difference <= tolerance
}

/// Processes alleles for a given locus and handle.
///
/// Coordinates the loading and processing of alleles for a family trio (father,
/// mother, and child). It calculates statistics, dropout probabilities, and de novo allele
/// information, and compiles the results into `AlleleResult` instances.
///
/// # Arguments
///
/// * `locus` - A reference to the `Locus` for which alleles are to be processed.
/// * `handle` - An `Arc` containing the `Handles` for the family members.
/// * `params` - Parameters such as the length of the clipping to be applied to alignments or the quantile used for parental allele frequency calculations.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// A result containing a vector of `AlleleResult` instances if successful, or an error if not.
pub fn process_alleles(
    locus: &Locus,
    handle: &mut TrioLocalData,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<Vec<AlleleResult>> {
    let father_alleles = load_alleles_handle('F', locus, &mut handle.father, params, aligner)?;
    let mother_alleles = load_alleles_handle('M', locus, &mut handle.mother, params, aligner)?;
    let child_alleles = load_alleles_handle('C', locus, &mut handle.child, params, aligner)?;

    if let Some(quick_mode) = &params.quick_mode {
        let should_skip = if quick_mode.is_zero() {
            check_field_equivalence(&father_alleles, &mother_alleles, &child_alleles, quick_mode)
        } else {
            check_field_similarity(&father_alleles, &mother_alleles, &child_alleles, quick_mode)
        };

        if should_skip {
            return Err(anyhow!("No significant difference at locus: {}", locus.id));
        }
    }

    // TODO: ongoing work
    // let matrix = TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles).unwrap();
    // let max_lhs = TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles)
    //     .and_then(|trinary_mat| snp::inheritance_prob(&trinary_mat))
    //     .map(|(_inherit_p, max_lh)| max_lh)
    //     .unwrap_or((-1.0, -1.0));
    // let max_lhs = [max_lhs.0, max_lhs.1];

    let mother_dropout_prob = math::get_dropout_prob(&mother_alleles);
    let father_dropout_prob = math::get_dropout_prob(&father_alleles);

    let father_reads = math::get_per_allele_reads(&father_alleles)
        .iter()
        .map(|a| a.to_string())
        .collect::<Vec<String>>()
        .join(",");
    let mother_reads = math::get_per_allele_reads(&mother_alleles)
        .iter()
        .map(|a| a.to_string())
        .collect::<Vec<String>>()
        .join(",");
    let child_reads = math::get_per_allele_reads(&child_alleles)
        .iter()
        .map(|a| a.to_string())
        .collect::<Vec<String>>()
        .join(",");

    let father_mc = join_allele_attribute(&father_alleles, |a| &a.motif_count);
    let mother_mc = join_allele_attribute(&mother_alleles, |a| &a.motif_count);
    let child_mc = join_allele_attribute(&child_alleles, |a| &a.motif_count);
    let father_al = join_allele_attribute(&father_alleles, |a| &a.allele_length);
    let mother_al = join_allele_attribute(&mother_alleles, |a| &a.allele_length);
    let child_al = join_allele_attribute(&child_alleles, |a| &a.allele_length);

    let mut out_vec = Vec::new();
    for dna in denovo::assess_denovo(
        &mother_alleles,
        &father_alleles,
        &child_alleles,
        params,
        aligner,
    ) {
        let allele_ratio = if dna.allele_coverage == 0 {
            0.0
        } else {
            dna.denovo_coverage as f64 / dna.allele_coverage as f64
        };
        let output = AlleleResult {
            trid: locus.id.clone(),
            genotype: dna.genotype,
            denovo_coverage: dna.denovo_coverage,
            allele_coverage: dna.allele_coverage,
            allele_ratio,
            child_coverage: dna.child_coverage,
            child_ratio: dna.denovo_coverage as f64 / dna.child_coverage as f64,
            mean_diff_father: dna.mean_diff_father,
            mean_diff_mother: dna.mean_diff_mother,
            father_dropout_prob,
            mother_dropout_prob,
            allele_origin: dna.allele_origin,
            denovo_status: dna.denovo_status,
            per_allele_reads_father: father_reads.to_owned(),
            per_allele_reads_mother: mother_reads.to_owned(),
            per_allele_reads_child: child_reads.to_owned(),
            father_dropout: father_alleles
                .get_naive_dropout(locus, &handle.father.karyotype)
                .to_string(),
            mother_dropout: mother_alleles
                .get_naive_dropout(locus, &handle.mother.karyotype)
                .to_string(),
            child_dropout: child_alleles
                .get_naive_dropout(locus, &handle.child.karyotype)
                .to_string(),
            index: dna.index,
            father_MC: father_mc.to_owned(),
            mother_MC: mother_mc.to_owned(),
            child_MC: child_mc.to_owned(),
            father_AL: father_al.to_owned(),
            mother_AL: mother_al.to_owned(),
            child_AL: child_al.to_owned(),
            father_overlap_coverage: dna
                .father_overlap_coverage
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<String>>()
                .join(","),
            mother_overlap_coverage: dna
                .mother_overlap_coverage
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<String>>()
                .join(","),
        };
        out_vec.push(output);
    }
    Ok(out_vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::ReadInfo;

    #[test]
    fn test_check_minimal_cost_inheritance() {
        let create_allele_set = |mc: Vec<&str>, al: Vec<&str>| AlleleSet {
            alleles: mc
                .into_iter()
                .zip(al)
                .map(|(mc, al)| Allele {
                    read_aligns: Vec::<(ReadInfo, i32)>::new(),
                    seq: Vec::new(),
                    motif_count: mc.to_string(),
                    allele_length: al.to_string(),
                    genotype: 0,
                    index: 0,
                })
                .collect(),
            hp_counts: [0, 0, 0],
        };

        let father_alleles = create_allele_set(vec!["10", "12"], vec!["100", "120"]);
        let mother_alleles = create_allele_set(vec!["11", "13"], vec!["110", "130"]);

        let child_alleles_match = create_allele_set(vec!["10000", "1100"], vec!["100", "130"]);
        assert!(check_field_equivalence(
            &father_alleles,
            &mother_alleles,
            &child_alleles_match,
            &QuickMode::AL(None)
        ));

        let child_alleles_no_match = create_allele_set(vec!["14", "15"], vec!["140", "150"]);
        assert!(!check_field_equivalence(
            &father_alleles,
            &mother_alleles,
            &child_alleles_no_match,
            &QuickMode::AL(None)
        ));

        let child_alleles_partial = create_allele_set(vec!["10", "14"], vec!["100", "140"]);
        assert!(!check_field_equivalence(
            &father_alleles,
            &mother_alleles,
            &child_alleles_partial,
            &QuickMode::AL(None)
        ));

        assert!(check_field_equivalence(
            &create_allele_set(vec!["1100"], vec!["120"]),
            &create_allele_set(vec!["1100"], vec!["100"]),
            &create_allele_set(vec!["1100"], vec!["100", "120"]),
            &QuickMode::AL(None)
        ));

        let child_alleles_small_diff = create_allele_set(vec!["10", "14"], vec!["90", "150"]);
        assert!(check_field_similarity(
            &father_alleles,
            &mother_alleles,
            &child_alleles_small_diff,
            &QuickMode::AL(Some(0.5))
        ));
    }
}
