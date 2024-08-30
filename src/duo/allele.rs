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
    handles::DuoLocalData,
    locus::Locus,
    math,
    util::{Params, QuickMode, Result},
};
use anyhow::anyhow;
use serde::Serialize;
use std::{cmp::Ordering, collections::HashSet};

/// Checks if the specified field values of two AlleleSet instances are equivalent by str comparison.
///
/// # Arguments
///
/// * `alleles_a` - A reference to the first AlleleSet.
/// * `alleles_b` - A reference to the second AlleleSet.
/// * `quick_mode` - A reference to the QuickMode enum specifying which field to compare and the tolerance.
///
/// # Returns
///
/// Returns true if the specified field values are equivalent, false otherwise.
pub fn check_allele_field_equivalence<'a>(
    alleles_a: &'a AlleleSet,
    alleles_b: &'a AlleleSet,
    quick_mode: &QuickMode,
) -> bool {
    let get_field_value = |allele: &'a Allele| -> &str {
        match quick_mode {
            QuickMode::AL(_) => &allele.allele_length,
            QuickMode::MC(_) => &allele.motif_count,
        }
    };

    let set_a: HashSet<&str> = alleles_a.iter().map(get_field_value).collect();
    let set_b: HashSet<&str> = alleles_b.iter().map(get_field_value).collect();

    set_a == set_b
}

fn check_field_similarity(
    alleles_a: &AlleleSet,
    alleles_b: &AlleleSet,
    quick_mode: &QuickMode,
) -> bool {
    let get_field_value = |allele: &Allele| -> i32 {
        match quick_mode {
            QuickMode::AL(_) => allele.allele_length.parse::<i32>().unwrap(),
            QuickMode::MC(_) => allele.motif_count.parse::<i32>().unwrap(),
        }
    };

    let tolerance = match quick_mode {
        QuickMode::AL(tol) | QuickMode::MC(tol) => tol.unwrap(),
    };

    let mut a_values: Vec<i32> = alleles_a.iter().map(get_field_value).collect();
    let mut b_values: Vec<i32> = alleles_b.iter().map(get_field_value).collect();

    fn is_similar(a: i32, b: i32, tolerance: f64) -> bool {
        let diff = (a - b).abs() as f64;
        let max_val = a.max(b) as f64;
        diff / max_val <= tolerance
    }

    if a_values.len() != b_values.len() {
        return false;
    }

    a_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    b_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    a_values
        .iter()
        .zip(b_values.iter())
        .all(|(&a, &b)| is_similar(a, b, tolerance))
}

#[derive(Debug, PartialEq, Serialize)]
#[allow(non_snake_case)]
pub struct AlleleResult {
    pub trid: String,
    pub genotype: usize,
    pub denovo_coverage: usize,
    pub allele_coverage: usize,
    #[serde(serialize_with = "serialize_with_precision")]
    pub allele_ratio: f64,
    pub a_coverage: usize,
    #[serde(serialize_with = "serialize_with_precision")]
    pub a_ratio: f64,
    #[serde(serialize_with = "serialize_with_precision")]
    pub mean_diff_b: f32,
    #[serde(serialize_with = "serialize_as_display")]
    pub per_allele_reads_a: String,
    pub per_allele_reads_b: String,
    pub a_dropout: String,
    pub b_dropout: String,
    pub index: usize,
    pub a_MC: String,
    pub b_MC: String,
    pub a_AL: String,
    pub b_AL: String,
    pub b_overlap_coverage: String,
}

pub fn process_alleles(
    locus: &Locus,
    handle: &mut DuoLocalData,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<Vec<AlleleResult>> {
    let a_alleles = load_alleles_handle('A', locus, &mut handle.sample1, params, aligner)?;
    let b_alleles = load_alleles_handle('B', locus, &mut handle.sample2, params, aligner)?;

    if let Some(quick_mode) = &params.quick_mode {
        let should_skip = if quick_mode.is_zero() {
            check_allele_field_equivalence(&a_alleles, &b_alleles, quick_mode)
        } else {
            check_field_similarity(&a_alleles, &b_alleles, quick_mode)
        };

        if should_skip {
            return Err(anyhow!("No significant difference at locus: {}", locus.id));
        }
    }

    let a_reads = math::get_per_allele_reads(&a_alleles)
        .iter()
        .map(|a| a.to_string())
        .collect::<Vec<String>>()
        .join(",");
    let b_reads = math::get_per_allele_reads(&b_alleles)
        .iter()
        .map(|a| a.to_string())
        .collect::<Vec<String>>()
        .join(",");

    let a_mc = join_allele_attribute(&a_alleles, |a| &a.motif_count);
    let b_mc = join_allele_attribute(&b_alleles, |a| &a.motif_count);
    let a_al = join_allele_attribute(&a_alleles, |a| &a.allele_length);
    let b_al = join_allele_attribute(&b_alleles, |a| &a.allele_length);

    let mut out_vec = Vec::new();
    for dna in denovo::assess_denovo(&a_alleles, &b_alleles, params, aligner) {
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
            a_coverage: dna.a_coverage,
            a_ratio: dna.denovo_coverage as f64 / dna.a_coverage as f64,
            mean_diff_b: dna.mean_diff_b,
            per_allele_reads_a: a_reads.to_owned(),
            per_allele_reads_b: b_reads.to_owned(),
            a_dropout: a_alleles
                .get_naive_dropout(locus, &handle.sample1.karyotype)
                .to_string(),
            b_dropout: b_alleles
                .get_naive_dropout(locus, &handle.sample2.karyotype)
                .to_string(),
            index: dna.index,
            a_MC: a_mc.to_owned(),
            b_MC: b_mc.to_owned(),
            a_AL: a_al.to_owned(),
            b_AL: b_al.to_owned(),
            b_overlap_coverage: dna
                .b_overlap_coverage
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
    use crate::allele::{Allele, AlleleSet};

    #[test]
    fn test_check_allele_field_equivalence() {
        let allele1 = Allele {
            read_aligns: vec![],
            seq: vec![],
            motif_count: "3".to_string(),
            allele_length: "15".to_string(),
            genotype: 0,
            index: 0,
        };

        let allele2 = Allele {
            read_aligns: vec![],
            seq: vec![],
            motif_count: "4".to_string(),
            allele_length: "20".to_string(),
            genotype: 1,
            index: 1,
        };

        let allele3 = Allele {
            read_aligns: vec![],
            seq: vec![],
            motif_count: "3".to_string(),
            allele_length: "15".to_string(),
            genotype: 2,
            index: 2,
        };

        let allele6 = Allele {
            read_aligns: vec![],
            seq: vec![],
            motif_count: "5".to_string(),
            allele_length: "250".to_string(),
            genotype: 3,
            index: 3,
        };

        let allele_set_a = AlleleSet {
            alleles: vec![allele1.clone(), allele2.clone()],
            hp_counts: [0, 0, 0],
        };

        let allele_set_b = AlleleSet {
            alleles: vec![allele2.clone(), allele1.clone()],
            hp_counts: [0, 0, 0],
        };

        let allele_set_c = AlleleSet {
            alleles: vec![allele1.clone(), allele2.clone()],
            hp_counts: [0, 0, 0],
        };

        let allele_set_d = AlleleSet {
            alleles: vec![allele1.clone()],
            hp_counts: [0, 0, 0],
        };

        let allele_set_e = AlleleSet {
            alleles: vec![allele1.clone(), allele3.clone()],
            hp_counts: [0, 0, 0],
        };

        let allele_set_f = AlleleSet {
            alleles: vec![allele1.clone(), allele6],
            hp_counts: [0, 0, 0],
        };

        // Allele length (AL) equivalence
        assert!(check_allele_field_equivalence(
            &allele_set_a,
            &allele_set_b,
            &QuickMode::AL(None)
        ));
        assert!(check_allele_field_equivalence(
            &allele_set_a,
            &allele_set_c,
            &QuickMode::AL(None)
        ));
        assert!(!check_allele_field_equivalence(
            &allele_set_a,
            &allele_set_d,
            &QuickMode::AL(None)
        ));
        assert!(!check_allele_field_equivalence(
            &allele_set_a,
            &allele_set_e,
            &QuickMode::AL(None)
        ));
        assert!(!check_allele_field_equivalence(
            &allele_set_a,
            &allele_set_f,
            &QuickMode::AL(None)
        ));
        assert!(!check_allele_field_equivalence(
            &allele_set_b,
            &allele_set_f,
            &QuickMode::AL(None)
        ));
        assert!(!check_allele_field_equivalence(
            &allele_set_c,
            &allele_set_f,
            &QuickMode::AL(None)
        ));

        assert!(check_field_similarity(
            &allele_set_a,
            &allele_set_e,
            &QuickMode::AL(Some(0.3))
        ));
    }
}
