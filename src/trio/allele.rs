//! Provides functionality for handling alleles and related operations.
//!
//! This module includes structures and functions to represent alleles, load them from VCF and BAM
//! files, perform read alignments, and process allele information for further analysis.
//!
use super::denovo;
use crate::{
    aligner::WFAligner,
    allele::{
        join_allele_attribute, load_alleles, serialize_as_display, serialize_with_precision,
        Allele, AlleleSet,
    },
    handles::TrioLocalData,
    locus::Locus,
    math,
    model::{AlleleOrigin, DenovoStatus, Params, QuickMode},
    util::Result,
};
use itertools::Itertools;
use serde::Serialize;
use std::{cmp::max, collections::HashSet};

/// Represents the result of allele processing, including various statistics and classifications.
#[derive(Debug, PartialEq, Serialize, Clone)]
#[allow(non_snake_case)]
pub struct AlleleResult {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub motifs: String,
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
    #[serde(skip)]
    pub read_ids: Option<Vec<String>>,
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

pub fn process_alleles(
    locus: &Locus,
    handle: &mut TrioLocalData,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<Vec<AlleleResult>> {
    let mut template_result = AlleleResult {
        chrom: locus.region.contig.to_string(),
        start: locus.region.start as usize,
        end: locus.region.end as usize,
        motifs: locus.motifs.join(","),
        trid: locus.id.clone(),
        genotype: 0,
        denovo_coverage: 0,
        allele_coverage: 0,
        allele_ratio: 0.0,
        child_coverage: 0,
        child_ratio: 0.0,
        mean_diff_father: 0.0,
        mean_diff_mother: 0.0,
        father_dropout_prob: 0.0,
        mother_dropout_prob: 0.0,
        allele_origin: AlleleOrigin::Unknown,
        denovo_status: DenovoStatus::Unknown,
        per_allele_reads_father: ".".to_string(),
        per_allele_reads_mother: ".".to_string(),
        per_allele_reads_child: ".".to_string(),
        father_dropout: ".".to_string(),
        mother_dropout: ".".to_string(),
        child_dropout: ".".to_string(),
        index: 0,
        father_MC: ".".to_string(),
        mother_MC: ".".to_string(),
        child_MC: ".".to_string(),
        father_AL: ".".to_string(),
        mother_AL: ".".to_string(),
        child_AL: ".".to_string(),
        father_overlap_coverage: ".".to_string(),
        mother_overlap_coverage: ".".to_string(),
        read_ids: None,
    };

    let father_alleles = load_alleles(locus, &mut handle.father, params, aligner);
    let mother_alleles = load_alleles(locus, &mut handle.mother, params, aligner);
    let child_alleles = load_alleles(locus, &mut handle.child, params, aligner);

    if let Ok(ref alleles) = father_alleles {
        template_result.father_dropout = alleles
            .get_naive_dropout(locus, &handle.father.karyotype)
            .to_string();
        template_result.father_dropout_prob = math::get_dropout_prob(alleles);
        template_result.per_allele_reads_father = math::get_per_allele_reads(alleles)
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<String>>()
            .join(",");
        template_result.father_MC = join_allele_attribute(alleles, |a| &a.motif_count);
        template_result.father_AL = join_allele_attribute(alleles, |a| &a.allele_length);
    }
    if let Ok(ref alleles) = mother_alleles {
        template_result.mother_dropout = alleles
            .get_naive_dropout(locus, &handle.mother.karyotype)
            .to_string();
        template_result.mother_dropout_prob = math::get_dropout_prob(alleles);
        template_result.per_allele_reads_mother = math::get_per_allele_reads(alleles)
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<String>>()
            .join(",");
        template_result.mother_MC = join_allele_attribute(alleles, |a| &a.motif_count);
        template_result.mother_AL = join_allele_attribute(alleles, |a| &a.allele_length);
    }
    if let Ok(ref alleles) = child_alleles {
        template_result.child_dropout = alleles
            .get_naive_dropout(locus, &handle.child.karyotype)
            .to_string();
        template_result.per_allele_reads_child = math::get_per_allele_reads(alleles)
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<String>>()
            .join(",");
        template_result.child_MC = join_allele_attribute(alleles, |a| &a.motif_count);
        template_result.child_AL = join_allele_attribute(alleles, |a| &a.allele_length);
    }

    let error_labels: Vec<&str> = [
        ("F", &father_alleles),
        ("M", &mother_alleles),
        ("C", &child_alleles),
    ]
    .iter()
    .filter_map(|&(label, res)| if res.is_err() { Some(label) } else { None })
    .collect();

    if !error_labels.is_empty() {
        log::warn!(
            "Skipping TRID={} missing genotyping in: {}",
            locus.id,
            error_labels.join(",")
        );
        return Ok(vec![template_result]);
    }

    let father_alleles = father_alleles.unwrap();
    let mother_alleles = mother_alleles.unwrap();
    let child_alleles = child_alleles.unwrap();

    if let Some(quick_mode) = &params.quick_mode {
        let should_skip = if quick_mode.is_zero() {
            check_field_equivalence(&father_alleles, &mother_alleles, &child_alleles, quick_mode)
        } else {
            check_field_similarity(&father_alleles, &mother_alleles, &child_alleles, quick_mode)
        };

        if should_skip {
            return Ok(vec![template_result]);
        }
    }

    let mut out_vec = Vec::new();
    for dna in denovo::assess_denovo(
        &mother_alleles,
        &father_alleles,
        &child_alleles,
        params,
        aligner,
    ) {
        let mut result = template_result.clone();
        result.genotype = dna.genotype;
        result.denovo_coverage = dna.denovo_coverage;
        result.allele_coverage = dna.allele_coverage;
        result.allele_ratio = if dna.allele_coverage == 0 {
            0.0
        } else {
            dna.denovo_coverage as f64 / dna.allele_coverage as f64
        };
        result.child_coverage = dna.child_coverage;
        result.child_ratio = dna.denovo_coverage as f64 / dna.child_coverage as f64;
        result.mean_diff_father = dna.mean_diff_father;
        result.mean_diff_mother = dna.mean_diff_mother;
        result.allele_origin = dna.allele_origin;
        result.denovo_status = dna.denovo_status;
        result.index = dna.index;
        result.father_overlap_coverage = dna
            .father_overlap_coverage
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<String>>()
            .join(",");
        result.mother_overlap_coverage = dna
            .mother_overlap_coverage
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<String>>()
            .join(",");
        result.read_ids = dna.read_ids;
        out_vec.push(result);
    }

    if out_vec.is_empty() {
        out_vec.push(template_result);
    }

    Ok(out_vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::TrgtRead;

    #[test]
    fn test_check_minimal_cost_inheritance() {
        let create_allele_set = |mc: Vec<&str>, al: Vec<&str>| AlleleSet {
            alleles: mc
                .into_iter()
                .zip(al)
                .map(|(mc, al)| Allele {
                    read_aligns: Vec::<(TrgtRead, i32)>::new(),
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
