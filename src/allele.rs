use crate::aligner::{AlignmentStatus, WFAligner};
use crate::denovo::{self, AlleleOrigin, DenovoStatus};
use crate::handles;
use crate::locus::Locus;
use crate::read::ReadInfo;
use crate::snp::{self, TrinaryMatrix};
use crate::stats;
use crate::util::Result;
use anyhow::anyhow;
use itertools::Itertools;
use log;
use noodles::{
    bam,
    bgzf::Reader,
    sam,
    vcf::{
        self,
        record::genotypes::{self, sample::value::Genotype},
        record::info::field,
        Record,
    },
};
use once_cell::sync::Lazy;
use serde::Serialize;
use std::{
    fs::File,
    sync::{Arc, Mutex},
};

#[derive(Debug, PartialEq)]
pub struct Allele {
    pub seq: Vec<u8>,
    pub genotype: usize,
    pub read_aligns: Vec<(ReadInfo, i32)>,
    pub motif_count: String,
    pub index: usize,
}

impl Allele {
    pub fn dummy(reads: Vec<ReadInfo>) -> Allele {
        Allele {
            seq: vec![],
            genotype: 0,
            read_aligns: reads.into_iter().map(|info| (info, 0i32)).collect_vec(),
            motif_count: String::from("0"),
            index: 0,
        }
    }
}

pub static FILTERS: &[&'static (dyn snp::ReadFilter + Sync)] =
    &[&snp::FilterByDist, &snp::FilterByFreq];

pub fn load_alleles(
    locus: &Locus,
    subhandle: &handles::SubHandle,
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Result<Vec<Allele>> {
    let (allele_seqs, motif_counts, genotype_indices) =
        get_allele_seqs(locus, &subhandle.vcf, &subhandle.vcf_header)?;
    let mut reads = get_reads(&subhandle.bam, &subhandle.bam_header, locus)?;

    snp::apply_read_filters(&mut reads, FILTERS);

    let reads_by_allele = assign_reads(&allele_seqs, reads, clip_len, aligner);
    let alleles = allele_seqs
        .into_iter()
        .enumerate()
        .map(|(index, allele_seq)| Allele {
            seq: allele_seq,
            genotype: genotype_indices[index],
            read_aligns: reads_by_allele[index].to_owned(),
            motif_count: motif_counts[index].to_owned(),
            index,
        })
        .collect();
    Ok(alleles)
}

fn assign_reads(
    alleles: &[Vec<u8>],
    reads: Vec<ReadInfo>,
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Vec<Vec<(ReadInfo, i32)>> {
    let mut reads_by_allele = vec![Vec::new(); alleles.len()];
    let mut index_flip: usize = 0;
    for read in reads {
        let mut max_score = None;
        let mut max_aligns = Vec::new();
        for (i, a) in alleles.iter().enumerate() {
            if let AlignmentStatus::StatusAlgCompleted = aligner.align_end_to_end(&read.bases, a) {
                let score = aligner.cigar_score_clipped(clip_len);
                match max_score {
                    None => {
                        max_score = Some(score);
                        max_aligns = vec![(i, score)];
                    }
                    Some(max) if score > max => {
                        max_score = Some(score);
                        max_aligns = vec![(i, score)];
                    }
                    Some(max) if score == max => {
                        max_aligns.push((i, score));
                    }
                    _ => (),
                }
            }
        }
        if !max_aligns.is_empty() {
            if max_aligns.len() > 1 {
                index_flip = (index_flip + 1) % max_aligns.len();
            }
            let (allele_index, align) = max_aligns[index_flip % max_aligns.len()];
            reads_by_allele[allele_index].push((read, align));
        }
    }
    reads_by_allele
}

pub fn get_reads(
    bam: &Arc<Mutex<bam::IndexedReader<Reader<File>>>>,
    bam_header: &Arc<sam::Header>,
    locus: &Locus,
) -> Result<Vec<ReadInfo>> {
    let mut bam = bam
        .lock()
        .map_err(|e| anyhow!("Failed to acquire lock: {}", e))?;
    let query = bam.query(bam_header, &locus.region)?;

    let reads: Vec<_> = query
        .map(|record| record.map_err(|e| e.into()).map(ReadInfo::new))
        .collect::<Result<Vec<_>>>()?;

    Ok(reads)
}

// TODO: improve, for now it just pulls out the genotype index and checks duplicates
fn is_homozygous(genotypes: &Genotype) -> bool {
    let alleles: Vec<_> = genotypes
        .iter()
        .filter_map(|allele| allele.position())
        .collect();
    alleles.into_iter().unique().count() == 1
}

static TRID_KEY: Lazy<field::Key> = Lazy::new(|| "TRID".parse().unwrap());
static MC_KEY: Lazy<genotypes::keys::Key> = Lazy::new(|| "MC".parse().unwrap());

fn get_allele_seqs(
    locus: &Locus,
    vcf: &Arc<Mutex<vcf::IndexedReader<File>>>,
    vcf_header: &Arc<vcf::Header>,
) -> Result<(Vec<Vec<u8>>, Vec<String>, Vec<usize>)> {
    let mut vcf = vcf
        .lock()
        .map_err(|_| anyhow!("Error locking Mutex for vcf::IndexedReader"))?;

    let query = vcf.query(vcf_header, &locus.region)?;
    let locus_id = &locus.id;

    if let Some(result) = query.into_iter().next() {
        let record = result?;
        let info = record.info();

        let trid = info.get(&*TRID_KEY).unwrap().unwrap().to_string();
        if &trid != locus_id {
            return Err(anyhow!("TRID={} missing", locus_id));
        }

        let format = record.genotypes().get_index(0).unwrap();
        let genotype = format.genotype().unwrap();
        if genotype.is_err() {
            return Err(anyhow!("TRID={} misses genotyping", locus_id));
        }
        let genotype = genotype.unwrap();

        let mc_field = format.get(&*MC_KEY).unwrap().unwrap().to_string();
        let motif_counts: Vec<String> = mc_field.split(',').map(ToString::to_string).collect();

        let (alleles, genotype_indices) = process_genotypes(&genotype, &record, locus)?;
        return Ok((alleles, motif_counts, genotype_indices));
    }
    Err(anyhow!("TRID={} missing", &locus.id))
}

fn process_genotypes(
    genotype: &Genotype,
    record: &Record,
    locus: &Locus,
) -> Result<(Vec<Vec<u8>>, Vec<usize>)> {
    let reference_bases = record.reference_bases().to_string().into_bytes();
    let alternate_bases = record
        .alternate_bases()
        .iter()
        .map(|base| base.to_string().into_bytes())
        .collect::<Vec<_>>();
    let mut seq = locus.left_flank.clone();
    let mut alleles = Vec::new();
    let mut genotype_indices = Vec::new();

    for allele in genotype.iter() {
        let allele_index: usize = allele
            .position()
            .ok_or_else(|| anyhow!("Allele position missing for TRID={}", locus.id))?;
        match allele_index {
            0 => seq.extend_from_slice(&reference_bases),
            _ => seq.extend_from_slice(
                alternate_bases
                    .get(allele_index - 1)
                    .ok_or_else(|| anyhow!("Invalid allele index for TRID={}", locus.id))?,
            ),
        }
        seq.extend_from_slice(&locus.right_flank);
        alleles.push(seq.clone());
        seq.truncate(locus.left_flank.len()); // reset the sequence for the next allele

        genotype_indices.push(allele_index);

        if genotype.len() > 1 && is_homozygous(genotype) {
            break;
        }
    }
    Ok((alleles, genotype_indices))
}

#[derive(Debug, PartialEq, Serialize)]
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
    pub index: usize,
    pub father_motif_counts: String,
    pub mother_motif_counts: String,
    pub child_motif_counts: String,
    #[serde(serialize_with = "serialize_with_precision")]
    pub maxlh: f64,
}

fn serialize_as_display<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.collect_str(value)
}

fn serialize_with_precision<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    let formatted_value = format!("{:.4}", value);
    serializer.serialize_str(&formatted_value)
}

fn load_alleles_handle(
    role: char,
    locus: &Locus,
    subhandle: &handles::SubHandle,
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Result<Vec<Allele>> {
    load_alleles(locus, subhandle, clip_len, aligner).map_err(|err| {
        log::warn!("Skipping {} in {}", err, role);
        err
    })
}

pub fn process_alleles(
    locus: &Locus,
    handle: Arc<handles::Handles>,
    clip_len: usize,
    parent_quantile: f64,
    aligner: &mut WFAligner,
) -> Result<Vec<AlleleResult>> {
    let father_alleles = load_alleles_handle('F', locus, &handle.father, clip_len, aligner)?;
    let mother_alleles = load_alleles_handle('M', locus, &handle.mother, clip_len, aligner)?;
    let child_alleles = load_alleles_handle('C', locus, &handle.child, clip_len, aligner)?;

    // TODO: ongoing work, the maximum likelihood is obtained naively
    let max_lhs = TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles)
        .and_then(|trinary_mat| snp::inheritance_prob(&trinary_mat))
        .map(|(_inherit_p, max_lh)| max_lh)
        .unwrap_or((-1.0, -1.0));
    let max_lhs = vec![max_lhs.0, max_lhs.1];

    let mother_dropout_prob = stats::get_dropout_prob(&mother_alleles);
    let father_dropout_prob = stats::get_dropout_prob(&father_alleles);

    let father_reads = stats::get_per_allele_reads(&father_alleles);
    let mother_reads = stats::get_per_allele_reads(&mother_alleles);
    let child_reads = stats::get_per_allele_reads(&child_alleles);

    let father_motifs = father_alleles
        .iter()
        .map(|a| a.motif_count.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let mother_motifs = mother_alleles
        .iter()
        .map(|a| a.motif_count.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let child_motifs = child_alleles
        .iter()
        .map(|a| a.motif_count.to_string())
        .collect::<Vec<_>>()
        .join(",");

    let mut out_vec = Vec::<AlleleResult>::new();
    for (i, dna) in denovo::assess_denovo(
        &mother_alleles,
        &father_alleles,
        &child_alleles,
        clip_len,
        parent_quantile,
        aligner,
    )
    .enumerate()
    {
        let output = AlleleResult {
            trid: locus.id.clone(),
            genotype: dna.genotype,
            denovo_coverage: dna.denovo_score,
            allele_coverage: dna.allele_coverage,
            allele_ratio: dna.denovo_score as f64 / dna.allele_coverage as f64,
            child_coverage: dna.child_coverage,
            child_ratio: dna.denovo_score as f64 / dna.child_coverage as f64,
            mean_diff_father: dna.mean_diff_father,
            mean_diff_mother: dna.mean_diff_mother,
            father_dropout_prob,
            mother_dropout_prob,
            allele_origin: dna.allele_origin,
            denovo_status: dna.denovo_status,
            per_allele_reads_father: father_reads
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(","),
            per_allele_reads_mother: mother_reads
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(","),
            per_allele_reads_child: child_reads
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(","),
            index: dna.index,
            father_motif_counts: father_motifs.clone(),
            mother_motif_counts: mother_motifs.clone(),
            child_motif_counts: child_motifs.clone(),
            maxlh: max_lhs[i],
        };
        out_vec.push(output);
    }
    Ok(out_vec)
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_is_homozygous() {
//         let mut header = Header::new();
//         let header_contig_line = r#"##contig=<ID=1,length=10>"#;
//         header.push_record(header_contig_line.as_bytes());
//         let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
//         header.push_record(header_gt_line.as_bytes());
//         header.push_sample("test_sample".as_bytes());
//         let vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
//         let mut record = vcf.empty_record();

//         let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Phased(0)];
//         record.push_genotypes(alleles).unwrap();
//         let genotypes = record.genotypes().unwrap().get(0);
//         // assert!(is_homozygous(&genotypes));

//         record.clear();

//         let alleles = &[GenotypeAllele::Unphased(2), GenotypeAllele::Phased(1)];
//         record.push_genotypes(alleles).unwrap();
//         let genotypes = record.genotypes().unwrap().get(0);
//         // assert!(!is_homozygous(&genotypes));
//     }
// }
