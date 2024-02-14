//! Provides functionality for handling alleles and related operations.
//!
//! This module includes structures and functions to represent alleles, load them from VCF and BAM
//! files, perform read alignments, and process allele information for further analysis.
//!
use crate::{
    aligner::{AlignmentStatus, WFAligner},
    denovo::{self, AlleleOrigin, DenovoStatus},
    handles::{self, Karyotype, Ploidy},
    locus::Locus,
    math,
    read::ReadInfo,
    snp::{self},
    util::{Params, Result},
};
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
    fmt,
    fs::File,
    ops::Index,
    sync::{Arc, Mutex},
};

#[derive(Debug, PartialEq)]
pub struct AlleleSet {
    pub alleles: Vec<Allele>,
    pub hp_counts: [usize; 3],
}

impl AlleleSet {
    pub fn iter(&self) -> std::slice::Iter<Allele> {
        self.alleles.iter()
    }

    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.len() == 0
    }

    pub fn get_naive_dropout(&self, locus: &Locus, karyotype: &Karyotype) -> DropoutType {
        const MIN_COVERAGE: usize = 2;

        let hp1_reads = self.hp_counts[0];
        let hp2_reads = self.hp_counts[1];
        let unphased_reads = self.hp_counts[2];

        let phased_reads = hp1_reads + hp2_reads;
        let total_reads = phased_reads + unphased_reads;

        let ploidy = karyotype.get_ploidy(locus.region.name()).unwrap();

        if total_reads < MIN_COVERAGE {
            return DropoutType::FullDropout;
        } else if (unphased_reads < MIN_COVERAGE)
            & (ploidy != Ploidy::One)
            & (hp1_reads < MIN_COVERAGE || hp2_reads < MIN_COVERAGE)
        {
            return DropoutType::HaplotypeDropout;
        }
        DropoutType::None
    }
}

#[derive(Debug, PartialEq, Serialize)]
pub enum DropoutType {
    FullDropout,
    HaplotypeDropout,
    None,
}

impl fmt::Display for DropoutType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DropoutType::FullDropout => write!(f, "FD"),
            DropoutType::HaplotypeDropout => write!(f, "HD"),
            DropoutType::None => write!(f, "N"),
        }
    }
}

impl IntoIterator for AlleleSet {
    type Item = Allele;
    type IntoIter = std::vec::IntoIter<Allele>;

    fn into_iter(self) -> Self::IntoIter {
        self.alleles.into_iter()
    }
}

impl<'a> IntoIterator for &'a AlleleSet {
    type Item = &'a Allele;
    type IntoIter = std::slice::Iter<'a, Allele>;

    fn into_iter(self) -> Self::IntoIter {
        self.alleles.iter()
    }
}

impl Index<usize> for AlleleSet {
    type Output = Allele;

    fn index(&self, index: usize) -> &Self::Output {
        &self.alleles[index]
    }
}

/// Allele representation
#[derive(Debug, PartialEq)]
pub struct Allele {
    /// The sequence of the allele
    pub seq: Vec<u8>,
    /// TRGT genotype index of the allele.
    pub genotype: usize,
    /// Vector of alignments of reads supporting this allele, storing the read and alignment score
    pub read_aligns: Vec<(ReadInfo, i32)>, // TODO: Do we need the score here?
    /// TRGT motif count in this allele
    pub motif_count: String,
    /// Index used to identify the allele
    pub index: usize,
}

/// Implementation of `Allele` struct.
///
/// Provides methods for creating and manipulating `Allele` instances.
impl Allele {
    /// Creates a dummy `Allele` with given reads and default values for other fields.
    ///
    /// # Arguments
    ///
    /// * `reads` - A vector of `ReadInfo` instances representing the reads.
    ///
    /// # Returns
    ///
    /// Returns an `Allele` instance with default values and provided reads.
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

/// Static array of filters to be applied to reads.
///
/// These filters are used to preprocess reads before assigning them to alleles. Filters can include
/// distance-based and frequency-based criteria to ensure the quality and relevance of the reads.
pub static FILTERS: &[&'static (dyn snp::ReadFilter + Sync)] =
    &[&snp::FilterByDist, &snp::FilterByFreq];

/// Calculates the haplotype counts from a vector of reads.
///
/// # Arguments
///
/// * `reads` - A vector of `ReadInfo` instances representing the reads.
///
/// # Returns
///
/// Returns an array containing the counts of reads for each haplotype.
pub fn calculate_hp_counts(reads: &[ReadInfo]) -> [usize; 3] {
    let mut hp_counts = [0; 3];
    for read in reads {
        match read.haplotype {
            Some(1) => hp_counts[0] += 1,
            Some(2) => hp_counts[1] += 1,
            None => hp_counts[2] += 1,
            _ => {}
        }
    }
    hp_counts
}

/// Loads alleles from VCF and BAM files for a given locus and subhandle.
///
/// Extracts allele sequences from a VCF file and retrieves corresponding reads from a
/// BAM file. It applies read filters and assigns reads to the closest matching allele based on
/// alignment.
///
/// # Arguments
///
/// * `locus` - A reference to the `Locus` for which alleles are to be loaded.
/// * `subhandle` - A reference to the `SubHandle` containing file handles.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// A result containing a vector of `Allele` instances if successful, or an error if not.
pub fn load_alleles(
    locus: &Locus,
    subhandle: &handles::SubHandle,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<AlleleSet> {
    let (allele_seqs, motif_counts, genotype_indices) =
        get_allele_seqs(locus, &subhandle.vcf, &subhandle.vcf_header)?;
    let mut reads = get_reads(&subhandle.bam, &subhandle.bam_header, locus)?;
    snp::apply_read_filters(&mut reads, FILTERS);

    let hp_counts = calculate_hp_counts(&reads);

    let reads_by_allele = if params.partition_by_alignment {
        assign_reads_by_alignment(&allele_seqs, reads, params.clip_len, aligner)
    } else {
        assign_reads_by_classification(&allele_seqs, reads)
    };

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

    Ok(AlleleSet { alleles, hp_counts })
}

/// Assigns reads to the closest matching allele based on alignment scores.
///
/// This function aligns reads to each allele sequence and assigns each read to the allele with the
/// highest alignment score. In case of ties, it distributes reads evenly.
///
/// # Arguments
///
/// * `alleles` - A slice of allele sequences.
/// * `reads` - A vector of `ReadInfo` instances representing the reads.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// Returns a vector of vectors, each containing tuples of `ReadInfo` and alignment scores.
fn assign_reads_by_alignment(
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

fn assign_reads_by_classification(
    alleles: &[Vec<u8>],
    reads: Vec<ReadInfo>,
) -> Vec<Vec<(ReadInfo, i32)>> {
    let mut reads_by_allele = vec![Vec::new(); alleles.len()];
    for read in reads {
        let al = read.classification.unwrap() as usize;
        reads_by_allele[al].push((read, 0));
    }
    reads_by_allele
}

/// Retrieves reads from a BAM file for a given locus.
///
/// This function queries a BAM file to extract reads that overlap with the genomic region of the
/// given locus. It then converts the BAM records into `ReadInfo` instances for further processing.
///
/// # Arguments
///
/// * `bam` - A reference to an `Arc<Mutex<bam::IndexedReader<Reader<File>>>>`.
/// * `bam_header` - A reference to an `Arc<sam::Header>`.
/// * `locus` - A reference to the `Locus` for which reads are to be retrieved.
///
/// # Returns
///
/// A result containing a vector of `ReadInfo` instances if successful, or an error if not.
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

    if reads.is_empty() {
        Err(anyhow!("no reads found for locus: {}", locus.id))
    } else {
        Ok(reads)
    }
}

// TODO: improve, for now it just pulls out the genotype index and checks duplicates
// fn is_homozygous(genotypes: &Genotype) -> bool {
//     let alleles: Vec<_> = genotypes
//         .iter()
//         .filter_map(|allele| allele.position())
//         .collect();
//     alleles.into_iter().unique().count() == 1
// }

// Define custom tags for lookup in TRGT VCFS
static TRID_KEY: Lazy<field::Key> = Lazy::new(|| "TRID".parse().unwrap());
static MC_KEY: Lazy<genotypes::keys::Key> = Lazy::new(|| "MC".parse().unwrap());

/// Extracts allele sequences and associated information from a VCF file for a given locus.
///
/// Reads a VCF file to obtain allele sequences, motif counts, and genotype indices
/// for a specific locus. It ensures that the alleles correspond to the locus ID and processes
/// genotype information to construct allele sequences.
///
/// # Arguments
///
/// * `locus` - A reference to the `Locus` for which allele sequences are to be extracted.
/// * `vcf` - A reference to an `Arc<Mutex<vcf::IndexedReader<File>>>`.
/// * `vcf_header` - A reference to an `Arc<vcf::Header>`.
///
/// # Returns
///
/// A result containing a tuple with allele sequences, motif counts, and genotype indices if successful, or an error if not.
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

/// Processes genotype information from a VCF record to generate allele sequences.
///
/// # Arguments
///
/// * `genotype` - A reference to the `Genotype` from the VCF record.
/// * `record` - A reference to the VCF `Record`.
/// * `locus` - A reference to the `Locus` associated with the alleles.
///
/// # Returns
///
/// A result containing a tuple with al0lele sequences and genotype indices if successful, or an error if not.
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
    }
    Ok((alleles, genotype_indices))
}

/// Represents the result of allele processing, including various statistics and classifications.
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
    pub father_dropout: String,
    pub mother_dropout: String,
    pub child_dropout: String,
    pub index: usize,
    pub father_motif_counts: String,
    pub mother_motif_counts: String,
    pub child_motif_counts: String,
    // #[serde(serialize_with = "serialize_with_precision")]
    // pub maxlh: f64,
}

/// Serialize a value as a display string.
///
/// Serialize values for JSON output, ensuring that they are represented
/// as human-readable strings.
///
/// # Arguments
///
/// * `value` - A reference to the value to serialize.
/// * `serializer` - The serializer to use.
///
/// # Returns
///
/// Returns a serialized string representation of the value.
fn serialize_as_display<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.collect_str(value)
}

/// Serialize a value with a specified precision as a string.
///
/// Formats numerical values with a fixed precision before serialization, ensuring
/// consistent representation in the output.
///
/// # Arguments
///
/// * `value` - A reference to the value to serialize.
/// * `serializer` - The serializer to use.
///
/// # Returns
///
/// Returns a serialized string representation of the value with precision.
fn serialize_with_precision<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    let formatted_value = format!("{:.4}", value);
    serializer.serialize_str(&formatted_value)
}

/// Loads alleles for a specific family member role (father, mother, or child).
///
/// Wrapper around `load_alleles` that handles errors and logs warnings
/// specific to the family member role being processed.
///
/// # Arguments
///
/// * `role` - A character representing the family member role ('F', 'M', or 'C').
/// * `locus` - A reference to the `Locus` for which alleles are to be loaded.
/// * `subhandle` - A reference to the `SubHandle` containing file handles.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
fn load_alleles_handle(
    role: char,
    locus: &Locus,
    subhandle: &handles::SubHandle,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<AlleleSet> {
    load_alleles(locus, subhandle, params, aligner).map_err(|err| {
        log::warn!("Skipping {} in {}", err, role);
        err
    })
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
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `parent_quantile` - The quantile used for parental allele frequency calculations.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// A result containing a vector of `AlleleResult` instances if successful, or an error if not.
pub fn process_alleles(
    locus: &Locus,
    handle: Arc<handles::Handles>,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<Vec<AlleleResult>> {
    let father_alleles = load_alleles_handle('F', locus, &handle.father, params, aligner)?;
    let mother_alleles = load_alleles_handle('M', locus, &handle.mother, params, aligner)?;
    let child_alleles = load_alleles_handle('C', locus, &handle.child, params, aligner)?;

    // TODO: ongoing work, the maximum likelihood is obtained naively
    // let max_lhs = TrinaryMatrix::new(&child_alleles, &father_alleles, &mother_alleles)
    //     .and_then(|trinary_mat| snp::inheritance_prob(&trinary_mat))
    //     .map(|(_inherit_p, max_lh)| max_lh)
    //     .unwrap_or((-1.0, -1.0));
    // let max_lhs = [max_lhs.0, max_lhs.1];

    let mother_dropout_prob = math::get_dropout_prob(&mother_alleles);
    let father_dropout_prob = math::get_dropout_prob(&father_alleles);

    let father_reads = math::get_per_allele_reads(&father_alleles);
    let mother_reads = math::get_per_allele_reads(&mother_alleles);
    let child_reads = math::get_per_allele_reads(&child_alleles);

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
    for (_i, dna) in denovo::assess_denovo(
        &mother_alleles,
        &father_alleles,
        &child_alleles,
        params,
        aligner,
    )
    .enumerate()
    {
        let output = AlleleResult {
            trid: locus.id.clone(),
            genotype: dna.genotype,
            denovo_coverage: dna.denovo_coverage,
            allele_coverage: dna.allele_coverage,
            allele_ratio: dna.denovo_coverage as f64 / dna.allele_coverage as f64,
            child_coverage: dna.child_coverage,
            child_ratio: dna.denovo_coverage as f64 / dna.child_coverage as f64,
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
            father_dropout: father_alleles
                .get_naive_dropout(&locus, &handle.father.karyotype)
                .to_string(),
            mother_dropout: mother_alleles
                .get_naive_dropout(&locus, &handle.mother.karyotype)
                .to_string(),
            child_dropout: child_alleles
                .get_naive_dropout(&locus, &handle.child.karyotype)
                .to_string(),
            index: dna.index,
            father_motif_counts: father_motifs.clone(),
            mother_motif_counts: mother_motifs.clone(),
            child_motif_counts: child_motifs.clone(),
            // maxlh: max_lhs[i],
        };
        out_vec.push(output);
    }
    Ok(out_vec)
}
