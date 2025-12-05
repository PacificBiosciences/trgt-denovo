use crate::{
    aligner::{AlignmentStatus, WFAligner},
    handles::{Karyotype, Ploidy, SampleLocalData},
    locus::Locus,
    model::Params,
    read::TrgtRead,
    snp::{self},
    util::Result,
};
use anyhow::anyhow;
use core::fmt;
use itertools::{izip, Itertools};
use rust_htslib::{
    bam::{self, Read as BamRead, Record},
    bcf::{
        self,
        header::HeaderView,
        record::GenotypeAllele::{PhasedMissing, UnphasedMissing},
        Read as VcfRead,
    },
};
use serde::Serialize;
use std::{ops::Index, str};

pub fn serialize_as_display<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.collect_str(value)
}

pub fn serialize_with_precision<
    T: std::fmt::Display + PartialOrd + Into<f64>,
    S: serde::Serializer,
>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    let mut formatted_value = format!("{:.1}", value);
    // Check if the value ends with ".0" after formatting with one decimal place.
    // If not, attempt to format with up to three decimal places without trailing zeros.
    if !formatted_value.ends_with(".0") {
        formatted_value = format!("{:.3}", value).trim_end_matches('0').to_string();
        // If the result ends with a decimal point, add a '0' to conform to the desired format.
        if formatted_value.ends_with('.') {
            formatted_value.push('0');
        }
    }
    serializer.serialize_str(&formatted_value)
}

/// Allele representation
#[derive(Debug, PartialEq, Clone)]
pub struct Allele {
    /// Vector of alignments of reads supporting this allele, storing the read and alignment score
    pub read_aligns: Vec<(TrgtRead, i32)>,
    /// The sequence of the allele
    pub seq: Vec<u8>,
    /// TRGT motif count in this allele
    pub motif_count: String,
    /// TRGT allele length in this allele
    pub allele_length: String,
    /// TRGT genotype value of the allele.
    pub genotype: usize,
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
    /// * `reads` - A vector of `TrgtRead` instances representing the reads.
    ///
    /// # Returns
    ///
    /// Returns an `Allele` instance with default values and provided reads.
    pub fn dummy(reads: Vec<TrgtRead>) -> Allele {
        Allele {
            seq: vec![],
            genotype: 0,
            read_aligns: reads.into_iter().map(|info| (info, 0i32)).collect_vec(),
            motif_count: String::from("0"),
            allele_length: String::from("0"),
            index: 0,
        }
    }
}

pub fn join_allele_attribute(
    alleles: &AlleleSet,
    attribute: impl Fn(&Allele) -> &String,
) -> String {
    alleles.iter().map(attribute).join(",")
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
pub fn calculate_hp_counts(reads: &[TrgtRead]) -> [usize; 3] {
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
/// * `params` - Parameters used downstream.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// A result containing a vector of `Allele` instances if successful, or an error if not.
pub fn load_alleles(
    locus: &Locus,
    subhandle: &mut SampleLocalData,
    params: &Params,
    aligner: &mut WFAligner,
) -> Result<AlleleSet> {
    let vcf_data = get_vcf_data(
        locus,
        &mut subhandle.vcf,
        &subhandle.vcf_header,
        subhandle.is_trgt_v1,
    )?;

    let mut reads = get_reads(&mut subhandle.bam, locus)?;
    snp::apply_read_filters(&mut reads, FILTERS);

    let hp_counts = calculate_hp_counts(&reads);

    let reads_by_allele = if params.partition_by_alignment {
        assign_reads_by_alignment(&vcf_data.allele_seqs, reads, params.clip_len, aligner)
    } else {
        assign_reads_by_classification(&vcf_data.allele_seqs, reads)
    };

    let alleles = izip!(
        vcf_data.allele_seqs,
        vcf_data.genotype_indices,
        reads_by_allele,
        vcf_data.motif_counts,
        vcf_data.allele_lengths,
    )
    .enumerate()
    .map(
        |(index, (allele_seq, genotype, read_aligns, motif_count, allele_length))| Allele {
            seq: allele_seq,
            genotype,
            read_aligns,
            motif_count,
            allele_length,
            index,
        },
    )
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
/// * `reads` - A vector of `TrgtRead` instances representing the reads.
/// * `clip_len` - The length of the clipping to be applied to alignments.
/// * `aligner` - A mutable reference to the `WFAligner` for performing alignments.
///
/// # Returns
///
/// Returns a vector of vectors, each containing tuples of `ReadInfo` and alignment scores.
fn assign_reads_by_alignment(
    alleles: &[Vec<u8>],
    reads: Vec<TrgtRead>,
    clip_len: usize,
    aligner: &mut WFAligner,
) -> Vec<Vec<(TrgtRead, i32)>> {
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
    reads: Vec<TrgtRead>,
) -> Vec<Vec<(TrgtRead, i32)>> {
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
/// * `bam` - A reference to an `bam::IndexedReader`.
/// * `locus` - A reference to the `Locus` for which reads are to be retrieved.
///
/// # Returns
///
/// A result containing a vector of `ReadInfo` instances if successful, or an error if not.
pub fn get_reads(bam: &mut bam::IndexedReader, locus: &Locus) -> Result<Vec<TrgtRead>> {
    bam.fetch((
        locus.region.contig.as_str(),
        locus.region.start,
        locus.region.end,
    ))
    .map_err(|e| anyhow!("Error querying in {} in BAM: {}", locus.id, e))?;

    let mut reads = Vec::new();
    let mut record = Record::new();
    while let Some(()) = bam.read(&mut record).transpose()? {
        if let Some(trgt_read) = TrgtRead::new(&record, locus)? {
            reads.push(trgt_read);
        }
    }

    if reads.is_empty() {
        Err(anyhow!("No reads found for locus: {}", locus.id))
    } else {
        Ok(reads)
    }
}

/// Extracts allele sequences and associated information from a VCF file for a given locus.
///
/// Reads a VCF file to obtain allele sequences, motif counts, and genotype indices
/// for a specific locus. It ensures that the alleles correspond to the locus ID and processes
/// genotype information to construct allele sequences.
///
/// # Arguments
///
/// * `locus` - A reference to the `Locus` for which allele sequences are to be extracted.
/// * `vcf` - A reference to an `bcf::IndexedReader`.
/// * `vcf_header` - A reference to an `bcf::HeaderView`.
/// * `is_trgt_v1` - Flag if TRGT v1.0 is used, if it is, skip the padding base in the allele sequences.
///
/// # Returns
///
/// A result containing VcfData if successful, or an error if not.
fn get_vcf_data(
    locus: &Locus,
    vcf: &mut bcf::IndexedReader,
    vcf_header: &HeaderView,
    is_trgt_v1: bool,
) -> Result<VcfData> {
    let locus_id = &locus.id;
    // TODO: make part of locus
    let rid = vcf_header
        .name2rid(locus.region.contig.as_bytes())
        .map_err(|_| anyhow!("Contig {} not found in VCF header", locus.region.contig))?;

    vcf.fetch(
        rid,
        locus.region.start as u64,
        Some(locus.region.end as u64),
    )
    .map_err(|e| anyhow!("Error querying {} in VCF: {}", locus_id, e))?;

    let mut record = vcf.empty_record();
    while let Some(()) = vcf.read(&mut record).transpose()? {
        let trid_matches = if let Ok(Some(trid_field)) = record.info(b"TRID").string() {
            let vcf_trid_parts = trid_field.iter().map(|s| str::from_utf8(s).ok());
            let locus_id_parts = locus_id.split(',').map(Some);
            vcf_trid_parts.eq(locus_id_parts)
        } else {
            false
        };

        if trid_matches {
            let (allele_seqs, genotype_indices) =
                process_genotypes(&record, locus, is_trgt_v1 as usize, locus_id)?;

            let motif_counts: Vec<String> = record
                .format(b"MC")
                .string()?
                .first()
                .and_then(|bytes| str::from_utf8(bytes).ok())
                .map(|s| s.split(',').map(str::to_string).collect())
                .ok_or_else(|| anyhow!("Missing or malformed MC field for TRID={}", locus_id))?;

            let allele_lengths: Vec<String> = record
                .format(b"AL")
                .integer()?
                .first()
                .map(|int_slice| int_slice.iter().map(|num| num.to_string()).collect())
                .ok_or_else(|| anyhow!("Missing or malformed AL field for TRID={}", locus_id))?;

            return Ok(VcfData {
                allele_seqs,
                motif_counts,
                allele_lengths,
                genotype_indices,
            });
        }
    }
    Err(anyhow!("TRID={} missing", locus_id))
}

fn process_genotypes(
    record: &bcf::Record,
    locus: &Locus,
    is_trgt_v1: usize,
    trid: &str,
) -> Result<(Vec<Vec<u8>>, Vec<usize>)> {
    let genotype = record
        .genotypes()
        .map_err(|_| anyhow!("Missing genotype for TRID={}", trid))?
        .get(0);
    if genotype[0] == UnphasedMissing || genotype[0] == PhasedMissing {
        return Err(anyhow!("Missing genotype for TRID={}", trid));
    }

    let vcf_allele_seqs: Vec<Vec<u8>> = record
        .alleles()
        .iter()
        .map(|allele| allele.iter().skip(is_trgt_v1).copied().collect())
        .collect();

    let genotype_indices: Vec<usize> = genotype
        .iter()
        .filter_map(|allele| allele.index())
        .map(|idx| idx as usize)
        .collect();

    Ok(build_allele_seqs(
        &vcf_allele_seqs,
        genotype_indices,
        &locus.left_flank,
        &locus.right_flank,
    ))
}

fn build_allele_seqs(
    vcf_allele_seqs: &[Vec<u8>],
    mut genotype_indices: Vec<usize>,
    left_flank: &[u8],
    right_flank: &[u8],
) -> (Vec<Vec<u8>>, Vec<usize>) {
    // Sort the indices to ensure phased genotypes are handled consistently.
    genotype_indices.sort_unstable();

    let mut seq = left_flank.to_vec();
    let mut alleles = Vec::new();

    for &allele_index in &genotype_indices {
        let allele_seq = &vcf_allele_seqs[allele_index];
        seq.extend_from_slice(allele_seq);
        seq.extend_from_slice(right_flank);
        alleles.push(seq.clone());
        seq.truncate(left_flank.len());
    }
    (alleles, genotype_indices)
}

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

        let ploidy = karyotype.get_ploidy(&locus.region.contig).unwrap();

        if total_reads < MIN_COVERAGE {
            DropoutType::FullDropout
        } else if (unphased_reads < MIN_COVERAGE)
            && (ploidy != Ploidy::One)
            && (hp1_reads < MIN_COVERAGE || hp2_reads < MIN_COVERAGE)
        {
            DropoutType::HaplotypeDropout
        } else {
            DropoutType::None
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

/// Represents data extracted from a VCF file for a given locus.
#[derive(Debug)]
pub struct VcfData {
    pub allele_seqs: Vec<Vec<u8>>,
    pub motif_counts: Vec<String>,
    pub allele_lengths: Vec<String>,
    pub genotype_indices: Vec<usize>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_allele_seqs_sorts_indices() {
        let vcf_allele_seqs = vec![b"REF".to_vec(), b"ALT1".to_vec(), b"ALT2".to_vec()];
        let left_flank = b"LF";
        let right_flank = b"RF";

        // 1|0 == 0|1
        let (alleles_1_0, indices_1_0) =
            build_allele_seqs(&vcf_allele_seqs, vec![1, 0], left_flank, right_flank);
        let (alleles_0_1, indices_0_1) =
            build_allele_seqs(&vcf_allele_seqs, vec![0, 1], left_flank, right_flank);

        assert_eq!(alleles_1_0, alleles_0_1);
        assert_eq!(indices_1_0, indices_0_1);
        assert_eq!(indices_1_0, vec![0, 1]);

        assert_eq!(alleles_1_0[0], b"LFREFRF".to_vec());
        assert_eq!(alleles_1_0[1], b"LFALT1RF".to_vec());
    }

    #[test]
    fn test_build_allele_seqs_with_homozygous() {
        let vcf_allele_seqs = vec![b"REF".to_vec(), b"ALT1".to_vec()];
        let left_flank = b"LF";
        let right_flank = b"RF";

        // 1/1 == two identical alleles
        let (alleles, indices) =
            build_allele_seqs(&vcf_allele_seqs, vec![1, 1], left_flank, right_flank);

        assert_eq!(indices, vec![1, 1]);
        assert_eq!(alleles.len(), 2);
        assert_eq!(alleles[0], alleles[1]);
        assert_eq!(alleles[0], b"LFALT1RF".to_vec());
    }

    #[test]
    fn test_build_allele_seqs_with_three_alleles() {
        let vcf_allele_seqs = vec![b"REF".to_vec(), b"ALT1".to_vec(), b"ALT2".to_vec()];
        let left_flank = b"LF";
        let right_flank = b"RF";

        // 2|0 == 0|2
        let (alleles_2_0, indices_2_0) =
            build_allele_seqs(&vcf_allele_seqs, vec![2, 0], left_flank, right_flank);
        let (alleles_0_2, indices_0_2) =
            build_allele_seqs(&vcf_allele_seqs, vec![0, 2], left_flank, right_flank);

        assert_eq!(alleles_2_0, alleles_0_2);
        assert_eq!(indices_2_0, indices_0_2);
        assert_eq!(indices_2_0, vec![0, 2]);
    }
}
