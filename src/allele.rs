use crate::{
    aligner::{AlignmentStatus, WFAligner},
    handles::{Karyotype, Ploidy, SampleLocalData},
    locus::Locus,
    read::ReadInfo,
    snp::{self},
    util::{Params, Result},
};
use anyhow::anyhow;
use core::fmt;
use itertools::Itertools;
use noodles::{
    bam,
    bgzf::Reader,
    sam,
    vcf::{
        self,
        record::{
            genotypes::{self, sample::value::Genotype},
            info::field::{self},
            Record,
        },
    },
};
use once_cell::sync::Lazy;
use serde::Serialize;
use std::{fs::File, ops::Index};

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
pub fn serialize_as_display<T: std::fmt::Display, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.collect_str(value)
}

/// Serialize a value with a specified precision as a string.
///
/// # Arguments
///
/// * `value` - A reference to the value to serialize.
/// * `serializer` - The serializer to use.
///
/// # Returns
///
/// Returns a serialized string representation of the value with precision.
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
    pub read_aligns: Vec<(ReadInfo, i32)>,
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
            allele_length: String::from("0"),
            index: 0,
        }
    }
}

pub fn join_allele_attribute(
    alleles: &AlleleSet,
    attribute: impl Fn(&Allele) -> &String,
) -> String {
    alleles
        .iter()
        .map(attribute)
        .map(|s| s.to_string())
        .collect::<Vec<String>>()
        .join(",")
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

    let mut reads = get_reads(&mut subhandle.bam, &subhandle.bam_header, locus)?;
    snp::apply_read_filters(&mut reads, FILTERS);

    let hp_counts = calculate_hp_counts(&reads);

    let reads_by_allele = if params.partition_by_alignment {
        assign_reads_by_alignment(&vcf_data.allele_seqs, reads, params.clip_len, aligner)
    } else {
        assign_reads_by_classification(&vcf_data.allele_seqs, reads)
    };

    let alleles = vcf_data
        .allele_seqs
        .into_iter()
        .enumerate()
        .map(|(index, allele_seq)| Allele {
            seq: allele_seq,
            genotype: vcf_data.genotype_indices[index],
            read_aligns: reads_by_allele[index].to_owned(),
            motif_count: vcf_data.motif_counts[index].to_owned(),
            allele_length: vcf_data.allele_lengths[index].to_owned(),
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
    bam: &mut bam::IndexedReader<Reader<File>>,
    bam_header: &sam::Header,
    locus: &Locus,
) -> Result<Vec<ReadInfo>> {
    let query = bam.query(bam_header, &locus.region)?;
    let reads: Vec<ReadInfo> = query
        .filter_map(|record| match record {
            Ok(r) => match ReadInfo::new(r, locus) {
                Ok(Some(read_info)) => Some(Ok(read_info)),
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            },
            Err(e) => Some(Err(e.into())),
        })
        .collect::<Result<Vec<ReadInfo>>>()?;

    if reads.is_empty() {
        Err(anyhow!("no reads found for locus: {}", locus.id))
    } else {
        Ok(reads)
    }
}

// Define custom tags for lookup in TRGT VCFs
static TRID_KEY: Lazy<field::Key> = Lazy::new(|| "TRID".parse().unwrap());
static MC_KEY: Lazy<genotypes::keys::Key> = Lazy::new(|| "MC".parse().unwrap());
static AL_KEY: Lazy<genotypes::keys::Key> = Lazy::new(|| "AL".parse().unwrap());

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
/// * `is_trgt_v1` - Flag if TRGT v1.0 is used, if it is, skip the padding base in the allele sequences.
///
/// # Returns
///
/// A result containing VcfData if successful, or an error if not.
fn get_vcf_data(
    locus: &Locus,
    vcf: &mut vcf::IndexedReader<File>,
    vcf_header: &vcf::Header,
    is_trgt_v1: bool,
) -> Result<VcfData> {
    let locus_id = &locus.id;
    let query = match vcf.query(vcf_header, &locus.region) {
        Ok(query) => query,
        Err(e) => {
            if e.kind() == std::io::ErrorKind::InvalidInput {
                return Err(anyhow!("TRID={} missing", locus_id));
            } else {
                return Err(anyhow!("Error querying {} in VCF: {}", locus_id, e));
            }
        }
    };

    for result in query {
        let record = result?;
        let info = record.info();
        let trid = match info.get(&*TRID_KEY) {
            Some(Some(field::Value::String(s))) => s.to_string(),
            _ => {
                log::trace!(
                    "Record missing TRID or TRID is not a string, skipping: {:?}",
                    record
                );
                continue;
            }
        };

        if &trid == locus_id {
            let format = record.genotypes().get_index(0).unwrap();

            let genotype = match format.genotype() {
                Some(Ok(genotype)) => genotype,
                _ => return Err(anyhow!("TRID={} misses genotyping", locus_id)),
            };

            let (allele_seqs, genotype_indices) =
                process_genotypes(&genotype, &record, locus, is_trgt_v1 as usize)?;

            let mc_field = format.get(&*MC_KEY).unwrap().unwrap().to_string();
            let motif_counts = mc_field.split(',').map(ToString::to_string).collect();

            let al_field = format.get(&*AL_KEY).unwrap().unwrap().to_string();
            let allele_lengths = al_field.split(',').map(ToString::to_string).collect();

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

/// Processes genotype information from a VCF record to generate allele sequences.
///
/// # Arguments
///
/// * `genotype` - A reference to the `Genotype` from the VCF record.
/// * `record` - A reference to the VCF `Record`.
/// * `locus` - A reference to the `Locus` associated with the alleles.
/// * `is_trgt_v1` - Flag if TRGT v1.0 is used, if it is, skip the padding base in the allele sequences.
///
/// # Returns
///
/// A result containing a tuple with allele sequences and genotype indices if successful, or an error if not.
fn process_genotypes(
    genotype: &Genotype,
    record: &Record,
    locus: &Locus,
    is_trgt_v1: usize,
) -> Result<(Vec<Vec<u8>>, Vec<usize>)> {
    let mut reference_bases = record.reference_bases().to_string().into_bytes();
    reference_bases = reference_bases.into_iter().skip(is_trgt_v1).collect();

    let alternate_bases = record
        .alternate_bases()
        .iter()
        .map(|base| {
            let mut base_bytes = base.to_string().into_bytes();
            base_bytes = base_bytes.into_iter().skip(is_trgt_v1).collect();
            base_bytes
        })
        .collect::<Vec<Vec<u8>>>();

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
