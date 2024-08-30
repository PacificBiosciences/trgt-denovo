//! Handles module for managing file handles and providing access to BAM and VCF files.
//!
//! This module defines structures and functions to manage file handles for BAM and VCF files,
//! allowing for concurrent access and operations on genomic data.

use crate::handles::vcf::header::record::value::Collection::Unstructured;
use crate::util::{self, Result};
use anyhow::anyhow;
use noodles::{bam, bgzf::Reader, sam, vcf};
use std::{fs::File, path::PathBuf};

#[derive(Debug, PartialEq, Clone)]
pub struct Karyotype {
    ploidy: PloidyInfo,
}

#[derive(Debug, PartialEq, Clone)]
enum PloidyInfo {
    PresetXX,
    PresetXY,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Ploidy {
    Zero,
    One,
    Two,
}

impl Karyotype {
    pub fn new(encoding: &str) -> Self {
        let ploidy = match encoding {
            "XX" => PloidyInfo::PresetXX,
            "XY" => PloidyInfo::PresetXY,
            _ => PloidyInfo::PresetXY,
        };
        Self { ploidy }
    }

    pub fn get_ploidy(&self, chrom: &str) -> Result<Ploidy> {
        match &self.ploidy {
            PloidyInfo::PresetXX => match chrom {
                "Y" | "chrY" => Ok(Ploidy::Zero),
                _ => Ok(Ploidy::Two),
            },
            PloidyInfo::PresetXY => match chrom {
                "X" | "chrX" | "Y" | "chrY" => Ok(Ploidy::One),
                _ => Ok(Ploidy::Two),
            },
        }
    }
}

/// Builds file paths for BAM and VCF files based on a given prefix.
///
/// # Arguments
///
/// * `prefix` - The prefix string used to construct file paths.
///
/// # Returns
///
/// A result containing a vector of `PathBuf` with the constructed file paths if successful,
/// or an error if the paths do not exist.
pub fn build_paths(prefix: &str) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for ext in &["spanning.sorted.bam", "sorted.vcf.gz"] {
        let mut path = PathBuf::from(prefix);
        if let Some(file_name) = path.file_name() {
            let file_name = file_name.to_string_lossy().to_string();
            path.set_file_name(format!("{}.{}", file_name, ext));
        }
        util::try_exists(&path)?;
        paths.push(path);
    }
    Ok(paths)
}

fn extract_karyotype(input: &str) -> Karyotype {
    const KARYOTYPE_PREFIX: &str = "--karyotype ";
    let karyotype_string = input
        .find(KARYOTYPE_PREFIX)
        .and_then(|start| {
            let remaining = &input[start + KARYOTYPE_PREFIX.len()..];
            remaining.split_whitespace().next()
        })
        .unwrap_or("XY")
        .to_string();
    Karyotype::new(&karyotype_string)
}

pub struct SampleLocalData {
    /// The indexed VCF reader.
    pub vcf: vcf::IndexedReader<File>,
    /// The VCF header information.
    pub vcf_header: vcf::Header,
    /// The indexed BAM reader.
    pub bam: bam::IndexedReader<Reader<File>>,
    /// The BAM header information.
    pub bam_header: sam::Header,
    /// Karyotype string extracted from TRGT BAM (TEMP)
    pub karyotype: Karyotype,
    /// Flag indicating if the TRGT version is 1.x.x or higher.
    pub is_trgt_v1: bool,
}

impl SampleLocalData {
    /// Initializes a new `LocalData` instance for a single family member.
    ///
    /// # Arguments
    ///
    /// * `prefix` - The file prefix for the family member's genomic data.
    ///
    /// # Returns
    ///
    /// A result containing the new `LocalData` instance if successful, or an error if the file paths cannot be parsed or file handles cannot be created.
    pub fn new(prefix: &str) -> Result<SampleLocalData> {
        let paths = build_paths(prefix)?;
        let paths_slice = paths.as_slice();

        if paths_slice.len() < 2 {
            return Err(anyhow!("Failed to parse paths"));
        }

        let bam_path = &paths_slice[0];
        let vcf_path = &paths_slice[1];

        let mut bam = match bam::indexed_reader::Builder::default().build_from_path(bam_path) {
            Ok(reader) => reader,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                return Err(anyhow!("Index not found for BAM: {}", bam_path.display()))
            }
            Err(e) => return Err(anyhow!("Failed to create BAM reader: {}", e)),
        };

        let bam_header = bam
            .read_header()
            .map_err(|e| anyhow!("Failed to read BAM header: {}", e))?;

        let mut vcf = match vcf::indexed_reader::Builder::default().build_from_path(vcf_path) {
            Ok(reader) => reader,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                return Err(anyhow!("Index not found for VCF: {}", vcf_path.display()))
            }
            Err(e) => return Err(anyhow!("Failed to create VCF reader: {}", e)),
        };

        let vcf_header = vcf
            .read_header()
            .map_err(|e| anyhow!("Failed to read VCF header: {}", e))?;

        let is_trgt_v1 = check_trgt_versions(&bam_header, &vcf_header)?;

        let karyotype = extract_karyotype(
            bam_header
                .programs()
                .get("trgt")
                .and_then(|trgt| trgt.command_line())
                .unwrap_or("XY"),
        );

        Ok(SampleLocalData {
            vcf,
            vcf_header,
            bam,
            bam_header,
            karyotype,
            is_trgt_v1,
        })
    }
}

/// Checks if the TRGT versions in the BAM and VCF headers match and returns a boolean flag indicating if the TRGT version is 1.x.x or higher.
///
/// # Arguments
///
/// * `bam_header` - The BAM file header.
/// * `vcf_header` - The VCF file header.
///
/// # Returns
///
/// A result containing a boolean flag if successful, or an error if the TRGT versions do not match.
fn check_trgt_versions(bam_header: &sam::Header, vcf_header: &vcf::Header) -> Result<bool> {
    let trgt_bam_version = bam_header
        .programs()
        .get("trgt")
        .and_then(|trgt| trgt.version())
        .unwrap_or("");

    let trgt_vcf_version = vcf_header
        .other_records()
        .get("trgtVersion")
        .and_then(|record| {
            if let Unstructured(values) = record {
                values.first().map(String::as_str)
            } else {
                None
            }
        })
        .unwrap_or("");

    if trgt_bam_version != trgt_vcf_version {
        return Err(anyhow!(
            "TRGT version mismatch: BAM version is '{}', VCF version is '{}'",
            trgt_bam_version,
            trgt_vcf_version
        ));
    }

    Ok(trgt_bam_version
        .split('.')
        .next()
        .and_then(|major| major.parse::<u8>().ok())
        .map_or(false, |major| major > 0))
}
/// A structure to hold local data for each family member.
pub struct TrioLocalData {
    /// Local data for the mother.
    pub mother: SampleLocalData,
    /// Local data for the father.
    pub father: SampleLocalData,
    /// Local data for the child.
    pub child: SampleLocalData,
}

impl TrioLocalData {
    /// Initializes a new `TrioLocalData` instance for a family trio.
    ///
    /// # Arguments
    ///
    /// * `mother_prefix` - The file prefix for the mother's genomic data.
    /// * `father_prefix` - The file prefix for the father's genomic data.
    /// * `child_prefix` - The file prefix for the child's genomic data.
    ///
    /// # Returns
    ///
    /// A result containing the new `TrioLocalData` instance if successful, or an error if any local data cannot be initialized.
    pub fn new(
        mother_prefix: &str,
        father_prefix: &str,
        child_prefix: &str,
    ) -> Result<TrioLocalData> {
        Ok(TrioLocalData {
            mother: SampleLocalData::new(mother_prefix)?,
            father: SampleLocalData::new(father_prefix)?,
            child: SampleLocalData::new(child_prefix)?,
        })
    }
}

/// A structure to hold local data for two samples.
pub struct DuoLocalData {
    /// Local data for the first sample.
    pub sample1: SampleLocalData,
    /// Local data for the second sample.
    pub sample2: SampleLocalData,
}

impl DuoLocalData {
    /// Initializes a new `DuoLocalData` instance for two samples.
    ///
    /// # Arguments
    ///
    /// * `sample1_prefix` - The file prefix for the first sample's genomic data.
    /// * `sample2_prefix` - The file prefix for the second sample's genomic data.
    ///
    /// # Returns
    ///
    /// A result containing the new `DuoLocalData` instance if successful, or an error if any local data cannot be initialized.
    pub fn new(sample1_prefix: &str, sample2_prefix: &str) -> Result<DuoLocalData> {
        Ok(DuoLocalData {
            sample1: SampleLocalData::new(sample1_prefix)?,
            sample2: SampleLocalData::new(sample2_prefix)?,
        })
    }
}
