//! Handles module for managing file handles and providing access to BAM and VCF files.
//!
//! This module defines structures and functions to manage file handles for BAM and VCF files,
//! allowing for concurrent access and operations on genomic data.

use crate::util::{self, Result};
use anyhow::anyhow;
use noodles::{bam, bgzf::Reader, sam, vcf};
use std::{
    fs::File,
    path::PathBuf,
    sync::{Arc, Mutex},
};

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

/// A collection of file handles for a trio of family members.
///
/// This struct contains `SubHandle` instances for the mother, father, and child,
/// providing access to their respective BAM and VCF files.
#[derive(Clone)]
pub struct Handles {
    /// The file handle for the mother's genomic data.
    pub mother: SubHandle,
    /// The file handle for the father's genomic data.
    pub father: SubHandle,
    /// The file handle for the child's genomic data.
    pub child: SubHandle,
}

/// Implements the creation of `Handles` instances.
impl Handles {
    /// Creates a new `Handles` instance for a trio of family members.
    ///
    /// # Arguments
    ///
    /// * `mother_prefix` - The file prefix for the mother's genomic data.
    /// * `father_prefix` - The file prefix for the father's genomic data.
    /// * `child_prefix` - The file prefix for the child's genomic data.
    ///
    /// # Returns
    ///
    /// A result containing the new `Handles` instance if successful, or an error if any file handle cannot be created.
    pub fn new(mother_prefix: &str, father_prefix: &str, child_prefix: &str) -> Result<Handles> {
        Ok(Handles {
            mother: SubHandle::new(mother_prefix)?,
            father: SubHandle::new(father_prefix)?,
            child: SubHandle::new(child_prefix)?,
        })
    }
}

/// A sub-handle for managing access to a single family member's BAM and VCF files.
///
/// This struct contains file handles and headers for both BAM and VCF files,
/// allowing for operations on genomic data.
#[derive(Clone)]
pub struct SubHandle {
    /// The indexed VCF reader wrapped in a thread-safe `Arc<Mutex<>>`.
    pub vcf: Arc<Mutex<vcf::IndexedReader<File>>>,
    /// The VCF header information.
    pub vcf_header: Arc<vcf::Header>,
    /// The indexed BAM reader wrapped in a thread-safe `Arc<Mutex<>>`.
    pub bam: Arc<Mutex<bam::IndexedReader<Reader<File>>>>,
    /// The BAM header information.
    pub bam_header: Arc<sam::Header>,
    /// Karyotype string extracted from TRGT BAM (TEMP)
    pub karyotype: Karyotype,
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

/// Implements the creation of `SubHandle` instances.
impl SubHandle {
    /// Creates a new `SubHandle` instance for a single family member.
    ///
    /// # Arguments
    ///
    /// * `prefix` - The file prefix for the family member's genomic data.
    ///
    /// # Returns
    ///
    /// A result containing the new `SubHandle` instance if successful, or an error if the file paths cannot be parsed or file handles cannot be created.
    pub fn new(prefix: &str) -> Result<SubHandle> {
        let paths = build_paths(prefix)?;
        let paths_slice = paths.as_slice();

        if paths_slice.len() < 2 {
            return Err(anyhow!("Failed to parse paths"));
        }

        let bam_path = &paths_slice[0];
        let vcf_path = &paths_slice[1];

        let mut bam = bam::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .map_err(|e| anyhow!("Failed to create bam reader: {}", e))?;

        let bam_header = bam
            .read_header()
            .map_err(|e| anyhow!("Failed to read bam header: {}", e))?;

        let karyotype = extract_karyotype(
            &bam_header
                .programs()
                .get("trgt")
                .and_then(|trgt| trgt.command_line())
                .unwrap_or("XY"),
        );

        let bam = Arc::new(Mutex::new(bam));
        let bam_header = Arc::new(bam_header);

        let mut vcf = vcf::indexed_reader::Builder::default()
            .build_from_path(vcf_path)
            .map_err(|e| anyhow!("Failed to create vcf reader: {}", e))?;

        let vcf_header = vcf
            .read_header()
            .map_err(|e| anyhow!("Failed to read vcf header: {}", e))?;

        let vcf = Arc::new(Mutex::new(vcf));
        let vcf_header = Arc::new(vcf_header);

        Ok(SubHandle {
            vcf,
            vcf_header,
            bam,
            bam_header,
            karyotype,
        })
    }
}
