//! Handles module for managing file handles and providing access to BAM and VCF files.
//!
//! This module defines structures and functions to manage file handles for BAM and VCF files,
//! allowing for concurrent access and operations on genomic data.

use crate::readers::open_vcf_reader;
use crate::util::{self, Result};
use anyhow::anyhow;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};
use std::path::{Path, PathBuf};

/// Represents the input source for a sample's VCF and BAM files.
#[derive(Debug, Clone)]
pub enum SampleInput<'a> {
    /// Prefix-based input: files are derived as `{prefix}.sorted.vcf.gz` and `{prefix}.spanning.sorted.bam`
    Prefix(&'a str),
    /// Explicit file paths for VCF and BAM
    Explicit { vcf: &'a Path, bam: &'a Path },
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

fn get_bam_pg_string(header: &bam::HeaderView, tag: &str) -> Option<String> {
    let header_text = String::from_utf8_lossy(header.as_bytes());
    for line in header_text.lines() {
        if !line.starts_with("@PG") {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        let is_trgt = fields.iter().any(|&f| f == "ID:trgt" || f == "PN:trgt");
        if is_trgt {
            if let Some(cl_field) = fields.iter().find(|&f| f.starts_with(tag)) {
                return Some(cl_field[3..].to_string());
            }
        }
    }
    None
}

fn get_vcf_trgt_version(header: &bcf::header::HeaderView) -> Option<String> {
    let mut trgt_version = None;
    for record in header.header_records().iter() {
        if let bcf::HeaderRecord::Generic { key, value } = record {
            if key == "trgtVersion" {
                trgt_version = Some(value.clone());
            }
        }
    }
    trgt_version
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
    pub vcf: bcf::IndexedReader,
    pub vcf_header: bcf::header::HeaderView,
    pub bam: bam::IndexedReader,
    pub karyotype: Karyotype,
    pub is_trgt_v1: bool,
}

impl SampleLocalData {
    pub fn new(input: SampleInput) -> Result<SampleLocalData> {
        let (bam_path, vcf_path) = match input {
            SampleInput::Prefix(prefix) => {
                let paths = build_paths(prefix)?;
                if paths.len() < 2 {
                    return Err(anyhow!("Failed to parse paths from prefix"));
                }
                (paths[0].clone(), paths[1].clone())
            }
            SampleInput::Explicit { vcf, bam } => (bam.to_path_buf(), vcf.to_path_buf()),
        };

        let bam = bam::IndexedReader::from_path(&bam_path).unwrap_or_else(|e| {
            panic!(
                "Failed to initialize BAM/CRAM reader for path {}: {}",
                bam_path.display(),
                e
            )
        });
        let bam_header = bam.header();

        let vcf = open_vcf_reader(&vcf_path)?;
        let vcf_header = vcf.header().to_owned();

        let vcf_version = get_vcf_trgt_version(&vcf_header).unwrap_or_default();
        let bam_version = get_bam_pg_string(bam_header, "VN:").unwrap_or_default();

        if bam_version != vcf_version {
            return Err(anyhow!(
                "TRGT version mismatch: BAM version is '{}', VCF version is '{}'",
                bam_version,
                vcf_version
            ));
        }
        let is_trgt_v1 = bam_version
            .split('.')
            .next()
            .and_then(|major| major.parse::<u8>().ok())
            .is_some_and(|major| major > 0);

        let bam_string = get_bam_pg_string(bam_header, "CL:").unwrap_or_default();
        let karyotype = extract_karyotype(&bam_string);

        Ok(SampleLocalData {
            vcf,
            vcf_header,
            bam,
            karyotype,
            is_trgt_v1,
        })
    }
}

pub struct TrioLocalData {
    pub mother: SampleLocalData,
    pub father: SampleLocalData,
    pub child: SampleLocalData,
}

impl TrioLocalData {
    pub fn new(
        mother_input: SampleInput,
        father_input: SampleInput,
        child_input: SampleInput,
    ) -> Result<TrioLocalData> {
        Ok(TrioLocalData {
            mother: SampleLocalData::new(mother_input)?,
            father: SampleLocalData::new(father_input)?,
            child: SampleLocalData::new(child_input)?,
        })
    }
}

pub struct DuoLocalData {
    pub sample1: SampleLocalData,
    pub sample2: SampleLocalData,
}

impl DuoLocalData {
    pub fn new(sample1_input: SampleInput, sample2_input: SampleInput) -> Result<DuoLocalData> {
        Ok(DuoLocalData {
            sample1: SampleLocalData::new(sample1_input)?,
            sample2: SampleLocalData::new(sample2_input)?,
        })
    }
}
