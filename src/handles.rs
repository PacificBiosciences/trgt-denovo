use crate::util::{self, Result};
use anyhow::anyhow;
use noodles::{bam, bgzf::Reader, sam, vcf};
use std::{
    fs::File,
    path::PathBuf,
    sync::{Arc, Mutex},
};

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

#[derive(Clone)]
pub struct Handles {
    pub mother: SubHandle,
    pub father: SubHandle,
    pub child: SubHandle,
}

impl Handles {
    pub fn new(mother_prefix: &str, father_prefix: &str, child_prefix: &str) -> Result<Handles> {
        Ok(Handles {
            mother: SubHandle::new(mother_prefix)?,
            father: SubHandle::new(father_prefix)?,
            child: SubHandle::new(child_prefix)?,
        })
    }
}

#[derive(Clone)]
pub struct SubHandle {
    pub vcf: Arc<Mutex<vcf::IndexedReader<File>>>,
    pub vcf_header: Arc<vcf::Header>,
    pub bam: Arc<Mutex<bam::IndexedReader<Reader<File>>>>,
    pub bam_header: Arc<sam::Header>,
}

impl SubHandle {
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
        })
    }
}
