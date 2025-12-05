use crate::util::Result;
use anyhow::anyhow;
use rust_htslib::{bcf, bgzf, faidx};
use std::{io::BufReader, path::Path};

pub fn open_vcf_reader(path: &Path) -> Result<bcf::IndexedReader> {
    let vcf = match bcf::IndexedReader::from_path(path) {
        Ok(vcf) => vcf,
        Err(e) => return Err(anyhow!("Failed to open VCF file {}: {}", path.display(), e)),
    };
    Ok(vcf)
}

pub fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let extension = path.extension().unwrap().to_str().unwrap();
    let fai_path = path.with_extension(extension.to_owned() + ".fai");
    if !fai_path.exists() {
        return Err(anyhow!(
            "Reference index file not found: {}. Create it using 'samtools faidx {}'",
            fai_path.display(),
            path.display()
        ));
    }
    faidx::Reader::from_path(path).map_err(|e| anyhow!(e))
}

pub type CatalogReader = BufReader<bgzf::Reader>;
const BUFFER_CAPACITY: usize = 128 * 1024;

pub fn open_catalog_reader(path: &Path) -> Result<CatalogReader> {
    let inner = bgzf::Reader::from_path(path)
        .map_err(|e| anyhow!("Failed to open catalog from {}: {}", path.display(), e))?;
    Ok(BufReader::with_capacity(BUFFER_CAPACITY, inner))
}
