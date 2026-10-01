use crate::util::Result;
use anyhow::{Context, anyhow};
use rust_htslib::{bcf, bgzf, faidx};
use std::{
    io::BufReader,
    path::{Path, PathBuf},
};

pub fn open_vcf_reader(path: &Path) -> Result<bcf::IndexedReader> {
    let vcf = match bcf::IndexedReader::from_path(path) {
        Ok(vcf) => vcf,
        Err(e) => return Err(anyhow!("Failed to open VCF file {}: {}", path.display(), e)),
    };
    Ok(vcf)
}

pub fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let mut index_name = path.as_os_str().to_owned();
    index_name.push(".fai");
    let fai_path = PathBuf::from(index_name);
    if !fai_path.exists() {
        return Err(anyhow!(
            "Reference index file not found: {}. Create it using 'samtools faidx {}'",
            fai_path.display(),
            path.display()
        ));
    }
    faidx::Reader::from_path(path)
        .with_context(|| format!("Failed to open reference genome: {}", path.display()))
}

pub type CatalogReader = BufReader<bgzf::Reader>;
const BUFFER_CAPACITY: usize = 128 * 1024;

pub fn open_catalog_reader(path: &Path) -> Result<CatalogReader> {
    let inner = bgzf::Reader::from_path(path)
        .map_err(|e| anyhow!("Failed to open catalog from {}: {}", path.display(), e))?;
    Ok(BufReader::with_capacity(BUFFER_CAPACITY, inner))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genome_reader_supports_extensionless_fasta() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("reference");
        std::fs::write(&genome, b">chr1\nACGT\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fai"),
            b"chr1\t4\t6\t4\t5\n",
        )
        .unwrap();

        let reader = open_genome_reader(&genome).unwrap();
        assert_eq!(reader.fetch_seq("chr1", 0, 3).unwrap(), b"ACGT");
    }
}
