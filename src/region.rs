use crate::util::Result;
use anyhow::anyhow;

#[derive(Debug, PartialEq)]
pub struct GenomicRegion {
    pub contig: String,
    pub start: u32,
    pub end: u32,
}

impl GenomicRegion {
    pub fn new(contig: impl Into<String>, start: u32, end: u32) -> Result<Self> {
        if start >= end {
            return Err(anyhow!("Invalid region: start {} >= end {}", start, end));
        }

        Ok(Self {
            contig: contig.into(),
            start,
            end,
        })
    }

    pub fn from_str_components(chrom: &str, start_str: &str, end_str: &str) -> Result<Self> {
        let start: u32 = start_str
            .parse()
            .map_err(|_| anyhow!("Invalid start position: {}", start_str))?;
        let end: u32 = end_str
            .parse()
            .map_err(|_| anyhow!("Invalid end position: {}", end_str))?;
        Self::new(chrom, start, end)
    }

    pub fn from_string(encoding: &str) -> Result<Self> {
        let elements: Vec<&str> = encoding.split(&[':', '-']).collect();

        if elements.len() != 3 {
            return Err(anyhow!("Invalid region encoding: {}", encoding));
        }

        let start: u32 = elements[1]
            .parse()
            .map_err(|_| anyhow!("Invalid region encoding: {}", encoding))?;
        let end: u32 = elements[2]
            .parse()
            .map_err(|_| anyhow!("Invalid region encoding: {}", encoding))?;

        Self::new(elements[0].to_string(), start, end)
    }

    pub fn intersect_position(&self, position: u32) -> bool {
        position >= self.start && position <= self.end
    }
}

#[cfg(test)]
mod tests {
    use super::GenomicRegion;
    #[test]
    fn init_region_from_valid_string_ok() {
        let region = GenomicRegion::from_string("chr1:100-200").unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_invalid_string_err() {
        assert_eq!(
            GenomicRegion::from_string("chr:1:100-200")
                .unwrap_err()
                .to_string(),
            "Invalid region encoding: chr:1:100-200".to_string()
        );
    }

    #[test]
    fn init_region_from_invalid_start_err() {
        assert_eq!(
            GenomicRegion::from_string("chr1:a-200")
                .unwrap_err()
                .to_string(),
            "Invalid region encoding: chr1:a-200".to_string()
        );
    }

    #[test]
    fn init_region_from_invalid_interval_err() {
        assert_eq!(
            GenomicRegion::from_string("chr1:200-100")
                .unwrap_err()
                .to_string(),
            "Invalid region: start 200 >= end 100".to_string()
        );
    }

    #[test]
    fn init_region_from_invalid_interval_new() {
        assert_eq!(
            GenomicRegion::new("chr1", 200, 100)
                .unwrap_err()
                .to_string(),
            "Invalid region: start 200 >= end 100".to_string()
        );
    }

    #[test]
    fn init_region_from_str_components_ok() {
        let region = GenomicRegion::from_str_components("chr1", "100", "200").unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_str_components_invalid_start_err() {
        assert_eq!(
            GenomicRegion::from_str_components("chr1", "a", "200")
                .unwrap_err()
                .to_string(),
            "Invalid start position: a".to_string()
        );
    }

    #[test]
    fn init_region_from_str_components_invalid_end_err() {
        assert_eq!(
            GenomicRegion::from_str_components("chr1", "100", "b")
                .unwrap_err()
                .to_string(),
            "Invalid end position: b".to_string()
        );
    }

    #[test]
    fn init_region_from_str_components_invalid_interval_err() {
        assert_eq!(
            GenomicRegion::from_str_components("chr1", "200", "100")
                .unwrap_err()
                .to_string(),
            "Invalid region: start 200 >= end 100".to_string()
        );
    }
}
