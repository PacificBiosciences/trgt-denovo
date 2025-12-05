//! Defines the `Locus` struct and associated functions for handling genomic loci.
//!
use crate::{
    readers::{open_catalog_reader, open_genome_reader, CatalogReader},
    region::GenomicRegion,
    util::Result,
};
use anyhow::{anyhow, Context};
use crossbeam_channel::Sender;
use rust_htslib::faidx;
use std::{collections::HashMap, io::BufRead, path::PathBuf};

/// Represents a genomic locus with its associated information.
///
/// A `Locus` contains the identifier, motifs, and flanking sequences
/// of a genomic region, as well as the region itself.
#[derive(Debug, PartialEq)]
pub struct Locus {
    /// A unique identifier for the locus.
    pub id: String,
    /// A list of motifs associated with the locus.
    pub motifs: Vec<String>,
    ///  The sequence of the left flanking region.
    pub left_flank: Vec<u8>,
    /// The sequence of the right flanking region.
    pub right_flank: Vec<u8>,
    /// The genomic region of the locus.
    pub region: GenomicRegion,
}

/// Decodes a single info field from a BED file into a name-value pair.
///
/// # Arguments
///
/// * `encoding` - A string slice representing the encoded info field.
///
/// # Returns
///
/// A result containing a tuple with the name and value of the info field if successful,
/// or an error if the field cannot be decoded.
fn decode_info_field(encoding: &str) -> Result<(&str, &str)> {
    let mut name_and_value = encoding.splitn(2, '=');
    let name = name_and_value
        .next()
        .ok_or_else(|| anyhow!("Invalid entry: {}", encoding))?;
    let value = name_and_value
        .next()
        .ok_or_else(|| anyhow!("Invalid entry: {}", encoding))?;
    Ok((name, value))
}

/// Decodes multiple info fields from a BED file into a hashmap.
///
/// # Arguments
///
/// * `info_fields` - A string slice representing the encoded info fields.
///
/// # Returns
///
/// A result containing a hashmap with the names and values of the info fields if successful,
/// or an error if any field cannot be decoded or if there are duplicate fields.
fn decode_fields(info_fields: &str) -> Result<HashMap<&str, String>> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding)?;
        if fields.insert(name, value.to_string()).is_some() {
            return Err(anyhow!("Duplicate field: {}", name));
        }
    }
    Ok(fields)
}

/// Represents a genomic locus with associated information.
impl Locus {
    /// Creates a new `Locus` instance from a BED record and genome data.
    ///
    /// # Arguments
    ///
    /// * `genome` - A mutable reference to an indexed FASTA reader.
    /// * `chrom_lookup` - A reference to a set containing chromosome names.
    /// * `line` - A reference to a BED record.
    /// * `flank_len` - The length of the flanking sequences to retrieve.
    ///
    /// # Returns
    ///
    /// A result containing the new `Locus` instance if successful, or an error if the
    /// chromosome is not found, required fields are missing, or flanks cannot be retrieved.
    pub fn new(
        genome: &faidx::Reader,
        chrom_lookup: &HashMap<String, u32>,
        line: &str,
        flank_len: usize,
    ) -> Result<Self> {
        const EXPECTED_FIELD_COUNT: usize = 4;
        let split_line: Vec<&str> = line.split_whitespace().collect();
        if split_line.len() != EXPECTED_FIELD_COUNT {
            return Err(anyhow!(
                "Expected {} fields in the format 'chrom start end info', found {}: {}",
                EXPECTED_FIELD_COUNT,
                split_line.len(),
                line
            ));
        }

        let (chrom, start, end, info_fields) = match &split_line[..] {
            [chrom, start, end, info_fields] => (*chrom, *start, *end, *info_fields),
            _ => unreachable!(),
        };

        if !chrom_lookup.contains_key(chrom) {
            return Err(anyhow!("Chromosome {} not found in reference", chrom));
        }

        let region = GenomicRegion::from_str_components(chrom, start, end)?;

        let fields = decode_fields(info_fields)?;

        let id = fields
            .get("ID")
            .ok_or_else(|| anyhow!("ID field missing"))?
            .to_string();
        let motifs = fields
            .get("MOTIFS")
            .ok_or_else(|| anyhow!("MOTIFS field missing"))?
            .split(',')
            .map(|s| s.to_string())
            .collect();

        let (left_flank, right_flank) = get_flanks(genome, &region, flank_len)?;

        Ok(Locus {
            id,
            motifs,
            left_flank,
            right_flank,
            region,
        })
    }
}

pub fn get_field(fields: &HashMap<&str, String>, key: &str) -> Result<String> {
    fields
        .get(key)
        .ok_or_else(|| anyhow!("{} field missing", key))
        .map(|s| s.to_string())
}

pub fn create_chrom_lookup(reader: &faidx::Reader) -> Result<HashMap<String, u32>> {
    let num_seqs = reader.n_seqs() as usize;
    let mut map = HashMap::with_capacity(num_seqs);
    for i in 0..num_seqs {
        let name = reader.seq_name(i as i32)?;
        let len = reader.fetch_seq_len(&name);
        let len_u32 = u32::try_from(len).map_err(|_| {
            anyhow!(
                "Sequence length for '{}' is negative and cannot be converted to u32",
                &name
            )
        })?;
        map.insert(name, len_u32);
    }
    Ok(map)
}

/// Retrieves a specific locus by its ID from a BED file using a genome reader.
///
/// # Arguments
///
/// * `genome_reader` - A mutable reference to an indexed FASTA reader.
/// * `catalog_reader` - A mutable reference to a BED reader.
/// * `tr_id` - The ID of the target repeat.
/// * `flank_len` - The length of the flanking sequences to retrieve.
///
/// # Returns
///
/// A result containing the `Locus` instance if found, or an error if the locus with the
/// specified ID cannot be found or processed.
pub fn get_locus(
    genome_reader: &faidx::Reader,
    catalog_reader: &mut CatalogReader,
    tr_id: &str,
    flank_len: usize,
) -> Result<Locus> {
    let chrom_lookup = create_chrom_lookup(genome_reader)?;
    let query = format!("ID={tr_id};");
    for (line_number, result_line) in catalog_reader.lines().enumerate() {
        let line =
            result_line.with_context(|| format!("Error reading BED line {}", line_number + 1))?;
        if line.contains(&query) {
            return Locus::new(genome_reader, &chrom_lookup, &line, flank_len)
                .with_context(|| format!("Error processing BED line {}", line_number + 1));
        }
    }
    Err(anyhow!("Unable to find locus {tr_id}"))
}

pub fn stream_loci_into_channel(
    repeats_path: PathBuf,
    genome_path: PathBuf,
    flank_len: usize,
    sender: Sender<Result<Locus>>,
) -> Result<()> {
    let catalog_reader = open_catalog_reader(&repeats_path)?;
    let genome_reader = open_genome_reader(&genome_path)?;
    let chrom_lookup = create_chrom_lookup(&genome_reader)?;

    for (line_number, result_line) in catalog_reader.lines().enumerate() {
        let line = match result_line {
            Ok(line) => line,
            Err(err) => {
                sender.send(Err(
                    anyhow!(err).context(format!("Error at BED line {}", line_number + 1))
                ))?;
                break;
            }
        };

        let locus_result = Locus::new(&genome_reader, &chrom_lookup, &line, flank_len);
        let locus = match locus_result {
            Ok(locus) => Ok(locus),
            Err(e) => {
                sender.send(Err(anyhow!("Error at BED line {}: {}", line_number + 1, e)))?;
                continue;
            }
        };

        sender
            .send(locus)
            .expect("Failed to send locus through channel");
    }
    Ok(())
}

fn get_flanks(
    genome: &faidx::Reader,
    region: &GenomicRegion,
    flank_len: usize,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let fetch_start = region.start as usize - flank_len;
    let fetch_end = region.end as usize + flank_len - 1;

    let full_seq = genome
        .fetch_seq(&region.contig, fetch_start, fetch_end)
        .map_err(|e| {
            anyhow!(
                "Error fetching sequence for region {}:{}-{}: {}",
                &region.contig,
                fetch_start,
                fetch_end,
                e
            )
        })?
        .to_ascii_uppercase();

    let left_flank = full_seq[..flank_len].to_vec();
    let right_flank = full_seq[full_seq.len() - flank_len..].to_vec();

    Ok((left_flank, right_flank))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn can_decode_info_field() {
        let decoded = decode_info_field("ID=AFF2").unwrap();
        assert_eq!("ID", decoded.0);
        assert_eq!("AFF2", decoded.1);
    }

    #[test]
    #[should_panic]
    fn panic_invalid_decode_info_field() {
        decode_info_field("ID:AFF2").unwrap();
    }
}
