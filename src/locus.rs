//! Defines the `Locus` struct and associated functions for handling genomic loci.
//!
use crate::util::Result;
use anyhow::{anyhow, Context};
use crossbeam_channel::Sender;
use noodles::{
    bed,
    core::{Position, Region},
    fasta::{self, io::BufReadSeek, record::Sequence},
};
use std::{
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::BufReader,
    path::PathBuf,
};

/// Represents a genomic locus with its associated information.
///
/// A `Locus` contains the identifier, structure, motifs, and flanking sequences
/// of a genomic region, as well as the region itself.
#[derive(Debug, PartialEq)]
pub struct Locus {
    /// A unique identifier for the locus.
    pub id: String,
    ///  The structure of the locus.
    pub struc: String,
    /// A list of motifs associated with the locus.
    pub motifs: Vec<String>,
    ///  The sequence of the left flanking region.
    pub left_flank: Vec<u8>,
    /// The sequence of the right flanking region.
    pub right_flank: Vec<u8>,
    /// The genomic region of the locus.
    pub region: Region,
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
        genome: &mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
        chrom_lookup: &HashSet<String>,
        line: &bed::Record<3>,
        flank_len: usize,
    ) -> Result<Self> {
        if !chrom_lookup.contains(line.reference_sequence_name()) {
            return Err(anyhow!(
                "Chromosome {} not found in reference",
                line.reference_sequence_name()
            ));
        }

        // -1 because bed is 0-based
        let region = Region::new(
            line.reference_sequence_name(),
            Position::new(usize::from(line.start_position()) - 1).unwrap()..=line.end_position(),
        );

        let fields = decode_fields(&line.optional_fields()[0])?;

        let id = fields
            .get("ID")
            .ok_or_else(|| anyhow!("ID field missing"))?
            .to_string();
        let struc = fields
            .get("STRUC")
            .ok_or_else(|| anyhow!("STRUC field missing"))?
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
            struc,
            motifs,
            left_flank,
            right_flank,
            region,
        })
    }
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
    genome_reader: &mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
    catalog_reader: &mut bed::Reader<BufReader<File>>,
    tr_id: &str,
    flank_len: usize,
) -> Result<Locus> {
    let chrom_lookup = HashSet::from_iter(
        genome_reader
            .index()
            .iter()
            .map(|entry| entry.name().to_string()),
    );

    let query = format!("ID={tr_id};");
    for (line_number, result) in catalog_reader.records::<3>().enumerate() {
        let line = result.with_context(|| format!("Error reading BED line {}", line_number + 1))?;
        if line.optional_fields().first().unwrap().contains(&query) {
            return Locus::new(genome_reader, &chrom_lookup, &line, flank_len)
                .with_context(|| format!("Error processing BED line {}", line_number + 1));
        }
    }
    Err(anyhow!("Unable to find locus {tr_id}"))
}

/// Returns an iterator over all loci in a BED file using a genome reader.
///
/// # Arguments
///
/// * `genome_reader` - A mutable reference to an indexed FASTA reader.
/// * `catalog_reader` - A mutable reference to a BED reader.
/// * `flank_len` - The length of the flanking sequences to retrieve for each locus.
///
/// # Returns
///
/// An iterator that yields results containing `Locus` instances or errors encountered
/// during processing.
pub fn get_loci<'a>(
    genome_reader: &'a mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
    catalog_reader: &'a mut bed::Reader<BufReader<File>>,
    flank_len: usize,
) -> impl Iterator<Item = Result<Locus>> + 'a {
    let chrom_lookup = HashSet::from_iter(
        genome_reader
            .index()
            .iter()
            .map(|entry| entry.name().to_string()),
    );

    catalog_reader
        .records::<3>()
        .enumerate()
        .map(move |(line_number, result_line)| {
            result_line
                .map_err(|err| {
                    anyhow!(err).context(format!("Error at BED line {}", line_number + 1))
                })
                .and_then(|line| {
                    Locus::new(genome_reader, &chrom_lookup, &line, flank_len)
                        .with_context(|| format!("Error processing BED line {}", line_number + 1))
                })
        })
}

pub fn stream_loci_into_channel(
    bed_filename: PathBuf,
    reference_filename: PathBuf,
    flank_len: usize,
    sender: Sender<Result<Locus>>,
) {
    let mut catalog_reader: bed::Reader<BufReader<File>> = fs::File::open(bed_filename.clone())
        .map(BufReader::new)
        .map(bed::Reader::new)
        .unwrap();

    let mut genome_reader = fasta::indexed_reader::Builder::default()
        .build_from_path(&reference_filename)
        .unwrap();

    let chrom_lookup = HashSet::from_iter(
        genome_reader
            .index()
            .iter()
            .map(|entry| entry.name().to_string()),
    );

    for (line_number, result_line) in catalog_reader.records::<3>().enumerate() {
        let line = match result_line {
            Ok(line) => line,
            Err(err) => {
                let error = anyhow!(err).context(format!("Error at BED line {}", line_number + 1));
                sender
                    .send(Err(error))
                    .expect("Failed to send error through channel");
                continue;
            }
        };

        let locus = Locus::new(&mut genome_reader, &chrom_lookup, &line, flank_len)
            .with_context(|| format!("Error processing BED line {}", line_number + 1));

        sender
            .send(locus)
            .expect("Failed to send locus through channel");
    }
}

/// Retrieves a flank sequence from the genome given a region and start/end positions.
///
/// # Arguments
///
/// * `genome` - A mutable reference to an indexed FASTA reader.
/// * `region` - A reference to the region of interest.
/// * `start` - The start position of the flank.
/// * `end` - The end position of the flank.
///
/// # Returns
///
/// A result containing the flank sequence if successful, or an error if the sequence
/// cannot be extracted.
fn get_flank(
    genome: &mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
    region: &Region,
    start: usize,
    end: usize,
) -> Result<Sequence> {
    let start_pos = Position::try_from(start + 1)?;
    let end_pos = Position::try_from(end)?;
    let query_region = Region::new(region.name(), start_pos..=end_pos);
    match genome.query(&query_region) {
        Ok(seq) => Ok(seq.sequence().to_owned()),
        Err(_) => Err(anyhow!("Unable to extract: {:?}", region)),
    }
}

/// Retrieves both left and right flank sequences for a given region from the genome.
///
/// # Arguments
///
/// * `genome` - A mutable reference to an indexed FASTA reader.
/// * `region` - A reference to the region of interest.
/// * `flank_len` - The length of the flanking sequences to retrieve.
///
/// # Returns
///
/// A result containing a tuple with the left and right flank sequences if successful,
/// or an error if the sequences cannot be extracted.
fn get_flanks(
    genome: &mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
    region: &Region,
    flank_len: usize,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let (lf_start, lf_end) = (
        region.interval().start().unwrap().get() - flank_len,
        region.interval().start().unwrap().get(),
    );
    let (rf_start, rf_end) = (
        region.interval().end().unwrap().get(),
        region.interval().end().unwrap().get() + flank_len,
    );

    let left_flank = get_flank(genome, region, lf_start, lf_end)?;
    let right_flank = get_flank(genome, region, rf_start, rf_end)?;

    Ok((
        left_flank.as_ref().to_ascii_uppercase(),
        right_flank.as_ref().to_ascii_uppercase(),
    ))
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
