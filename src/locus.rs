use crate::util::Result;
use anyhow::{anyhow, Context};
use noodles::{
    bed,
    core::{Position, Region},
    fasta::{self, io::BufReadSeek, record::Sequence},
};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufReader,
};

#[derive(Debug, PartialEq)]
pub struct Locus {
    pub id: String,
    pub struc: String,
    pub motifs: Vec<String>,
    pub left_flank: Vec<u8>,
    pub right_flank: Vec<u8>,
    pub region: Region,
}

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

impl Locus {
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
        if line.optional_fields().get(0).unwrap().contains(&query) {
            return Locus::new(genome_reader, &chrom_lookup, &line, flank_len)
                .with_context(|| format!("Error processing BED line {}", line_number + 1));
        }
    }
    Err(anyhow!("Unable to find locus {tr_id}"))
}

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
        .filter_map(move |(line_number, result_line)| match result_line {
            Ok(line) => Some(Ok((line, line_number))),
            Err(err) => Some(Err(
                anyhow!(err).context(format!("Error at BED line {}", line_number + 1))
            )),
        })
        .map(move |result| {
            result.and_then(|(line, line_number)| {
                Locus::new(genome_reader, &chrom_lookup, &line, flank_len)
                    .with_context(|| format!("Error processing BED line {}", line_number + 1))
            })
        })
}

fn get_flank(
    genome: &mut fasta::IndexedReader<Box<dyn BufReadSeek>>,
    region: &Region,
    start: usize,
    end: usize,
) -> Result<Sequence> {
    let start_pos = Position::try_from(start + 1).unwrap();
    let end_pos = Position::try_from(end).unwrap();
    let query_region = Region::new(region.name(), start_pos..=end_pos);
    match genome.query(&query_region) {
        Ok(seq) => Ok(seq.sequence().to_owned()),
        Err(_) => Err(anyhow!("Unable to extract: {:?}", region)),
    }
}

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
