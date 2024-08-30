//! Module for representing and building read information from alignment records.
//!
//! This module defines structures and functions for extracting and representing
//! various pieces of information from read alignment records, such as bases,
//! classification, and mismatch offsets

use noodles::sam::{
    alignment::Record,
    record::data::field::{tag, value::Array, Value},
};
use once_cell::sync::Lazy;

/// Represents a single read from an alignment record.
#[derive(Debug, PartialEq, Clone)]
pub struct ReadInfo {
    /// The sequence of bases for the read.
    pub bases: Box<[u8]>,
    /// The TRGT classification of the read.
    pub classification: Option<u8>,
    /// HP tag
    pub haplotype: Option<u8>,
    /// The start offset of the read relative to the reference and the locus start position.
    pub start_offset: Option<i32>,
    /// The end offset of the read relative to the reference and the locus end position.
    pub end_offset: Option<i32>,
    /// A list of offsets where SNPs occur between the original HiFi read and the reference relative to the locus start and end positions.
    pub mismatch_offsets: Option<Vec<i32>>,
    /// Read alignment start position
    pub pos: usize,
}

// Define custom tags for lookup in TRGT BAMlets
static AL_KEY: Lazy<tag::Tag> = Lazy::new(|| "AL".parse().unwrap());
static HP_KEY: Lazy<tag::Tag> = Lazy::new(|| "HP".parse().unwrap());
static MO_KEY: Lazy<tag::Tag> = Lazy::new(|| "MO".parse().unwrap());
static SO_KEY: Lazy<tag::Tag> = Lazy::new(|| "SO".parse().unwrap());
static EO_KEY: Lazy<tag::Tag> = Lazy::new(|| "EO".parse().unwrap());

/// Constructs a `ReadInfo` from a SAM record.
impl ReadInfo {
    /// Creates a new `ReadInfo` instance from a given SAM record.
    ///
    /// # Arguments
    ///
    /// * `record` - The SAM record to extract read information from.
    ///
    /// # Returns
    ///
    /// A `ReadInfo` instance containing the extracted information.
    pub fn new(record: Record) -> Self {
        let bases = record
            .sequence()
            .as_ref()
            .iter()
            .map(|base| u8::from(*base))
            .collect::<Vec<u8>>()
            .into_boxed_slice();

        let data = record.data();

        let classification = data.get(&*AL_KEY).and_then(|value| match value {
            Value::UInt8(v) => Some(*v),
            Value::Int32(v) => Some(*v as u8),
            _ => None,
        });

        let haplotype = data.get(&*HP_KEY).and_then(|value| match value {
            Value::UInt8(v) => Some(*v),
            _ => None,
        });

        let start_offset = data.get(&*SO_KEY).and_then(|value| match value {
            Value::Int32(v) => Some(*v),
            _ => None,
        });

        let end_offset = data.get(&*EO_KEY).and_then(|value| match value {
            Value::Int32(v) => Some(*v),
            _ => None,
        });

        let mismatch_offsets = data
            .get(&*MO_KEY)
            .and_then(|value| value.as_array())
            .and_then(|array| match array {
                Array::Int32(vec) => Some(vec.clone()),
                _ => None,
            });

        Self {
            bases,
            classification,
            haplotype,
            start_offset,
            end_offset,
            mismatch_offsets,
            pos: record.alignment_start().unwrap().get(),
        }
    }
}

/// Builder for constructing `ReadInfo` instances with optional fields.
#[derive(Debug, PartialEq, Clone)]
pub struct ReadInfoBuilder {
    /// The sequence of bases for the read.
    bases: Box<[u8]>,
    /// The TRGT classification of the read.
    classification: Option<u8>,
    /// HP tag
    haplotype: Option<u8>,
    /// The start offset of the read relative to the reference and the locus start position.
    start_offset: Option<i32>,
    /// The end offset of the read relative to the reference and the locus end position.
    end_offset: Option<i32>,
    /// A list of offsets where SNPs occur between the original HiFi read and the reference relative to the locus start and end positions.
    mismatch_offsets: Option<Vec<i32>>,
    /// Alignment start position
    pos: usize,
}

/// Provides an interface for building `ReadInfo` instances.
impl ReadInfoBuilder {
    /// Creates a new `ReadInfoBuilder` instance with default values.
    ///
    /// # Returns
    ///
    /// A new instance of `ReadInfoBuilder`.
    pub fn new() -> Self {
        Self {
            bases: Box::new([0u8; 10]),
            classification: Some(0),
            haplotype: None,
            start_offset: Some(0),
            end_offset: Some(0),
            mismatch_offsets: Some(vec![0, 0, 0, 0, 0]),
            pos: 0,
        }
    }

    /// Sets the bases for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `bases` - A sequence of bases represented as a byte array.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_bases<T: Into<Box<[u8]>>>(mut self, bases: T) -> Self {
        self.bases = bases.into();
        self
    }

    /// Sets the classification for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `classification` - An optional classification code as a byte.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_classification(mut self, classification: Option<u8>) -> Self {
        self.classification = classification;
        self
    }

    /// Sets the haplotype for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `haplotype` - An optional haplotype code as a byte.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_haplotype(mut self, haplotype: Option<u8>) -> Self {
        self.haplotype = haplotype;
        self
    }

    /// Sets the start offset for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `start_offset` - An optional start offset as an integer.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_start_offset(mut self, start_offset: Option<i32>) -> Self {
        self.start_offset = start_offset;
        self
    }

    /// Sets the end offset for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `end_offset` - An optional end offset as an integer.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_end_offset(mut self, end_offset: Option<i32>) -> Self {
        self.end_offset = end_offset;
        self
    }

    /// Sets the mismatch offsets for the `ReadInfo` being built.
    ///
    /// # Arguments
    ///
    /// * `mismatch_offsets` - An optional vector of mismatch offsets as integers.
    ///
    /// # Returns
    ///
    /// This method returns the builder itself.
    pub fn with_mismatch_offsets(mut self, mismatch_offsets: Option<Vec<i32>>) -> Self {
        self.mismatch_offsets = mismatch_offsets;
        self
    }

    pub fn with_pos(mut self, pos: usize) -> Self {
        self.pos = pos;
        self
    }

    /// Builds a `ReadInfo` instance from the builder.
    ///
    /// # Returns
    ///
    /// A result containing the `ReadInfo` instance if successful, or an error message if not.
    pub fn build(self) -> Result<ReadInfo, &'static str> {
        Ok(ReadInfo {
            bases: self.bases,
            classification: self.classification,
            haplotype: self.haplotype,
            start_offset: self.start_offset,
            end_offset: self.end_offset,
            mismatch_offsets: self.mismatch_offsets,
            pos: self.pos,
        })
    }
}

/// Provides default values for building a `ReadInfo` instance.
///
/// This implementation sets all optional fields to `None` and initializes `bases`
/// with a default value of 10 zeroed bytes.
impl Default for ReadInfoBuilder {
    /// Returns a `ReadInfoBuilder` instance with default values.
    fn default() -> Self {
        Self {
            bases: Box::new([0u8; 10]),
            classification: None,
            haplotype: None,
            start_offset: None,
            end_offset: None,
            mismatch_offsets: None,
            pos: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_info_builder() {
        let dummy_read_info = ReadInfoBuilder::default()
            .with_bases(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            .with_classification(Some(1))
            .with_start_offset(Some(100))
            .with_end_offset(Some(200))
            .with_mismatch_offsets(Some(vec![1, 2, 3, 4, 5]))
            .build()
            .unwrap();

        assert_eq!(
            dummy_read_info.bases,
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10].into()
        );
        assert_eq!(dummy_read_info.classification, Some(1));
        assert_eq!(dummy_read_info.start_offset, Some(100));
        assert_eq!(dummy_read_info.end_offset, Some(200));
        assert_eq!(dummy_read_info.mismatch_offsets, Some(vec![1, 2, 3, 4, 5]));
    }

    #[test]
    fn test_read_info_builder_defaults() {
        let default_read_info = ReadInfoBuilder::default().build().unwrap();

        assert_eq!(
            default_read_info.bases,
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0].into()
        );
        assert_eq!(default_read_info.classification, None);
        assert_eq!(default_read_info.start_offset, None);
        assert_eq!(default_read_info.end_offset, None);
        assert_eq!(default_read_info.mismatch_offsets, None);
    }
}
