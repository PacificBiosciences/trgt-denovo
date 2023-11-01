use noodles::sam::alignment::Record;
use noodles::sam::record::data::field::value::Array;
use noodles::sam::record::data::field::{tag, Value};
use once_cell::sync::Lazy;

#[derive(Debug, PartialEq, Clone)]
pub struct FlankingReadInfo {
    is_left_flank: bool,
}

#[derive(Debug, PartialEq, Clone)]
pub struct ReadInfo {
    pub bases: Box<[u8]>,
    pub classification: Option<u8>,
    pub start_offset: Option<i32>,
    pub end_offset: Option<i32>,
    pub mismatch_offsets: Option<Vec<i32>>,
    pub flank_info: Option<FlankingReadInfo>,
}

static AL_KEY: Lazy<tag::Tag> = Lazy::new(|| "AL".parse().unwrap());
static MO_KEY: Lazy<tag::Tag> = Lazy::new(|| "MO".parse().unwrap());
static SO_KEY: Lazy<tag::Tag> = Lazy::new(|| "SO".parse().unwrap());
static EO_KEY: Lazy<tag::Tag> = Lazy::new(|| "EO".parse().unwrap());

impl ReadInfo {
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
            start_offset,
            end_offset,
            mismatch_offsets,
            flank_info: None,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]

pub struct ReadInfoBuilder {
    bases: Box<[u8]>,
    classification: Option<u8>,
    start_offset: Option<i32>,
    end_offset: Option<i32>,
    mismatch_offsets: Option<Vec<i32>>,
    flank_info: Option<FlankingReadInfo>,
}

impl ReadInfoBuilder {
    pub fn new() -> Self {
        Self {
            bases: Box::new([0u8; 10]),
            classification: Some(0),
            start_offset: Some(0),
            end_offset: Some(0),
            mismatch_offsets: Some(vec![0, 0, 0, 0, 0]),
            flank_info: None,
        }
    }

    pub fn with_bases<T: Into<Box<[u8]>>>(mut self, bases: T) -> Self {
        self.bases = bases.into();
        self
    }

    pub fn with_classification(mut self, classification: Option<u8>) -> Self {
        self.classification = classification;
        self
    }

    pub fn with_start_offset(mut self, start_offset: Option<i32>) -> Self {
        self.start_offset = start_offset;
        self
    }

    pub fn with_end_offset(mut self, end_offset: Option<i32>) -> Self {
        self.end_offset = end_offset;
        self
    }

    pub fn with_mismatch_offsets(mut self, mismatch_offsets: Option<Vec<i32>>) -> Self {
        self.mismatch_offsets = mismatch_offsets;
        self
    }

    pub fn with_flank_info(mut self, flank_info: Option<FlankingReadInfo>) -> Self {
        self.flank_info = flank_info;
        self
    }

    pub fn build(self) -> Result<ReadInfo, &'static str> {
        Ok(ReadInfo {
            bases: self.bases,
            classification: self.classification,
            start_offset: self.start_offset,
            end_offset: self.end_offset,
            mismatch_offsets: self.mismatch_offsets,
            flank_info: self.flank_info,
        })
    }
}

impl Default for ReadInfoBuilder {
    fn default() -> Self {
        Self {
            bases: Box::new([0u8; 10]),
            classification: None,
            start_offset: None,
            end_offset: None,
            mismatch_offsets: None,
            flank_info: None,
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
            .with_flank_info(Some(FlankingReadInfo {
                is_left_flank: true,
            }))
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
        assert_eq!(
            dummy_read_info.flank_info,
            Some(FlankingReadInfo {
                is_left_flank: true
            })
        );
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
        assert_eq!(default_read_info.flank_info, None);
    }
}
