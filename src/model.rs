use serde::Serialize;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct AlnScoring {
    pub mismatch: i32,
    pub gap_opening1: i32,
    pub gap_extension1: i32,
    pub gap_opening2: i32,
    pub gap_extension2: i32,
}

impl Default for AlnScoring {
    fn default() -> Self {
        Self {
            mismatch: 8,
            gap_opening1: 4,
            gap_extension1: 2,
            gap_opening2: 24,
            gap_extension2: 1,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum QuickMode {
    AL(Option<f64>),
    MC(Option<f64>),
}

impl QuickMode {
    pub fn is_zero(&self) -> bool {
        matches!(self, QuickMode::AL(None) | QuickMode::MC(None))
    }
}

pub struct Params {
    pub clip_len: usize,
    pub parent_quantile: f64,
    pub partition_by_alignment: bool,
    pub quick_mode: Option<QuickMode>,
}

/// Enumerates the types of de novo events that can occur.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoType {
    Expansion,
    Contraction,
    Substitution,
    Unclear,
}

/// Represents the status of a de novo event, indicating whether it is de novo and its type.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum DenovoStatus {
    Denovo(DenovoType),
    NotDenovo,
    Unknown,
}

impl std::fmt::Display for DenovoType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DenovoType::Expansion => write!(f, "+"),
            DenovoType::Contraction => write!(f, "-"),
            DenovoType::Substitution => write!(f, "="),
            DenovoType::Unclear => write!(f, "?"),
        }
    }
}

impl std::fmt::Display for DenovoStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DenovoStatus::Denovo(denovo_type) => write!(f, "Y:{}", denovo_type),
            DenovoStatus::NotDenovo => write!(f, "X"),
            DenovoStatus::Unknown => write!(f, "."),
        }
    }
}

/// Represents the number of alleles inherited from a parent.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum AlleleNum {
    One,
    Two,
    Unclear,
}

/// Represents the origin of an allele, specifying which parent it was inherited from.
#[derive(Debug, PartialEq, Clone, Serialize)]
pub enum AlleleOrigin {
    Father { allele: AlleleNum },
    Mother { allele: AlleleNum },
    Unclear,
    Unknown,
}

impl AlleleOrigin {
    pub fn new(value: usize) -> Option<Self> {
        match value {
            0 => Some(AlleleOrigin::Mother {
                allele: AlleleNum::One,
            }),
            1 => Some(AlleleOrigin::Mother {
                allele: AlleleNum::Two,
            }),
            2 => Some(AlleleOrigin::Father {
                allele: AlleleNum::One,
            }),
            3 => Some(AlleleOrigin::Father {
                allele: AlleleNum::Two,
            }),
            _ => None,
        }
    }
}

impl fmt::Display for AlleleNum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlleleNum::One => write!(f, "1"),
            AlleleNum::Two => write!(f, "2"),
            AlleleNum::Unclear => write!(f, "?"),
        }
    }
}

impl fmt::Display for AlleleOrigin {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlleleOrigin::Father { allele } => write!(f, "F:{}", allele),
            AlleleOrigin::Mother { allele } => write!(f, "M:{}", allele),
            AlleleOrigin::Unclear => write!(f, "?"),
            AlleleOrigin::Unknown => write!(f, "."),
        }
    }
}
