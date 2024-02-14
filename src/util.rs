//! Utility functions and types for error handling and file path validation.
//!
//! This module provides common utility functions used throughout the program,
//! including custom result types and error handling mechanisms.

use anyhow::anyhow;
use log;
use std::path::Path;

/// Custom result type for error handling throughout the program.
pub type Result<T> = anyhow::Result<T>;

/// Logs the provided error and exits the program.
///
/// # Arguments
///
/// * `err` - The error to log before exiting.
pub fn handle_error_and_exit(err: anyhow::Error) -> ! {
    log::error!("{:#}", err);
    std::process::exit(1);
}

/// Checks if the provided file path exists.
///
/// # Arguments
///
/// * `path` - The file path to check.
///
/// # Returns
///
/// A result indicating success if the path exists, or an error if it does not.
pub fn try_exists(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(anyhow!("Path does not exist: {}", path.display()));
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub struct AlnScoring {
    pub mismatch: i32,
    pub gap_opening1: i32,
    pub gap_extension1: i32,
    pub gap_opening2: i32,
    pub gap_extension2: i32,
}

#[derive(Debug)]
pub struct Params {
    pub clip_len: usize,
    pub parent_quantile: f64,
    pub partition_by_alignment: bool,
}
