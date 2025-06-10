use crate::util::{AlnScoring, QuickMode, Result};
use anyhow::anyhow;
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use log::{Level, LevelFilter, Record};
use once_cell::sync::Lazy;
use owo_colors::{OwoColorize, Style};
use std::{
    collections::HashMap,
    io::Write,
    ops::Deref,
    path::{Path, PathBuf},
};

/// Full version string including the crate version and git description.
///
/// This version string is used in the command-line interface to provide detailed version information.
/// It includes the crate version from Cargo.toml and additional build information such as the git commit hash.
/// # Examples
/// * `0.1.0-1ba958a-dirty` - while on a dirty branch
/// * `0.1.0-1ba958a` - with a fresh commit
pub static FULL_VERSION: Lazy<String> = Lazy::new(|| {
    let git_describe = env!("VERGEN_GIT_DESCRIBE");
    if git_describe.is_empty() {
        env!("CARGO_PKG_VERSION").to_string()
    } else {
        format!("{}-{}", env!("CARGO_PKG_VERSION"), git_describe)
    }
});

#[derive(Parser, Debug)]
#[command(name="trgt-denovo",
          author="Tom Mokveld <tmokveld@pacificbiosciences.com>\nEgor Dolzhenko <edolzhenko@pacificbiosciences.com>", 
          version=&**FULL_VERSION,
          about="Tandem repeat de novo caller",
          long_about = None,
          after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
          This program comes with ABSOLUTELY NO WARRANTY; it is intended for
          Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
    #[arg(
        short = 'v',
        long = "verbose",
        action = ArgAction::Count,
        global = true
    )]
    pub verbosity: u8,

    /// Silence all output
    #[arg(
        long = "quiet",
        action = ArgAction::SetTrue,
        global = true,
        conflicts_with = "verbosity",
    )]
    pub quiet: bool,
}

impl Command {
    pub fn name(&self) -> &'static str {
        match self {
            Command::Trio(_) => "trio",
            Command::Duo(_) => "duo",
        }
    }
}

#[derive(Parser, Debug, Clone)]
pub struct SharedArgs {
    /// Path to reference genome FASTA
    #[arg(
        required = true,
        short = 'r',
        long = "reference",
        value_name = "FASTA",
        value_parser = check_file_exists
    )]
    pub reference_filename: PathBuf,

    /// BED file with repeat coordinates
    #[arg(
        required = true,
        short = 'b',
        long = "bed",
        value_name = "BED",
        value_parser = check_file_exists
    )]
    pub bed_filename: PathBuf,

    /// Output tsv path
    #[arg(
        required = true,
        short = 'o',
        long = "out",
        value_name = "TSV",
        value_parser = check_prefix_path
    )]
    pub output_path: String,

    /// TRID of a specific repeat to analyze, should be in the BED file (note: this is assumed to be unique where the first match will be analyzed
    #[arg(long = "trid", value_name = "TRID")]
    pub trid: Option<String>,

    #[arg(
        short = '@',
        value_name = "THREADS",
        default_value = "1",
        value_parser = threads_in_range
    )]
    pub num_threads: usize,

    /// Amount of additional flanking sequence that should be used during alignment
    #[arg(
        help_heading("Advanced"),
        long = "flank-len",
        value_name = "FLANK_LEN",
        default_value = "50"
    )]
    pub flank_len: usize,

    /// Score alignments without stripping the flanks
    #[arg(help_heading("Advanced"), long = "no-clip-aln", value_name = "CLIP")]
    pub no_clip_aln: bool,

    /// Quantile of alignments scores to determine the threshold
    #[arg(
        help_heading("Advanced"),
        long = "p-quantile",
        value_name = "QUANTILE",
        default_value = "1.0",
        value_parser = parse_quantile
    )]
    pub p_quantile: f64,

    /// Scoring function for 2-piece gap affine alignment (non-negative values): mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2
    #[arg(
        help_heading("Advanced"),
        long = "aln-scoring",
        value_name = "SCORING",
        default_value = "8,4,2,24,1",
        value_parser = scoring_from_string
    )]
    pub aln_scoring: AlnScoring,

    /// Within-sample partitioning using alignment rather than the TRGT BAMlet AL field
    #[arg(
        help_heading("Advanced"),
        long = "partition-by-aln",
        value_name = "ALN"
    )]
    pub partition_by_alignment: bool,

    /// Only test loci that differ by a certain fraction in allele length. Format: <field>,<fraction> (e.g. AL,0.1 or AL)
    #[arg(
        help_heading("Advanced"),
        long = "quick",
        value_parser = parse_quick_option
    )]
    pub quick: Option<QuickMode>,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    Trio(TrioArgs),
    Duo(DuoArgs),
}

impl Deref for TrioArgs {
    type Target = SharedArgs;

    fn deref(&self) -> &Self::Target {
        &self.shared
    }
}

impl Deref for DuoArgs {
    type Target = SharedArgs;

    fn deref(&self) -> &Self::Target {
        &self.shared
    }
}

#[derive(Parser, Debug, Clone)]
#[command(group(ArgGroup::new("trio")))]
#[command(arg_required_else_help(true))]
pub struct TrioArgs {
    #[command(flatten)]
    pub shared: SharedArgs,

    /// Prefix of mother VCF and spanning reads BAM files
    #[arg(
        required = true,
        short = 'm',
        long = "mother",
        value_name = "PREFIX",
        value_parser = check_prefix_path
    )]
    pub mother_prefix: String,

    /// Prefix of father VCF and spanning reads BAM files
    #[arg(
        required = true,
        short = 'f',
        long = "father",
        value_name = "PREFIX",
        value_parser = check_prefix_path
    )]
    pub father_prefix: String,

    /// Prefix of child VCF and spanning reads BAM files
    #[arg(
        required = true,
        short = 'c',
        long = "child",
        value_name = "PREFIX",
        value_parser = check_prefix_path
    )]
    pub child_prefix: String,
}

#[derive(Parser, Debug, Clone)]
#[command(group(ArgGroup::new("duo")))]
#[command(arg_required_else_help(true))]
pub struct DuoArgs {
    #[command(flatten)]
    pub shared: SharedArgs,

    /// Prefix of first sample VCF and spanning reads BAM files
    #[arg(
        required = true,
        short = '1',
        long = "sample-a",
        value_name = "PREFIX",
        value_parser = check_prefix_path
    )]
    pub a_prefix: String,

    /// Prefix of second sample VCF and spanning reads BAM files
    #[arg(
        required = true,
        short = '2',
        long = "sample-b",
        value_name = "PREFIX",
        value_parser = check_prefix_path
    )]
    pub b_prefix: String,
}

/// Initializes the verbosity level for logging based on the command-line arguments.
///
/// Sets up the logger with a specific verbosity level that is determined
/// by the number of occurrences of the `-v` or `--verbose` flag in the command-line arguments.
///
/// # Arguments
///
/// * `args` - A reference to the parsed command-line arguments.
pub fn init_verbose(args: &Cli) {
    let filter_level: LevelFilter = if args.quiet {
        LevelFilter::Off
    } else {
        match args.verbosity {
            0 => LevelFilter::Info,  // -v
            1 => LevelFilter::Debug, // -vv
            _ => LevelFilter::Trace, // -vvv or more
        }
    };

    env_logger::Builder::from_default_env()
        .format(format_log)
        .filter_level(filter_level)
        .init();
}

static LEVEL_STYLES: Lazy<HashMap<Level, (&'static str, Style)>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert(Level::Error, ("ERROR", Style::new().red()));
    m.insert(Level::Warn, ("WARN", Style::new().yellow()));
    m.insert(Level::Info, ("INFO", Style::new().green()));
    m.insert(Level::Debug, ("DEBUG", Style::new().blue()));
    m.insert(Level::Trace, ("TRACE", Style::new().magenta()));
    m
});

fn format_log(buf: &mut env_logger::fmt::Formatter, record: &Record) -> std::io::Result<()> {
    let (level_text, style) = LEVEL_STYLES.get(&record.level()).unwrap();
    let level_str = level_text.style(*style).to_string();
    writeln!(
        buf,
        "{} [{}] - {}",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
        level_str,
        record.args()
    )
}

/// Checks if the provided path prefix exists.
///
/// Validates that the path prefix provided as an argument exists in the file system.
/// It is used to ensure that the file paths constructed using this prefix will be valid.
///
/// # Arguments
///
/// * `s` - A string slice representing the path prefix to check.
///
/// # Returns
///
/// Returns a `Result<String>` which is Ok if the path prefix exists, or an Err with a descriptive message if not.
fn check_prefix_path(s: &str) -> Result<String> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(anyhow!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(s.to_string())
}

/// Validates that the provided string represents a valid number of threads.
///
/// Checks if the string argument can be parsed into a non-zero positive integer
/// that represents the number of threads to use. It ensures that the number of threads is within
/// a valid range.
///
/// # Arguments
///
/// * `s` - A string slice representing the number of threads.
///
/// # Returns
///
/// Returns a `Result<usize>` which is Ok if the number is valid, or an Err with a descriptive message if not.
fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{}` is not a valid thread number", s))?;
    if thread == 0 {
        return Err(anyhow!("Number of threads must be >= 1"));
    }
    Ok(thread)
}

/// Parses a string into a floating-point number representing a quantile.
///
/// Attempts to parse a string into a `f64` that represents a quantile value.
/// It validates that the value is within the range [0.0, 1.0].
///
/// # Arguments
///
/// * `s` - A string slice representing the quantile to parse.
///
/// # Returns
///
/// Returns a `Result<f64>` which is Ok if the value is within the valid range, or an Err with a descriptive message if not.
fn parse_quantile(s: &str) -> Result<f64> {
    let value = s
        .parse::<f64>()
        .map_err(|e| anyhow!("Could not parse float: {}", e))?;
    if !(0.0..=1.0).contains(&value) {
        Err(anyhow!(
            "The value must be between 0.0 and 1.0, got: {}",
            value
        ))
    } else {
        Ok(value)
    }
}

/// Checks if the provided file path exists.
///
/// Validates that the file path provided as an argument exists in the file system.
/// It is used to ensure that the file paths provided for input files are valid before attempting to process them.
///
/// # Arguments
///
/// * `s` - A string slice representing the file path to check.
///
/// # Returns
///
/// Returns a `Result<PathBuf>` which is Ok if the file exists, or an Err with a descriptive message if not.
fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        return Err(anyhow!("File does not exist: {}", path.display()));
    }
    Ok(path.to_path_buf())
}

/// Parses a string of comma-separated alignment scoring parameters into an `AlnScoring` struct.
///
/// # Arguments
///
/// * `s` - A string slice containing the scoring parameters in the format "mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2".
///
/// # Returns
///
/// Returns a `Result<AlnScoring>` which is Ok if the string is correctly formatted and contains valid values,
/// or an Err with a descriptive message if the input is invalid.
/// ```
fn scoring_from_string(s: &str) -> Result<AlnScoring> {
    const NUM_EXPECTED_VALUES: usize = 5;
    let values: Vec<i32> = s.split(',').filter_map(|x| x.parse().ok()).collect();
    if values.len() != NUM_EXPECTED_VALUES {
        return Err(anyhow!(
            "Expected {} comma-separated integers values in scoring. Got {} -> {}",
            NUM_EXPECTED_VALUES,
            values.len(),
            s
        ));
    }

    let (mismatch, gap_opening1, gap_extension1, gap_opening2, gap_extension2) =
        (values[0], values[1], values[2], values[3], values[4]);

    if mismatch <= 0
        || gap_opening1 < 0
        || gap_extension1 <= 0
        || gap_opening2 < 0
        || gap_extension2 <= 0
    {
        return Err(anyhow!(
            "Invalid penalties. Got (mismatch={}, gap_opening1={}, gap_extension1={}, gap_opening2={},
gap_extension2={}) -> (mismatch>0, gap_opening1>=0, gap_extension1>0, gap_opening2>=0,
gap_extension2>0)",
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2
        ));
    }

    Ok(AlnScoring {
        mismatch,
        gap_opening1,
        gap_extension1,
        gap_opening2,
        gap_extension2,
    })
}

fn parse_quick_option(s: &str) -> Result<QuickMode> {
    let parts: Vec<&str> = s.split(',').collect();

    let (field, fraction) = match parts.len() {
        1 => (parts[0].to_uppercase(), None),
        2 => {
            let frac = parts[1]
                .parse()
                .map_err(|_| anyhow!("Invalid fraction value"))?;
            if !(0.0..=1.0).contains(&frac) {
                return Err(anyhow!("Fraction must be between 0.0 and 1.0"));
            }
            (
                parts[0].to_uppercase(),
                if frac == 0.0 { None } else { Some(frac) },
            )
        }
        _ => {
            return Err(anyhow!(
                "Invalid quick option format. Expected <field> or <field>,<fraction>"
            ))
        }
    };

    match field.as_str() {
        "AL" => Ok(QuickMode::AL(fraction)),
        // "MC" => Ok(QuickMode::MC(fraction)),
        _ => Err(anyhow!("Invalid field for quick option. Must be AL")),
    }
}
