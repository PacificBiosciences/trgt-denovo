use crate::util::QuickMode;
use crate::util::{AlnScoring, Result};
use anyhow::anyhow;
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use env_logger::fmt::Color;
use log::{Level, LevelFilter};
use once_cell::sync::Lazy;
use std::ops::Deref;
use std::{
    io::Write,
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

    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = ArgAction::Count)]
    pub verbosity: u8,
}

#[derive(Parser, Debug, Clone)]
pub struct SharedArgs {
    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub reference_filename: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "bed")]
    #[clap(help = "BED file with repeat coordinates")]
    #[clap(value_name = "BED")]
    #[arg(value_parser = check_file_exists)]
    pub bed_filename: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "out")]
    #[clap(help = "Output tsv path")]
    #[clap(value_name = "TSV")]
    #[arg(value_parser = check_prefix_path)]
    pub output_path: String,

    #[clap(long = "trid")]
    #[clap(
        help = "TRID of a specific repeat to analyze, should be in the BED file (note: this is assumed to be unique where the first match will be analyzed"
    )]
    #[clap(value_name = "TRID")]
    pub trid: Option<String>,

    #[clap(short = '@')]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    #[arg(value_parser = threads_in_range)]
    pub num_threads: usize,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(help = "Amount of additional flanking sequence that should be used during alignment")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(default_value = "50")]
    pub flank_len: usize,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "no-clip-aln")]
    #[clap(help = "Score alignments without stripping the flanks")]
    #[clap(value_name = "CLIP")]
    pub no_clip_aln: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "p-quantile")]
    #[clap(help = "Quantile of alignments scores to determine the threshold")]
    #[clap(value_name = "QUANTILE")]
    #[clap(default_value = "1.0")]
    #[arg(value_parser = parse_quantile)]
    pub p_quantile: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "aln-scoring")]
    #[clap(value_name = "SCORING")]
    #[clap(
        help = "Scoring function for 2-piece gap affine alignment (non-negative values): mismatch,gap_opening1,gap_extension1,gap_opening2,gap_extension2"
    )]
    #[clap(default_value = "8,4,2,24,1")]
    #[arg(value_parser = scoring_from_string)]
    pub aln_scoring: AlnScoring,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "partition-by-aln")]
    #[clap(
        help = "Within-sample partitioning using alignment rather than the TRGT BAMlet AL field"
    )]
    #[clap(value_name = "ALN")]
    pub partition_by_alignment: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "quick")]
    #[clap(
        help = "Only test loci that differ by a certain fraction in allele length. Format: <field>,<fraction> (e.g. AL,0.1 or AL)"
    )]
    #[clap(value_parser = parse_quick_option)]
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

    #[clap(required = true)]
    #[clap(short = 'm')]
    #[clap(long = "mother")]
    #[clap(help = "Prefix of mother VCF and spanning reads BAM files")]
    #[clap(value_name = "PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub mother_prefix: String,

    #[clap(required = true)]
    #[clap(short = 'f')]
    #[clap(long = "father")]
    #[clap(help = "Prefix of father VCF and spanning reads BAM files")]
    #[arg(value_parser = check_prefix_path)]
    #[clap(value_name = "PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub father_prefix: String,

    #[clap(required = true)]
    #[clap(short = 'c')]
    #[clap(long = "child")]
    #[clap(help = "Prefix of child VCF and spanning reads BAM files")]
    #[clap(value_name = "PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub child_prefix: String,
}

#[derive(Parser, Debug, Clone)]
#[command(group(ArgGroup::new("duo")))]
#[command(arg_required_else_help(true))]
pub struct DuoArgs {
    #[command(flatten)]
    pub shared: SharedArgs,

    #[clap(required = true)]
    #[clap(short = '1')]
    #[clap(long = "sample-a")]
    #[clap(help = "Prefix of first sample VCF and spanning reads BAM files")]
    #[clap(value_name = "PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub a_prefix: String,

    #[clap(required = true)]
    #[clap(short = '2')]
    #[clap(long = "sample-b")]
    #[clap(help = "Prefix of second sample VCF and spanning reads BAM files")]
    #[clap(value_name = "PREFIX")]
    #[arg(value_parser = check_prefix_path)]
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
    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::from_default_env()
        .format(|buf, record| {
            let level = record.level();
            let mut style = buf.style();
            match record.level() {
                Level::Error => style.set_color(Color::Red),
                Level::Warn => style.set_color(Color::Yellow),
                Level::Info => style.set_color(Color::Green),
                Level::Debug => style.set_color(Color::Blue),
                Level::Trace => style.set_color(Color::Cyan),
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                style.value(level),
                record.args()
            )
        })
        .filter_level(filter_level)
        .init();
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
