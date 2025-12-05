use crate::{
    model::{AlnScoring, QuickMode},
    util::Result,
};
use anyhow::anyhow;
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use log::{Level, LevelFilter};
use owo_colors::{
    colors::{Blue, Green, Magenta, Red, Yellow},
    OwoColorize, Stream, Style,
};
use std::{
    io::Write,
    ops::Deref,
    path::{Path, PathBuf},
};

#[cfg(has_git_describe)]
pub const FULL_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "-", env!("VERGEN_GIT_DESCRIBE"));

#[cfg(not(has_git_describe))]
pub const FULL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(name="trgt-denovo",
          author="Tom Mokveld <tmokveld@pacificbiosciences.com>\nEgor Dolzhenko <edolzhenko@pacificbiosciences.com>", 
          version=FULL_VERSION,
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

    /// Amount of additional flanking sequence that should be used during alignment (this should match --output-flank-len in TRGT)
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
#[command(arg_required_else_help(true))]
#[command(group(ArgGroup::new("mother_input").required(true).args(["mother_prefix", "mother_vcf"])))]
#[command(group(ArgGroup::new("father_input").required(true).args(["father_prefix", "father_vcf"])))]
#[command(group(ArgGroup::new("child_input").required(true).args(["child_prefix", "child_vcf"])))]
pub struct TrioArgs {
    #[command(flatten)]
    pub shared: SharedArgs,

    /// Prefix of mother VCF and spanning reads BAM files
    #[arg(
        short = 'm',
        long = "mother",
        value_name = "PREFIX",
        value_parser = check_prefix_path,
        conflicts_with_all = ["mother_vcf", "mother_bam"]
    )]
    pub mother_prefix: Option<String>,

    /// Path to mother VCF file (use with --mother-bam instead of --mother)
    #[arg(
        long = "mother-vcf",
        value_name = "VCF",
        value_parser = check_file_exists,
        requires = "mother_bam"
    )]
    pub mother_vcf: Option<PathBuf>,

    /// Path to mother BAM file (use with --mother-vcf instead of --mother)
    #[arg(
        long = "mother-bam",
        value_name = "BAM",
        value_parser = check_file_exists,
        requires = "mother_vcf"
    )]
    pub mother_bam: Option<PathBuf>,

    /// Prefix of father VCF and spanning reads BAM files
    #[arg(
        short = 'f',
        long = "father",
        value_name = "PREFIX",
        value_parser = check_prefix_path,
        conflicts_with_all = ["father_vcf", "father_bam"]
    )]
    pub father_prefix: Option<String>,

    /// Path to father VCF file (use with --father-bam instead of --father)
    #[arg(
        long = "father-vcf",
        value_name = "VCF",
        value_parser = check_file_exists,
        requires = "father_bam"
    )]
    pub father_vcf: Option<PathBuf>,

    /// Path to father BAM file (use with --father-vcf instead of --father)
    #[arg(
        long = "father-bam",
        value_name = "BAM",
        value_parser = check_file_exists,
        requires = "father_vcf"
    )]
    pub father_bam: Option<PathBuf>,

    /// Prefix of child VCF and spanning reads BAM files
    #[arg(
        short = 'c',
        long = "child",
        value_name = "PREFIX",
        value_parser = check_prefix_path,
        conflicts_with_all = ["child_vcf", "child_bam"]
    )]
    pub child_prefix: Option<String>,

    /// Path to child VCF file (use with --child-bam instead of --child)
    #[arg(
        long = "child-vcf",
        value_name = "VCF",
        value_parser = check_file_exists,
        requires = "child_bam"
    )]
    pub child_vcf: Option<PathBuf>,

    /// Path to child BAM file (use with --child-vcf instead of --child)
    #[arg(
        long = "child-bam",
        value_name = "BAM",
        value_parser = check_file_exists,
        requires = "child_vcf"
    )]
    pub child_bam: Option<PathBuf>,
}

#[derive(Parser, Debug, Clone)]
#[command(arg_required_else_help(true))]
#[command(group(ArgGroup::new("sample_a_input").required(true).args(["a_prefix", "a_vcf"])))]
#[command(group(ArgGroup::new("sample_b_input").required(true).args(["b_prefix", "b_vcf"])))]
pub struct DuoArgs {
    #[command(flatten)]
    pub shared: SharedArgs,

    /// Prefix of first sample VCF and spanning reads BAM files
    #[arg(
        short = '1',
        long = "sample-a",
        value_name = "PREFIX",
        value_parser = check_prefix_path,
        conflicts_with_all = ["a_vcf", "a_bam"]
    )]
    pub a_prefix: Option<String>,

    /// Path to first sample VCF file (use with --sample-a-bam instead of --sample-a)
    #[arg(
        long = "sample-a-vcf",
        value_name = "VCF",
        value_parser = check_file_exists,
        requires = "a_bam"
    )]
    pub a_vcf: Option<PathBuf>,

    /// Path to first sample BAM file (use with --sample-a-vcf instead of --sample-a)
    #[arg(
        long = "sample-a-bam",
        value_name = "BAM",
        value_parser = check_file_exists,
        requires = "a_vcf"
    )]
    pub a_bam: Option<PathBuf>,

    /// Prefix of second sample VCF and spanning reads BAM files
    #[arg(
        short = '2',
        long = "sample-b",
        value_name = "PREFIX",
        value_parser = check_prefix_path,
        conflicts_with_all = ["b_vcf", "b_bam"]
    )]
    pub b_prefix: Option<String>,

    /// Path to second sample VCF file (use with --sample-b-bam instead of --sample-b)
    #[arg(
        long = "sample-b-vcf",
        value_name = "VCF",
        value_parser = check_file_exists,
        requires = "b_bam"
    )]
    pub b_vcf: Option<PathBuf>,

    /// Path to second sample BAM file (use with --sample-b-vcf instead of --sample-b)
    #[arg(
        long = "sample-b-bam",
        value_name = "BAM",
        value_parser = check_file_exists,
        requires = "b_vcf"
    )]
    pub b_bam: Option<PathBuf>,
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

#[inline(always)]
fn level_style(level: Level) -> (&'static str, Style) {
    match level {
        Level::Error => ("ERROR", Style::new().fg::<Red>().bold()),
        Level::Warn => ("WARN", Style::new().fg::<Yellow>()),
        Level::Info => ("INFO", Style::new().fg::<Green>()),
        Level::Debug => ("DEBUG", Style::new().fg::<Blue>()),
        Level::Trace => ("TRACE", Style::new().fg::<Magenta>()),
    }
}

fn format_log(buf: &mut env_logger::fmt::Formatter, record: &log::Record) -> std::io::Result<()> {
    let (label, style) = level_style(record.level());
    let ts = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
    let painted_label = label.if_supports_color(Stream::Stderr, |t| style.style(t));
    writeln!(buf, "{ts} [{}] - {}", painted_label, record.args())
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Helper struct to create temporary test files
    struct TestFiles {
        reference: NamedTempFile,
        bed: NamedTempFile,
        mother_vcf: NamedTempFile,
        mother_bam: NamedTempFile,
        father_vcf: NamedTempFile,
        father_bam: NamedTempFile,
        child_vcf: NamedTempFile,
        child_bam: NamedTempFile,
    }

    impl TestFiles {
        fn new() -> Self {
            let mut reference = NamedTempFile::new().unwrap();
            writeln!(reference, ">chr1\nACGT").unwrap();
            let mut bed = NamedTempFile::new().unwrap();
            writeln!(bed, "chr1\t100\t200\ttest_trid").unwrap();
            Self {
                reference,
                bed,
                mother_vcf: NamedTempFile::new().unwrap(),
                mother_bam: NamedTempFile::new().unwrap(),
                father_vcf: NamedTempFile::new().unwrap(),
                father_bam: NamedTempFile::new().unwrap(),
                child_vcf: NamedTempFile::new().unwrap(),
                child_bam: NamedTempFile::new().unwrap(),
            }
        }
    }

    #[test]
    fn test_trio_with_prefix_args() {
        let files = TestFiles::new();
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "trio",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "-m",
            "mother_prefix",
            "-f",
            "father_prefix",
            "-c",
            "child_prefix",
        ]);

        assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
        if let Command::Trio(args) = result.unwrap().command {
            assert_eq!(args.mother_prefix, Some("mother_prefix".to_string()));
            assert_eq!(args.father_prefix, Some("father_prefix".to_string()));
            assert_eq!(args.child_prefix, Some("child_prefix".to_string()));
            assert!(args.mother_vcf.is_none());
            assert!(args.father_vcf.is_none());
            assert!(args.child_vcf.is_none());
        } else {
            panic!("Expected Trio command");
        }
    }

    #[test]
    fn test_trio_with_explicit_files() {
        let files = TestFiles::new();
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "trio",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "--mother-vcf",
            files.mother_vcf.path().to_str().unwrap(),
            "--mother-bam",
            files.mother_bam.path().to_str().unwrap(),
            "--father-vcf",
            files.father_vcf.path().to_str().unwrap(),
            "--father-bam",
            files.father_bam.path().to_str().unwrap(),
            "--child-vcf",
            files.child_vcf.path().to_str().unwrap(),
            "--child-bam",
            files.child_bam.path().to_str().unwrap(),
        ]);

        assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
        if let Command::Trio(args) = result.unwrap().command {
            assert!(args.mother_prefix.is_none());
            assert!(args.father_prefix.is_none());
            assert!(args.child_prefix.is_none());
            assert!(args.mother_vcf.is_some());
            assert!(args.mother_bam.is_some());
            assert!(args.father_vcf.is_some());
            assert!(args.father_bam.is_some());
            assert!(args.child_vcf.is_some());
            assert!(args.child_bam.is_some());
        } else {
            panic!("Expected Trio command");
        }
    }

    #[test]
    fn test_trio_prefix_and_explicit_conflict() {
        let files = TestFiles::new();
        // Providing both --mother and --mother-vcf should fail
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "trio",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "-m",
            "mother_prefix",
            "--mother-vcf",
            files.mother_vcf.path().to_str().unwrap(),
            "--mother-bam",
            files.mother_bam.path().to_str().unwrap(),
            "-f",
            "father_prefix",
            "-c",
            "child_prefix",
        ]);

        assert!(
            result.is_err(),
            "Should fail when both prefix and explicit files are provided"
        );
    }

    #[test]
    fn test_trio_explicit_requires_both_vcf_and_bam() {
        let files = TestFiles::new();
        // Providing --mother-vcf without --mother-bam should fail
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "trio",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "--mother-vcf",
            files.mother_vcf.path().to_str().unwrap(),
            // Missing --mother-bam
            "-f",
            "father_prefix",
            "-c",
            "child_prefix",
        ]);

        assert!(
            result.is_err(),
            "Should fail when --mother-vcf is provided without --mother-bam"
        );
    }

    #[test]
    fn test_trio_mixed_prefix_and_explicit() {
        let files = TestFiles::new();
        // Using prefix for some samples and explicit for others should work
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "trio",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "-m",
            "mother_prefix", // prefix for mother
            "--father-vcf",
            files.father_vcf.path().to_str().unwrap(),
            "--father-bam",
            files.father_bam.path().to_str().unwrap(), // explicit for father
            "-c",
            "child_prefix", // prefix for child
        ]);

        assert!(
            result.is_ok(),
            "Should allow mixing prefix and explicit across different samples: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_duo_with_prefix_args() {
        let files = TestFiles::new();
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "duo",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "-1",
            "sample_a_prefix",
            "-2",
            "sample_b_prefix",
        ]);

        assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
        if let Command::Duo(args) = result.unwrap().command {
            assert_eq!(args.a_prefix, Some("sample_a_prefix".to_string()));
            assert_eq!(args.b_prefix, Some("sample_b_prefix".to_string()));
            assert!(args.a_vcf.is_none());
            assert!(args.b_vcf.is_none());
        } else {
            panic!("Expected Duo command");
        }
    }

    #[test]
    fn test_duo_with_explicit_files() {
        let files = TestFiles::new();
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "duo",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "--sample-a-vcf",
            files.mother_vcf.path().to_str().unwrap(),
            "--sample-a-bam",
            files.mother_bam.path().to_str().unwrap(),
            "--sample-b-vcf",
            files.father_vcf.path().to_str().unwrap(),
            "--sample-b-bam",
            files.father_bam.path().to_str().unwrap(),
        ]);

        assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
        if let Command::Duo(args) = result.unwrap().command {
            assert!(args.a_prefix.is_none());
            assert!(args.b_prefix.is_none());
            assert!(args.a_vcf.is_some());
            assert!(args.a_bam.is_some());
            assert!(args.b_vcf.is_some());
            assert!(args.b_bam.is_some());
        } else {
            panic!("Expected Duo command");
        }
    }

    #[test]
    fn test_duo_explicit_requires_both_vcf_and_bam() {
        let files = TestFiles::new();
        // Providing --sample-a-vcf without --sample-a-bam should fail
        let result = Cli::try_parse_from([
            "trgt-denovo",
            "duo",
            "-r",
            files.reference.path().to_str().unwrap(),
            "-b",
            files.bed.path().to_str().unwrap(),
            "-o",
            "output.tsv",
            "--sample-a-vcf",
            files.mother_vcf.path().to_str().unwrap(),
            // Missing --sample-a-bam
            "-2",
            "sample_b_prefix",
        ]);

        assert!(
            result.is_err(),
            "Should fail when --sample-a-vcf is provided without --sample-a-bam"
        );
    }
}
