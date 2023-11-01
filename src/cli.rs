use crate::util::Result;
use anyhow::anyhow;
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use env_logger::fmt::Color;
use log::{Level, LevelFilter};
use once_cell::sync::Lazy;
use std::{
    io::Write,
    path::{Path, PathBuf},
};

pub static FULL_VERSION: Lazy<String> = Lazy::new(|| {
    format!(
        "{}-{}",
        env!("CARGO_PKG_VERSION"),
        env!("VERGEN_GIT_DESCRIBE")
    )
});

#[derive(Parser)]
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

#[derive(Subcommand)]
pub enum Command {
    Trio(TrioArgs),
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("trio")))]
#[command(arg_required_else_help(true))]
pub struct TrioArgs {
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

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "out")]
    #[clap(help = "Output csv path")]
    #[clap(value_name = "CSV")]
    #[arg(value_parser = check_prefix_path)]
    pub output_path: String,

    #[clap(long = "trid")]
    #[clap(help = "TRID of a specific repeat to analyze, should be in the BED file")]
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
    #[clap(long = "parental-quantile")]
    #[clap(
        help = "Quantile of alignments scores to determine the parental threshold (default is strict and takes only the top score)"
    )]
    #[clap(value_name = "QUANTILE")]
    #[clap(default_value = "1.0")]
    #[arg(value_parser = parse_quantile)]
    pub parent_quantile: f64,
}

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

fn check_prefix_path(s: &str) -> Result<String> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(anyhow!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(s.to_string())
}

fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{}` is not a valid thread number", s))?;
    if thread <= 0 {
        return Err(anyhow!("Number of threads must be >= 1"));
    }
    Ok(thread)
}

fn parse_quantile(s: &str) -> Result<f64> {
    let value = s
        .parse::<f64>()
        .map_err(|e| anyhow!("Could not parse float: {}", e))?;
    if value < 0.0 || value > 1.0 {
        Err(anyhow!(
            "The value must be between 0.0 and 1.0, got: {}",
            value
        ))
    } else {
        Ok(value)
    }
}

fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        return Err(anyhow!("File does not exist: {}", path.display()));
    }
    Ok(path.to_path_buf())
}
