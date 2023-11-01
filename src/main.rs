use clap::Parser;
use trgt_denovo::{
    cli::{init_verbose, Cli, Command},
    commands::trio,
    util::Result,
};

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    match cli.command {
        Command::Trio(args) => trio(args)?,
    }
    Ok(())
}
