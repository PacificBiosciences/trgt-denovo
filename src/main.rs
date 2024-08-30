use clap::Parser;
use trgt_denovo::{
    cli::{init_verbose, Cli, Command},
    commands::{duo, trio},
    util::Result,
};

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    log::trace!("CLI options set: {:?}", cli);
    match cli.command {
        Command::Trio(args) => trio(args)?,
        Command::Duo(args) => duo(args)?,
    }
    Ok(())
}
