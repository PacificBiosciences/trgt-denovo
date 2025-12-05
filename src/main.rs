use std::time::Instant;

use clap::Parser;
use trgt_denovo::{
    cli::{init_verbose, Cli, Command, FULL_VERSION},
    commands::{duo, trio},
    util::Result,
};

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        FULL_VERSION,
        cli.command.name()
    );
    let start_timer = Instant::now();
    log::trace!("CLI options set: {:?}", cli);
    match cli.command {
        Command::Trio(args) => trio(args)?,
        Command::Duo(args) => duo(args)?,
    }
    log::info!("Total execution time: {:.2?}", start_timer.elapsed());
    log::info!("{} end", env!("CARGO_PKG_NAME"));

    Ok(())
}
