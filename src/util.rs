use anyhow::anyhow;
use log;
use std::path::Path;

pub type Result<T> = anyhow::Result<T>;

pub fn handle_error_and_exit(err: anyhow::Error) -> ! {
    log::error!("{:#}", err);
    std::process::exit(1);
}

pub fn try_exists(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(anyhow!("Path does not exist: {}", path.display()));
    }
    Ok(())
}
