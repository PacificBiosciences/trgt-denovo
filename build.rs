use std::error::Error;
use vergen_gitcl::{Emitter, GitclBuilder};

fn main() -> Result<(), Box<dyn Error>> {
    println!("cargo::rustc-check-cfg=cfg(has_git_describe)");

    let mut git = GitclBuilder::default();
    git.describe(false, true, Some("ThisPatternShouldNotMatchAnythingEver"));

    let emitter_res = Emitter::default()
        .fail_on_error()
        .quiet()
        .add_instructions(&git.build()?)
        .and_then(|e| e.emit());

    if emitter_res.is_ok() {
        println!("cargo:rustc-cfg=has_git_describe");
    }

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");
    println!("cargo:rerun-if-changed=.git/refs/tags");

    Ok(())
}
