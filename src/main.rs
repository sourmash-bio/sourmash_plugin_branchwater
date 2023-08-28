use std::path::PathBuf;

use anyhow::Result;

pub mod utils;
pub mod manysearch;

use manysearch::manysearch;

use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// List of queries (one sig path per line in the file)
    #[clap(parse(from_os_str))]
    querylist: PathBuf,

    /// List of signatures to search
    #[clap(parse(from_os_str))]
    siglist: PathBuf,

    /// ksize
    #[clap(short, long, default_value = "31")]
    ksize: u8,

    /// threshold
    #[clap(short, long, default_value = "0.85")]
    threshold: f64,

    /// scaled
    #[clap(short, long, default_value = "1000")]
    scaled: usize,

    /// The path for output
    #[clap(parse(from_os_str), short, long)]
    output: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let opts = Cli::parse();

    manysearch(
        opts.querylist,
        opts.siglist,
        opts.threshold,
        opts.ksize,
        opts.scaled,
        opts.output,
    )?;

    Ok(())
}
