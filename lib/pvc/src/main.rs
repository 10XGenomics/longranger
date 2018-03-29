// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#![allow(dead_code)]

// Utils
extern crate docopt;
extern crate fern;
extern crate regex;
extern crate csv;
extern crate time;
extern crate probability;
extern crate backtrace;
extern crate itertools;

#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate serde;

#[macro_use(s)]
extern crate ndarray;

#[macro_use]
extern crate log;

// Bio& BAM
extern crate bio;
extern crate rust_htslib;
extern crate debruijn;

pub mod call;
pub mod locus;
pub mod event;
pub mod detector;
pub mod hmm;
pub mod bam_event;
pub mod asm_caller;
pub mod validate;

use docopt::Docopt;
use backtrace::Backtrace;
use std::thread;
use std::panic;
use std::io;

const USAGE: &'static str = "
Phased/POA Variant Caller

Usage:
  pvc [options] call-one <fasta> <bam> <locus>
  pvc [options] call-bed -o <out> <fasta> <bam> <bed> [<which>]
  pvc [options] bam-svs <out> <bam>
  pvc [options] cands <bam> <locus> <out>
  pvc [options] validate-one <bam> <fasta> <locus>
  pvc [options] validate-bedpe <bam> <fasta> <bed>
  pvc (-h | --help)
  pvc --version

Options:
  --min-size=<m>       Mininum event size
  --min-kmer-obs=<k>   Minimum number of kmer observations
  --het-read-prob=<p>  Deletion detection parameter
  --asm                Validate events from an assembly (allows secondary mappings & uses larger windows)
  -h --help            Show this screen.
  --version            Show version.
  --trace              Trace logging
  -d --debug           Debug logging
  --coverage-json=<s>  Json file containing coverage distribution. If provided, regions with excessive coverage are filtered out
";

const VERSION: &'static str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Deserialize, Clone)]
pub struct Args {
    flag_min_size: Option<usize>,
    flag_min_kmer_obs: Option<usize>,
    flag_het_read_prob: Option<f64>,
    flag_asm: bool,
    flag_coverage_json: Option<String>,

    flag_trace: bool,
    flag_debug: bool,

    arg_locus: Option<String>,
    arg_bed: String,
    arg_fasta: String,
    arg_bam: String,
    arg_out: Option<String>,
    arg_which: Option<usize>,

    cmd_call_one: bool,
    cmd_call_bed: bool,
    cmd_bam_svs: bool,
    cmd_cands: bool,
    cmd_validate_one: bool,
    cmd_validate_bedpe: bool,
}

fn setup_log(level: log::LogLevelFilter) {

    let base_config = fern::Dispatch::new().level(level);

    let logger_config = fern::Dispatch::new()
        .format(|out, msg, record| {
                    out.finish(format_args!("[{}][{}] {}",
                                            time::now().strftime("%H:%M:%S").unwrap(),
                                            record.level(),
                                            msg))
                })
        .chain(fern::log_file("output.log").expect("couldn't open log file"))
        .chain(io::stdout());

    let cfg = base_config.chain(logger_config).apply();

    if let Err(e) = cfg {
        panic!("Failed to initialize global logger: {}", e);
    }
}


fn main() {

    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);

    let level = if args.flag_trace {
        log::LogLevelFilter::Trace
    } else if args.flag_debug {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    };

    setup_log(level);
    info!("pvc version: {}", VERSION);

    // Setup panic hook. If a stage panics, we'll shutdown cleanly to martian
    let p = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        let backtrace = Backtrace::new();

        let thread = thread::current();
        let thread = thread.name().unwrap_or("unnamed");

        let msg = match info.payload().downcast_ref::<&'static str>() {
            Some(s) => *s,
            None => {
                match info.payload().downcast_ref::<String>() {
                    Some(s) => &**s,
                    None => "Box<Any>",
                }
            }
        };

        let msg = match info.location() {
            Some(location) => {
                format!("thread '{}' panicked at '{}': {}:{}{:?}",
                        thread,
                        msg,
                        location.file(),
                        location.line(),
                        backtrace)
            }
            None => format!("thread '{}' panicked at '{}'{:?}", thread, msg, backtrace),
        };

        println!("{}", msg);
        p(info);
    }));


    if args.cmd_call_one {
        println!("locus: '{:?}'", args.arg_locus);
        call::call_one(&args);
    }

    if args.cmd_call_bed {
        println!("bed: {}", args.arg_bed);
        call::call_bed(&args);
    }

    if args.cmd_bam_svs {
        bam_event::go(&args);
    }

    if args.cmd_cands {
        detector::go(&args);
    }

    if args.cmd_validate_one {
        validate::go_one(&args);
    }

    if args.cmd_validate_bedpe {
        validate::go_bedpe(&args);
    }
}
