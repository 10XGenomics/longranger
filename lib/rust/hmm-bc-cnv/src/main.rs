// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
extern crate docopt;
extern crate csv;
extern crate rust_htslib;
extern crate ndarray;
extern crate probability;
extern crate statrs;
extern crate itertools;
extern crate special;
extern crate assert;

#[macro_use]
extern crate serde_derive;
extern crate serde;


use docopt::Docopt;

pub mod hmm;
pub mod hmm_cnv;
pub mod hmm_cnv_asread;
pub mod hmm_cnv_read;
pub mod fragment;
pub mod math;
pub mod support;
pub mod bed;

use hmm_cnv_read::*;
use hmm_cnv::*;
use hmm_cnv_asread::*;



const USAGE: &'static str = "
HMM Barcode CNV caller
Usage:
  hmm-bc-cnv read   <bam> <badregion> <outcnv> [options]
  hmm-bc-cnv asread <fragment> <bam> <badregion> <outcnv> [options]
  hmm-bc-cnv bc     <fragment> <bam> <badregion> <outcnv> [options]
  hmm-bc-cnv (-h | --help)
  hmm-bc-cnv --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --binsize=<binsize>       Size of each bin [default: 10000].
  --numstates=<numstates>   Number of copy number states [default: 11].
  --fragver=<fragver>       Version of fragments.h5 file. 0 for the latest, 1 for the oldest and 2 for the second verison. [default: 0].
  --probchange=<probchange>     Probability that a copy number will be changed. [default: 1.0e-4].
  --minprob=<minprob>       Minimum state or trans probability. [default: 3.3e-4].
  --covbed=<covbed>         Bin level coverage file in bedgraph format. [default: barcode_coverage.bedgraph]
  --covbin=<covbin>         Final bin size for coverage output. [default: 10000]
  --primary-contigs=<file>  File containing list of primary contigs. Calls will be restricted to these contigs.
";




#[derive(Debug, Deserialize)]
struct Args {
    arg_fragment:   Option<String>, // the input fragments csv file
    arg_bam:        Option<String>, // input bam file
    arg_outcnv:     String, // output cnv file
    arg_badregion:  String,
    flag_binsize:   i32,
    flag_numstates: usize,
    flag_fragver:   i32,
    flag_probchange: f64,
    flag_minprob:   f64,
    flag_primary_contigs: Option<String>,
    cmd_read:       bool,
    cmd_asread:     bool,
    cmd_bc:         bool,
    flag_covbed:    String,
    flag_covbin:    usize,
}



fn main() {

    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.deserialize())
                            .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);

    let max_num_run = 1;

    if args.cmd_read {
        let mut cnv = CNVRead::new(args.arg_bam.unwrap(), args.flag_binsize, args.flag_numstates,
            args.flag_probchange, 2.0, args.arg_badregion, args.flag_primary_contigs);
        cnv.multiple_run(max_num_run);
        cnv.output_cnvs(args.arg_outcnv);
    } else if args.cmd_bc {
        let mut cnv = CNV::new(args.arg_fragment.unwrap(), args.arg_bam.unwrap(), args.flag_binsize,
            args.flag_numstates, args.flag_fragver, args.flag_probchange, 2.0, args.arg_badregion, args.flag_primary_contigs,
            args.flag_minprob);
        cnv.multiple_run(max_num_run);
        cnv.output_cnvs(args.arg_outcnv);
    } else if args.cmd_asread {
        let zoomoutfactor = if args.flag_covbin <= args.flag_binsize as usize {1usize} else { args.flag_covbin / args.flag_binsize as usize };
        println!("setup");
        let mut cnv = CNVAsread::new(args.arg_fragment.unwrap(), args.arg_bam.unwrap(), args.flag_binsize,
            args.flag_numstates, args.flag_fragver, args.flag_probchange, 2.0, args.arg_badregion, args.flag_primary_contigs,
            args.flag_minprob);
        println!("mr");
        cnv.multiple_run(max_num_run);
        println!("output");
        cnv.output_cnvs(args.arg_outcnv);
        cnv.write_out_bc_cov(args.flag_covbed, zoomoutfactor, args.flag_covbin);
    }
}
