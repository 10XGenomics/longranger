// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
extern crate rust_htslib;
extern crate docopt;
extern crate regex;
extern crate csv;
extern crate itertools;
extern crate bincode;

#[cfg(test)]
extern crate tempfile;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;

mod aux;
mod ec_las_fg;
mod bc_correction;
mod bed;
mod correction;
mod fragment;
mod fx;
mod pval;
mod thread_iterator;

use docopt::Docopt;
use ec_las_fg::BLF;
use aux::*;

const USAGE: &'static str = "
Report Singe Partition

Usage:
  report_single_partition <bcwl> <bccount> <bam> <contigs> <left7> <right9> <lectable> <rectable> <outbam>  [--rldist=<rldist>] [--pval=<pval>]  [--freq=<freq>] [--out=<out>] [--mapq=<mapq>] [--maxerr=<maxerr>] [--nofragmerg] [--mappingbc] [--pseudo=<pseudo>] [--target=<target>] [--centromere=<centromere>] [--mapqbc=<mapqbc>] [--readbuffsz=<readbuffsz>] [--allowindel]
  report_single_partition (-h | --help)
  report_single_partition --version

Options:
  -h --help             Show this screen.
  --version             Show version.
  --rldist=<rldist>     Read Link Distance [default: 60000].
  --mapq=<mapq>         Mapq threshold below which reads will not serve as fragment boundaries [default: 30].
  --pval=<pval>         Pval threshold for fragment merging [default: 1e-3].
  --freq=<freq>         Frequency threshold for fragment merging [default: 2].
  --out=<out>           Output file for the final fragments [default: fragments.csv].
  --maxerr=<maxerr>     Max number of errors allowed in UMI correction [default: 2].
  --nofragmerg          Skipping fragment merging step.
  --mappingbc           Use mapping information for barcode correction.
  --pseudo=<pseudo>     Psudo count added to each barcode count [default: 1.0].
  --target=<target>     Target file for WES [default: ].
  --centromere=<centromere>     Centromere file [default: ].
  --mapqbc=<mapqbc>     Mapping quality above which mapping information will be used for barcode correction [default: 30].
  --readbuffsz=<readbuffsz>     Size of the read buffer [default: 100000].
  --allowindel          Whether indel is allowed in barcode error correction.
";

#[derive(Debug, Deserialize)]
struct Args {
    arg_bccount: String,
    arg_bcwl: String,
    arg_left7: String, // left 7-mer code
    arg_right9: String, // right 9-mer code
    arg_lectable: String,
    arg_rectable: String,
    arg_bam: String,
    arg_contigs: String,
    arg_outbam: String,
    flag_rldist: i32,
    flag_mapq: u8,
    flag_pval: f64,
    flag_freq: i32,
    flag_out: String,
    flag_maxerr: i32,
    flag_pseudo: f64,
    flag_nofragmerg: bool,
    flag_mappingbc: bool,
    flag_centromere: String,
    flag_target: String,
    flag_mapqbc: u8,
    flag_readbuffsz: usize,
    flag_allowindel: bool,
}


fn main() {

    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.deserialize())
                            .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);

    get_barcode_segment(&args.arg_bcwl, &args.arg_left7, &args.arg_right9);
    get_ec_table(&args.arg_left7, &args.arg_lectable, false, args.flag_allowindel);
    get_ec_table(&args.arg_right9, &args.arg_rectable, false, args.flag_allowindel);

    let mut blf = BLF::new(args.flag_out,
        &args.flag_centromere, args.arg_bcwl, args.arg_lectable, args.arg_rectable,
        args.arg_bccount, 0.01, 0.95,
        args.flag_maxerr, args.flag_pseudo, args.flag_rldist, args.flag_mapq,
        args.flag_mapqbc, args.flag_readbuffsz, args.flag_mappingbc);


    blf.find_raw_info(args.arg_bam, args.arg_outbam, args.arg_contigs, &args.flag_target, args.flag_nofragmerg, args.flag_pval, args.flag_freq);
    blf.print_stats();

}
