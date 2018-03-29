// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::fs::File;
use std::cmp::min;
use std::ascii::AsciiExt;
use std::cmp::max;
use csv;
use std::str;
use std::str::FromStr;

use rust_htslib::bam;
use rust_htslib::bam::record::{Record, Aux, Cigar};
use bio::io::{fasta, bed};

use locus::{Locus, CoverageStat};
use event;
use asm_caller;
use Args;




pub fn get_haplotype(rec: &Record) -> Option<u8> {
    match rec.aux(b"HP") {
        Some(Aux::Integer(hp)) => Some(hp as u8),
        _ => None,
    }
}

pub fn get_haplotype_qual(rec: &Record) -> Option<u8> {
    match rec.aux(b"PC") {
        Some(Aux::Integer(q)) => Some(q as u8),
        _ => None,
    }
}

fn chrom_len(chrom: &str, fa: &mut fasta::IndexedReader<File>) -> u64 {
    for s in fa.index.sequences() {
        if s.name == chrom {
            return s.len;
        }
    }

    0
}

pub fn read_locus(fa: &mut fasta::IndexedReader<File>,
                  loc: &Locus,
                  pad_left: u32,
                  pad_right: u32)
                  -> (Vec<u8>, usize) {
    let mut seq = Vec::new();

    let new_start = max(0, loc.start as i32 - pad_left as i32) as u64;
    let new_end = u64::from(min(loc.end + pad_right, chrom_len(&loc.chrom, fa) as u32));

    fa.read(&loc.chrom, new_start, new_end, &mut seq).unwrap();
    assert!(new_end - new_start <= seq.len() as u64);

    let slc = seq.as_mut_slice();
    let new_slc = slc.to_ascii_uppercase();
    (new_slc.into_iter().collect(), new_start as usize)
}

pub fn call_one(args: &Args) -> (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>) {
    let mut fa = fasta::IndexedReader::from_file(&args.arg_fasta).unwrap();
    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");

    let locus = FromStr::from_str(&args.arg_locus.clone().unwrap()).unwrap();
    let call = asm_caller::call_locus(&mut bam, &mut fa, &locus, args);
    info!("Got call: {:?}", call);
    call
}


pub fn call_bed(args: &Args) -> (usize, usize) {

    let mut fa = fasta::IndexedReader::from_file(&args.arg_fasta).expect("error opening fa");
    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");
    let mut bed = bed::Reader::from_file(&args.arg_bed).expect("couldn't open BED file");

    let which = args.arg_which;

    let mut calls = Vec::new();

    let mut num_tries = 0;
    let mut num_calls = 0;

    // Compute the coverage mean and stdev of coverage if coverage json is specified
    let coverage_stat = if let Some(ref path) = args.flag_coverage_json {
        Some(CoverageStat::from_coverage_json(path))
    } else {
        None
    };

    for (idx, r) in bed.records().enumerate() {

        if which.is_some() && which != Some(idx) {
            continue;
        }

        let rec = r.expect("BED reading error");
        info!("Processing locus: {}:{}-{}, line:{}",
              rec.chrom(),
              rec.start(),
              rec.end(),
              idx);

        if rec.end() - rec.start() > 4000 {
            info!("skipping long region: {}:{}-{}",
                  rec.chrom(),
                  rec.start(),
                  rec.end());
            continue;
        }

        if !fa.index.sequences().iter().any(|s| s.name == rec.chrom()) {
            info!("skipping non-existent chrom: {}:{}-{}",
                  rec.chrom(),
                  rec.start(),
                  rec.end());
            continue;
        }

        let locus = Locus {
            chrom: rec.chrom().to_string(),
            start: rec.start() as u32,
            end: rec.end() as u32,
        };

        // If the coverage-json file argument is specified, bail out if the coverge is
        // excessive in the locus. 
        if let Some(ref stat) = coverage_stat {
            if locus.has_excessive_coverage(&mut bam, stat) {
                info!("skipping locus with excessive coverage: {}:{}-{}",
                  locus.chrom,
                  locus.start,
                  locus.end);
                continue;
            }
        }

        let call = asm_caller::call_locus(&mut bam, &mut fa, &locus, args);

        let msg = if call == (None, None) {
            "(No call)"
        } else {
            "(Got call)"
        };

        info!("Result at {:?}. {}: {:?}", idx, msg, call);

        let evs = event::event_to_bed_pe(rec.chrom(), call);
        if !evs.is_empty() {
            num_calls += 1;
        }
        num_tries += 1;
        calls.extend(evs);
    }

    info!("Made {} calls in {} candidates", num_calls, num_tries);

    // Write calls to BEDPE
    if let Some(f) = args.arg_out.clone() {
        let mut w = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(f)
            .unwrap();
        for c in calls {
            w.serialize(c).unwrap();
        }
    }

    (num_calls, num_tries)
}


pub fn len(cig: &Cigar) -> u32 {
    match *cig {
        Cigar::Match(l) => l,
        Cigar::Ins(l) => l,
        Cigar::Del(l) => l,
        Cigar::RefSkip(l) => l,
        Cigar::SoftClip(l) => l,
        Cigar::HardClip(l) => l,
        Cigar::Pad(l) => l,
        Cigar::Equal(l) => l,
        Cigar::Diff(l) => l,
        Cigar::Back(l) => l,
    }
}

/// Find the position in the read of `ref_pos`.
pub fn find_ref_pos_in_read(ref_pos: u32, rec: &Record) -> usize {

    let mut read_pos = 0;
    let mut tpl_pos = rec.pos() as u32;

    if tpl_pos > ref_pos + 300 {
        return 0;
    }

    // Avoid wrap-around at start of contig
    //if tpl_pos < max(0, ref_pos as isize - 300) as u32 {
    //    return rec.seq().len();
    //}

    for cig_elem in &rec.cigar() {

        match *cig_elem {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                read_pos += l;
                tpl_pos += l;
            }

            Cigar::Ins(l) => read_pos += l,

            Cigar::Del(l) => tpl_pos += l,
            Cigar::SoftClip(l) => read_pos += l,
            _ => (),
        }

        if tpl_pos > ref_pos {
            return max(0_i32, (read_pos as i32) - (tpl_pos - ref_pos) as i32) as usize;
        }
    }
    rec.seq().len()
}


/// Return sequence of rec that aligns within the requested locus
#[inline(never)]
pub fn get_seq_bounded(rec: &Record, locus: &Locus) -> Vec<u8> {
    let r_start = find_ref_pos_in_read(locus.start, rec);
    let mut r_end = find_ref_pos_in_read(locus.end, rec);

    // FIXME: this is a hack to prevent the interval from going negative
    // there is an issue in find_ref_pos_in_read
    r_end = max(r_end, r_start);

    let rec_seq = rec.seq();
    let mut seq = Vec::with_capacity(r_end - r_start);
    for i in r_start..r_end {
        seq.push(rec_seq[i]);
    }

    seq
}

#[inline(never)]
pub fn get_qual_bounded(rec: &Record, locus: &Locus) -> Vec<u8> {
    let r_start = find_ref_pos_in_read(locus.start, rec);
    let r_end = find_ref_pos_in_read(locus.end, rec);

    let rec_qual = rec.qual();
    let mut qual = Vec::with_capacity(r_end - r_start);
    for i in r_start..r_end {
        qual.push(rec_qual[i]);
    }

    qual
}

pub fn print_seq(s: &Vec<u8>) -> String {
    let mut bs = String::new();
    for b in s.iter() {
        bs.push(*b as char);
    }
    bs
}

#[cfg(test)]
mod test {
    use Args;

    extern crate simple_logger;
    pub fn setup_test_logger() {
        simple_logger::init().unwrap();
    }

    /// Test for a bug which caused the coordinates of an event near the
    /// start of a contig to underflow & generate a call at 2^32 - x
    #[test]
    pub fn underflow_test() {
        setup_test_logger();

        let args = Args {
            flag_trace: true,
            flag_debug: true,
            flag_coverage_json: None,
            arg_locus: Some("chr10:0-1500".to_string()),
            arg_bed: "".to_string(),
            arg_fasta: "test/wrap_test.fa".to_string(),
            arg_bam: "test/wrap_test.bam".to_string(),
            arg_out: None,
            arg_which: None,
            flag_min_size: None,
            flag_asm: true,
            cmd_call_one: true,
            cmd_call_bed: false,
            cmd_bam_svs: false,
            cmd_cands: false,
            cmd_validate_one: false,
            cmd_validate_bedpe: false,
            flag_het_read_prob: None,
            flag_min_kmer_obs: Some(3),
        };

        let res = super::call_one(&args);
        println!("{:?}", res);
        if res.1.is_some() {
            assert!((res.1.unwrap().0).1 < 10000);
        }
    }
}
