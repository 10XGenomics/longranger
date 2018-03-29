// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
/*
use rust_htslib::bam::{self, Read, IndexedReader};
use std::path::{PathBuf};

use call::{self, read_locus, get_haplotype, get_seq_bounded};
use locus::Locus;
use poa::{Poa};
use event;
use detector;

use bio::io::{fasta};
use std::fs::File;
use std::cmp::min;
use std::str;
use Args;

const COV_WINDOW: u32 = 200;
const FLANKING_COV: f32 = 4.0;
const MIN_OVERLAP: u32 = 25;
const MIN_COVERAGE: u32 = 4;
const MISMATCH: i32 = 5;

pub fn call_locus(bam: &mut bam::IndexedReader, fa: &mut fasta::IndexedReader<File>, locus: &Locus, args: &Args) ->
 (Option<(usize, usize)>, Option<(usize, usize)>) {
     let debug = args.flag_debug;

    let (ref_seq, ref_start_pos) = read_locus(fa, locus, 5, 5);
    let mut haps = build_hap_poa(bam, locus, args);

    let mut calls = [None, None];

    let start_loc = Locus {
        chrom: locus.chrom.clone(),
        start: locus.start,
        end: locus.start + (COV_WINDOW as u32)
    };

    let stop_loc = Locus {
        chrom: locus.chrom.clone(),
        start: locus.end - (COV_WINDOW as u32),
        end: locus.end
    };

    let cov_start = detector::mean_phased_coverage(bam, &start_loc);
    let cov_end  = detector::mean_phased_coverage(bam, &stop_loc);

    for (idx, hap) in haps.iter_mut().enumerate() {
        let hap_seq = hap.consensus();
        debug!("seq len: {}, ref len: {}", hap_seq.len(), ref_seq.len());
        let ev = event::find_deletion(&hap_seq, &ref_seq, ref_start_pos, args.flag_min_size);

        if cov_start[idx] > FLANKING_COV && cov_end[idx] > FLANKING_COV {
            calls[idx] = ev;
        } else {
            debug!("rejected call: {:?}", ev);
            println!("rejected call: {:?}", ev);
        }
    }

    if debug {
        haps[0].to_dot(PathBuf::from("h1.dot"));
        haps[1].to_dot(PathBuf::from("h2.dot"));
    }

    (calls[0], calls[1])
}


const MAX_READS: usize = 2000;

#[inline(never)]
pub fn build_hap_poa(bam: &mut IndexedReader, locus: &Locus, args: &Args) -> [Poa; 2] {
    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    let mut num_reads = 0;

    let hap1 = Poa::new_config(-MISMATCH as i16, MIN_OVERLAP as i16, MIN_COVERAGE as i16);
    let hap2 = Poa::new_config(-MISMATCH as i16, MIN_OVERLAP as i16, MIN_COVERAGE as i16);
    let mut haps = [hap1, hap2];
    let mut non_phased = Vec::new();

    // Add phased reads to haplotype buckets
    bam.seek(tid, locus.start, locus.end).ok().expect("Error seeking BAM file.");
    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();

        // Skip poor mappings & secondary alignments
        if rec.mapq() < 10 || rec.is_secondary() { continue; }

        match call::get_haplotype(&rec) {
            Some(hap) => {
                let mut h = haps.get_mut((hap - 1) as usize).unwrap();
                let chars = call::get_seq_bounded(&rec, locus);

                if chars.len() < 40 {
                    continue;
                }

                if h.num_reads < 2 {
                    h.add_read(&chars);
                } else {
                    let aln = h.try_align_read(&chars);

                    // Require a reasonable match to the available gra
                    if h.accept_aln(&aln) {
                        h.commit_alignment(&chars, aln);
                        let name = format!("hap{}_r{}.dot", hap, num_reads);
                        println!("added: {}, name:{}", str::from_utf8(rec.qname()).unwrap(), name);
                        h.to_dot(PathBuf::from(name));
                    } else {
                        println!("dropped: {}, score:{}", str::from_utf8(rec.qname()).unwrap(), aln.best_alignments.get()[0].score);
                        h.rejected_phased_reads += 1;
                    }
                }
            }
            None => non_phased.push(rec),
        }

        num_reads += 1;
        if num_reads > MAX_READS {
            break
        }
    }

    // Go back through and add unphased reads to the best-fitting bucket
    for rec in non_phased
    {
        match get_haplotype(&rec)
        {
            Some(_) => (),
            None => {
                let chars = get_seq_bounded(&rec, locus);

                //trace!("p: {:?}, sec: {:?}, sup: {:?}, cig: {:?}", rec.pos(), rec.is_secondary(), rec.is_supplementary(), rec.cigar());
                //trace!("{}", print_seq(&chars));

                if chars.len() < 40 {
                    continue;
                }

                let aln1 = {
                    let h1 = haps.get_mut(0).unwrap();
                    h1.try_align_read(&chars) };

                let aln2 = {
                    let h2 = haps.get_mut(1).unwrap();
                    h2.try_align_read(&chars) };

                let s1 = aln1.best_alignments.get()[0].score;
                let s2 = aln2.best_alignments.get()[0].score;

                if s1 > s2 && s1 > 24 {
                    haps.get_mut(0).unwrap().commit_alignment(&chars, aln1);
                }

                if  s1 < s2 && s2 > 24 {
                    haps.get_mut(1).unwrap().commit_alignment(&chars, aln2);
                }
            },
        }

        num_reads += 1;
        if num_reads > MAX_READS {
            break
        }
    }

    info!("Used reads: {}", num_reads);
    info!("Rejected reads: {},{}", haps[0].rejected_phased_reads, haps[1].rejected_phased_reads);
    haps
}


#[inline(never)]
pub fn build_single_poa(bam: &mut IndexedReader, locus: &Locus) -> Poa {
    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    let mut num_reads = 0;

    let mut hap = Poa::new();

    // Add phased reads to haplotype buckets
    bam.seek(tid, locus.start, locus.end).ok().expect("Error seeking BAM file.");
    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();

        if rec.mapq() < 5 || rec.is_secondary() { continue; }

        let chars = get_seq_bounded(&rec, locus);

        if chars.len() < 40 {
            continue;
        }

        if hap.num_reads < 2 {
            hap.add_read(&chars);
        } else {
            let aln = hap.try_align_read(&chars);
            let score = aln.best_alignments.get()[0].score;

            if score > min(24, (hap.num_nodes() << 1) as i16) {
                hap.commit_alignment(&chars, aln);
            }
        }

        num_reads += 1;
        if num_reads > MAX_READS {
            break
        }
    }

    info!("Used reads: {}", num_reads);
    hap
}
*/
