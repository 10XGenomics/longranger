// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use rust_htslib::bam::{self, Read, IndexedReader};
use rust_htslib::bam::record::{Record, Cigar};
use serde::Serialize;

use ndarray::{Array, Ix};
use call::get_haplotype;
use locus::Locus;
use std::cmp::{max, min};
use std::str::FromStr;
use hmm::{HmmModel, viterbi};
use probability::distribution::{Binomial, Discrete};
use std::ops::Range;
use event::BedpeRow;

use std::path::Path;

use csv;
use Args;

struct PosRecord {
    chrom: String,
    pos: usize,
    hap1_coverage: u32,
    hap2_coverage: u32,
    unphased_coverage: u32,
    total_coverage: u32,

    clip_left: u32,
    clip_right: u32,

    is_sum: u32,
    is_sum_sq: u32,
    is_n: u32,
}

pub fn write_table<P: AsRef<Path>, T: Serialize>(data: Vec<T>, path: P) {
    let mut w = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();

    for d in data {
        w.serialize(d).unwrap();
    }
}


pub fn go(_args: &Args) {
    let args = _args.clone();

    let locus = FromStr::from_str(&args.arg_locus.clone().unwrap()).unwrap();
    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");
    let candidates = find_del_candidates(&mut bam, &locus, args.flag_het_read_prob.unwrap_or(0.97));

    let cb = candidates
        .iter()
        .map(|range| BedpeRow::simple(&locus.chrom, range))
        .collect();
    write_table(cb, args.arg_out.unwrap());
}

#[derive(PartialEq, Eq, Clone, Copy, Debug, Ord, PartialOrd)]
pub enum Edit {
    Match,
    Mismatch,
    Insertion,
    Deletion,
    Clip,
}

impl Edit {
    pub fn extends_read(self) -> bool {
        self == Edit::Match || self == Edit::Mismatch || self == Edit::Insertion
    }

    pub fn extends_tpl(self) -> bool {
        self == Edit::Match || self == Edit::Mismatch || self == Edit::Deletion
    }

    pub fn is_indel(self) -> bool {
        self == Edit::Insertion || self == Edit::Deletion
    }
}

#[derive(Clone,Copy,Debug, PartialEq, Eq, Ord, PartialOrd)]
pub struct AlignedSegment {
    tpl_pos: u32,
    read_pos: u32,
    length: u32,
    edit: Edit,
}

impl AlignedSegment {
    pub fn read_end(&self) -> u32 {
        self.read_pos +
        if self.edit.extends_read() {
            self.length
        } else {
            0
        }
    }

    pub fn tpl_end(&self) -> u32 {
        self.tpl_pos +
        if self.edit.extends_tpl() {
            self.length
        } else {
            0
        }
    }

    pub fn read_interval(&self) -> (u32, u32) {
        (self.read_pos, self.read_pos + self.length)
    }

    pub fn tpl_interval(&self) -> (u32, u32) {
        (self.tpl_pos, self.tpl_pos + self.length)
    }

    pub fn is_indel(&self) -> bool {
        self.edit == Edit::Insertion || self.edit == Edit::Deletion
    }

    pub fn cmp_edit(&self, other: &AlignedSegment) -> bool {
        self.tpl_pos == other.tpl_pos && self.length == other.length && self.edit == other.edit
    }
}


/// Convert a BAM alignment into a sequence of `AlignedSegments`
pub fn aligned_pairs(rec: &Record) -> Vec<AlignedSegment> {

    let mut read_pos = 0;
    let mut tpl_pos = rec.pos() as u32;
    let mut segs = Vec::new();

    for cig_elem in &rec.cigar() {
        match *cig_elem {
            Cigar::Match(l) => {
                let seg = AlignedSegment {
                    tpl_pos: tpl_pos,
                    read_pos: read_pos,
                    length: l,
                    edit: Edit::Match,
                };
                segs.push(seg);
                tpl_pos += l;
                read_pos += l;
            }

            Cigar::Ins(l) => {
                let seg = AlignedSegment {
                    tpl_pos: tpl_pos,
                    read_pos: read_pos,
                    length: l,
                    edit: Edit::Insertion,
                };
                segs.push(seg);
                read_pos += l;
            }

            Cigar::Del(l) => {
                let seg = AlignedSegment {
                    tpl_pos: tpl_pos,
                    read_pos: read_pos,
                    length: l,
                    edit: Edit::Deletion,
                };
                segs.push(seg);
                tpl_pos += l;
            }

            Cigar::SoftClip(l) |
            Cigar::HardClip(l) => read_pos += l,
            _ => (),
        }
    }

    segs
}

pub fn incr_cov(events: &Vec<AlignedSegment>, array: &mut Array<u32, Ix>, offset: u32) {
    for ev in events {
        if ev.edit == Edit::Match {
            let (tstart, tend) = ev.tpl_interval();

            let start = min(tstart.saturating_sub(offset) as isize,
                            array.shape()[0] as isize);
            let end = min(tend.saturating_sub(offset) as isize,
                          array.shape()[0] as isize);

            if end > 0 || start > array.shape()[0] as isize {
                let mut slc = array.slice_mut(s![start..end]);
                slc += 1;
            }
        }
    }
}

struct PhasedCovSummary {
    h1: Array<u32, Ix>,
    h2: Array<u32, Ix>,
    unphased: Array<u32, Ix>,
    total: Array<u32, Ix>,
    het_read_prob: f64,
    min_prob: f64,
    trans_prob: f64,
}

impl HmmModel for PhasedCovSummary {
    fn num_positions(&self) -> usize {
        self.h1.shape()[0]
    }

    fn num_states(&self) -> usize {
        2
    }

    fn prior(&self, state: usize) -> f64 {
        if state == 0 { 0.99 } else { 0.01 }
    }

    fn score_state(&self, pos: usize, state: usize) -> f64 {
        let n = (self.h1[pos] + self.h2[pos]) as usize;
        let x = max(self.h1[pos], self.h2[pos]) as usize;

        let p = if state == 0 { 0.65 } else { self.het_read_prob };
        Binomial::new(n, p).mass(x).max(self.min_prob)
    }

    fn score_trans(&self, _: usize, prev_state: usize, next_state: usize) -> f64 {
        if prev_state != next_state {
            self.trans_prob
        } else {
            1.0 - self.trans_prob
        }
    }
}

pub fn find_intervals<T, I: Iterator<Item = T>, F: Fn(T) -> bool>(iter: I,
                                                                  f: F)
                                                                  -> Vec<Range<u32>> {

    let mut segs: Vec<Range<u32>> = Vec::new();

    let mut active_interval = None;
    for (_idx, v) in iter.enumerate() {
        let idx = _idx as u32;
        let is_true = f(v);

        active_interval = match (is_true, active_interval) {
            (false, None) => None,
            (false, Some(interval)) => {
                segs.push(interval);
                None
            }
            (true, None) => Some(idx..idx + 1),
            (true, Some(rng)) => Some(rng.start..idx + 1),
        }
    }

    match active_interval {
        Some(interval) => segs.push(interval),
        _ => (),
    }

    segs
}

pub fn coverage_summary(bam: &mut IndexedReader,
                        locus: &Locus)
                        -> (Array<u32, Ix>, Array<u32, Ix>, Array<u32, Ix>, Array<u32, Ix>) {

    let locus_size = (locus.end - locus.start) as usize;
    let offset = locus.start;

    let mut hap1 = Array::from_elem((locus_size), 0);
    let mut hap2 = Array::from_elem((locus_size), 0);
    let mut unphased = Array::from_elem((locus_size), 0);
    let mut total = Array::from_elem((locus_size), 0);

    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    bam.fetch(tid, locus.start, locus.end - 1)
        .expect("Error fetching BAM file.");

    for _rec in bam.records() {
        let rec = _rec.unwrap();

        if rec.mapq() < 20 {
            continue;
        }

        let events = aligned_pairs(&rec);

        incr_cov(&events, &mut total, offset);

        match get_haplotype(&rec) {
            None => incr_cov(&events, &mut unphased, offset),
            Some(1) => incr_cov(&events, &mut hap1, offset),
            Some(2) => incr_cov(&events, &mut hap2, offset),
            _ => unreachable!("unrecognized haplotype"),
        }
    }

    (hap1, hap2, unphased, total)
}

const MIN_GAP: u32 = 200;

pub fn find_del_candidates(bam: &mut IndexedReader,
                           locus: &Locus,
                           het_read_prob: f64)
                           -> Vec<Range<u32>> {
    let (h1, h2, un, total) = coverage_summary(bam, locus);

    let cov_model = PhasedCovSummary {
        h1: h1,
        h2: h2,
        unphased: un,
        total: total.clone(),
        het_read_prob: het_read_prob,
        min_prob: (-10.0 as f64).exp(),
        trans_prob: 1e-5,
    };


    let vit = viterbi(cov_model);
    let het_dels: Vec<_> = find_intervals(vit.iter(), |x| *x == 1)
        .into_iter()
        .map(|int| (int.start + locus.start..int.end + locus.start))
        .collect();
    println!("het dels: {:?}", het_dels);

    let min_cov = 2;
    let hom_dels: Vec<_> = find_intervals(total.iter(), |x| *x < min_cov)
        .into_iter()
        .map(|int| (int.start + locus.start..int.end + locus.start))
        .collect();
    println!("hom dels: {:?}", hom_dels);

    let mut dels: Vec<Range<u32>> = het_dels.into_iter().chain(hom_dels).collect();
    dels.sort_by_key(|x| x.start);


    let mut merged_dels = Vec::new();

    if dels.len() > 0 {
        let mut last_locus = dels[0].clone();
        for locus in dels.into_iter().skip(1) {
            if locus.end - last_locus.start > MIN_GAP {
                merged_dels.push(last_locus.clone());
                last_locus = locus;
            } else {
                last_locus = min(last_locus.start, locus.start)..max(last_locus.end, locus.end);
            }
        }
        merged_dels.push(last_locus);
    }

    println!("all dels: {:?}", merged_dels);
    merged_dels
}






#[inline(never)]
pub fn mean_phased_coverage(bam: &mut IndexedReader, locus: &Locus) -> [f32; 3] {

    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");

    let locus_size = (locus.end - locus.start) as usize;
    let mut hap1 = 0.0;
    let mut hap2 = 0.0;
    let mut unphased = 0.0;

    for _pi in bam.pileup() {
        let pi = _pi.unwrap();

        for al in pi.alignments() {
            let rec = al.record();
            match get_haplotype(&rec) {
                None => unphased += 1.0,
                Some(1) => hap1 += 1.0,
                Some(2) => hap2 += 1.0,
                _ => unreachable!("unrecognized haplotype"),
            }
        }
    }

    [hap1 / (locus_size as f32),
     hap2 / (locus_size as f32),
     unphased / (locus_size as f32)]
}




pub fn clip_map(bam: &IndexedReader, locus: &Locus, bandwidth: i32) -> Array<u32, (Ix, Ix)> {

    let locus_size = (locus.end - locus.start) as usize;
    let mut clips = Array::from_elem((locus_size, 2), 0);

    {
        let mut add = |clip, col| match clip {
            Some((len, tpl_pos)) if len > 5 => {

                let start = max(0, tpl_pos - bandwidth - locus.start as i32) as isize;
                let end = min((locus_size - 1) as usize,
                              (tpl_pos + bandwidth) as usize - locus.start as usize) as
                          isize;

                let mut view = clips.slice_mut(s![start..end, col..col + 1]);
                view += 1;
            }
            _ => (),
        };

        for _rec in bam.records() {
            let rec = _rec.ok().unwrap();

            let (start, end) = clip_positions(&rec);
            add(start, 0);
            add(end, 1);
        }
    }

    clips
}
/*
/// Summarize insert sizes over locus
pub fn insert_size_map(bam: &IndexedReader, locus: &Locus, bandwidth: i32) -> (Array<f32, Ix>, Array<f32, Ix>, Array<f32, Ix>) {

    let locus_size = (locus.end - locus.start) as usize;
    let l_start = locus.start as i32;

    let mut sum_dist = Array::from_elem((locus_size), 0_f32);
    let mut sum_sq_dist = Array::from_elem((locus_size), 0_f32);
    let mut weight = Array::from_elem((locus_size), 0_f32);

    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();

        if !rec.is_first_in_pair() || rec.insert_size() > 50000 {
            continue
        }

        let pos = rec.pos();
        let insert_size = min(1500, rec.insert_size());
        let center = pos + insert_size / 2;

        let start = pos as isize;
        let end =  min((locus_size-1) as i32, pos + insert_size) as isize;

        let mut v = sum_dist.slice_mut(s![start..end]);
        v += insert_size as f32;

        let mut v = sum_sq_dist.slice_mut(s![start..end]);
        v += (insert_size * insert_size) as f32;

        let mut v = weight.slice_mut(s![start..end]);
        v += 1.0;
    }

    (sum_dist, sum_sq_dist, weight)
}
*/


pub fn clip_positions(rec: &Record) -> (Option<(i32, i32)>, Option<(i32, i32)>) {
    let mut start_clip = None;
    let mut end_clip = None;
    let mut have_match = false;

    let mut tpl_pos = rec.pos() as i32;

    for cig_elem in &rec.cigar() {
        match *cig_elem {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                tpl_pos += l as i32;
                have_match = true;
            }

            Cigar::Del(l) => tpl_pos += l as i32,

            Cigar::SoftClip(l) => {
                if !have_match {
                    start_clip = Some((l as i32, tpl_pos));
                } else {
                    end_clip = Some((l as i32, tpl_pos));
                }
            }
            _ => (),
        }
    }

    (start_clip, end_clip)
}


/*
pub fn collect_stats(bam: &IndexedReader, locus: &Locus) {
    let (h1, h2, un, total) = coverage_map(bam, locus);
    let clips = clip_map(bam, locus, 10);
    let (is_sum, is_sum_sq, mut is_n) = insert_size_map(bam, locus, 100);

    let mean_cov = total.mean(Axis(0));


    // Insert size z-score
    is_n += 1.0;
    let mean_is = &is_sum / &is_n;
    let overall_mean_is = mean_is.mean(Axis(0));
    let std_is = (is_sum_sq - (&is_sum * &is_sum)) / (&is_n);
    let insert_length_z = (mean_is - overall_mean_is) / std_is;


    // Het deletion emission scores
    //let p_no_del =
    //let p_del =

}
*/
