// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::banded;
use bio::alignment::pairwise::Scoring;
use hmm::{self, HmmModel};
use std::ops::Range;

#[derive(Copy, Clone, Eq, PartialEq)]
enum EvSt {
    Start,
    LeftMatch,
    DelMatch,
    RightMatch,
}

struct AlnHmm {
    p_trans: f64,
    obs: Vec<AlignmentOperation>,
}

impl HmmModel for AlnHmm {
    fn num_positions(&self) -> usize {
        self.obs.len()
    }

    fn num_states(&self) -> usize {
        3
    }

    // Must start in first state
    fn prior(&self, state: usize) -> f64 {
        if state == 0 { 1.0 } else { 0.0 }
    }

    fn score_state(&self, pos: usize, state: usize) -> f64 {
        let obs = self.obs[pos];

        let is_del = obs == AlignmentOperation::Del;

        if state == 0 || state == 2 {
            if is_del { 0.1 } else { 0.9 }
        } else {
            if is_del { 0.9 } else { 0.1 }
        }
    }

    // Only single-step forward transitions are allowed
    fn score_trans(&self, _: usize, prev_state: usize, state: usize) -> f64 {
        if prev_state == 2 && state == 2 {
            // Not normalized, but we want to keep this the same as 0 -> 0
            return 1.0 - self.p_trans;
        }

        if prev_state == state {
            1.0 - self.p_trans
        } else if state == prev_state + 1 {
            self.p_trans
        } else {
            0.0
        }
    }
}

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("unrecognized"),
    }
}

pub fn rc_seq(vec: &Vec<u8>) -> Vec<u8> {
    let mut res = Vec::new();

    for b in vec.iter().rev() {
        res.push(complement(*b));
    }

    res
}


pub fn find_deletion(query: &Vec<u8>,
                     target: &Vec<u8>,
                     ref_pos: usize,
                     min_size: Option<usize>)
                     -> Option<(usize, usize)> {
    // Need to make sure that the sparse DP gets the right match score
    let score = Scoring::from_scores(-50, -1, 20, -80);
    let mut aligner = banded::Aligner::with_scoring(score, 12, 15);

    let alignment1 = aligner.local(query, target);

    let query_rc = rc_seq(query);
    let alignment2 = aligner.local(&query_rc, target);

    let (alignment, query) = if alignment1.score > alignment2.score {
        (alignment1, query)
    } else {
        (alignment2, &query_rc)
    };

    if query.len() < 5 || alignment.x_aln_len() < 5 {
        return None;
    }

    let min_sz = min_size.unwrap_or(40);
    debug!("aln:\n{}", alignment.pretty(query, target));

    // Search for a deletion > 25bp in the query seq
    // Do some filtering for 'clean' deletions:
    // 1. Require 15bp of matches on each side of event
    // 2.

    let mut event = None;
    let mut tpl_pos = ref_pos + alignment.ystart;

    let hmm_model = AlnHmm {
        p_trans: 0.01,
        obs: alignment.operations.clone(),
    };
    let traceback = hmm::viterbi(hmm_model);
    //println!("{:?}", traceback);
    //println!("aln len: {:?}", alignment.operations.len());

    for (aln, tb) in alignment.operations.iter().zip(traceback) {

        event = match (event, tb) {
            (None, 1) => Some((tpl_pos, tpl_pos + 1)),
            (Some((start, _)), 1) => Some((start, tpl_pos + 1)),
            (_, _) => event,
        };

        match *aln {
            AlignmentOperation::Match |
            AlignmentOperation::Subst |
            AlignmentOperation::Del => tpl_pos += 1,
            _ => (),
        }
    }

    debug!("{:?}", event);

    // Don't emit event shorter than min_sz
    match event {
        Some((start, stop)) if stop - start > min_sz => Some((stop - start, start)),
        _ => None,
    }
}

pub enum Zyg {
    Hap1,
    Hap2,
    Hom,
}


pub fn make_event(chrom: &str,
                  ev: (usize, usize),
                  qual: usize,
                  zyg: Zyg,
                  name: Option<String>,
                  extra_info: Option<String>)
                  -> BedpeRow {
    let filters = ".";
    let name = name.unwrap_or(".".to_string());

    let haps = match zyg {
        Zyg::Hap1 => "ZS=HET;HAPS=0,0",
        Zyg::Hap2 => "ZS=HET;HAP=1,1",
        Zyg::Hom => "ZS=HOM;HAP=.,.",
    };
    let mut info = format!("{};TYPE=DEL;SOURCE=LOCAL_ASM", haps);
    if let Some(s) = extra_info {
        info.push(';');
        info.push_str(&s)
    };

    let (ev_len, ev_start) = ev;
    let ev_end = ev_start + ev_len;

    BedpeRow {
        chrom1: chrom.to_string(),
        start1: ev_start,
        end1: ev_start + 1,
        chrom2: chrom.to_string(),
        start2: ev_end,
        end2: ev_end + 1,
        name: name.to_string(),
        qual: qual as u8,
        strand1: '.',
        strand2: '.',
        filters: filters.to_string(),
        info: info,
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct BedpeRow {
    pub chrom1: String,
    pub start1: usize,
    pub end1: usize,
    pub chrom2: String,
    pub start2: usize,
    pub end2: usize,
    pub name: String,
    qual: u8,
    strand1: char,
    strand2: char,
    filters: String,
    info: String,
}

impl BedpeRow {
    pub fn simple(chrom: &str, range: &Range<u32>) -> BedpeRow {
        BedpeRow {
            chrom1: chrom.to_string(),
            start1: range.start as usize,
            end1: (range.start + 1) as usize,
            chrom2: chrom.to_string(),
            start2: range.end as usize,
            end2: (range.end + 1) as usize,
            name: ".".to_string(),
            qual: 100,
            strand1: '.',
            strand2: '.',
            filters: ".".to_string(),
            info: ".".to_string(),
        }
    }

    pub fn size(&self) -> usize {
        self.start2 - self.start1
    }
}


pub fn event_to_bed_pe(chr: &str,
                       del: (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>))
                       -> Vec<BedpeRow> {
    match del {
        (None, None) => vec![],
        (Some((e1, s1)), Some((e2, s2))) if e1 == e2 => {
            vec![make_event(chr, e1, (s1 + s2), Zyg::Hom, None, None)]
        }
        (Some((e1, s1)), Some((e2, s2))) => {
            vec![make_event(chr, e1, s1, Zyg::Hap1, None, None),
                 make_event(chr, e2, s2, Zyg::Hap2, None, None)]
        }
        (Some((e, s)), None) => vec![make_event(chr, e, s, Zyg::Hap1, None, None)],
        (None, Some((e, s))) => vec![make_event(chr, e, s, Zyg::Hap2, None, None)],
    }
}
