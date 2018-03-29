// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::pileup::Indel;
use locus::Locus;
use std::str;
//use std::str::FromStr;
use event::{BedpeRow, make_event, Zyg};
use csv;
use Args;



pub fn go(_args: &Args) {
    let args = _args.clone();

    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");
    let chroms: Vec<String> = bam.header
        .target_names()
        .iter()
        .map(|n| str::from_utf8(n).unwrap().to_string())
        .collect();
    let mut dels = Vec::new();

    for idx in 0..bam.header.target_count() as usize {
        let chrom = chroms[idx].clone();
        println!("Processing: {}", chrom);

        let chrom = chroms[idx].clone();
        let loc = Locus {
            chrom: chrom,
            start: 0,
            end: bam.header.target_len(idx as u32).unwrap(),
        };
        call_dels(&mut bam, &loc, &mut dels);
        println!("Got events: {}", dels.len());
    }

    //let locus = FromStr::from_str(&args.arg_locus.unwrap()).unwrap();
    //call_dels(&mut bam, &locus);

    // Write calls to BEDPE
    if let Some(f) = args.arg_out.clone() {
        let mut w = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(f)
            .unwrap();
        for c in dels {
            if c.size() > 35 {
                w.serialize(c).unwrap();
            }
        }
    }
}


pub fn call_dels(bam: &mut bam::IndexedReader, locus: &Locus, dels: &mut Vec<BedpeRow>) {

    let chr = &locus.chrom;
    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");

    for _pi in bam.pileup() {
        let pi = _pi.unwrap();
        let pos = pi.pos() as usize;
        let depth = pi.depth();
        //println!("pos: {}, depth: {}", pos, depth);

        //if depth != 2 { continue; }
        let depth_info = format!("DEPTH={}", depth);

        let mut al1 = Indel::None;
        let mut al2 = Indel::None;

        for (idx, al) in pi.alignments().enumerate() {
            let indel = al.indel();
            if idx == 0 {
                al1 = indel;
            } else if idx == 1 {
                al2 = indel;
            }
        }

        let name = format!("sn{}", dels.len());

        match (al1, al2) {
            (Indel::Del(e1), Indel::Del(e2)) if e1 == e2 => {
                dels.push(make_event(chr,
                                     (e1 as usize, pos),
                                     100,
                                     Zyg::Hom,
                                     Some(name),
                                     Some(depth_info)))
            }
            (Indel::Del(e1), Indel::Del(e2)) => {
                dels.push(make_event(chr,
                                     (e1 as usize, pos),
                                     100,
                                     Zyg::Hap1,
                                     Some(name.clone()),
                                     Some(depth_info.clone())));
                let name2 = format!("sn{}", dels.len());
                dels.push(make_event(chr,
                                     (e2 as usize, pos),
                                     100,
                                     Zyg::Hap2,
                                     Some(name2),
                                     Some(depth_info)));
            }
            (Indel::Del(e), Indel::None) => {
                let zyg = if depth == 2 { Zyg::Hap1 } else { Zyg::Hom };
                dels.push(make_event(chr,
                                     (e as usize, pos),
                                     100,
                                     zyg,
                                     Some(name),
                                     Some(depth_info)))
            }
            (Indel::None, Indel::Del(e)) => {
                dels.push(make_event(chr,
                                     (e as usize, pos),
                                     100,
                                     Zyg::Hap2,
                                     Some(name),
                                     Some(depth_info)))
            }
            _ => (),
        }
    }
}
