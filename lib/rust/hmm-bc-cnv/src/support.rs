// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
// support functions
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use csv;
use fragment::*;
use probability::prelude::*;
use std::path::Path;
use serde::de::DeserializeOwned;
use std::iter::FromIterator;

#[derive(Debug)]
pub struct GenomeInfo {
    pub primary_contigs: HashSet<String>,
    pub chroms: Vec<String>,
    pub num_chrom: usize,
    pub chroms_to_len: HashMap<String, i32>,
    pub chroms_to_idx: HashMap<String, usize>,
    pub num_bins: Vec<usize>,
    pub max_num_bin: usize,
    pub total_bp_primary: f64,
    pub total_bin_primary: i32,
}

pub struct FragmentInfo {
    pub frag_enter: Vec<Vec<i32>>,
    pub frag_leave: Vec<Vec<i32>>,
    pub readcount: Vec<Vec<f32>>,
    pub dist_per_read:  f64,
    pub scale:          f64,
    pub mol_per_hap:    f64,
    pub prob_terminal:  f64,
    pub total_bp:       f64,
    pub total_count_asread: f64,
}


pub struct ReadcountInfo {
    pub readcount: Vec<Vec<i32>>,
    pub total_count: f64,
}


pub fn read_table<T: DeserializeOwned, P: AsRef<Path>>(filename: P, has_headers: bool, delim: u8) -> Vec<T> {

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(has_headers)
        .delimiter(delim)
        .from_path(filename)
        .expect("couldn't load table file");

    let mut tbl = Vec::new();
    for r in rdr.deserialize() {
        tbl.push(r.unwrap());
    }

    tbl
}


impl GenomeInfo {
    pub fn new(bam_file: String, bin_size: i32, primary_contigs_file: Option<String>) -> GenomeInfo {

        let bam_reader = bam::Reader::from_path(bam_file).ok().expect("Error opening BAM file");
        let chroms: Vec<String> = bam_reader.header().target_names().clone().into_iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap()).collect();
        let chrom_lens: Vec<usize> = (0..chroms.len())
            .map(|x| bam_reader.header().target_len(x as u32).unwrap() as usize).collect();


        let primary_contigs = match primary_contigs_file {
            Some(f) => read_table(f, false, b','),
            None => chroms.clone(),
        };
        let primary_contig_set = HashSet::from_iter(primary_contigs.into_iter());


        let mut chroms_to_len: HashMap<String, i32> = HashMap::new();
        let mut chroms_to_idx: HashMap<String, usize> = HashMap::new();
        let mut new_chroms: Vec<String> = Vec::new();
        let mut num_bins: Vec<usize> = Vec::new();
        let mut max_num_bin = 0;
        let mut num_chrom = 0;
        for i in 0..chroms.len() {
            if primary_contig_set.contains(&chroms[i]) {
                chroms_to_len.insert(chroms[i].clone(), chrom_lens[i] as i32);
                chroms_to_idx.insert(chroms[i].clone(), num_chrom as usize);
                new_chroms.push(chroms[i].clone());
                num_chrom += 1;
                let nbin = (chrom_lens[i]/(bin_size as usize)) + 1;
                num_bins.push(nbin);
                if nbin > max_num_bin {max_num_bin = nbin;}
            }
        }

        let total_bp_primary = chroms_to_len.values().fold(0.0f64, |acc, &x| acc+x as f64);
        let total_bin_primary: i32 = num_bins.iter().fold(0, |acc, &x| acc+x as i32);

        GenomeInfo {
            primary_contigs:    primary_contig_set,
            chroms:             new_chroms,
            num_chrom:          num_chrom,
            chroms_to_len:      chroms_to_len,
            chroms_to_idx:      chroms_to_idx,
            num_bins:           num_bins,
            max_num_bin:        max_num_bin,
            total_bp_primary:   total_bp_primary,
            total_bin_primary:  total_bin_primary,
        }
    }

    pub fn get_bin_arrays<T: Clone>(&self, value: T) -> Vec<Vec<T>> {
        let mut res = Vec::new();
        for sz in self.num_bins.iter() {
            let v = vec![value.clone(); *sz];
            res.push(v);
        }

        res
    }
}

impl FragmentInfo {
    pub fn new(frag_file: String, gi: &GenomeInfo, bin_size: i32, fragver: i32) -> FragmentInfo {

        let fragments_original = match fragver {
            0 => read_table(frag_file, true, b','),
            1 => read_table::<FragmentV1, _>(frag_file, true, b',').into_iter().map(|x| x.to_simple()).collect(),
            2 => read_table::<FragmentV2, _>(frag_file, true, b',').into_iter().map(|x| x.to_simple()).collect(),
            _ => panic!("Wrong fragments.h5 format choice. Valid versions are 0, 1 and 2."),
        };

        let fragments: Vec<Fragment> = fragments_original.into_iter().
            filter(|x| x.num_reads >=2 && gi.chroms_to_idx.contains_key(&x.chrom)).collect();

        let total_bp: f64 = fragments.iter().fold(0.0f64, |acc, x| acc + x.obs_len as f64);

        let total_reads: f64 = fragments.iter().fold(0.0f64, |acc, x| acc + x.num_reads as f64 - 1.0f64);
        let dist_per_read = total_bp / total_reads;
        let scale = total_bp / fragments.len() as f64;
        let mol_per_hap = total_bp / gi.total_bp_primary/ 2f64;
        let expon_distr = Exponential::new(1.0/scale);
        let prob_terminal = expon_distr.distribution(bin_size as f64);
        println!("{} {} {} {} {}", fragments.len(), dist_per_read, scale, mol_per_hap, prob_terminal);

        let mut f_enter: Vec<Vec<i32>> = gi.get_bin_arrays(0i32);
        let mut f_leave: Vec<Vec<i32>> = gi.get_bin_arrays(0i32);
        let mut readcount: Vec<Vec<f32>> = gi.get_bin_arrays(0.0f32);
        let mut total_count_asread = 0.0f64;

        for f in fragments {
            if !gi.chroms_to_idx.contains_key(&f.chrom) {continue;}
            let chrom_idx = gi.chroms_to_idx[&f.chrom] as usize;

            let enter_bin = (f.start_pos / bin_size) as usize;
            let leave_bin = (f.end_pos / bin_size) as usize;

            let enter_frac = ((enter_bin as i32 + 1) * bin_size - f.start_pos) as f32 / bin_size as f32;
            let leave_frac = (f.end_pos - leave_bin as i32 * bin_size) as f32 / bin_size as f32;

            if enter_bin == leave_bin {
                readcount[chrom_idx][enter_bin] += enter_frac - (1.0 - leave_frac);
            } else {

                for i in enter_bin+1..leave_bin { 
                    readcount[chrom_idx][i] += 1.0; 
                }
                readcount[chrom_idx][enter_bin] += enter_frac;
                readcount[chrom_idx][leave_bin] += leave_frac;
            }

            f_enter[chrom_idx][enter_bin] +=1;
            f_leave[chrom_idx][leave_bin] +=1;

            total_count_asread += (f.end_pos - f.start_pos) as f64 / (bin_size as f64);
        }

        FragmentInfo {
            frag_enter:     f_enter,
            frag_leave:     f_leave,
            readcount:      readcount,
            dist_per_read:  dist_per_read,
            scale:          scale,
            mol_per_hap:    mol_per_hap,
            prob_terminal:  prob_terminal,
            total_bp:       total_bp,
            total_count_asread: total_count_asread,
        }
    }
}

impl ReadcountInfo {
    pub fn new(bam_file: String, gi: &GenomeInfo, bin_size: i32) -> ReadcountInfo {
        let bam_reader = bam::Reader::from_path(bam_file).ok().expect("Error opening BAM file");
        let mut readcount: Vec<Vec<i32>> = vec![vec![0; gi.max_num_bin]; gi.num_chrom];
        let chroms: Vec<String> = bam_reader.header().target_names().clone().into_iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap()).collect();
        let mut total_count = 0.0f64;
        for _rec in bam_reader.records() {
            let rec = _rec.unwrap();
            if rec.tid() < 0 {continue;}
            let c = &chroms[rec.tid() as usize];
            if gi.primary_contigs.contains(c) {
                let p = rec.pos();
                if p >= 0 {
                    readcount[gi.chroms_to_idx[c]][rec.pos() as usize / bin_size as usize] +=1;
                }
                total_count += 1.0;
            }
        }

        ReadcountInfo {
            readcount:      readcount,
            total_count:    total_count,
        }
    }
}


// pub struct FragmentInfoAsread {
//     pub readcount:      Vec<Vec<i32>>,
//     pub total_count:    f64,
// }
//
// impl FragmentInfoAsread {
//     pub fn new(frag_file: String, gi: &GenomeInfo, bin_size: i32, fragver: i32) -> FragmentInfoAsread {
//
//         let fragments_original = match fragver {
//             0 => csv::Reader::from_file(frag_file)
//                 .expect("Fail to read in fragments.csv file")
//                 .decode().collect::<csv::Result<Vec<Fragment>>>().unwrap(),
//             1 => csv::Reader::from_file(frag_file)
//                 .expect("Fail to read in fragments.csv file")
//                 .decode().collect::<csv::Result<Vec<FragmentV1>>>().unwrap()
//                 .into_iter().map(|x| x.to_simple()).collect(),
//             2 => csv::Reader::from_file(frag_file)
//                 .expect("Fail to read in fragments.csv file")
//                 .decode().collect::<csv::Result<Vec<FragmentV2>>>().unwrap()
//                 .into_iter().map(|x| x.to_simple()).collect(),
//             _ => panic!("Wrong fragments.h5 format choice. Valid versions are 0, 1 and 2."),
//         };
//
//         let fragments: Vec<Fragment> = fragments_original.into_iter().
//             filter(|x| x.num_reads >=2 && gi.primary_contigs.contains_key(&x.chrom)).collect();
//         let mut total_count = 0.0f64;
//         let mut readcount: Vec<Vec<i32>> = vec![vec![0; gi.max_num_bin]; gi.num_chrom];
//
//         for f in fragments {
//             if !gi.chroms_to_idx.contains_key(&f.chrom) {continue;}
//             let chrom_idx = gi.chroms_to_idx[&f.chrom] as usize;
//             let enter_bin = (f.start_pos / bin_size) as usize;
//             let leave_bin = (f.end_pos / bin_size) as usize;
//             total_count += (leave_bin - enter_bin +1) as f64;
//             for i in enter_bin..leave_bin+1 { readcount[chrom_idx][i] +=1; }
//         }
//
//         FragmentInfoAsread {
//             readcount:      readcount,
//             total_count:    total_count,
//         }
//     }
// }
