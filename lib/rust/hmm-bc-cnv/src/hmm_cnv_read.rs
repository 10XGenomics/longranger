// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use hmm::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use math::*;
use support::*;
use bed::*;

pub struct CNVChromHMMRead {
    pub readcount:      Vec<i32>,
    pub is_good_bin:    Vec<bool>,
    pub copy_number:    Vec<u16>,
    pub num_bin:        usize,
    pub num_states:     usize,
    pub ploidy:         f64,
    pub read_per_hap:   f64,
    pub log_prob_trans: f64,
    pub log_prob_remain:f64,
}

pub struct CNVRead {
    pub chrom_cnvs:     Vec<CNVChromHMMRead>,
    pub copy_number:    Vec<Vec<u16>>,
    pub num_states:     usize,
    pub num_chrom:      usize,
    pub max_num_bin:    usize,
    pub num_bins:       Vec<usize>,
    pub bin_size:       i32,
    pub total_count:    f64,
    pub total_bin_primary:i32,
    pub chrom_names:    Vec<String>,
    pub ploidy:         f64,
}

impl HmmModel for CNVChromHMMRead {
    #![allow(unused_variables)]
    fn num_positions(&self) -> usize { self.num_bin }

    fn num_states(&self) -> usize {self.num_states}

    fn prior(&self, state:usize) -> f64 {0.0f64}

    fn score_state(&self, pos: usize, state: usize) -> f64 {
        if !self.is_good_bin[pos] {return 0.0f64;}
        let expected: f64 = (self.read_per_hap * state as f64).max(0.01f64);
        loglikelihood(expected, self.readcount[pos])
    }

    fn score_trans(&self, pos:usize, prev_state: usize, state: usize) -> f64 {
        if prev_state == state { self.log_prob_remain } else {self.log_prob_trans}
    }
}

impl CNVRead {
    pub fn new(bam_file: String, bin_size: i32, num_states: usize, prob_status_change: f64,
        initial_ploidy: f64, bad_regions_file: String, primary_contigs: Option<String>)
        -> CNVRead {
        println!("Getting genome inforamtion");
        let gi = GenomeInfo::new(bam_file.clone(), bin_size, primary_contigs);
        println!("Getting reads inforamtion");
        let rc = ReadcountInfo::new(bam_file, &gi, bin_size);
        println!("Getting bad regions inforamtion");
        let bed_info = BadRegion::new(bad_regions_file, gi.chroms_to_len.clone(), bin_size);
        let mut is_good_bin: Vec<Vec<bool>> = Vec::new();
        for _ in 0..gi.num_chrom {is_good_bin.push(vec![true;0]);}
        for (chr, status) in bed_info.get_all_bin_status() {
            if !gi.chroms_to_idx.contains_key(&chr) {continue;}
            is_good_bin[gi.chroms_to_idx[&chr]]=status;
        }
        // for i in 0..gi.num_chrom {
        //     println!("\n\nchromosome {}", i);
        //     for j in 0..gi.num_bins[i] {
        //         if j%50 ==0 {print!("\n{}:\t",i*100)};
        //         print!{"{} ",is_good_bin[i][j]};
        //     }
        // }

        let mut chrom_cnvs: Vec<CNVChromHMMRead> = Vec::with_capacity(gi.num_chrom);
        for i in 0..gi.num_chrom {
            chrom_cnvs.push(CNVChromHMMRead{
                readcount:      rc.readcount[i].clone(),
                copy_number:    vec![0u16;gi.num_bins[i]],
                is_good_bin:    is_good_bin[i].clone(),
                num_bin:        gi.num_bins[i],
                num_states:     num_states,
                ploidy:         initial_ploidy,
                read_per_hap:   rc.total_count / gi.total_bin_primary as f64/ initial_ploidy,
                log_prob_trans: prob_status_change.ln(),
                log_prob_remain:(1.0f64-(num_states-1) as f64 * prob_status_change).ln(),
            });
        }

        println!("read per hap {}", rc.total_count / gi.total_bin_primary as f64/ initial_ploidy);

        CNVRead{
            chrom_cnvs:     chrom_cnvs,
            copy_number:    vec![vec![0u16; gi.max_num_bin]; gi.num_chrom],
            num_states:     num_states,
            num_chrom:      gi.num_chrom,
            max_num_bin:    gi.max_num_bin,
            num_bins:       gi.num_bins,
            bin_size:       bin_size,
            total_count:    rc.total_count,
            total_bin_primary:gi.total_bin_primary,
            chrom_names:    gi.chroms,
            ploidy:         initial_ploidy,
        }
    }

    pub fn hmm_inference(&mut self) -> f64 {
        let mut sum_ploidy: i32 = 0;
        let mut sum_num_bin: i32 = 0;
        let mut num_non_zero: i32 = 0;
        self.num_bins.iter().fold(0, |acc, &x| acc +x) as i32;
        for i in 0..self.num_chrom {
            println!("Processing chromosome {}", self.chrom_names[i]);
            self.copy_number[i] = self.chrom_cnvs[i].viterbi();
            for j in 0..self.chrom_cnvs[i].num_bin {
                if self.chrom_cnvs[i].is_good_bin[j] && (self.copy_number[i][j] < self.num_states as u16 -1) {
                    sum_num_bin += 1;
                    sum_ploidy += self.copy_number[i][j] as i32;
                    if self.copy_number[i][j] > 0 {num_non_zero += 1;}
                }
            }
        }
        println!("{} {} {}", sum_ploidy, num_non_zero, sum_num_bin);
        //println!("{} {}", sum_ploidy, sum_num_bin);
        //let num_non_zero_bin = self.copy_number[0].iter().fold(0, |acc, &x| if x > 0 {acc+1} else {acc});
        //println!{"# non zero bins {}", num_non_zero_bin};
        let ploidy: f64 = sum_ploidy as f64 / sum_num_bin as f64;
        println!("ploidy {}", ploidy);
        let read_per_hap = self.total_count / self.total_bin_primary as f64/ ploidy;
        // for i in 0..self.chrom_cnvs[0].readcount.len() {
        //     if i%20 ==0 {print!("\n{}:\t",i)};
        //     print!{"{}->{} ",self.chrom_cnvs[0].readcount[i], self.copy_number[0][i]};
        // }
        // println!("");
        for i in 0..self.num_chrom {
            self.chrom_cnvs[i].ploidy = ploidy;
            self.chrom_cnvs[i].read_per_hap = read_per_hap;
        }
        ploidy
    }

    pub fn multiple_run(&mut self, times: i32) {
        let mut last_ploidy = 2.0f64;
        for _ in 0..times {
            let ploidy = self.hmm_inference();
            if (ploidy - last_ploidy).abs() < 1e-6 {break;}
            last_ploidy = ploidy;
        }
    }

    pub fn output_cnvs(&self, outfile:String) {
        let mut cnv_out = BufWriter::new(File::create(outfile).expect("Error in creating output cnv file."));

        for i in 0..self.num_chrom {
            let mut bin_start = 0usize;
            let mut is_cnv = false;
            let mut copy_number = 0;
            for j in 0..self.num_bins[i] {
                if is_cnv {
                    if self.copy_number[i][j]!=self.copy_number[i][j-1] {
                        // close the previous cnv ending at j-1
                        write!(cnv_out, "{}\t{}\t{}\t{}\n", self.chrom_names[i],
                            bin_start*self.bin_size as usize, j*self.bin_size as usize, copy_number)
                            .unwrap();
                        is_cnv = false;
                    } else if j==self.num_bins[i]-1 {
                        // close the terminal cnv
                        write!(cnv_out, "{}\t{}\t{}\t{}\n", self.chrom_names[i],
                            bin_start*self.bin_size as usize, (j+1)*self.bin_size as usize, copy_number)
                            .unwrap();
                        is_cnv = false;
                    }
                }

                if !is_cnv && self.copy_number[i][j]!=2 {
                    is_cnv = true;
                    bin_start = j;
                    copy_number = self.copy_number[i][j];
                }
            }
        }
    }
}
