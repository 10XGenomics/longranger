// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use hmm::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use math::*;
use support::*;
use bed::*;

pub struct CNVChromHMM {
    pub frag_enter:     Vec<i32>,
    pub frag_leave:     Vec<i32>,
    pub readcount:      Vec<f32>,
    pub is_good_bin:    Vec<bool>,
    pub copy_number:    Vec<u16>,
    pub num_bin:        usize,
    pub num_states:     usize,
    pub scale:          f64,
    pub mol_per_hap:    f64,
    pub prob_terminal:  f64,
    pub prob_status_change: f64,
    pub min_log_prob:   f64,
}

pub struct CNV {
    pub chrom_cnvs:     Vec<CNVChromHMM>,
    pub copy_number:    Vec<Vec<u16>>,
    pub num_states:     usize,
    pub num_chrom:      usize,
    pub max_num_bin:    usize,
    pub num_bins:       Vec<usize>,
    pub bin_size:       i32,
    pub dist_per_read:  f64,
    pub scale:          f64,
    pub mol_per_hap:    f64,
    pub prob_terminal:  f64,
    pub total_count_asread: f64,
    pub total_good_bins:    f64,
    pub chrom_names:    Vec<String>,
}


impl HmmModel for CNVChromHMM {
    #![allow(unused_variables)]
    fn num_positions(&self) -> usize { self.num_bin }

    fn num_states(&self) -> usize {self.num_states}

    fn prior(&self, state:usize) -> f64 {0.0f64}

    fn score_state(&self, pos: usize, state: usize) -> f64 {0.0f64}

    fn score_trans(&self, pos:usize, prev_state: usize, state: usize) -> f64 {
        let mut exp_die: f64 = self.mol_per_hap* state as f64 * self.prob_terminal;
        let mut exp_create: f64 = exp_die;
        if state < prev_state {
            exp_die += self.mol_per_hap * (prev_state - state) as f64;
        } else if state > prev_state {
            exp_create += self.mol_per_hap * (state - prev_state) as f64;
        }

        let mut score = if self.is_good_bin[pos] {
            loglikelihood(exp_die, self.frag_leave[pos]) +
            loglikelihood(exp_create, self.frag_enter[pos])
        } else {0.0f64};
        score = score.max(self.min_log_prob);


        let cp_transition_prob = if state != prev_state {
            self.prob_status_change.ln()
        } else {
            (1.0f64 - self.prob_status_change * (self.num_states-1) as f64).ln()
        };

        // if self.is_good_bin[pos] && prev_state == 2 && state==2 {
        //     println!("{}:{}->{}\tending:{:.3} vs {}\tbeginning:{:.3} vs {}\t{}\t{}", pos, prev_state,
        //         state, exp_die, self.frag_leave[pos], exp_create, self.frag_enter[pos],
        //         cp_transition_prob, score);
        // }

        if self.prob_status_change < 0.0 { score } else {score + cp_transition_prob}

    }

}

impl CNV {
    pub fn new(frag_file: String, bam_file: String, bin_size: i32, num_states: usize,
        fragver: i32, prob_status_change: f64, initial_ploidy: f64, bad_regions_file: String,
        primary_contigs: Option<String>, min_prob: f64) -> CNV {

        let gi = GenomeInfo::new(bam_file, bin_size, primary_contigs);
        let fi = FragmentInfo::new(frag_file, &gi, bin_size, fragver);
        let bed_info = BadRegion::new(bad_regions_file, gi.chroms_to_len.clone(), bin_size);
        let mut is_good_bin: Vec<Vec<bool>> = Vec::new();
        for _ in 0..gi.num_chrom {is_good_bin.push(vec![true;0]);}
        let mut total_count_asread = 0.0f64;
        let mut total_good_bins = 0.0f64;
        for (chr, status) in bed_info.get_all_bin_status() {
            if !gi.chroms_to_idx.contains_key(&chr) {continue;}
            let idx = gi.chroms_to_idx[&chr];
            is_good_bin[idx]=status;
            for i in 0..gi.num_bins[idx] {
                //if is_good_bin[idx][i] {
                total_count_asread += fi.readcount[idx][i] as f64;
                total_good_bins += 1.0;
                //}
            }
        }


        let mol_per_hap =  total_count_asread / total_good_bins / initial_ploidy;
        println!("mol_per_hap {} {} {} {}", mol_per_hap, total_count_asread, total_good_bins, initial_ploidy);
        // for i in 0..gi.num_chrom {
        //     println!("\n\nchromosome {}", i);
        //     for j in 0..gi.num_bins[i] {
        //         if j%50 ==0 {print!("\n{}:\t",j*100)};
        //         print!{"{} ",if is_good_bin[i][j] {1} else {0}};
        //     }
        // }
        // println!("");

        let mut chrom_cnvs: Vec<CNVChromHMM> = Vec::with_capacity(gi.num_chrom);
        for i in 0..gi.num_chrom {
            chrom_cnvs.push(CNVChromHMM{
                frag_enter:     fi.frag_enter[i].clone(),
                frag_leave:     fi.frag_leave[i].clone(),
                readcount:      fi.readcount[i].clone(),
                copy_number:    vec![0u16;gi.num_bins[i]],
                is_good_bin:    is_good_bin[i].clone(),
                num_bin:        gi.num_bins[i],
                num_states:     num_states,
                scale:          fi.scale,
                mol_per_hap:    mol_per_hap,
                prob_terminal:  fi.prob_terminal,
                prob_status_change: prob_status_change,
                min_log_prob:   min_prob.ln(),
            });
        }

        CNV{
            chrom_cnvs:     chrom_cnvs,
            copy_number:    vec![vec![0u16; gi.max_num_bin]; gi.num_chrom],
            num_states:     num_states,
            num_chrom:      gi.num_chrom,
            max_num_bin:    gi.max_num_bin,
            num_bins:       gi.num_bins.clone(),
            bin_size:       bin_size,
            dist_per_read:  fi.dist_per_read,
            scale:          fi.scale,
            mol_per_hap:    mol_per_hap,
            prob_terminal:  fi.prob_terminal,
            total_count_asread: total_count_asread,
            total_good_bins:    total_good_bins,
            chrom_names:    gi.chroms.clone(),
        }
    }

    pub fn hmm_inference(&mut self) -> f64 {
        let mut sum_ploidy: i32 = 0;
        let mut sum_num_bin: i32 = 0;
        let mut num_non_zero: i32 = 0;
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

        let ploidy: f64 = sum_ploidy as f64 / sum_num_bin as f64;
        println!("ploidy {}", ploidy);
        self.mol_per_hap = self.total_count_asread / self.total_good_bins / ploidy;
        for i in 0..self.num_chrom { self.chrom_cnvs[i].mol_per_hap = self.mol_per_hap;}
        for i in 0..self.chrom_cnvs[0].frag_enter.len() {
            if i%20 ==0 {print!("\n{}:\t",i)};
            if !self.chrom_cnvs[0].is_good_bin[i] {continue;}
            print!{"{:3}:{:3}:{:3}->{}  ",self.chrom_cnvs[0].frag_enter[i], self.chrom_cnvs[0].frag_leave[i],
                self.chrom_cnvs[0].readcount[i], self.copy_number[0][i]};}
        println!("");
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

                //if i==0 {
                //    println!("{} {} {} {} {} {} {} {}", "chr1", j*self.bin_size as usize, (j+1)*self.bin_size as usize,
                //        self.copy_number[0][j], self.chrom_cnvs[0].frag_enter[j], self.chrom_cnvs[0].frag_leave[j],
                //        is_cnv, bin_start);
                //}
            }
        }
    }
}
