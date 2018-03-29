// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std; 
use hmm::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use math::*;
use support::*;
use bed::*;
use itertools::Itertools;

pub struct CNVChromHMMAsread {
    pub is_good_bin:    Vec<bool>,
    pub readcount:      Vec<f32>,
    pub copy_number:    Vec<u16>,
    pub num_bin:        usize,
    pub num_states:     usize,
    pub read_per_hap:   f64,
    pub log_prob_trans: f64,
    pub log_prob_remain:f64,
    pub min_log_prob:   f64,
}

pub struct CNVAsread {
    pub chrom_cnvs:     Vec<CNVChromHMMAsread>,
    pub copy_number:    Vec<Vec<u16>>,
    pub num_states:     usize,
    pub num_chrom:      usize,
    pub max_num_bin:    usize,
    pub num_bins:       Vec<usize>,
    pub bin_size:       i32,
    pub total_count:    f64,
    pub total_good_bins:    f64,
    pub chrom_names:    Vec<String>,
    pub ploidy:         f64,
}


impl HmmModel for CNVChromHMMAsread {
    #![allow(unused_variables)]
    fn num_positions(&self) -> usize { self.num_bin }

    fn num_states(&self) -> usize {self.num_states}

    fn prior(&self, state:usize) -> f64 {0.0f64}

    fn score_state(&self, pos: usize, state: usize) -> f64 {
        let expected: f64 = self.read_per_hap * (state as f64).max(0.05f64);
        //loglikelihood(expected, self.readcount[pos].round() as i32).max(self.min_log_prob)
        negative_binomial_loglikelihood(expected, 1.1, self.readcount[pos].round() as i32).max(self.min_log_prob)
    }

    fn score_trans(&self, pos:usize, prev_state: usize, state: usize) -> f64 {
        if prev_state == state { self.log_prob_remain } else {self.log_prob_trans}
    }
}

impl CNVAsread {
    pub fn new(frag_file: String, bam_file: String, bin_size: i32, num_states: usize,
        fragver: i32, prob_status_change: f64, initial_ploidy: f64, bad_regions_file: String, 
        primary_contigs: Option<String>, min_prob: f64) -> CNVAsread {

        let gi = GenomeInfo::new(bam_file, bin_size, primary_contigs);
        let fi = FragmentInfo::new(frag_file, &gi, bin_size, fragver);
        let bed_info = BadRegion::new(bad_regions_file, gi.chroms_to_len.clone(), bin_size);

        let mut is_good_bin: Vec<Vec<bool>> = Vec::new();
        let mut total_count_asread = 0.0f64;
        let mut total_good_bins = 0.0f64;
        for _ in 0..gi.num_chrom {is_good_bin.push(vec![true;0]);}

        for (chr, status) in bed_info.get_all_bin_status() {
            if !gi.chroms_to_idx.contains_key(&chr) {continue;}
            let idx = gi.chroms_to_idx[&chr];
            is_good_bin[idx]=status;
            for i in 0..gi.num_bins[idx] {
                if is_good_bin[idx][i] {
                    total_count_asread += fi.readcount[idx][i] as f64;
                    total_good_bins += 1.0;
                }
            }
        }

        let read_per_hap =  total_count_asread / total_good_bins / initial_ploidy;
        println!("mol_per_hap {} {} {} {}", read_per_hap, total_count_asread, total_good_bins, initial_ploidy);

        let mut chrom_cnvs: Vec<CNVChromHMMAsread> = Vec::with_capacity(gi.num_chrom);
        for i in 0..gi.num_chrom {
            chrom_cnvs.push(CNVChromHMMAsread{
                copy_number:    vec![0u16;gi.num_bins[i]],
                is_good_bin:    is_good_bin[i].clone(),
                num_bin:        gi.num_bins[i],
                readcount:      fi.readcount[i].clone(),
                num_states:     num_states,
                read_per_hap:   read_per_hap,
                log_prob_trans:  prob_status_change.ln(),
                log_prob_remain: (1.0f64-(num_states-1) as f64 * prob_status_change).ln(),
                min_log_prob:   min_prob.ln(),
            });
        }

        CNVAsread{
            chrom_cnvs:     chrom_cnvs,
            copy_number:    vec![vec![0u16; gi.max_num_bin]; gi.num_chrom],
            num_states:     num_states,
            num_chrom:      gi.num_chrom,
            max_num_bin:    gi.max_num_bin,
            num_bins:       gi.num_bins.clone(),
            bin_size:       bin_size,
            total_count:    total_count_asread,
            total_good_bins:    total_good_bins,
            chrom_names:    gi.chroms.clone(),
            ploidy:         initial_ploidy,
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
                if self.chrom_cnvs[i].is_good_bin[j] && (self.copy_number[i][j] < self.num_states as u16 - 1) {
                    sum_num_bin += 1;
                    sum_ploidy += self.copy_number[i][j] as i32;
                    if self.copy_number[i][j] > 0 {num_non_zero += 1;}
                }
            }
        }
        println!("{} {} {}", sum_ploidy, num_non_zero, sum_num_bin);

        let ploidy: f64 = sum_ploidy as f64 / sum_num_bin as f64;
        println!("ploidy {}", ploidy);
        let read_per_hap = self.total_count / self.total_good_bins / ploidy;
        for i in 0..self.num_chrom {self.chrom_cnvs[i].read_per_hap = read_per_hap;}
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
            // get forward-backward inference
            let mut fb_inference = self.chrom_cnvs[i].get_inference(2);
            fb_inference.compute();

            let mut bin_start = 0usize;
            let mut is_cnv = false;
            for j in 0..self.num_bins[i] {
                if is_cnv {
                    if self.copy_number[i][j] != self.copy_number[i][j-1] || j == self.num_bins[i] - 1 {
                        let end_pos = 
                            if j == self.num_bins[i]-1 {
                                j + 1
                            } else {
                                j
                            };

                        // close the previous cnv ending at j-1
                        // get p-value
                        let (cp, _, p_val) = fb_inference.interval_significance(&self.chrom_cnvs[i], bin_start, end_pos);

                        let rc = &self.chrom_cnvs[i].readcount[bin_start..end_pos];
                        let rc_mean: f64 = rc.iter().map(|x| *x as f64).sum::<f64>() / (rc.len() as f64);
                        let rc_max = rc.iter().cloned().fold(std::f32::MIN, |m,v| if m > v { m } else { v });
                        let rc_min = rc.iter().cloned().fold(std::f32::MAX, |m,v| if m < v { m } else { v });

                        write!(cnv_out, "{}\t{}\t{}\t{}\t{:.4e}\t{:.3}\t{:.3}\t{:.3}\n", self.chrom_names[i],
                            bin_start*self.bin_size as usize, end_pos*self.bin_size as usize, cp, p_val, rc_mean, rc_max, rc_min)
                            .unwrap();
                        is_cnv = false;
                    }
                }

                if !is_cnv && self.copy_number[i][j]!=2 {
                    is_cnv = true;
                    bin_start = j;
                }

            }
        }
    }

    pub fn write_out_bc_cov(&self, bedgraph_file: String, zoomout_factor: usize, final_bin: usize) {
        let mut cov_out = BufWriter::new(File::create(bedgraph_file).expect("Error in creating output coverage file."));
        write!(cov_out, "#chrom\tstart\tend\tcov\n").unwrap();

        for (i, chrom_instance) in self.chrom_cnvs.iter().enumerate() {
            let ch = format!("chr{}", i+1);
            for (j, chk) in chrom_instance.readcount.clone().into_iter().take(chrom_instance.num_bin).chunks(zoomout_factor).into_iter().enumerate() {
                let chk: Vec<i32> = chk.map(|x| x.round() as i32).collect();
                let chk_sum: i32 = chk.iter().sum();
                let mean_cov: usize = chk_sum as usize / chk.len();
                write!(cov_out, "{}\t{}\t{}\t{}\n", ch, j*final_bin, (j+1)*final_bin, mean_cov).unwrap();
            }
        }
    }
}
