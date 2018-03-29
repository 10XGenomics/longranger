// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//     1. barcode error correction
//     2. lariat adversarial rescue
//     3. fragment generation and merging

use std::collections::VecDeque;
use std::cmp::{min, max};
use std::fs::File;
use std::collections::HashSet;
use std::io::{BufReader, BufRead};
use std::f64::consts::PI;

use itertools::Itertools;
use csv;

use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Cigar;

use pval::poisson_pval as ext_poission_pval;
use bc_correction::{BarcodeValidator, CorrectionRst};
use bed::*;
use fragment::{Fragment, BarcodeState};
use correction::*;

// using FxHashMap is critical to genrate reproducible results
// the regular HashMap will lead to slightly different fragmentation results and 
// slightly different mapping values.
use ::fx::*;
use ::thread_iterator::*;

const BASES: [u8;5] = [b'N',b'A',b'T',b'G',b'C'];

fn log_factorial(x: f64) -> f64 {
    x*x.ln() - x + 0.5 * (2.0*PI*x).ln()
}

fn fast_poisson_pval(freq: i32, lambda:f64) -> f64 {
    let max_res: f64 = 0.01;
    assert!(lambda > 0.0, "lambda should be positive");
    assert!(freq > 0, "frequency should be positive");

    let mut x: f64 = freq as f64;
    let mut ratio: f64 = x / lambda;
    assert!(ratio > 5.0, "fast pvalue calculation only applies to freq >> lambda");
    println!("\tfreq {} lambda {} ratio {}", freq, lambda, ratio);
    let mut cur: f64 = (-lambda + x*lambda.ln() - log_factorial(x)).exp();
    println!("factorial of  {} is {} ", x, (log_factorial(x)).exp());

    let mut sum = cur;
    loop {
        let res: f64 = 1.0 as f64 / (ratio - 1.0 as f64);
        if res < max_res { break; }
        x += 1.0;
        ratio = x / lambda;
        cur = cur * lambda / x;
        sum += cur;
    }

    println!("pvalue of Poisson({}, {}) is  {}", freq, lambda, sum);

    sum
}

fn poisson_pval(freq: i32, lambda: f64) -> f64 {
    if freq >= 50 {
        if freq as f64 / lambda > 5.0 {
            println!("\tfreq {} lambda {} ratio {}", freq, lambda, freq as f64 / lambda);
            fast_poisson_pval(freq, lambda)
        } else  {
            1.0f64
        }
    } else {
        ext_poission_pval(freq, lambda)
    }
}

pub struct BLF {  // Barcode correction, Lariat adversarial rescue, Fragment generation
    frag_out:                 csv::Writer<File>,
    centromere:               Bed,
    barcode_validator:        BarcodeValidator,
    barcode_states:           FxHashMap<usize, BarcodeState>,
    buff_size:                usize,
    pos_seperation:           i32,
    read_buff:                VecDeque<bam::Record>,
    variant_cnt:              FxHashMap<LocSeq, i32>,
    rldist:                   i32, // read link distant, default 60,0000
    mapq:                     u8, // threshold for a read to be used to delimit a fragment
    mapping_q_thr:            u8, // threshold for when mapping based error correction should apply
    using_mapping:             bool,
    num_reads:                i32,
    num_reads_correct:        i32,
    num_reads_uncorrectable:  i32,
    num_reads_low_prob:       i32,
    num_reads_no_rx_qx:       i32,
    num_reads_corrected:      i32,
    num_reads_newly_corrected:i32,
    num_reads_corrected_different: i32,
    cur_pos:                  i32,
    chrom:                    String,
    tid:                      i32,
}

impl BLF {
    pub fn new(frag_out_file: String,
        centromere_file: &str, barcode_whitelist_file: String, barcode_left_ec_table: String,
        barcode_right_ec_table: String, barcode_count_file: String, max_expected_barcode_errors: f64,
        bc_confidence_threshold: f64, max_num_err: i32, psudo_count: f64,
        rldist: i32, mapq: u8, mapping_q_thr: u8, read_buff_sz: usize, using_mapping: bool)
        -> BLF {
        let frag_wrt = csv::WriterBuilder::new().has_headers(false).from_path(frag_out_file).expect("Fail to create the fragment output file.");
        let bcv = BarcodeValidator::new(barcode_whitelist_file, barcode_left_ec_table,
            barcode_right_ec_table, barcode_count_file, max_expected_barcode_errors,
            bc_confidence_threshold, max_num_err, psudo_count);
        let centro_bed = match centromere_file {
            "" => Bed::new(),
            _ =>  Bed::from_simple_bed(&SimpleBed::from_tsv(String::from(centromere_file))),
        };
        //let barcode_states: Vec<BarcodeState>= (0..bcv.whitelist.len()).map(|x| BarcodeState::new(x)).collect();
        let barcode_states = FxHashMap();
        let read_buff: VecDeque<bam::Record> = VecDeque::with_capacity(3000);
        let variant_cnt: FxHashMap<LocSeq, i32> = FxHashMap();

        BLF{
            frag_out:                 frag_wrt,
            centromere:               centro_bed,
            barcode_validator:        bcv,
            barcode_states:           barcode_states,
            buff_size:                read_buff_sz,
            pos_seperation:           300,
            read_buff:                read_buff,
            variant_cnt:              variant_cnt,
            rldist:                   rldist,
            mapq:                     mapq,
            mapping_q_thr:            mapping_q_thr,
            using_mapping:            using_mapping,
            num_reads:                0,
            num_reads_correct:        0,
            num_reads_uncorrectable:  0,
            num_reads_low_prob:       0,
            num_reads_no_rx_qx:       0,
            num_reads_corrected:      0,
            num_reads_newly_corrected:0,
            num_reads_corrected_different: 0,
            cur_pos:                  0,
            chrom:                    String::new(),
            tid:                      -1,
        }
    }


    pub fn time_to_act(&self) -> bool {
        self.read_buff.len()> 0 && self.cur_pos - &self.read_buff[0].pos() > self.pos_seperation ||
            self.read_buff.len() >= self.buff_size
    }

    pub fn find_raw_info(&mut self, bam_file: String, outbam_file: String, contigs: String,
        target_file: &str, nofragmerg: bool, pval: f64, freq: i32) {
        let target_bed= match target_file {
            "" => Bed::new(),
            _ =>  Bed::from_simple_bed(&SimpleBed::from_tsv(String::from(target_file))),
        };
        let bam_reader = bam::Reader::from_path(&bam_file).ok().expect("Error opening BAM file");
        let mut bam_writer = bam::Writer::from_path(&outbam_file,
            &bam::header::Header::from_template(&bam_reader.header())).
            ok().expect("Error in creating the output BAM file");
        bam_writer.set_threads(2).ok();
        // load in primary contigs
        let file = File::open(contigs).expect("no such file");
        let buf = BufReader::new(file);
        let primary_contigs: HashSet<String> = buf.lines().map(|l| l.expect("Could not parse line")).collect();

        let chroms: Vec<String> = bam_reader.header().target_names().clone().into_iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap()).collect();
        let chrom_lens: Vec<usize> = (0..chroms.len())
            .map(|x| bam_reader.header().target_len(x as u32).unwrap() as usize).collect();

        let iter = ThreadProxyIterator::bam_file(bam_file.clone(), 64);

        for (tid, tid_iter) in iter.map(|r| r.unwrap()).group_by(|r| r.tid()).into_iter() {
            println!("tid={}",tid);

            //self.barcode_states = (0..self.barcode_validator.whitelist.len()).map(|x| BarcodeState::new(x)).collect();
            self.barcode_states = FxHashMap();
            self.read_buff = VecDeque::with_capacity(3000);
            self.variant_cnt = FxHashMap();
            self.tid = tid;
            self.chrom = if tid <0 {"unmapped".to_string()} else {chroms[tid as usize].clone()};
            let target_len = if tid < 0 {0usize} else {chrom_lens[tid as usize]};
            let to_skip: bool = tid == -1 || !primary_contigs.contains(&chroms[tid as usize]);

            println!("is to_skip fragment generation and LAR {} >>", to_skip);

            for rec in tid_iter {
                self.num_reads +=1;
                if !to_skip && !rec.is_duplicate() {
                    for v in Variant::get_variant(&rec, "AC").into_iter() {
                        *(self.variant_cnt.entry(v.get_loc_seq()).or_insert(0)) += 1;
                    }
                }
                self.cur_pos = rec.pos();
                self.read_buff.push_back(rec); // push and move rec into read_buff.
                if self.time_to_act() { self.handle_bc_fragment_las(&mut bam_writer, to_skip); }
            }

            self.cur_pos += 10000; // to clean up all remaining reads
            self.handle_bc_fragment_las(&mut bam_writer, to_skip);

            // fragment generations and outputing
            self.finalize_fragment();

            let fragments = if nofragmerg {
                self.abstract_raw_fragments()
            } else {
                unimplemented!("fragment merging is not implemented");
                //self.get_final_fragments(target_len, pval, freq)
            };

            println!("number of fragments for this chromosome: {}", fragments.len());

            if target_bed.exists {
                for f in fragments.into_iter()
                    .filter(|f| target_bed.is_overlap(&chroms[f.tid as usize], f.start, f.end)) {
                    let result = self.frag_out.serialize(f);
                    assert!(result.is_ok());
                }
            } else {
                for f in fragments {
                    let result = self.frag_out.serialize(f);
                    assert!(result.is_ok());
                }
            }

        }
    }

    pub fn handle_bc_fragment_las(&mut self, bam_writer: &mut bam::Writer, to_skip: bool) {
        while self.time_to_act() {
            let mut rec = self.read_buff.pop_front().expect("error in reading read_buff");
            // update_updatable
            // LAR
            let mut mapq: u8 = rec.mapq();
            let aux_value = bam::record::Aux::Integer(rec.mapq() as i32);
            rec.push_aux("OM".as_bytes(), &aux_value);
            if !to_skip && !rec.is_duplicate() && (!self.centromere.exists ||
                !self.centromere.is_overlap(&self.chrom, rec.pos(), rec.pos()+100)) {
                let mut v_best: FxHashMap<BaseRange, Variant> = FxHashMap();
                for v in Variant::get_variant(&rec, "AC").into_iter() { v_best.insert(v.get_base_range(), v); }
                for v in Variant::get_variant(&rec, "XC").into_iter() { v_best.remove(&v.base_range); }
                for (_, v) in v_best.into_iter() {
                    let count = self.variant_cnt[&v.loc_seq];
                    if count > 1 {
                        let penalty: f64 = match v.base_range.length { 1 => -2.0, _ => -3.0};
                        let modified_penalty: f64 = penalty/count as f64;
                        let as_tag: f64 = rec.aux(String::from("AS").as_bytes())
                            .expect("error in getting AS")
                            .integer() as f64;
                        let calculated_xs = as_tag - (mapq as f64 / 10.0);
                        let mapq_f = (as_tag - penalty + modified_penalty - calculated_xs) * 10.0;
                        let mapq0 = 60f64.min(mapq_f).round() as u8;
                        if mapq0 > mapq { mapq = mapq0;}
                    }
                }
            }
            rec.remove_aux("AC".as_bytes());
            rec.remove_aux("XC".as_bytes());
            rec.set_mapq(mapq);  

            // To fix -- need to handle multiple GEM groups
            // barcode correction
            //let bc = self.correct_barcode_with_ec_table(&rec);
            //if let CorrectionRst::NEWBC(_bc) = bc {
            //    rec.push_aux("BX".as_bytes(), &bam::record::Aux::String(_bc.as_bytes()))
            //};

            // fragment generation
            if let Some(bx_id) = self.get_bx_id(&rec) {
                    let bcs = self.barcode_states.entry(bx_id).or_insert(BarcodeState::new(bx_id));
                    bcs.add_read(&rec, self.rldist, self.mapq);
            }

            bam_writer.write(&rec).expect("Error in writing reads");
        }

        let mut key_to_del: Vec<LocSeq> = Vec::new();
        for (k, _) in self.variant_cnt.iter() {
            if self.cur_pos - k.pos > self.pos_seperation+100 {
                key_to_del.push(k.clone());
            }
        }
        for k in key_to_del.iter() {
            self.variant_cnt.remove(k);
        }
    }

    fn get_bx_id(&self, rec: &bam::Record) -> Option<usize> {
        if let Some(bx) = get_tag(&rec, "BX") {
            let parts = bx.split('-').collect::<Vec<&str>>();
            let barcode_idx = *self.barcode_validator.barcode_2_idx.get(parts[0]).expect("Corrected barcode is not in whitelist");

            let gg = parts[1].parse::<usize>().unwrap();
            let offset = (gg as usize - 1) * self.barcode_validator.whitelist.len();
            Some(barcode_idx + offset)

        } else {
            None
        }        
    }

    fn get_gem_group(&self, rec: &bam::Record) -> Option<u8> {
        if let Some(rg) = get_tag(&rec, "RG") {
            let rg_comps = rg.split(':').collect::<Vec<&str>>();
            if rg_comps.len() < 3 {
                return None
            }
            
            let gem_group = rg_comps[2].parse::<u8>();
            return match gem_group {
                Ok(g) => Some(g),
                _ => None,
            }
        }

        None
    }

    // FIXME! handle GEM groups
    pub fn correct_barcode_with_ec_table(&mut self, rec: &bam::Record)
        -> CorrectionRst {
        let (barcode, qual_bytes) = match (rec.aux("RX".as_bytes()), rec.aux("QX".as_bytes())) {
            (Some(_raw), Some(_qual)) =>
                 (String::from_utf8(_raw.string().to_vec()).expect("Fail to convert RX into String."),
                 _qual.string().to_vec()),
            _ => {
                self.num_reads_no_rx_qx +=1;
                return CorrectionRst::NA
            },
        };

        if self.barcode_validator.barcode_2_idx.get(&barcode).is_some()
        {
            self.num_reads_correct += 1;
            return CorrectionRst::CORRECT;
        }

        // get the left 7-mer
        let left_7: String = barcode.chars().take(7).collect();
        let left_7_choices = match self.barcode_validator.ec_table_left.get(&left_7) {
            Some(s) => { s },
            None => {
                self.num_reads_uncorrectable += 1;
                return CorrectionRst::UNCORRECTABLE;
            }
        };

        let mut likelihoodies: FxHashMap<String, f64> = FxHashMap();

        let mut min_dist = 9999;
        let range = (rec.pos(), rec.cigar().end_pos().unwrap_or(rec.pos()));
        let mapq = rec.mapq();

        for l in &left_7_choices.info {
            let mut qv0 = 1.0f64;
            for p in &l.pos {qv0 *= probability(qual_bytes[*p]);}

            let (r9s, pos_shift) = match l.mtype {
                MutType::MUTATION | MutType::NA  => {
                    let right_9: String = barcode.chars().skip(7).take(9).collect();
                    (vec![right_9], 0isize)
                },
                MutType::DELETION => {
                    let right_9: String = barcode.chars().skip(6).take(9).collect();
                    (vec![right_9], -1isize)
                },
                MutType::INSERTION => {
                    let right_9: String = barcode.chars().skip(8).take(8).collect();
                    let mut r9s: Vec<String> = Vec::new();
                    for i in 1..BASES.len() {
                        let mut r: String = right_9.clone();
                        r.push(BASES[i] as char);
                        r9s.push(r);
                    }
                    (r9s, 1isize)
                },
            };

            for right_9 in r9s {
                let right_9_choices = match self.barcode_validator.ec_table_right.get(&right_9) {
                    Some(s) => { s },
                    None => { continue;}
                };
                // skip if total number of errors exceeds max allowed
                if left_7_choices.min_dist + right_9_choices.min_dist > self.barcode_validator.max_num_err { continue;}
                for r in &right_9_choices.info {
                    let mut qv = qv0;
                    // keep track the total number of errors
                    // only keep the candidates with the smallest number of errors
                    let d = left_7_choices.min_dist + right_9_choices.min_dist;
                    if d < min_dist {
                        min_dist = d;
                        likelihoodies.clear();
                    }
                    for p in &r.pos {
                        let idx: usize = ((*p+7) as isize +pos_shift) as usize;
                        if idx >= barcode.len() {continue;}
                        qv *= probability(qual_bytes[idx]);
                    }
                    let candidate: String = l.code.clone() + &r.code.clone();
                    if let Some(idx) = self.barcode_validator.barcode_2_idx.get(&candidate) {
                        let likelihood_frag: f64 = if self.using_mapping && mapq > self.mapping_q_thr {
                            let bcs = self.barcode_states.entry(*idx).or_insert(BarcodeState::new(*idx));
                            frag_dist_2_likelihood(bcs.distance_from_closest_fragment(range).0)
                        } else {1.0};

                        let gg = self.get_gem_group(rec).unwrap();
                        let likelihood = self.barcode_validator.barcode_rates[&gg][*idx]*(0.0005 as f64).
                            max(qv)*likelihood_frag;
                        let mut v = likelihoodies.entry(candidate.clone()).or_insert(likelihood);
                        if likelihood > *v {*v = likelihood;}
                    }
                }
            }
        }

        if likelihoodies.len() == 0 {
            self.num_reads_uncorrectable += 1;
            return CorrectionRst::UNCORRECTABLE;
        }

        let mut max_prob = 0.0f64;
        let mut best_code: String = "".to_string();
        //let mut lvalues = Vec::new();
        let mut total_likelihood = 0.0;
        for (code, l) in &likelihoodies {
            total_likelihood += *l;
            if *l > max_prob {
                //lvalues.push(*l);
                max_prob = *l;
                best_code = code.clone();
            }
        }

        max_prob /= total_likelihood;
        if best_code != "" && max_prob > self.barcode_validator.bc_confidence_threshold {
            // NOTE: the following code will be needed when fragment generation proceeds barcode correction
            // adjust fragmentation when a barcode got corrected.
            // if self.using_mapping && mapq > self.mapping_q_thr {
            //     let idx = self.barcode_validator.barcode_2_idx.get(&best_code).unwrap();
            //     if rec.mtid() == rec.tid() && (rec.pos()-rec.mpos()).abs() < self.rldist {
            //         self.barcode_states[*idx].adjust_fragment(range, self.rldist);
            //     }
            // }

            self.num_reads_corrected += 1;
            match get_tag(rec, "BX") {
                Some(bx) => {
                    if bx != best_code.clone()+"-1" {
                        self.num_reads_corrected_different +=1;
                        CorrectionRst::DIFFBC(best_code+"-1")
                    } else {
                        CorrectionRst::CORRECT
                    }
                },
                None => {
                    self.num_reads_newly_corrected += 1;
                    CorrectionRst::NEWBC(best_code+"-1")
                },
            }
        } else {
            self.num_reads_low_prob += 1;
            CorrectionRst::LOWPROB
        }
    }

    pub fn finalize_fragment(&mut self) {
        for bs in &mut self.barcode_states.values_mut() {
            bs.finalize_molecule();
        }
    }


    pub fn conditional_finalize_fragment(&mut self, min_pos:i32) {
        for bs in &mut self.barcode_states.values_mut() {
            bs.condidtional_finalize_molecule(min_pos);
        }
    }

    pub fn abstract_raw_fragments(&self) -> Vec<Fragment> {
        let mut fragments : Vec<Fragment> = Vec::new();
        for bc_state in self.barcode_states.values() {
            for f in bc_state.get_fragments() { 
                fragments.push(f) 
            }
        }
        fragments.sort();
        fragments
    }


    /*
   pub fn get_final_fragments(&self, chr_length: usize, pval_thr: f64, freq_thr: i32) -> Vec<Fragment> {

       let expt_same_bc_pairs: f64 = self.barcode_validator.barcode_rates[1].iter().
           fold(0.0, |sum, &x| sum+(x as f64).powf(2.0));
        println!("expected number of pairs of fragements with the same barcode {}", expt_same_bc_pairs);

        let bin: usize = 10000;
        let max_length: usize = chr_length/bin + 1;
        let max_right: usize = 200000/bin + 1;

        println!("max_length = {}\tmax_right = {}", max_length, max_right);

        let mut collision_freq_2d: Vec<Vec<i32>> = vec![vec![0; max_right]; max_length];
        let mut start_endpoints: Vec<i32> = vec![0; max_length];
        let mut end_endpoints: Vec<i32> = vec![0; max_length];

        println!("after array initialization");


        for bc_state in &self.barcode_states {
            let fragments: Vec<Fragment> = bc_state.get_fragments();
            for f in &fragments {
                let left_bin = f.start as usize /bin;
                let right_bin = f.end as usize /bin;
                start_endpoints[left_bin] += 1;
                end_endpoints[right_bin] += 1;
            }


            let frag_len = fragments.len();
            for i in 0..frag_len {
                let right_bin_i = fragments[i].end as usize/ bin;
                for j in (i+1)..frag_len {
                    let left_bin_j = fragments[j].start as usize/ bin;
                    let gap_size = left_bin_j - right_bin_i ;
                    if gap_size < max_right {
                        collision_freq_2d[right_bin_i][gap_size] +=1;
                    }

                }
            }
        }

        println!("finishing recording fragment endpoint and gap informations");


        for i in 0..max_length {
            for j in 0..max_right{
                if collision_freq_2d[i][j] >= freq_thr {
                    let ttl_num_pair = end_endpoints[i]*start_endpoints[i+j];
                    let expt_cnt = (ttl_num_pair as f64) * expt_same_bc_pairs;
                    if expt_cnt < collision_freq_2d[i][j] as f64 {
                        println!("{}-{} {} > {}", i, j, collision_freq_2d[i][j], expt_cnt);
                        let pval = poisson_pval(collision_freq_2d[i][j],expt_cnt);
                        if pval <= pval_thr {
                            collision_freq_2d[i][j] = 1;
                        } else {
                            collision_freq_2d[i][j] = 0;
                        }
                    }
                }
            }
        }
        println!("updated the colision array");

        let mut final_fragments: Vec<Fragment> = Vec::new();
        for bc_state in &self.barcode_states {
            let fragments = bc_state.get_fragments();
            let frag_len = fragments.len();
            let mut is_to_merge: Vec<bool> = vec![false; frag_len];
            for i in 0..frag_len {
                let right_bin_i = fragments[i].end as usize/ bin;
                for j in (i+1)..frag_len {
                    let left_bin_j = fragments[j].start as usize/ bin;
                    let right_bin_j = fragments[j].end as usize/ bin;
                    let gap_size = left_bin_j - right_bin_i ;
                    if gap_size >= max_right {continue;}
                    if collision_freq_2d[right_bin_j][gap_size] == 1 {
                        for k in i..j { is_to_merge[k]=true;}
                        if j!=i+1 {println!("Merge non adjacent fragments {} {}", i, j);}
                    }
                }
            }

            let mut i = 0;
            while i < frag_len {
                let mut frag = fragments[i].clone();
                while is_to_merge[i] {
                    frag.merge(&fragments[i+1]);
                    i+=1;
                }
                final_fragments.push(frag);
                i+=1;
            }
        }

        final_fragments.sort();
        final_fragments
    }
    */

    pub fn print_stats(&self) {
         let num_reads_with_error = self.num_reads - self.num_reads_correct;

         println!{"\n"};
         println!{"Total number of reads                       {}", self.num_reads};
         println!{"Total number of reads already correct       {}", self.num_reads_correct};
         println!{"Total number of reads with errors           {}", self.num_reads - self.num_reads_correct};
         println!{"percent of reads with errors                {}", num_reads_with_error as f64 / self.num_reads as f64};
         println!{"Total number of reads uncorrectable         {}", self.num_reads_uncorrectable};
         println!{"percent reads uncorrectable                 {}", self.num_reads_uncorrectable as f64 / self.num_reads as f64};
         println!{"Total number of reads no barcode            {}", self.num_reads_no_rx_qx};
         println!{"Total number of reads with barcode          {}", self.num_reads - num_reads_with_error + self.num_reads_corrected};
         println!{"percent reads with barcode                  {}", 1.0 - (num_reads_with_error - self.num_reads_corrected) as f64 / self.num_reads as f64};
         println!{"Total number of reads low probability       {}", self.num_reads_low_prob};
         println!{"percent reads low probability               {}", self.num_reads_low_prob as f64 / self.num_reads as f64};
         println!{"Total number of reads corrected             {}", self.num_reads_corrected};
         println!{"percent reads corrected                     {}", self.num_reads_corrected as f64 / self.num_reads as f64};
         println!{"Total number of reads newly corrected       {}", self.num_reads_newly_corrected};
         println!{"percent reads newly corrected               {}", self.num_reads_newly_corrected as f64 / self.num_reads as f64};
         println!{"Total number of reads corrected differently {}", self.num_reads_corrected_different};
         println!{"\n"};
    }
}


pub fn probability(qual:u8) -> f64 {
    (10.0 as f64).powf(-(qual as i32 - 33) as f64/10.0)
}

pub fn frag_dist_2_likelihood(dist:i32) -> f64 {
    const MEAN: f64 = 1e4; // likelihood drops to 0.01 at 50 KB
    if dist <= 0 { 1.0f64 }
    else { (-dist as f64 / MEAN).exp()}
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub struct BaseRange {
    pos_read:  i32,
    length:    i32,
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub struct LocSeq {
    pos:        i32,
    seq:        Vec<u8>,
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub struct Variant {
    base_range: BaseRange,
    loc_seq:    LocSeq,
}


impl Variant {
    pub fn get_variant(rec: &bam::Record, tag: &str) -> Vec<Variant> {
        let mut rst: Vec<Variant> = Vec::new();
        if let Some(info) = get_tag(rec, tag) {
            for opt in info.split(';') {
                let parts: Vec<&str> = opt.split(',').collect();
                if parts.len() >= 3 {
                    let p = parts[0].parse().unwrap();
                    let p_r: i32 = parts[1].parse().unwrap();
                    let l: i32 = parts[2].parse().unwrap();


                    let bam_seq = rec.seq();
                    let mut seq = Vec::with_capacity(max(l,0i32) as usize);
                    
                    //let mut notes: String = "".to_string();
                    if rec.is_reverse() {

                        // If the read is hard-clipped, adjust where we get the bases form.
                        let mut last_cigar = rec.cigar()[0].clone();
                        for c in rec.cigar().into_iter() {
                            last_cigar = c.clone();
                        }

                        let read_offset = 
                            match last_cigar {
                                Cigar::HardClip(v) => v as usize,
                                _ => 0
                             };

                        // Where the variant is in the clipped read sequence
                        let p_clipped_r = p_r - (read_offset as i32);

                        if p_clipped_r < 0 {
                            continue;
                        }
                        
                        let read_len = bam_seq.len() as i32;
                        for i in (max(0, (read_len - p_clipped_r - max(l,0)) as usize)) .. (read_len - p_clipped_r) as usize {
                            seq.push(bam_seq[i]);
                        }

                        if p_clipped_r as usize >= bam_seq.len() {
                            println!("bad tags: {}", String::from_utf8_lossy(rec.qname()));
                        }

                    } else {
                        // If the read is hard-clipped, adjust where we get the bases form.
                        let read_offset = 
                            match rec.cigar()[0] {
                                Cigar::HardClip(v) => v as usize,
                                _ => 0
                             };

                        // Where the variant is in the clipped read sequence
                        let p_clipped_r = p_r - (read_offset as i32);
                        if p_clipped_r < 0 {
                            continue;
                        }

                        for i in (p_clipped_r as usize) .. min(bam_seq.len(), (p_clipped_r + max(l, 0)) as usize) {
                            seq.push(bam_seq[i]);
                        }

                        if p_clipped_r as usize >= bam_seq.len() {
                            println!("bad tags: {}", String::from_utf8_lossy(rec.qname()));
                        }

                    }

                    rst.push(
                        Variant{
                            base_range: BaseRange { pos_read: p_r, length: l },
                            loc_seq:    LocSeq { pos: p, seq: seq }
                        }
                    );
                }
            }
        }
        rst
    }

    #[inline]
    pub fn get_base_range(&self) -> BaseRange {self.base_range.clone()}

    #[inline]
    pub fn get_loc_seq(&self) -> LocSeq { self.loc_seq.clone() }
}


pub fn get_tag(r: &bam::record::Record, tag: &str) -> Option<String> {
    let _barcode = r.aux(tag.as_bytes());
    match _barcode {
           Some(bc) => {
               let a:String = String::from_utf8(bc.string().to_vec()).expect("Fail to convert bc bytes to string");
               Some(a.clone())
           },
           None => None
    }
}
