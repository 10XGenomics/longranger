// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use rust_htslib::bam;
use std::fmt;
use std::fmt::Write;
use bc_correction::BarcodeValidator;
use std::cmp::max;
use std::collections::HashSet;
use std::hash::{Hash, SipHasher, Hasher};

#[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Deserialize, Serialize)]
pub struct Fragment {
    pub tid: i32,
    pub start: i32,
    pub end: i32,
    pub barcode_idx: usize,
    pub num_reads: i32, // Number of read pairs
}

impl fmt::Display for Fragment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Fragment:: {} {}-{} bc:{} num reads: {}", self.tid, self.start, self.end, self.barcode_idx, self.num_reads)
    }
}

impl Fragment {
    pub fn to_string(&self, barcode_validator: &BarcodeValidator) -> String {
        let mut output = String::new();
        write!(&mut output, "Fragment:: {} {}-{} bc:{} num reads: {}",
        self.tid, self.start, self.end, barcode_validator.whitelist[self.barcode_idx],
        self.num_reads).expect("Fail to write out this Fragment.");
        output
    }

    pub fn merge(&mut self, other: &Fragment) {
        if self.end >= other.start && self.start <= other.end {
            panic!("Two adjacent fragments belong to the same barcode overlap  {} {}", self, other);
        }
        if self.start < other.start {
            self.end = other.end;
        } else {
            self.start = other.start;
        }
        self.num_reads += other.num_reads;
    }
}

pub struct BarcodeState {
  pub last_pos : i32,
  pub last_pos_all_read: i32,
  pub first_pos : i32,
  pub tid: i32,
  pub barcode_idx : usize,
  pub frag_starts : Vec<i32>,
  pub frag_ends : Vec<i32>,
  pub frag_num_reads : Vec<i32>,
  pub qname_set: Option<HashSet<u32>>,
  pub low_mapq_qname_set: Option<HashSet<u32>>,
}

impl BarcodeState {

  pub fn new(bc_idx:usize) -> BarcodeState {
    BarcodeState{
        last_pos:0,
        last_pos_all_read: 0,
        first_pos:0,
        tid: -1,
        barcode_idx:bc_idx,
        frag_starts : Vec::new(),
        frag_ends : Vec::new(),
        frag_num_reads : Vec::new(),
        qname_set: None,
        low_mapq_qname_set: None,
    }
  }

  // Update the barcode state in with the new read.
  // If this read should actually create a new fragment, we should return the old fragment and set the state
  // to the new fragments
  pub fn add_read(&mut self, record: &bam::record::Record, rldist: i32, mapq: u8)
  {
      let mut s = SipHasher::new();
      record.qname().hash(&mut s);
      let hash_key: u32 = s.finish() as u32;

      if self.qname_set.is_some()  && (record.tid() != self.tid ||
        record.pos() > self.last_pos_all_read + rldist) {
          self.frag_starts.push(self.first_pos);
          self.frag_ends.push(self.last_pos);
          if let Some(ref q) = self.qname_set {
             self.frag_num_reads.push(q.len() as i32);
          }
          self.qname_set = None;
          self.low_mapq_qname_set = None;
      }

      // whether and how to add the read info
      let end_pos = record.cigar().end_pos().unwrap();

      if record.mapq() < mapq {
          if let Some(ref mut q) = self.low_mapq_qname_set{
              q.insert(hash_key);
              self.last_pos_all_read = end_pos;
          }
      } else {
          if self.qname_set.is_none() {
              if record.mapq() >= mapq {
                  self.first_pos = record.pos();
                  self.last_pos = end_pos;
                  self.last_pos_all_read = self.last_pos;
                  self.tid = record.tid();
                  self.qname_set = Some(HashSet::new());
                  self.low_mapq_qname_set = Some(HashSet::new());
              }
          }

          self.last_pos_all_read = end_pos;
          self.last_pos = self.last_pos_all_read;
          if let Some(ref mut q) = self.qname_set {
              if let Some(ref mut l) = self.low_mapq_qname_set {
                  q.insert(hash_key);
                  for k in l.drain() {q.insert(k);}
              }
          }
      }
  }


   pub fn finalize_molecule(&mut self)
   {
     if self.qname_set.is_some() {
         self.frag_starts.push(self.first_pos);
         self.frag_ends.push(self.last_pos);
         if let Some(ref q) = self.qname_set {
             self.frag_num_reads.push(q.len() as i32);
         }
         self.qname_set = None;
         self.low_mapq_qname_set = None;
     }
   }

   pub fn condidtional_finalize_molecule(&mut self, min_pos: i32) {
       if self.qname_set.is_some() && self.last_pos < min_pos {
            self.frag_starts.push(self.first_pos);
            self.frag_ends.push(self.last_pos);
            if let Some(ref q) = self.qname_set {
                 self.frag_num_reads.push(q.len() as i32);
            }
            self.qname_set = None;
            self.low_mapq_qname_set = None;
       }
   }

   pub fn get_fragments(&self) -> Vec<Fragment> {
       let mut fragments : Vec<Fragment> = Vec::new();
       for i in 0..self.frag_starts.len() {
           fragments.push(Fragment{barcode_idx : self.barcode_idx, tid:self.tid,
               start: self.frag_starts[i] as i32,
                end: self.frag_ends[i] as i32,
                num_reads: self.frag_num_reads[i]}
           );
       }
       fragments
   }


  pub fn distance_from_closest_fragment(&self, pos:(i32, i32)) -> (i32, Option<usize>) {
      let fragment_idx = match self.frag_ends.binary_search(&pos.1) {
          Ok(idx) => idx,
          Err(idx) => idx
      };
      const MAXDIST: i32 = 999999999;
      if fragment_idx == 0 || fragment_idx == self.frag_ends.len() {return (MAXDIST, None);}
      let dist_pre = pos.0 - self.frag_ends[fragment_idx-1];
      let dist_next = self.frag_starts[fragment_idx]-pos.1;
      if dist_pre < dist_next { (max(0, dist_pre), Some(fragment_idx-1)) }
      else {(max(0, dist_next), Some(fragment_idx))}
  }

  pub fn adjust_fragment(&mut self, pos:(i32,i32), rldist: i32) {
      let (dist, frag_idx) = self.distance_from_closest_fragment(pos);
      let frag_idx: usize = frag_idx.unwrap();
      if dist > 0 {
        let mut is_to_test = false;
        let mut test_idx1: usize = 0;
        let mut test_idx2: usize = 0;
        if pos.0 < self.frag_starts[frag_idx] &&  self.frag_starts[frag_idx] - pos.0 <= rldist {
            self.frag_starts[frag_idx] = pos.0;
            if frag_idx > 0 {
                is_to_test = true;
                test_idx1  = frag_idx-1;
                test_idx2 = frag_idx;
            }
            self.frag_num_reads[frag_idx] += 1;
        } else if pos.1 > self.frag_ends[frag_idx] &&  pos.0 - self.frag_ends[frag_idx] <= rldist {
            self.frag_ends[frag_idx] = pos.1;
            self.frag_num_reads[frag_idx] += 1;
            if frag_idx < self.frag_starts.len()-1 {
                is_to_test = true;
                test_idx1  = frag_idx;
                test_idx2 = frag_idx+1;
            }
        }

        if is_to_test {
            if self.frag_starts[test_idx2] - self.frag_ends[test_idx1] <= rldist {
                self.frag_ends[test_idx1] = self.frag_ends[test_idx2];
                self.frag_num_reads[test_idx1] += self.frag_num_reads[test_idx2];
                self.frag_starts.remove(test_idx2);
                self.frag_ends.remove(test_idx2);
                self.frag_num_reads.remove(test_idx2);
            }
        }
      }
  }
}
