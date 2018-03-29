// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use csv;
use std::collections::HashMap;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct Bed {
    pub chrom:      String,
    pub start_pos:  i32,
    pub end_pos:    i32,
}

pub struct BadRegion {
    pub regions:    HashMap<String, Vec<Bed>>,
    pub bin_size:   f64,
    pub chr_sizes:  HashMap<String, i32>,
}

impl BadRegion {
    pub fn new(bad_regions_file: String, chr_sizes: HashMap<String, i32>, bin_size: i32) -> BadRegion {
        println!("{}", bad_regions_file);
        //let full_regions = csv::Reader::from_file(chr_size_file).expect("Fail to read in bed file")
        //    .has_headers(false).delimiter(b'\t')
        //    .decode().collect::<csv::Result<Vec<Bed>>>().unwrap();
        //let chr_sizes = full_regions.into_iter().map(|x| (x.chrom, x.end_pos))
        //    .collect::<HashMap<String,i32>>();


        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(bad_regions_file)
            .expect("Fail to read in bed file");


        let beds = rdr.deserialize().map(|x| x.unwrap()).collect::<Vec<Bed>>();
        let mut regions: HashMap<String, Vec<Bed>> = HashMap::new();
        for b in beds {
            if !chr_sizes.contains_key(&b.chrom) {continue;}
            (*regions.entry(b.chrom.clone()).or_insert_with(|| Vec::new())).push(b);
        }
        BadRegion {
            regions:    regions,
            bin_size:   bin_size as f64,
            chr_sizes:  chr_sizes,
        }
    }

    pub fn get_bin_status(&self, chrom: String) -> Option<Vec<bool>> {
        if !self.chr_sizes.contains_key(&chrom) {return None;}
        let num_bins = (self.chr_sizes[&chrom] as usize)/(self.bin_size as usize) + 1;
        let mut status = vec![true; num_bins];
        for b in self.regions[&chrom].clone() {
            let s = (b.start_pos as f64 / self.bin_size).floor() as usize;
            let e = (b.end_pos as f64 / self.bin_size).ceil() as usize;
            for i in s..e+1 {status[i]=false;}
        }
        Some(status)
    }

    pub fn get_all_bin_status(&self) -> HashMap<String, Vec<bool>> {
        let mut rst: HashMap<String, Vec<bool>> = HashMap::new();
        for (chr, size) in &self.chr_sizes {
            let num_bins = (*size as usize / self.bin_size as usize) + 1;
            let mut status = vec![true; num_bins];
            if self.regions.contains_key(chr) {
                for b in self.regions[chr].clone() {
                    let s = (b.start_pos as f64 / self.bin_size).floor() as usize;
                    let e = (b.end_pos as f64 / self.bin_size).ceil() as usize;
                    for i in s..e {status[i]=false;}
                }
            }
            rst.insert(chr.clone(), status);
        }
        rst
    }
}
