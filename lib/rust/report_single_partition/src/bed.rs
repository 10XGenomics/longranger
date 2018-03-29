// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use csv;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead};


#[derive(Serialize, Deserialize)]
pub struct SimpleBedRecord {
    chrom: String,
    start: i32,
    end:   i32,
}

pub struct SimpleBed {
    regions: Vec<SimpleBedRecord>,
}

impl SimpleBed {
    pub fn new() -> SimpleBed {
        SimpleBed { regions: Vec::new()}
    }

    pub fn to_tsv(self, out_file: String) {
        let mut wrt = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(out_file)
            .expect("Fail to create the output bed file.");
            
        for f in &self.regions {
            let result = wrt.serialize(f);
            assert!(result.is_ok());
        }
    }

    pub fn from_tsv(in_file: String) -> SimpleBed {
        println!("In reading file {}", in_file);
        let mut rdr = csv::ReaderBuilder::new()
                        .delimiter(b'\t')
                        .has_headers(false)
                        .from_path(in_file)
                        .expect("Fail to open the input bed file.");

        SimpleBed { 
            regions: rdr.deserialize().map(|x| x.unwrap()).collect() 
        }
    }

    pub fn load_centromere_file(in_file: String) -> SimpleBed {
        let mut regions: Vec<SimpleBedRecord> = Vec::new();
        let file = File::open(in_file).expect("no such file");
        let buf = BufReader::new(file);
        for _l in buf.lines() {
            let l = _l.expect("Error in reading centromere file.");
            if l.starts_with("#") {continue;}
            let tokens: Vec<&str> = l.split('\t').collect();
            if tokens.len() < 4 {continue;}
            regions.push(SimpleBedRecord{chrom: String::from(tokens[1]), start: tokens[2].parse().unwrap(),
                end: tokens[3].parse().unwrap()});
        }
        SimpleBed{regions: regions}
    }
}


pub type Chrom = String;
pub struct Pos {
    starts: Vec<i32>,
    ends:   Vec<i32>,
}

impl Pos {
    pub fn new() -> Pos {
        Pos{ starts: Vec::new(), ends: Vec::new()}
    }

    pub fn add(&mut self, s: i32, e: i32) {
        self.starts.push(s);
        self.ends.push(e);
    }
}

pub struct Bed {
    regions:  HashMap<Chrom, Pos>,
    pub exists:   bool,
}

impl Bed {
    pub fn new() -> Bed {
        Bed {regions: HashMap::new(), exists: false}
    }

    pub fn add_region(&mut self, c: Chrom, s: i32, e: i32) {
        let v: &mut Pos = self.regions.entry(c).or_insert(Pos::new());
        v.add(s, e);
        self.exists = true;
    }

    pub fn from_simple_bed(sb: &SimpleBed) -> Bed {
        let mut bed = Bed::new();
        for sbr in &sb.regions { bed.add_region(sbr.chrom.clone(), sbr.start, sbr.end);}
        bed
    }

    pub fn to_simple_bed(&self) -> SimpleBed {
        let mut sb = SimpleBed::new();
        for (c, ref pos) in self.regions.iter() {
            for i in 0..pos.starts.len() {
                sb.regions.push(SimpleBedRecord {chrom: c.clone(), start: pos.starts[i], end: pos.ends[i]});
            }
        }
        sb
    }

    pub fn is_overlap(&self, c: &Chrom, s: i32, e: i32) -> bool {
        if let Some(ref pos) = self.regions.get(c) {
            let idx = match pos.starts.binary_search(&s) {
                Ok(i) => i,
                Err(i) => i,
            };
            return (idx < pos.starts.len() && e > pos.starts[idx]) || (idx > 0 && s < pos.ends[idx-1]);
        }
        false
    }
}
