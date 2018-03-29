// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#[derive(Serialize, Deserialize, Debug)]
pub enum MutType {
    MUTATION,
    DELETION,
    INSERTION,
    NA,
}


#[derive(Serialize, Deserialize, Debug)]
pub struct CorrectionEntry {
    pub code: String,
    pub pos: Vec<usize>,
    pub mtype: MutType,
}

impl CorrectionEntry {
    pub fn new(c: String, p:Vec<usize>, t:MutType) -> CorrectionEntry {
        CorrectionEntry{code:c, pos:p, mtype:t}
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Correction {
    pub min_dist: i32,
    pub info: Vec<CorrectionEntry>,
}

impl Correction {
    pub fn new() -> Correction{
        Correction{min_dist:9999, info: Vec::new()}
    }

    pub fn add(& mut self, ce:CorrectionEntry, d:i32) {
        if d == self.min_dist {
            self.info.push(ce);
        } else if d < self.min_dist {
            self.min_dist = d;
            self.info.clear();
            self.info.push(ce);
        }
    }
}
