// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::BufRead;
use serde_json;
use std::io::prelude::Write;
use correction::{CorrectionEntry, Correction, MutType};

/// get_ec_table takes a list of codes (either the left 7-mer or the right 9-mer) to generate
/// error correction table, wwhich is a HashMap of n-mer, Correction
pub fn get_ec_table(code: &str, ectable: &str, allow_2_errors: bool, allow_indel: bool) {
    let bases = "NATGC".as_bytes().to_vec();

    let file = File::open(code).expect("no code file");
    let buf = BufReader::new(file);
    let codes: Vec<String> = buf.lines().map(|l| l.expect("Could not parse line")).collect();
    let mut ec_table: HashMap<String, Correction> = HashMap::new();

    for c in &codes {
        // barcode with 0 errors
        (*ec_table.entry(c.clone()).or_insert_with(|| Correction::new()))
            .add(CorrectionEntry::new(c.clone(),vec![],MutType::NA),0);

        let code_bytes = c.as_bytes().to_vec();
        for i in 0..code_bytes.len() {
            let mut new_mer = code_bytes.clone();
            for n in &bases{
                if code_bytes[i] == *n {continue;}
                new_mer[i] = *n;
                let new_code: String = String::from_utf8(new_mer.clone()).
                    expect("Fail to convert a mutated code into a candidate barcode.");
                // barcode with 1 errors
                (*ec_table.entry(new_code.clone()).or_insert_with(|| Correction::new())).
                    add(CorrectionEntry::new(c.clone(),vec![i], MutType::MUTATION),1);

            }
        }

        // barcode with 2 errors
        if allow_2_errors {
            for i in 0..code_bytes.len() {
                for j in (i+1)..code_bytes.len() {
                    let mut new_mer = code_bytes.clone();
                    for n in &bases{
                        for m in &bases{
                            if code_bytes[i] == *n || code_bytes[j] ==*m { continue; }
                            new_mer[i] = *n;
                            new_mer[j] = *m;
                            let new_code: String = String::from_utf8(new_mer.clone()).
                                expect("Fail to convert a mutated code into a candidate barcode.");
                            (*ec_table.entry(new_code.clone()).or_insert_with(|| Correction::new())).
                                add(CorrectionEntry::new(c.clone(),vec![i,j], MutType::MUTATION),2);
                        }

                    }


                }
            }
        }

        if allow_indel {
            // deletion of 1 base
            for i in 0..code_bytes.len()-1 {
                let mut new_mer = code_bytes.clone();
                new_mer.remove(i);
                new_mer.push(bases[0]);
                for n in &bases {
                        new_mer[&code_bytes.len()-1] = *n;
                        let new_code: String = String::from_utf8(new_mer.clone()).
                            expect("Fail to convert a mutated code into a candidate barcode.");
                        (*ec_table.entry(new_code.clone()).or_insert_with(|| Correction::new())).
                            add(CorrectionEntry::new(c.clone(),vec![i], MutType::DELETION),1);
                }
            }

            // single base addition
            for i in 0..code_bytes.len()-1 {
                for k in 1..bases.len() {
                    let mut new_mer = code_bytes.clone();
                    new_mer.insert(i, bases[k]);
                    new_mer.pop();
                    let new_code: String = String::from_utf8(new_mer).
                        expect("Fail to convert a mutated code into a candidate barcode.");
                    (*ec_table.entry(new_code.clone()).or_insert_with(|| Correction::new())).
                        add(CorrectionEntry::new(c.clone(),vec![i], MutType::INSERTION),1);
                }
            }
        }
    }

    let mut d_min = 9999;
    let mut d_max = 0;
    let mut ttl_d = 0f64;
    let mut ttl_pair = 0;
    for i in 0..(codes.len()-1) {
        let code_i = codes[i].as_bytes();
        for j in (i+1)..(codes.len()) {
            let code_j = codes[j].as_bytes();
            let mut d = 0i32;
            for k in 0..code_i.len() {
                if code_i[k] != code_j[k] {d+=1;}
            }
            if d < d_min {d_min = d;}
            if d > d_max {d_max = d;}
            ttl_pair += 1;
            ttl_d += d as f64;
        }
    }
    println!("min d {}\tmax d {}\tmean d {}", d_min, d_max, ttl_d/(ttl_pair as f64));

    let mut ttl_case = 0;
    let mut ttl_uniq_case = 0;
    let mut ttl_dup = 0f64;
    let mut max_dup = 0;
    for (_, v) in &ec_table {
        if v.min_dist == 0 {continue;}
        ttl_case += 1;
        if v.info.len()==1 {ttl_uniq_case+=1;}
        ttl_dup += v.info.len() as f64;
        if v.info.len() > max_dup {max_dup = v.info.len();}
    }

    println!("total {}\nuniq {}\nmean dup {}\nmax duplicates {}\ntotal count {}", ttl_case, ttl_uniq_case,
        ttl_dup/(ttl_case as f64), max_dup, ttl_dup);

    let mut f_out = File::create(ectable).
        expect("error in creating the output json error correction table.");
    f_out.write_all(
        serde_json::to_string(&ec_table)
            .expect("Error in serializing the ec table into the json string.")
            .as_bytes()
        )
        .expect("Error in writing out the json file.");
}

/// take in the barcode white list file and output the lists of
/// left 7-mers and right 9-mers
pub fn get_barcode_segment(bcwl: &str, left_file: &str, right_file: &str) {
    let file = File::open(bcwl).expect("no barcode whitelist file");
    let buf = BufReader::new(file);
    let whitelist: Vec<String> = buf.lines().map(|l| l.expect("Could not parse line")).collect();

    let mut left_7_ec: HashMap<String, i32> = HashMap::new();
    let mut right_9_ec: HashMap<String, i32> = HashMap::new();
    for bc in &whitelist{
        let left_7: String = bc.chars().take(7).collect();
        let right_9: String = bc.chars().skip(7).take(9).collect();
        *left_7_ec.entry(left_7).or_insert(0) += 1;
        *right_9_ec.entry(right_9).or_insert(0) += 1;
    }
    println!("Component sizes: left {}\tright {}\n", left_7_ec.len(),right_9_ec.len());

    let mut f_left_7 = BufWriter::new(File::create(left_file)
        .expect("Error in creating the left 7-mer file"));
    for c in left_7_ec.keys() {
        writeln!(f_left_7, "{}", c).expect("Error in writing left 7-mers");
    }

    let mut f_right_9 = BufWriter::new(File::create(right_file)
        .expect("Error in creating the right 9-mer file"));
    for c in right_9_ec.keys() {
        writeln!(f_right_9, "{}", c).expect("Error in writing right 9-mers");
    }
}
