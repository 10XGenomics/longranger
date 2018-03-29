// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::str::FromStr;
use regex::Regex;
use std::fmt;
use serde_json;
use serde_json::Value;
use std::cmp::{max, min};
use std::fs::File;
use std::path::Path;
use rust_htslib::bam::{Read, IndexedReader};

#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

impl Locus {
    pub fn from_string(s: &str) -> Locus {
        FromStr::from_str(s).ok().unwrap()
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}



#[derive(Debug)]
pub struct LocusParseError;

fn remove_commas(s: &str) -> String {
    let ss = s.to_string();
    ss.replace(",", "")
}


impl FromStr for Locus {
    type Err = LocusParseError;

    fn from_str(s: &str) -> Result<Locus, LocusParseError> {
        let re = Regex::new(r"^(.*):([0-9,]+)(-|..)([0-9,]+)$").unwrap();
        let cap = re.captures(s);

        if cap.is_none() {
            return Result::Err(LocusParseError {});
        }

        let cap = cap.unwrap();

        let start_s = remove_commas(cap.get(2).unwrap().as_str());
        let end_s = remove_commas(cap.get(4).unwrap().as_str());

        Ok(Locus {
               chrom: cap.get(1).unwrap().as_str().to_string(),
               start: FromStr::from_str(&start_s).unwrap(),
               end: FromStr::from_str(&end_s).unwrap(),
           })
    }
}


pub struct CoverageStat {
    mean: f64,
    stdev: f64,
}

impl CoverageStat {
    pub fn from_coverage_json<P: AsRef<Path>>(path: P) -> Self {

        // Consider only the top percentile distribution while computing the stats
        let percentile = 0.999;

        // Open the file in read-only mode.
        let file = File::open(path).expect("Could not open file");

        // Read the JSON contents of the file as an instance of `User`.
        let data: Value = serde_json::from_reader(file)
            .expect("Error parsing the json file");

        // TODO: The layer could be abstracted and passed as an optional argument
        let map = data["confident_summary_depth_info_deduped"]
            .as_object()
            .unwrap();

        let mut coverage_dist = Vec::new();
        let mut sum_counts = 0u64;
        for (key, value) in map.iter() {
            let coverage = key.parse::<usize>().unwrap();
            let count = value.as_u64().unwrap();
            coverage_dist.push((coverage, count));
            sum_counts += count;
        }
        coverage_dist.sort();

        let cumulative_cutoff = (percentile * (sum_counts as f64)) as u64;
        // Find the top coverage (according to the percentile)
        let mut top_coverage = 0;
        let mut cumulative_counts = 0u64;
        for &(coverage, count) in &coverage_dist {
            if cumulative_counts > cumulative_cutoff {
                top_coverage = coverage;
                break;
            }
            cumulative_counts += count;
        }

        // Bin all the coverage count > top_coverage into top_coverage
        let mut new_coverage_dist = vec![0u64; top_coverage + 1];
        for &(coverage, count) in coverage_dist.iter() {
            new_coverage_dist[min(coverage, top_coverage)] += count;
        }

        // Find the mean and stdev
        let reductions = new_coverage_dist
            .iter()
            .enumerate()
            .fold((0u64, 0u64), |sum, (cov, &count)| {
                (sum.0 + (cov as u64) * count, sum.1 + (cov as u64) * (cov as u64) * count)
            });

        let mean = (reductions.0 as f64) / (sum_counts as f64);
        let stdev = ((reductions.1 as f64) / (sum_counts as f64) - mean * mean).sqrt();

        CoverageStat { mean, stdev }
    }
}


impl Locus {
    /// Check if the mean coverage within any window of size `window_size` within the locus is
    /// more that 10 standard deviations away from the mean coverage. The coverage statistics
    /// are computed from the input `coverage_json` file.
    pub fn has_excessive_coverage(&self, bam: &mut IndexedReader,
                                  coverage_stat: &CoverageStat)
                                  -> bool {

        let window_size = 200usize;
        let cov_cutoff = coverage_stat.mean + 10.0f64 * coverage_stat.stdev;

        let tid = bam.header.tid(self.chrom.as_bytes()).unwrap();
        // We wanto to fetch records from the region which is window_size
        // before the start and after the end
        let start = self.start.saturating_sub(window_size as u32);
        let end = self.end + window_size as u32;

        // A circular array to keep track of coverage within a window
        let mut window_coverage = vec![0; window_size];
        let mut total_coverage = 0;
        let mut max_window_coverage = 0;

        bam.fetch(tid, start, end)
            .expect("Error fetching BAM file.");
        // pileup over all covered sites
        for p in bam.pileup() {
            let pileup = p.unwrap();
            let idx = (pileup.pos() as usize - start as usize) % window_size;
            let mut current_coverage = 0;
            // Don't use pileup.depth() as it does not skip poor mappings.
            // Instead iterate through the records and count coverage.
            for alignment in pileup.alignments() {
                let record = alignment.record();
                // Skip poor mappings & secondary alignments
                if record.mapq() < 10 || record.is_secondary() {
                    continue;
                }
                current_coverage += 1;
            }
            total_coverage += current_coverage;
            total_coverage -= window_coverage[idx];
            window_coverage[idx] = current_coverage;
            max_window_coverage = max(max_window_coverage, total_coverage);
        }

        let max_avg_coverage = f64::from(max_window_coverage) / (window_size as f64);
        println!(" Max avg coverage = {}", max_avg_coverage);
        
        max_avg_coverage > cov_cutoff

    }
}
