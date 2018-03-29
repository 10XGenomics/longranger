//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command to generate barcoded MSP substrings from FASTQ data
//!
//! For each input read, compute the MSP substrings. Each substring is at most 2K-P long, and
//! is stored in a fixed-size struct Bsp, along with the barcode id, read id, and read position,
//! and an Exts bitfield indicating the neighboring bases (if any).

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::collections::{HashMap};
use std::env;
use std::marker::PhantomData;
use rayon::prelude::*;

use utils::BcIndexer;
use shardio::{ShardWriteManager, ShardSender, ShardDef};
use multifastq::{MultiFastqIter, InputRead, Fastq};
use debruijn::msp::simple_scan;
use std::mem;

use debruijn::dna_string::DnaString;
use debruijn::{Kmer, Vmer, KmerIter, base_to_bits};
use debruijn::kmer::{self};
use debruijn::vmer::{Lmer};

use utils;

pub type Kmer1 = kmer::Kmer64;
pub type Pmer = kmer::IntKmer<u16>;
pub type Lmer3 = Lmer<Kmer1, [u64; 3]>;
const P: usize = 8;
const MEM_ALLOC_GB: usize = 8;

use martian::{Json, JsonDict, MartianStage, Error};
use serde::Serialize;
use serde_json::{self, Value};
use serde_json::map::Map;
//pub type JsonDict = Map<String, Value>;
//pub type Json = Value;

pub trait Smer<K: Kmer>: Vmer<K> + Serialize + Sync + Send {}

impl Smer<Kmer1> for Lmer3 where Lmer<Kmer1, [u64; 3]>: Vmer<Kmer1> {}

#[derive(Debug, Serialize, Deserialize)]
pub struct Bsp<K: Kmer, V: Smer<K>> {
    pub sequence: V,
    pub bc: u32,
    k: PhantomData<K>,
}

impl<K: Kmer, V: Smer<K>> Bsp<K, V> {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
 }

struct BspShard;

impl<K: Kmer,V: Smer<K>> ShardDef<Bsp<K,V>> for BspShard {
    fn get_shard(v: &Bsp<K,V>) -> usize {
        v.bc as usize
    }
}


/// Entry point for an MSP run. bc_wl is the path to the BC whitelist file, which will be used to
/// convert BC sequences to ints. fastq is set of pseudo-FASTQ files created by the _SORT_FASTQ_BY_BARCODE
/// pipeline. Permutation is the permutation of p-mers to use, serialized into the give path.
/// The best permutation to use is determined by msp_get_pmer_permutation below.
#[inline(never)]
pub fn main_msp(trim_min_qual: u8, bc_wl: &Path, fastq: Vec<PathBuf>, permutation: &Path, out_path: &Path) {

    let bc_idx = BcIndexer::new(bc_wl);
    let permutation : &Vec<usize> = &utils::read_obj(permutation).expect("couldn't read permutation");

    let mem = (MEM_ALLOC_GB - 1) * (10 as usize).pow(9);
    let nshards = 8192;
    let bsp_size = mem::size_of::<Bsp<Kmer1, Lmer3>>();
    let buf_per_shard = mem / nshards / bsp_size * 2;       // the *2 here is an awful fudge
    info!("MSP mem config: mem: {}, bsp_size: {}, Shards: {}, Buffered items per shard: {}", mem,
            bsp_size,  nshards, buf_per_shard);

    let shard = ShardWriteManager::new(out_path, buf_per_shard, nshards, 2);
    let sender = shard.get_sender();

    let nreads: u64 = 
        fastq.into_par_iter().map_with(sender, |sender, fastq_fn| {
            let fastq = MultiFastqIter::new(&fastq_fn, &bc_idx);
            process_fastq_set::<_, Kmer1, Lmer3>(trim_min_qual, Pmer::k(), &permutation, fastq, sender)
        }).sum();

    info!( "MSP read count: {}", nreads);
}


/// Load a sampling of reads and determine the distribution of pmers.  We will redefine the order of
/// pmers to make the most common by the lexicographically highest, which will gives more uniform
/// chunk sizes.
pub fn msp_get_pmer_permutation(fastq: Vec<PathBuf>, bc_wl: &Path, permutation_path: &Path) {
    let bc_idx = &BcIndexer::new(bc_wl);
    let fq_src_small = fastq.iter().flat_map(|p| MultiFastqIter::new(&p, &bc_idx).take(25000));
    let pmer_counts = count_pmers_fastq(fq_src_small);
    let permutation = pmer_inverse_freq_permutation(&pmer_counts);

    info!("Writing pmer permutation: {:?}", permutation_path);
    utils::write_obj(&permutation, permutation_path).expect("write failed");
}


/// Do the MSP splitting of sequence seq, and pack the MSP partitions into Bsp structs,
/// and send these struct out to disk.
#[inline(never)]
fn msp_read<K: Kmer + Send + Sync, V: Smer<K>>
           (p: usize,
            partition: u32,
            permutation: &Vec<usize>,
            _read: u16,
            seq: &[u8],
            shard_sender: &mut ShardSender<Bsp<K,V>, BspShard>) {

    // Can't do anything with strings shorter than K
    if seq.len() < K::k() { return; }

    let msp_parts = simple_scan(K::k(), p, seq, permutation);

    for (_, _, start, len) in msp_parts {

        assert!(len >= K::k());
        let sequence = V::from_slice(&seq[start .. start + len]);
        let b = Bsp { sequence, bc: partition, k: PhantomData };

        if (b.len() as usize) < K::k() || (b.len() as usize) != b.sequence.len()
        {
            println!("bad bsp: {:?}", b);
            panic!("bad bsp: {:?}", b);
        }

        shard_sender.send(b);
    }
}

/// Find the longest prefix of read such that all bases in the final kmer have qv >= min_qual
/// Should implement https://github.com/10XDev/supernova/blob/master/lib/assembly/src/paths/long/BuildReadQGraph48.cc#L70
fn find_trim_len(fq: &Fastq, min_qual: u8, k: usize) -> usize {

    let qvs = fq.qual.clone().into_bytes();
    let mut good = 0;

    for i in (0..qvs.len()).rev() {
        if qvs[i]-33 < min_qual {
            good = 0;
        } else {
            good += 1;
            if good == k {
                return i + k;
            }
        }
    }

    return 0;
}

fn process_fastq_set<T,K,V>(min_qual: u8,
                        p: usize,
                        permutation: &Vec<usize>,
                        input: T,
                        shard_sender: &mut ShardSender<Bsp<K,V>, BspShard>) -> u64
    where T: IntoIterator<Item = InputRead>, K:Kmer + Sync + Send, V:Smer<K>
{
    let mut iter = input.into_iter();
    let mut nreads : u64 = 0;
    loop {
        match iter.next() {
            Some(inp) => {
                // R1 is numbered with an even number, R2 is the next value
                let trim_len1 = find_trim_len(&inp.r1, min_qual, K::k());
                let r1b: Vec<u8> = inp.r1.read.as_bytes().iter().cloned().map(base_to_bits).collect();
                let r1_slice = &r1b[..trim_len1];
                msp_read::<K,V>(
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id,
                         r1_slice,
                         shard_sender);

                let trim_len2 = find_trim_len(&inp.r2, min_qual, K::k());
                let r2b: Vec<u8> = inp.r2.read.as_bytes().iter().cloned().map(base_to_bits).collect();
                let r2_slice = &r2b[..trim_len2];
                msp_read::<K,V>(
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id + 1,
                         r2_slice,
                         shard_sender);
                nreads += 2;
            }
            None => break,
        }
    }

    return nreads;
}

pub fn enumerate_pmers() -> Vec<Pmer> {
    let mut pmers: Vec<Pmer> = Vec::new();
    let p = Pmer::k();

    for i in 0..(1 << 2 * p) {
        let pmer = Pmer::from_u64(i << (64 - 2 * p) | p as u64);
        pmers.push(pmer);
    }

    pmers
}

/// Return a vector of length 1<<(2*p) ranks, one for each lexicographically
/// sorted pmer. The highest-ranked pmer is the least frequent in the counts HashMap.
fn pmer_inverse_freq_permutation(counts: &HashMap<Pmer, usize>) -> Vec<usize> {
    let mut freq_vec: Vec<(usize, Pmer)> = Vec::new();

    // Make a vector of (count, pmer) pairs
    for pmer in enumerate_pmers() {
        match counts.get(&pmer) {
            Some(count) => freq_vec.push((*count, pmer)),
            None => freq_vec.push((0, pmer)),
        }
    }

    // Sort to get the least common pmers at the top
    freq_vec.sort();

    // Tag each pmer with its rank, and sort back into pmer order
    let mut perm: Vec<(usize, usize)> = 
        freq_vec.into_iter().
        enumerate().
        map(|(rank, (_, pmer))| (pmer.to_u64() as usize, rank)).
        collect();

    perm.sort();

    // Return the vector of ranks the pmers
    perm.into_iter().map(|(_, rank)| rank).collect()
}


fn count_pmers_fastq<T: Iterator<Item = InputRead>>(mut input: T) -> HashMap<Pmer, usize> {
    let mut counter = HashMap::new();
    loop {
        match input.next() {
            Some(inp) => {

                // R1 is numbered with an even number, R2 is the next value
                let r1 = DnaString::from_dna_string(&inp.r1.read);
                let iter1: KmerIter<Pmer, _> = r1.iter_kmers();
                for pmer in iter1 {
                    let min_pmer = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

                let r2 = DnaString::from_dna_string(&inp.r2.read);
                let iter2: KmerIter<Pmer, _> = r2.iter_kmers();
                for pmer in iter2 {
                    let min_pmer = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

            }
            None => break,
        }
    }

    counter
}

pub struct MspMartian;

#[derive(Deserialize)]
struct Main {
    chunk: Vec<PathBuf>,
    whitelist: PathBuf,
    permutation: PathBuf,
    trim_min_qual: u8,
}

impl MartianStage for MspMartian {
    fn split(&self, args: JsonDict) -> Result<JsonDict, Error> {

        let fastqs : Vec<String> = serde_json::from_value(args["fastqs"].clone()).expect("invalid fastqs arg");


        let _whitelist = args["barcode_whitelist"].as_str().unwrap();
        let whitelist  = Path::new(_whitelist);

        let mut permutation_file = env::current_dir().unwrap();
        permutation_file.push("permutation.perm");

        // Compute the permutation
        let paths = fastqs.clone().iter().map(PathBuf::from).collect();
        msp_get_pmer_permutation(paths, &whitelist, &permutation_file);

        let mut chunks : Vec<Json> = Vec::new();
        for in_files in fastqs.chunks(8)
        {
            let chunk = json!({
                "chunk": in_files,
                "permutation": permutation_file,
                "__mem_gb": MEM_ALLOC_GB as f64,
                "__threads": 4
            });
            chunks.push(chunk);
        }

        let chunk_array = Value::Array(chunks);
        let mut cc =  Map::new();
        cc.insert("chunks".to_string(), chunk_array);
        Ok(cc)
    }


    fn main(&self, _args: JsonDict, outs: JsonDict) -> Result<JsonDict, Error> {
        let args: Main = serde_json::from_value(json!(_args)).unwrap();
        {
            let _out_fn = outs["chunks"].as_str().unwrap();
            let out_fn = Path::new(_out_fn);

            main_msp(args.trim_min_qual, &args.whitelist, args.chunk, &args.permutation, &out_fn);
        }
        Ok(outs)
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> Result<JsonDict, Error> {

        let mut chunk_vec: Vec<Json> = Vec::new();

        for c in chunk_outs {
            let chunk = c.get("chunks").unwrap();
            chunk_vec.push(chunk.clone());
        }

        let mut obj = Map::new();
        obj.insert("chunks".to_string(), json!(chunk_vec));
        Ok(obj)
    }
}


#[cfg(test)]
mod tests {
    use multifastq::Fastq;

    #[test]
    fn test_qv_trim_read() {
        let seq = "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC".to_string().into_bytes();
        let quals = "FFFFFFFFIFFFFFFFFFFIIIFFBFIFFFFIFBFFIBFIFBBFFIFFIFFFFFFFFFFFBBBBBBBBBB07BB7BB<BBBBBBBBBB".to_string().into_bytes();
        
        for k in 8..65 {
            for i in 0..quals.len() {
                // put in a bad qv and make sure we get the right length
                let mut myquals = quals.clone();
                myquals[i] = 34 as u8;

                let fq = Fastq {
                    read: String::from_utf8(seq.clone()).unwrap(),
                    qual: String::from_utf8(myquals).unwrap(),
                };

                let trim_length = super::find_trim_len(&fq, 10, k);

                let expected_length = if i < quals.len() - k { quals.len() } else if i >= k { i } else { 0 };
                assert_eq!(expected_length, trim_length);
            }
        }
    }
}
