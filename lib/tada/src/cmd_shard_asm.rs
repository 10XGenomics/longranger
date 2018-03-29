//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command for assembling MSP shard Kmers into short DeBruijn graph snippets
//!
//! Enumerates all Kmers observed in the MSP shard and the barcodes observed on each.  Filters
//! the kmers by an arbitrary rule.  Then labels each Kmer with all the observed extensions from it.
//! Note, the extension in this phase will be a super-set of to Kmers that pass the filtering.
//! This will be resolved later.


use std::path::{Path, PathBuf};
use std::collections::{HashMap};
use std::fs;
use std::hash::{Hash, BuildHasher};
use std::env;
use std::mem;
use std::cmp::max;
use rayon::prelude::*;
use rayon;

use shardio::{ShardWriteManager, ShardReaderSet, ShardDef};
use debruijn::Kmer;
use debruijn::{self, Exts};
use debruijn::paths::DebruijnGraph;
use debruijn::compression::{SimpleCompress, PathCompression};
use cmd_msp::{Kmer1, Lmer3, Bsp};

use martian::{Json, JsonDict, MartianStage, Error};
use serde_json::{self, Value};
use serde_json::map::Map;


pub type DbG = DebruijnGraph<Kmer1, u16>;

pub fn est_gb_hashmap<K: Eq + Hash, V, S: BuildHasher>(hash_map: HashMap<K,V,S>) -> f32 {
    11.0/10.0 * ((hash_map.capacity() * (mem::size_of::<K>() + mem::size_of::<V>() + mem::size_of::<u64>())) as f32) / 1e9
}

pub fn est_gb_vec<V>(vec: Vec<V>) -> f32 {
    (vec.capacity() * ( mem::size_of::<V>())) as f32 / 1e9
}

struct GraphShard;
impl<K:Kmer, D> ShardDef<DebruijnGraph<K, D>> for GraphShard {
    fn get_shard(_: &DebruijnGraph<K, D>) -> usize {
        0
    }
} 

pub fn main_shard_asm(min_kmer_obs: usize, chunk_id: usize, total_chunks: usize, shard_chunks: Vec<PathBuf>, sedge_asm_out: &Path, _sedge_bcs_out: &Path, num_threads: usize) {

    let reader = ShardReaderSet::open(&shard_chunks);
    let my_shards: Vec<(usize, usize)> = (0..reader.num_shards()).filter(|x| x % total_chunks == chunk_id).enumerate().collect();
    let my_total_shards = my_shards.len();

    if num_threads > 0 {
        println!("setup rayon with {} threads", num_threads);
        let cfg = rayon::Configuration::new().num_threads(num_threads);
        rayon::initialize(cfg).unwrap();
    }

    let graph_shards: ShardWriteManager<DbG, GraphShard> = ShardWriteManager::new(sedge_asm_out, 2048, 1, 1);
    let sender = graph_shards.get_sender();

    // process each shard, and collect results
    let mut results = Vec::new();
    my_shards.into_par_iter().map_with(sender, |sender, (shard_count, shard_id)| {
        info!("Processing shard {} of {}", shard_count, my_total_shards);
        let mut bsps: Vec<Bsp<Kmer1, Lmer3>> = Vec::new();

        reader.read_shard(shard_id, &mut bsps);
        let seqs = bsps.iter().map(|x| (x.sequence, Exts::empty(), x.bc)).collect();

        // FIXME - define 3 read / 2 BC filter here
        let filter = debruijn::filter::CountFilter::new(min_kmer_obs);
        let (valid_kmers, _) = debruijn::filter::filter_kmers(&seqs, filter, false);
        info!("shard_asm: {} valid kmers", valid_kmers.len());

        let cmp = SimpleCompress::new(|a: u16, b: &u16| { max(a, *b) });
        let path_comp: PathCompression<Kmer1, u16, _> = PathCompression::new(false, cmp);
        let graph = path_comp.build_nodes(&valid_kmers).finish();

        sender.send(graph);

        info!("Shard:{}. BPSs:{}, valid_kmers:{}", shard_count, est_gb_vec(bsps), est_gb_vec(valid_kmers));
        1
    }).collect_into(&mut results);
}


pub struct ShardAsmMartian;

impl MartianStage for ShardAsmMartian {
    fn split(&self, _: JsonDict) -> Result<JsonDict, Error> {

        let total_chunks = 256;

        let mut chunks = Vec::new();
        for chunk_id in 0..total_chunks
        {
            let chunk = json!({
                "chunk_id": chunk_id,
                "total_chunks": total_chunks,
                "__threads": 4,
                "__mem_gb": 36.0
            });
            chunks.push(chunk);
        }

        let chunk_array = Value::Array(chunks);
        let mut cc = Map::new();
        cc.insert("chunks".to_string(), chunk_array);
        Ok(cc)
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> Result<JsonDict, Error> {

        let min_kmer_obs = args["min_kmer_obs"].as_u64().unwrap() as usize;
        let shard_chunks: Vec<PathBuf> = serde_json::from_value(args["chunks"].clone()).unwrap();

        let chunk_id = args["chunk_id"].as_u64().unwrap() as usize;
        let total_chunks = args["total_chunks"].as_u64().unwrap() as usize;
        let threads = args["__threads"].as_u64().unwrap() as usize;

        let sedge_asm = Path::new(outs.get("sedge_asm").unwrap().as_str().unwrap());
        let sedge_bcs = Path::new(outs.get("sedge_bcs").unwrap().as_str().unwrap());

        main_shard_asm(min_kmer_obs, chunk_id, total_chunks, shard_chunks, sedge_asm, sedge_bcs, threads);
        Ok(outs.clone())
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> Result<JsonDict, Error> {

        let mut sedge_asms: Vec<Json> = Vec::new();
        let mut sedge_bcss : Vec<Json> = Vec::new();

        let pwd = env::current_dir().unwrap();

        for (idx, c) in chunk_outs.into_iter().enumerate() {

            //let _old_name =
            let old_name = c["sedge_asm"].as_str().unwrap();
            let mut new_name = pwd.clone();
            new_name.push(format!("chunk{}.sedge_asm", idx));
            fs::rename(old_name, new_name.clone()).unwrap();
            sedge_asms.push(Value::String(new_name.to_str().unwrap().to_string()));

            let old_name = c["sedge_bcs"].as_str().unwrap();
            let mut new_name = pwd.clone();
            new_name.push(format!("chunk{}.sedge_bcs", idx));
            fs::rename(old_name, new_name.clone()).unwrap();
            sedge_bcss.push(Value::String(new_name.to_str().unwrap().to_string()));
        }

        let mut obj = Map::new();
        obj.insert("sedge_asm".to_string(), Value::Array(sedge_asms));
        obj.insert("sedge_bcs".to_string(), Value::Array(sedge_bcss));
        Ok(obj)
    }
}