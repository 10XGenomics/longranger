//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use utils;
use std::path::Path;

use debruijn::Mer;
use debruijn::paths::{DebruijnGraph};
use csv;

use cmd_msp::Kmer1;

#[derive(Serialize)]
struct GraphNode {
    id: usize,
    len: usize,
    //num_bcs: usize,
    exts_left: u8,
    exts_right: u8,
    sequence: String,
}


pub fn main_graph_stats(graph: &Path, stats: &Path)
{
    let graph : DebruijnGraph<Kmer1, Vec<u32>> = utils::read_obj(graph).expect("read graph");
    //let bcs : MultiVec<u32> = utils::read_obj(bcs).expect("read bcs");

    //println!("bcs: start_pos: {:?}, items: {:?}, vec_len: {:?}", &bcs.start_pos[0..10], &bcs.items[0..40], &bcs.vec_len[0..10]);

    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_path(stats).expect("open csv");

    for i in 0 .. graph.len()
    {
        let e = graph.get_node(i);

        let record = GraphNode {
            id: e.node_id,
            len: e.sequence().len(),
            //num_bcs: bcs.vec_len[e.id] as usize,
            exts_left: e.l_edges().len() as u8,
            exts_right: e.r_edges().len() as u8,
            sequence: e.sequence().to_string(),
        };

        wtr.serialize(record).expect("csv write error");
    }
}

pub fn main_write_gfa(graph: &Path, gfa_out: &Path)
{
    let graph : DebruijnGraph<Kmer1, Vec<u32>> = utils::read_obj(graph).expect("read graph");
    graph.to_gfa(gfa_out);
}

/*
pub fn write_graph_fasta(_graph: &Path, fasta: &Path) {

    let graph: TempGraph = utils::read_obj(_graph).expect("can't read graph");
    let edge_db = debruijn::EdgeDb::new(&graph);

    let mut _wtr = File::create(fasta).unwrap();
    let mut wtr = BufWriter::new(_wtr);

    for edge_id in 0 .. graph.start_pos.len() {
        let e = debruijn::make_vedge(&graph, &edge_db, edge_id);
        writeln!(wtr, ">{}", e.id).unwrap();
        writeln!(wtr, "{}", e.sequence.to_dna_string()).unwrap();
    }
}



pub fn write_graph_bcs_matrix(_graph: &Path, _bcs: &Path, mm_file: &Path) {
    let graph : TempGraph = utils::read_obj(_graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(_bcs).expect("read bcs");

    let max_bc_id = *bcs.items.iter().max().unwrap();
    let num_edges = graph.start_pos.len();

    let mut _wtr = File::create(mm_file).unwrap();
    let mut wtr = BufWriter::new(_wtr);

    // Main header
    writeln!(wtr, "%%MatrixMarket {} {} {} {}", "matrix", "coordinate", "pattern", "general").unwrap();
    
    // Comment
    writeln!(wtr, "% edge -> bc matrix for graph: {:?}", _graph).unwrap();

    // size line
    writeln!(wtr, "{} {} {}", num_edges, max_bc_id+1, bcs.items.len()).unwrap();

    for edge_id in 0 .. graph.start_pos.len() {
        let bc_slice = bcs.get_slice(edge_id);

        for bc_id in bc_slice {
            writeln!(wtr, "{} {}", edge_id + 1, bc_id + 1).unwrap();
        }
    }
}



pub fn write_basevector_graph(_graph: &Path, bv_path: &Path) {

    let graph: TempGraph = utils::read_obj(_graph).expect("can't read graph");

    let mut _wtr = File::create(bv_path).unwrap();
    let mut wtr = BufWriter::new(_wtr);
    graph.write_to_sn_format(&mut wtr);
}
*/