// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
/*

// FIXME - update to current petgraph/daggy, or drop

use std::cmp::{max, min, Ordering};
use std::iter;
use std::collections::HashMap;
use daggy::Dag;
use daggy::walker::Walker;
use petgraph;
use petgraph::graph::{DefIndex, NodeIndex};
use std::fs::File;
use std::fmt;
use std::io::Write;
use std::path::PathBuf;
use std::collections::HashSet;

/// Traceback pointers for POA alignment. This generalizes the standard SW
/// traceback by incudling the POA node id that you came from
#[derive(Debug)]
pub enum Path {
    Start,
    Match(NodeIndex),
    Mismatch(NodeIndex),
    Deletion(NodeIndex),
    Insertion,
}

/// Pointer to a position in the POA alignment -- for tracking
/// the best alignment position for traceback
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AlnPtr {
    pub score: i16,
    read_pos: u16,
    tpl_pos: NodeIndex,
}

/// All the data about a POA aligment. Allows the client to
/// inspect the alignment before deciding to include the read in the POA
/// Useful for the case of unbarcoded reads where we want to screen
/// a read against both haplotypes
pub struct CandidateAln {
    score_cols: Vec<Vec<(i16, Path)>>,
    node_lookup: HashMap<NodeIndex<DefIndex>, usize>,
    pub best_alignments: BestAlns,
    topo_sort: Vec<NodeIndex>,
}

impl CandidateAln {
    pub fn score() -> u16 {
        0
    }

    pub fn is_unique() -> bool {
        true
    }
}

impl Ord for AlnPtr {
    fn cmp(&self, other: &AlnPtr) -> Ordering {
        let me = (-self.score, self.read_pos, self.tpl_pos);
        let other = (-other.score, other.read_pos, other.tpl_pos);
        me.cmp(&other)
    }
}

impl PartialOrd for AlnPtr {
    fn partial_cmp(&self, other: &AlnPtr) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Little helper class for keeping track of the 'best'
/// objects seen so far.  Kinda like a top-k priority queue.
#[derive(Debug)]
pub struct TrackBest<T: Ord + PartialOrd + Copy> {
    max_n: usize,
    heap: Vec<T>,
    min_val: Option<T>,
}

impl<T: Ord + PartialOrd + Copy> TrackBest<T> {
    pub fn new(max_n: usize) -> TrackBest<T> {
        TrackBest {
            max_n: max_n,
            heap: Vec::with_capacity(max_n + 1),
            min_val: None,
        }
    }

    pub fn add(&mut self, v: T) {
        let min_val = match self.min_val {
            Some(min_val) => min(min_val, v),
            None => v,
        };

        self.min_val = Some(min_val);

        if self.heap.len() < self.max_n || v <= min_val {
            self.heap.push(v);
        }

        if self.heap.len() >= self.max_n * 2 {
            self.heap.sort();
            self.heap.truncate(self.max_n)
        }
    }

    pub fn get(&self) -> Vec<T> {
        let mut r = self.heap.clone();
        r.sort();
        r.truncate(self.max_n);
        r
    }

    pub fn min(&mut self) -> T {
        self.heap.sort();
        self.heap[0]
    }
}

pub type BestAlns = TrackBest<AlnPtr>;

/// POA nodes are labelled with a single base pair, ASCII encoded
struct NodeLabel(u8);

/// Display a node as it's ASCII char, rather than a u8
impl fmt::Display for NodeLabel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"{}", self.0 as char)
    }
}

/// Each node in the POA tracks the DNA base,
/// the set of reads aligned to this base, and
/// the number of reads that 'span' this node
struct PoaNode {
    base: NodeLabel,
    reads: Vec<u16>,
    spanning_reads: u16,
}

impl PoaNode {
    pub fn empty(base: u8) -> PoaNode {
        PoaNode {
            base: NodeLabel(base),
            reads: Vec::new(),
            spanning_reads: 0,
        }
    }

    pub fn init(base: u8, read_id: u16) -> PoaNode {
        PoaNode {
            base: NodeLabel(base),
            reads: vec![read_id],
            spanning_reads: 0,
        }
    }

    pub fn add_read(&mut self, read_id: u16) {
        self.reads.push(read_id);
    }
}

impl fmt::Debug for PoaNode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
               "[{} - {}, sp: {}]",
               self.base,
               self.reads.len(),
               self.spanning_reads)
    }
}


pub struct PoaSettings {

}

/// A POA is a DAG of POA nodes.  Also included for convenience
/// is a pointer to the initial start and end nodes of the POA,
/// the total number of reads in the POA, and the alignment scoring
/// parameters.
pub struct Poa {
    graph: Dag<PoaNode, ()>,
    start_node: NodeIndex,
    end_node: NodeIndex,
    mtch: i16,
    mismatch: i16,
    del: i16,
    ins: i16,
    min_overlap_score: i16,
    min_coverage: i16,
    low_cov_penalty: i16,
    pub num_reads: u16,
    pub rejected_phased_reads: u16,
}

impl Poa {
    pub fn new() -> Poa {
        let mut g = Dag::new();
        let start_node = g.add_node(PoaNode::empty(b'0'));
        let end_node = g.add_node(PoaNode::empty(b'1'));
        Poa {
            graph: g,
            start_node: start_node,
            end_node: end_node,
            mtch: 1,
            mismatch: -4,
            del: -5,
            ins: -5,
            min_overlap_score: 25,
            min_coverage: 2,
            low_cov_penalty: 2,
            num_reads: 0,
            rejected_phased_reads: 0,

        }
    }

    pub fn new_config(mismatch: i16, min_overlap_score: i16, min_coverage: i16) -> Poa {
        let mut g = Dag::new();
        let start_node = g.add_node(PoaNode::empty(b'0'));
        let end_node = g.add_node(PoaNode::empty(b'1'));
        Poa {
            graph: g,
            start_node: start_node,
            end_node: end_node,
            mtch: 1,
            mismatch: mismatch,
            del: -5,
            ins: -5,
            min_overlap_score: min_overlap_score,
            min_coverage: min_coverage,
            low_cov_penalty: 2,
            num_reads: 0,
            rejected_phased_reads: 0,
        }
    }

    pub fn num_reads(&self) -> usize
    {
        self.num_reads as usize
    }

    pub fn num_nodes(&self) -> usize
    {
        self.graph.node_count()
    }

    /// Create an initial POA graph from a single read
    pub fn add_first_read(&mut self, sequence: &Vec<u8>) {
        let gg = &mut self.graph;
        let mut last_node = self.start_node;

        for i in 0..sequence.len() {
            let mut node_weight = PoaNode::init(sequence[i], 0);
            node_weight.spanning_reads = 1;
            let new_node = gg.add_node(node_weight);
            let new_edge_res = gg.add_edge(last_node, new_node, ());
            assert!(new_edge_res.is_ok());
            last_node = new_node;
        }

        let new_edge_res = gg.add_edge(last_node, self.end_node, ());
        assert!(new_edge_res.is_ok());
    }

    /// Align a read to the graph, then traceback the alignment and add new nodes as required.
    pub fn add_read(&mut self, sequence: &Vec<u8>) {
        if self.num_reads == 0 {
            self.add_first_read(sequence);
            self.num_reads += 1;
            return;
        }

        let mut aln = self.try_align_read(sequence);
        self.thread_alignment(sequence, &aln.score_cols, &aln.node_lookup, &mut aln.best_alignments);
        self.num_reads += 1;
    }


    /// Add a candidate alignment to the POA.
    pub fn commit_alignment(&mut self, sequence: &Vec<u8>, mut aln: CandidateAln) {
        self.thread_alignment(sequence, &aln.score_cols, &aln.node_lookup, &mut aln.best_alignments);
        self.num_reads += 1;
    }

    /// Determine the consensus sequence implied by the POA
    #[inline(never)]
    pub fn consensus(&mut self) -> Vec<u8> {

        let topo = petgraph::algo::toposort(self.graph.graph());
        let mut score_vec: Vec<(i32, NodeIndex)> = Vec::with_capacity(topo.len());
        let mut ni_lookup: HashMap<NodeIndex<DefIndex>, usize> = HashMap::with_capacity(topo.len());

        for (idx, ni) in iter::Iterator::enumerate(topo.iter().cloned()) {
            ni_lookup.insert(ni, idx);
        }

        // Mark each node with a measure of 'coverage':
        // the number of reads that could have supported the node
        self.mark_aln_span(&topo, &ni_lookup);

        // Track overall best path
        let mut best_score = 0;
        let mut best_end_node = self.start_node;


        // Get the average coverage of nodes with >= 1 coverage
        let mut occ_weight = 0;
        let mut occ_nodes = 0;
        for ni in topo.iter().cloned() {
            let w = self.graph.node_weight(ni).unwrap();
            if w.reads.len() > 1 {
                occ_weight += w.reads.len();
                occ_nodes += 1;
            }
        }
        let avg_cov = (occ_weight / max(occ_nodes, 1)) as i16;

        //let min_reads = max(self.min_coverage, min(avg_cov-2, self.num_reads as i16 / 3)) as i32;
        let min_reads = max(self.min_coverage, min(self.min_coverage + 3, avg_cov / 2)) as i32;
        println!("min reads: {}, avg cov: {}", min_reads, avg_cov);

        // Viterbi recursion over DAG
        for ni in topo.iter().cloned() {

            let mut parents = self.graph.parents(ni);

            // Score the node -- score goes negative if we have less than 3 read,
            // or if less than half the reads go through this node
            let node = self.graph.node_weight(ni).unwrap();
            let mut cov_score = node.reads.len() as i32 - min_reads;
            if cov_score < 0 {
                // Extra penalty for going below the coverage threshold
                cov_score = cov_score * (self.low_cov_penalty as i32);
            }

            let cons_score = 2 * node.reads.len() as i32 - (node.spanning_reads as i32);
            let node_score = min(cov_score, cons_score);

            // Best local traceback
            let mut max_score = node_score;
            let mut max_parent: NodeIndex = self.start_node;

            loop {
                match parents.next_node(&self.graph) {
                    Some(parent_node) => {

                        let prev_score = score_vec[ni_lookup[&parent_node]].0;
                        let new_score = prev_score + node_score;
                        if new_score > max_score {
                            max_score = new_score;
                            max_parent = parent_node;
                        }
                    }
                    None => break,
                }
            }

            score_vec.push((max_score, max_parent));

            if max_score > best_score {
                best_score = max_score;
                best_end_node = ni;
            }
        }

        // Traceback max-path through DAG to get consensus sequence
        let mut rev_consensus = Vec::new();
        let mut state = best_end_node;

        if best_score > 0 {
            loop {
                let node_weight = self.graph.node_weight(state).unwrap();
                rev_consensus.push(node_weight.base.0);

                // Go to previous state
                state = score_vec[ni_lookup[&state]].1;

                // Pointer to 0 state means we are done
                if state == self.start_node {
                    break;
                }
            }
        }

        rev_consensus.reverse();
        rev_consensus
    }

    /// Take an alignment between the POA and a read, and add the new read to the POA
    pub fn thread_alignment(&mut self,
                            seq: &Vec<u8>,
                            cols: &Vec<Vec<(i16, Path)>>,
                            lu: &HashMap<NodeIndex<DefIndex>, usize>,
                            best_alns: &mut BestAlns)
                            -> (NodeIndex, NodeIndex) {


        let best_aln = best_alns.get()[0];
        let mut downstream_node = best_aln.tpl_pos;
        let mut current_node = best_aln.tpl_pos;
        let mut end_node: Option<NodeIndex> = None;

        //let mut read = best_aln.read_pos as usize;
        let mut read = seq.len();
        let read_end = best_aln.read_pos;

        let mut gg = &mut self.graph;
        let read_id = self.num_reads;

        if read > (read_end as usize)
        {
            downstream_node = self.end_node;
        }

        while read > (read_end as usize) {
            let new_upstream_node = gg.add_node(PoaNode::init(seq[read - 1], read_id));
            let new_edge_res = gg.add_edge(new_upstream_node, downstream_node, ());
            assert!(new_edge_res.is_ok());

            read = read - 1;
            downstream_node = new_upstream_node;
        }


        loop {
            let path = cols[lu[&current_node]].get(read).unwrap();

            match path.1 {
                Path::Start => break,

                // Target node exists, make sure an edge to downstream node exists
                Path::Match(upstream_node) => {
                    let edge = gg.find_edge(upstream_node, downstream_node);
                    match edge {
                        None => {
                            let _ = gg.add_edge(upstream_node, downstream_node, ());
                        }
                        _ => (),
                    };

                    gg.node_weight_mut(upstream_node).unwrap().add_read(read_id);

                    read = read - 1;
                    downstream_node = upstream_node;
                    current_node = upstream_node;
                }

                // Introduce new node representing mismatched base
                Path::Mismatch(upstream_node) => {
                    let new_upstream_node = gg.add_node(PoaNode::init(seq[read - 1], read_id));
                    let new_edge_res = gg.add_edge(new_upstream_node, downstream_node, ());
                    assert!(new_edge_res.is_ok());

                    read = read - 1;
                    downstream_node = new_upstream_node;
                    current_node = upstream_node;
                }

                // Jump over the current node
                Path::Deletion(upstream_node) => {
                    current_node = upstream_node;
                }

                // Introduce inserted node
                Path::Insertion => {
                    let new_upstream_node = gg.add_node(PoaNode::init(seq[read - 1], read_id));
                    let new_edge_res = gg.add_edge(new_upstream_node, downstream_node, ());
                    assert!(new_edge_res.is_ok());

                    read = read - 1;
                    downstream_node = new_upstream_node;
                }
            }

            end_node = end_node.or(Some(downstream_node))
        }

        (current_node, end_node.unwrap())
    }

    /// Mark each node with the number of reads 'spanning' it.  A read 'spans' a
    /// node if it is aligned both a predeccesor and a successor of the node.
    #[inline(never)]
    fn mark_aln_span(&mut self, topo: &Vec<NodeIndex>, lu: &HashMap<NodeIndex, usize>) {
        let mut fwd: Vec<HashSet<u16>> = Vec::new();
        let mut rev: Vec<HashSet<u16>> = Vec::new();
        let n = topo.len();

        for ni in topo.iter() {
            let mut set = HashSet::new();
            set.extend(self.graph.node_weight(*ni).unwrap().reads.iter());

            let mut parents = self.graph.parents(*ni);
            loop {
                match parents.next_node(&self.graph) {
                    Some(parent_node) => set.extend(fwd[lu[&parent_node]].iter()),
                    None => break,
                }
            }
            fwd.push(set);
        }

        for ni in topo.iter().rev() {
            let mut children = self.graph.children(*ni);

            let mut set = HashSet::new();
            set.extend(self.graph.node_weight(*ni).unwrap().reads.iter());

            loop {
                match children.next_node(&self.graph) {
                    Some(child_node) => set.extend(rev[n - 1 - lu[&child_node]].iter()),
                    None => break,
                }
            }

            rev.push(set);
        }

        for (idx, ni) in topo.iter().enumerate() {
            let node = self.graph.node_weight_mut(*ni).unwrap();
            node.spanning_reads = fwd[idx].intersection(&rev[n - 1 - idx]).fold(0, |c, _| c + 1);
        }
    }

    /// Align a read to the POA and return a CandidateAln, which can be used to traceback
    /// the alignment in order to add the read to the POA
    #[inline(never)]
    pub fn try_align_read(&self, seq: &Vec<u8>) -> CandidateAln {
        // Topo sort -- vector of NodeIndex
        let topo = petgraph::algo::toposort(self.graph.graph());
        let mut cols: Vec<Vec<(i16, Path)>> = Vec::with_capacity(topo.len());
        let mut ni_lookup: HashMap<NodeIndex<DefIndex>, usize> = HashMap::with_capacity(topo.len());

        for (idx, ni) in iter::Iterator::enumerate(topo.iter().cloned()) {
            ni_lookup.insert(ni, idx);
        }

        let mut best_alns: BestAlns = TrackBest::new(10);

        for ni in topo.iter().cloned() {
            let mut score_col: Vec<(i16, Path)> = Vec::with_capacity(seq.len() + 1);
            for _ in 0..seq.len() + 1 {
                score_col.push((-32000, Path::Start));
            }

            let mut parents = self.graph.parents(ni);
            let terminal_col = ni == self.start_node || ni == self.end_node;

            loop {
                match parents.next_node(&self.graph) {
                    Some(parent_node) => {
                        let tpl_base = self.graph.node_weight(parent_node).unwrap().base.0;
                        let prev_col = cols.get(ni_lookup[&parent_node]).unwrap();
                        self.aln_prev_col(&seq,
                                          tpl_base,
                                          ni,
                                          &mut score_col,
                                          parent_node,
                                          prev_col,
                                          &mut best_alns);
                    }
                    None => break,
                }
            }

            self.aln_cur_col(ni, &mut score_col, &mut best_alns, terminal_col);
            cols.push(score_col);
        }

        CandidateAln {
            score_cols: cols,
            node_lookup: ni_lookup,
            best_alignments: best_alns,
            topo_sort: topo,
        }
    }

    /// Update score_col of the alignment with the insertion moves
    #[inline(never)]
    fn aln_cur_col(&self,
                       cur_node: NodeIndex,
                       score_col: &mut Vec<(i16, Path)>,
                       best_alns: &mut BestAlns,
                       terminal_col: bool) {
        score_col[0] = (0, Path::Start);

        for i in 1..score_col.len() {
            let ins_score = score_col[i - 1].0 +
                            if terminal_col {
                0
            } else {
                self.ins
            };
            let ins_move = Path::Insertion;

            if ins_score > score_col[i].0 {
                score_col[i] = (ins_score, ins_move);
            }

            //if i == score_col.len() - 1
            if i > 20 || i == score_col.len() - 1
            {
                let aln_ptr = AlnPtr {
                    score: score_col[i].0,
                    read_pos: i as u16,
                    tpl_pos: cur_node,
                };
                best_alns.add(aln_ptr);
            }
        }
    }

    /// Update score_col with the match and deletion moves coming from prev_col
    #[inline(never)]
    fn aln_prev_col(&self,
                        seq: &Vec<u8>,
                        tpl_base: u8,
                        cur_node: NodeIndex,
                        score_col: &mut Vec<(i16, Path)>,
                        prev_idx: NodeIndex,
                        prev_col: &Vec<(i16, Path)>,
                        best_alns: &mut BestAlns) {
        score_col[0] = (0, Path::Start);

        for i in 1..score_col.len() {
            let del_score = prev_col[i].0 + self.del;
            let del_move = Path::Deletion(prev_idx);

            if del_score > score_col[i].0 {
                score_col[i] = (del_score, del_move);
            }

            let (match_score, match_move) = if seq[i - 1] == tpl_base {
                let match_score = prev_col[i - 1].0 + self.mtch;
                let match_move = Path::Match(prev_idx);
                (match_score, match_move)
            } else {
                let mismatch_score = prev_col[i - 1].0 + self.mismatch;
                let mismatch_move = Path::Mismatch(prev_idx);
                (mismatch_score, mismatch_move)
            };

            if match_score > score_col[i].0 {
                score_col[i] = (match_score, match_move);
            }

            //if i == score_col.len() - 1
            if i > 20 || i == score_col.len() - 1
            {
                let aln_ptr = AlnPtr {
                    score: score_col[i].0,
                    read_pos: i as u16,
                    tpl_pos: cur_node,
                };
                best_alns.add(aln_ptr);
            }
        }
    }

    /// Default rule for deciding whether to accept the alignment
    pub fn accept_aln(&self, aln: &CandidateAln) -> bool {
        let score = aln.best_alignments.get()[0].score;
        score > min(self.min_overlap_score, (self.num_nodes() / 2) as i16)
    }

    /// Write the POA to a dot file.
    pub fn to_dot(&self, path: PathBuf) {
        let dot = petgraph::dot::Dot::new(self.graph.graph());
        let text = format!("{:?}", dot);
        let mut f = File::create(path).expect("couldn't open file");
        f.write(text.as_bytes()).expect("io error");
    }
}

pub fn to_byte_vec(s: &str) -> Vec<u8> {
    let ss = s.to_string().replace(" ",  "");
    ss.bytes().collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn add_same() {
        let mut poa = Poa::new();
        let s1 = to_byte_vec("ACGATACAGATAG");
        poa.add_read(&s1);
        poa.add_read(&s1);

        poa.to_dot(PathBuf::from("test-same.dot"));

        assert_eq!(s1.len(), poa.graph.node_count() - 2);
    }


    #[test]
    fn add_overlapping() {
        let mut poa = Poa::new();
        let s1 = to_byte_vec("ACGATACAGATAG");
        let s2 = to_byte_vec(   "ATACAGATAGATC");
        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.to_dot(PathBuf::from("test-ovl.dot"));

        assert_eq!(s1.len() + 3, poa.graph.node_count() - 2);
    }


    #[test]
    fn add_snp() {
        let mut poa = Poa::new();
        let s1 = to_byte_vec("ACGATACAGATAG");
        let s2 = to_byte_vec(   "ATAGAGATAGATC");
        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.to_dot(PathBuf::from("test-snp.dot"));

        assert_eq!(s1.len() + 4, poa.graph.node_count() - 2);
    }


    #[test]
    fn add_del() {
        let mut poa = Poa::new();
        let s1 = to_byte_vec("ACGACGTTACAGATAG");
        let s2 = to_byte_vec("   ACGTTACAG TAGATC");
        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.to_dot(PathBuf::from("test-del.dot"));

        assert_eq!(s1.len() + 3, poa.graph.node_count() - 2);
        assert_eq!(s1.len() + 3 + 1 + 1 + 1, poa.graph.edge_count());
    }


    #[test]
    fn add_ins() {
        let mut poa = Poa::new();
        let s1 = to_byte_vec("ACGACGTTACA TATAG");
        let s2 = to_byte_vec("   ACGTTACAGTATAGATC");
        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.to_dot(PathBuf::from("test-ins.dot"));

        assert_eq!(s1.len() + 6, poa.graph.node_count());
        assert_eq!(s1.len() + 7, poa.graph.edge_count());
    }


    #[test]
    fn consensus_simple() {
        let mut poa = Poa::new();

        let s1 = to_byte_vec("ACGACATAGGTCGTTACAACGATATAGAGAAGTCATCGT");
        poa.add_read(&s1);
        poa.add_read(&s1);
        poa.add_read(&s1);
        assert_eq!(poa.consensus(), s1);
    }


    #[test]
    fn consensus_snp() {
        let mut poa = Poa::new();

        let s1 = to_byte_vec("ACGACATAGGTCGTTACAACGATATAGAGAAGTCATCGTAAAA");
        let s2 = to_byte_vec("ACGACATAGGTCGTTACAGCGATATAGAGAAGTCATCGT");
        //                                      *

        poa.add_read(&s1);
        poa.add_read(&s2);
        poa.add_read(&s2);
        assert_eq!(poa.consensus(), s2);
    }


    #[test]
    fn consensus_del() {
        let mut poa = Poa::new();

        let s1 = to_byte_vec("ACGACATAGGTCGTTAC  CGATATAGAGAAGTCATCGTAAAA");
        let s2 = to_byte_vec("ACGACATAGGTCGTTACAGCGATATAGAGAAGTCATCGT");
        //                                     **

        poa.add_read(&s1);
        poa.add_read(&s2);
        poa.add_read(&s2);
        assert_eq!(poa.consensus(), s2);
    }

    #[test]
    fn consensus_ins() {
        let mut poa = Poa::new();

        let s1 = to_byte_vec("GGGACGACATAGGTCGTTACAAAAAACGATATAGAGAAGTCATCGTAAAA");
        let s2 = to_byte_vec("   ACGACATAGGTCGTTAC    AGCGATATAGAGAAGTCATCGT");
        //                                        *****

        poa.add_read(&s1);
        poa.add_read(&s2);
        poa.add_read(&s2);
        assert_eq!(poa.consensus(), s2);
    }

    #[test]
    fn subtle_ins() {
        let mut poa = Poa::new();

        let s1 = to_byte_vec("GGGACGACATAGGTCGTTAC AAAAAACGATATAGAGAAGTCATCGTAAAA");
        let s2 = to_byte_vec("GGGACGACATAGGTCGTTACAAAAAAACGATATAGAGAAGTCATCGTAAAA");

        // Add 2 reads with insertion, 3 reads without
        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.add_read(&s1);
        poa.add_read(&s2);

        poa.add_read(&s1);

        assert_eq!(poa.consensus(), s1);
    }


    /*
    #[test]
    fn add_extend_add_mm() {
        let mut poa = Poa::new();

        // 3rd read should align to 2nd read until the last 4 bp
        let s1 = to_byte_vec("GGGACGACATAGGTCGTTACAAAAAACGATATAGAGAAGTCATCGTAAAA");
        let s2 = to_byte_vec("          GGTCGTTACAAAAAAACGATATAGAGAAGTCATCGTAAAATCGTACAGATATG");
        let s3 = to_byte_vec("          GGTCGTTACAAAAAAACGATATAGAGAAGTCATCGTAAAATCGTACAGCCGCAC");

        poa.add_read(&s1);
        poa.add_read(&s2);
        poa.add_read(&s3);
        poa.to_dot(PathBuf::from("test-extend.dot"));

        assert_eq!(poa.consensus(), s1);
    }
    */
}

*/
