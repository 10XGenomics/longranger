extern crate libc;
use self::libc::{c_int, size_t, c_char};

use self::libc::int32_t;
use self::libc::int64_t;
use std::ffi::{CString};
use std::fmt;
use std::ptr;
use std::io::Result;
use std::io::{Error, ErrorKind};

use bwt;


#[repr(C)]
struct bntann1_t {
    offset: i64,
    len: i64,
    n_ambs: i64,
    gi: u32,
    name: *const c_char,
    anno: *const c_char,
}

#[repr(C)]
struct bntamb1_t {
    offset: i64,
    len: i32,
    amb: u8,
}

#[repr(C)]
struct bntseq_t {
    l_pac: i64,
    n_seqs: i32,
    seed: u32,
    anns: *const bntann1_t, // n_seqs elements
    n_holes: i32,
    ambs: *const bntamb1_t, // n_holes elements
    fp_pac: *const c_int,     // FILE pointer
}

//impl bntseq_t {
//    fn get_ann(&self, i : int32_t) -> bntann1_t
//    {
//        unsafe { &self.anns[i] }
//    }
//}

impl fmt::Debug for bntseq_t
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.pad("bnt_seq")
        //f.debug_struct("bntseq_t")
        //.pad("l_pac", &self.l_pac)
        //.field("n_seqs", &self.n_seqs)
        //.field("seed", &self.seed)
        //.finish()
    }
}



#[repr(C)]
#[derive(Debug)]
struct bwaidx_t<'a> {
    bwt: &'a bwt::bwt_t,        // FM-index
    bns: &'a bntseq_t,     // information on the reference sequences
    pac: &'a u8,           // the actual 2-bit encoded reference sequences with 'N' converted to a random base
}


// Wrap bwa types
#[repr(C)]
struct mem_seed_t {
	rbeg: int64_t,
	qbeg: int32_t,
    len:  int32_t,
}

#[repr(C)]
struct mem_chain_t {
    n: c_int,
    m: c_int,
    pos: int64_t,
    seeds: *mut mem_seed_t,
}

#[repr(C)]
struct mem_chain_v {
    n : size_t,
    m : size_t,
    a : *mut mem_chain_t
}

#[repr(C)]
struct mem_opt_t {
    a: c_int,                // match score
    b: c_int,                // mismath penalty
    q: c_int,                // gap open penalty
    r: c_int,                // gap extension penalty
	// int a, b, q, r;      // match score, mismatch penalty and gap open/extension penalty. A gap of size k costs q+k*r

    pen_unpaired: c_int,     // phred-scaled penalty for unpaired reads
    pen_clip5: c_int,        // 5' clipping penalty
    pen_clip3: c_int,        // 3' clipping penalty

	w: c_int,                // band width
    zdrop: c_int,            // Z-dropoff

	thresh: c_int,               // output score threshold; only affecting output
	flag : c_int,           // see MEM_F_* macros
	min_seed_len: c_int,     // minimum seed length
	split_factor: c_int,      // split into a seed if MEM is longer than min_seed_len*split_factor
	split_width: c_int,       // split into a seed if its occurence is smaller than this value
	max_occ: c_int,            // skip a seed if its occurence is larger than this value
    max_chain_gap: c_int,     // do not chain seed if it is max_chain_gap-bp away from the closest seed
	n_threads: c_int,         // number of threads
	chunk_size: c_int,        // process chunk_size-bp sequences in a batch
    mask_level: f32,          // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
    chain_drop_ratio: f32,    // drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
	mask_level_redun: f32,
	mapq_coef_len: f32,
	mapq_coef_fac: c_int,
	max_ins: c_int,           // when estimating insert size distribution, skip pairs with insert longer than this value
	max_matesw: c_int,        // perform maximally max_matesw rounds of mate-SW for each end
	mat: [i8; 25],              // scoring matrix; mat[0] == 0 if unset
}

#[repr(C)]
struct mem_alnreg_t {
    rb: int64_t,            // [rb,re): reference sequence in the alignment
    re: int64_t,
    qb: c_int,              // [qb,qe): query sequence in the alignment
    qe: c_int,
    score: c_int,           // best local SW score
	truesc: c_int,          // actual score corresponding to the aligned region; possibly smaller than $score
    sub: c_int,             // 2nd best SW score
	csub: c_int,            // SW score of a tandem hit
	sub_n: c_int,           // approximate number of suboptimal hits
	w: c_int,               // actual band width used in extension
    seedcov: c_int,         // length of regions coverged by seeds
	secondary: c_int,       // index of the parent hit shadowing the current hit; <0 if primary
	hash: u64,
}

#[repr(C)]
struct mem_alnreg_v {
    n: size_t,
    m: size_t,
    a: *mut mem_alnreg_t,
}

#[repr(C)]
struct mem_pestat_t {
    low: c_int,         // lower and upper bounds within which a read pair is considered to be properly paired
    high: c_int,
    failed: c_int,      // non-zero if the orientation is not supported by sufficient data
    avg: f64,           // mean and stddev of the insert size distribution
    std: f64,
}

#[repr(C)]
struct mem_aln_t { // This struct is only used for the convenience of API.
    pos: int64_t,  // forward strand 5'-end mapping position
    rid: c_int,    // reference sequence index in bntseq_t; <0 for unmapped
	flag: c_int,   // extra flag

	//uint32_t is_rev:1, mapq:8, NM:23; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
	field: u32,
    n_cigar: c_int, // number of CIGAR operations
	cigar: *mut u32,    // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
    score: c_int,
    sub: c_int,
}

fn show_idx(idx: bwaidx_t){
    let nseqs = idx.bns.n_seqs;
    println!("contigs: {}", nseqs);
    for i in 0..nseqs   {
        let ann = unsafe { ptr::read(idx.bns.anns.offset(i as isize)) };
        println!("len: {}, offset: {}", ann.len, ann.offset);
        //let name = unsafe { CStr::from_ptr(ann.name).to_bytes() };
        //println!("name: {:?}", name);
        println!("name: {:?}", ann.name);
    }
}

//fn align_read(opts: *const mem_opt_t, bwaidx: bwaidx_t, sequence: &str)
//{
//    mem_chain(opts,&(bwaidx.bwt), bwaidx.bnsseq.l_pac, sequence.len )
//}

#[link(name="bwa")]
extern {
    // ASCII -> nucleotide code mapping table
    static nst_nt4_table: [c_char; 256];

    fn mem_chain(opt: *const mem_opt_t, bwt: *const bwt::bwt_t, l_pac: int64_t, len: c_int, seq: *const u8) -> mem_chain_v;

    /// Sort chains by the total size of their matches
    /// filter out chains that overlap better scoring chains
    fn mem_chain_flt(opt: *const mem_opt_t, n_chn: c_int, chains: *mut mem_chain_t) -> c_int;

    /// Print chains
    fn mem_print_chain(bns: *const bntseq_t, chn: *const mem_chain_v);

    fn mem_opt_init() -> *mut mem_opt_t;

    fn bwa_idx_load_bwt(hint: *const u8) -> *const bwt::bwt_t;

    fn bwa_idx_load(hint: *const libc::c_char, which: c_int) -> *const bwaidx_t;
}

fn load_idx<'a>(idx_file: &str) -> Result<bwaidx_t<'a>> {
    let bwaidx = unsafe {
        let filename = CString::new(idx_file).unwrap();
        let r = bwa_idx_load(filename.as_ptr(), 7);
        if r == ptr::null()
        {
            Err(Error::new(ErrorKind::NotFound, idx_file))
        }
        else
        {
            Ok(ptr::read(r))
        }
    };

    println!("idx: {:?}", bwaidx);
    bwaidx
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use std::path::{Path};
    use std::ptr;

    #[test]
    fn it_works()
    {
        println!("cd: {:?}", env::current_dir().unwrap());
    }

    #[test]
    fn load_test_index(){
        println!("cd: {:?}", env::current_dir().unwrap());
        let test_file = "test_data/test-ref.fasta";
        //let bwaidx = super::load_idx(test_file).unwrap();
        //println!("idx: {:?}", bwaidx);
        //show_idx(bwaidx);
        //super::show_idx(bwaidx);
        //panic!("ASdf")
    }

    #[test]
    fn test_align(){
        //let opts = super::mem_opt_init();


    }
}
