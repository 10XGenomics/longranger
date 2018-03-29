

extern crate libc;
use self::libc::{c_int, size_t};
use std::fmt;

#[allow(non_camel_case_types)]
type bwtint_t = u64;

#[repr(C)]
pub struct bwt_t {
	primary : u64, // S^{-1}(0), or the primary index of BWT
	l2: [u64; 5],  // C(), cumulative count
	seq_len: u64,  // sequence length
	bwt_size: u64, // size of bwt, about seq_len/4
	bwt: *mut u32, // BWT
	// occurance array, separated to two parts
    cnt_table: [u32; 256],
	// suffix array
	sa_intv: c_int,
	n_sa: u64,
    sa: *mut u64,
}

impl fmt::Debug for bwt_t {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		f.pad("bwt_t");
		f.write_fmt(format_args!("l2: {:?}", &self.l2));
		f.write_fmt(format_args!("seq_len: {}", &self.seq_len));
		f.write_fmt(format_args!("bwt_size: {}", &self.bwt_size))
		//f.debug_struct("bwt_t")
		//.field("l2", &self.l2)
		//.field("seq_len", &self.seq_len)
		//.field("bwt_size", &self.bwt_size)
		//.finish()
    }
}

#[repr(C)]
pub struct bwtintv_t {
    x: [bwtint_t; 3],
    info: bwtint_t,
}

#[repr(C)]
pub struct bwtintv_v {
    n: size_t,
    m: size_t,
    a: *mut bwtintv_t,
}
