// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct Fragment {
    pub chrom:      String,
    pub start_pos:  i32,
    pub end_pos:    i32,
    pub obs_len:    i32,
    pub est_len:    i32,
    pub bc:         String,
    pub num_reads:  i32,
    pub molecule_id:u32,
}

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct FragmentV1 {
    pub bc:         String,
    pub bc_est_len: i32,
    pub bc_mean_reads_per_fragment: f32,
    pub bc_num_reads:i32,
    pub chrom:      String,
    pub end_pos:    i32,
    pub est_len:    i32,
    pub num_reads:  i32,
    pub num_reads_se: i32,
    pub obs_len:    i32,
    pub start_pos:  i32,
    pub molecule_id:u32,
}

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct FragmentV2 {
    pub bc:         String,
    pub bc_est_len: i32,
    pub bc_mean_reads_per_fragment: f32,
    pub bc_num_reads:i32,
    pub bc_num_unmapped_reads: i32,
    pub chrom:      String,
    pub end_pos:    i32,
    pub est_len:    i32,
    pub num_reads:  i32,
    pub num_reads_se: i32,
    pub obs_len:    i32,
    pub start_pos:  i32,
    pub molecule_id:u32,
}

impl FragmentV1 {
    pub fn to_simple(self) -> Fragment {
        Fragment {
            chrom:      self.chrom,
            start_pos:  self.start_pos,
            end_pos:    self.end_pos,
            obs_len:    self.obs_len,
            est_len:    self.est_len,
            bc:         self.bc,
            num_reads:  self.num_reads,
            molecule_id: self.molecule_id,
        }
    }
}

impl FragmentV2 {
    pub fn to_simple(self) -> Fragment {
        Fragment {
            chrom:      self.chrom,
            start_pos:  self.start_pos,
            end_pos:    self.end_pos,
            obs_len:    self.obs_len,
            est_len:    self.est_len,
            bc:         self.bc,
            num_reads:  self.num_reads,
            molecule_id: self.molecule_id,
        }
    }
}
