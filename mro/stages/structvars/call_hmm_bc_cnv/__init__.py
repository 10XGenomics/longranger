#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import tenkit.hdf5 as tk_hdf5
import tenkit.reference as tk_reference
import subprocess
import os
import numpy as np
import longranger.sv.io as tk_sv_io
import longranger.genomic_tracks as lr_gt
from math import log10
import tenkit.bio_io as tk_io


__MRO__ = """
stage CALL_HMM_BC_CNV(
    in  string  subcommand,
    in  bam     possorted_bam,
    in  h5      fragments,
    in  int     fragment_version,
    in  int     bin_size,
    in  float   status_change_penalty,
    in  string  blacklist,
    in  float   min_prob,
    in  bool    allow_bin_size_adj,
    in  string  reference_path,
    in  float   max_pval,
    in  int     minimal_cnv_size,
    out int     final_bin_size,
    out bed     bc_cnv,
    out bedpe   bc_large_cnv,
    src py     "stages/structvars/call_hmm_bc_cnv",
) split using (
    in  int     chunk_id,
)
"""

NUM_FRAGMENT_BASES_20486_128G = 560196981191.0

def split(args):
    mem_in_gb = 16
    threads = 4
    chunk_defs = [{'chunk_id': 0, '__mem_gb': mem_in_gb,'__threads':threads}]
    return {'chunks': chunk_defs}


def main(args, outs):

    if args.fragments is None:
        outs.bc_cnv = None
        outs.bc_large_cnv = None
        return

    rust_env = os.environ.copy()
    rust_env["RUST_BACKTRACE"] = "1"
    final_blacklist = lr_gt.get_genomic_track(args.blacklist, "terminal_cnv", args.reference_path, "default_blacklist.bed")
    if final_blacklist is None:
        final_blacklist = args.possorted_bam+"_tmp"
        open(final_blacklist,'w').close()

    if args.subcommand == "bc" and args.fragments:
        frag_csv = outs.bc_cnv+".csv"
        bin_size, frag_version = convert_fragments_to_csv(args.fragments, frag_csv, args.bin_size, args.allow_bin_size_adj)
        cnv_args = ['hmm-bc-cnv', args.subcommand, frag_csv, args.possorted_bam, final_blacklist, outs.bc_cnv, "--fragver", str(frag_version), "--binsize", str(bin_size), "--probchange", str(args.status_change_penalty), "--minprob", str(args.min_prob)]
    elif args.subcommand == "read":
        cnv_args = ['hmm-bc-cnv', args.subcommand, args.possorted_bam, final_blacklist, outs.bc_cnv, "--binsize", str(args.bin_size), "--probchange", str(args.status_change_penalty)]
    elif args.subcommand == "asread":
        frag_csv = outs.bc_cnv+".csv"
        bin_size, frag_version = convert_fragments_to_csv(args.fragments, frag_csv, args.bin_size, args.allow_bin_size_adj)
        cnv_args = ['hmm-bc-cnv', args.subcommand, frag_csv, args.possorted_bam, final_blacklist, outs.bc_cnv, "--fragver", str(frag_version), "--binsize", str(bin_size), "--probchange", str(args.status_change_penalty), "--minprob", str(args.min_prob)]

    print cnv_args
    subprocess.check_call(cnv_args, env=rust_env)
    outs.final_bin_size = bin_size

    chroms = []
    starts1 = []
    end1 =[]
    starts2 = []
    end2 = []
    info_strs = []
    quals = []

    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)

    spikes = tk_io.get_target_regions(open(args.spikes))
    with open(outs.bc_cnv) as fin:
        for line in fin.readlines():
            if line.startswith('#') or line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track'): continue
            infos = line.strip().split("\t")
            cp = int(infos[3])
            ch = infos[0]
            s = int(infos[1])
            e = int(infos[2])

            # Some basic filtering
            if primary_contigs and ch not in primary_contigs:
                continue

            if cp == 2 or (e-s) < args.minimal_cnv_size:
                continue

            if cp > 2:
                if ch not in spikes: continue
                overlaps = spikes[ch].overlapping_regions(max(0, s-bin_size), e+bin_size)
                ln = len(overlaps)
                if ln > 0 and \
                    overlap(s-bin_size, s+bin_size, overlaps[0][0], overlaps[0][1]) and \
                    overlap(e-bin_size, e+bin_size, overlaps[ln-1][0], overlaps[ln-1][1]):
                        continue;

            chroms.append(infos[0])
            starts1.append(s)
            end1.append(s+1)
            starts2.append(e)
            end2.append(e+1)
            pval = float(infos[4])
            #if pval > args.max_pval:
            #    continue
            if pval < 1e-100:
                qual = 1000
            else:
                qual = int(-log10(pval)*10)
            quals.append(qual)
            if cp > 2:
                info_strs.append('TYPE=DUP;COPY=%d' % cp)
            elif cp < 2:
                info_strs.append('TYPE=DEL;COPY=%d' % cp)

    sv_df = tk_sv_io.create_sv_df(chroms, starts1, end1, chroms, starts2, end2, np.arange(len(chroms)), quals, info_strs = info_strs)
    tk_sv_io.write_sv_df_to_bedpe(sv_df, outs.bc_large_cnv)



def join(args, outs, chunk_defs, chunk_outs):
    co = chunk_outs[0]

    if co.bc_cnv is None:
        outs.bc_cnv = None
        outs.bc_large_cnv = None
        return
    else:
        os.rename(co.bc_cnv, outs.bc_cnv)
        os.rename(co.bc_large_cnv, outs.bc_large_cnv)


def convert_fragments_to_csv(h5_file, frag_csv, oiringal_bin_size, allow_bin_size_adj):
    data = tk_hdf5.read_data_frame(h5_file)
    data = data[data["num_reads"]>1]
    data.to_csv(frag_csv, index=False)
    num_cols = data.shape[1]
    frag_version = 0
    if num_cols == 13: frag_version = 2
    if num_cols == 12: frag_version = 1

    total_frag_bases = data["obs_len"].sum()
    if allow_bin_size_adj and total_frag_bases < NUM_FRAGMENT_BASES_20486_128G:
       return (min(100000, get_opt_bin_size(oiringal_bin_size*1.0/total_frag_bases*NUM_FRAGMENT_BASES_20486_128G)), frag_version)
    else:
        return (oiringal_bin_size, frag_version)


def get_opt_bin_size(raw_adj_bin_size):
    level = 1
    while level < 3.0e9:
        level2 = level*10
        if raw_adj_bin_size > level2:
            level = level2
            continue
        else: # level < raw_adj_bin_size < level2
            second_pos = (raw_adj_bin_size % level) / level
            base = int(raw_adj_bin_size)/level
            if second_pos < 0.5: return base*level
            else: return (base+1)*level

def overlap(s1, e1, s2, e2):
    return not ((e1 < s2) or (s1 > e2))
