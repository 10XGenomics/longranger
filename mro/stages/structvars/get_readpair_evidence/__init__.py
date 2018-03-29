#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import numpy.random
import tenkit.pandas as pd
import martian
import json

import tenkit.bam as tk_bam
import longranger.sv.utils as tk_sv_utils
import longranger.sv.io as tk_sv_io
import longranger.sv.readpairs as tk_readpairs

__MRO__ = """
stage GET_READPAIR_EVIDENCE(
    in  bam   input,
    in  bedpe sv_variants,
    in  int   break_extend,
    in  int   min_mapq,
    in  int   min_reads_to_call,
    in  int   min_lr_to_call,
    in  float rp_lr_multiplier,
    in  json  insert_sizes,
    in  json  basic_summary,
    in  bool  best_only,
    out bedpe sv_variants,
    src py    "stages/structvars/get_readpair_evidence",
) split using (
    in  int   start_idx,
    in  int   stop_idx,
)
"""

MAX_READPAIRS = 500

def split(args):
    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_variants)

    nsvs = sv_df.shape[0]
    nbreaks_per_chunk = max(20, int(np.ceil(nsvs / 100.0))) # avoid overchunking
    nchunks = int(np.ceil(nsvs / float(nbreaks_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nbreaks_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nbreaks_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':16})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0}]
    return {'chunks': chunk_defs}


def main(args, outs):

    bedpe_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_variants)
    bedpe_df = tk_sv_utils.get_dataframe_loc(bedpe_df, list(range(int(args.start_idx), int(args.stop_idx))))

    max_insert, ins_logsf_fun = tk_sv_utils.get_insert_size_info(args.insert_sizes)
    if max_insert is None:
        martian.throw('No Q60 reads')

    with open(args.basic_summary, 'r') as f:
        summary = json.load(f)

    chimera_rate_del = summary['far_chimera_rate']
    chimera_rate_inv = summary['far_chimera_rate'] + summary['same_dir_chimera_rate']
    chimera_rate_trans = summary['far_chimera_rate']
    chimera_rate_dup = summary['far_chimera_rate'] + summary['outward_dir_chimera_rate']

    chimera_rates = {tk_readpairs.DEL_STR:chimera_rate_del,
                     tk_readpairs.INV_STR:chimera_rate_inv,
                     tk_readpairs.TDUP_STR:chimera_rate_dup,
                     tk_readpairs.TRANS_FF_STR:chimera_rate_trans,
                     tk_readpairs.TRANS_FR_STR:chimera_rate_trans,
                     tk_readpairs.TRANS_RR_STR:chimera_rate_trans,
                     tk_readpairs.TRANS_RF_STR:chimera_rate_trans}

    in_bam = tk_bam.create_bam_infile(args.input)

    out_quals = []
    out_infos = []
    out_chroms1 = []
    out_starts1 = []
    out_stops1 = []
    out_chroms2 = []
    out_starts2 = []
    out_stops2 = []

    for i, (_, row) in enumerate(bedpe_df.iterrows()):
        in_svt = tk_sv_io.get_sv_type(row.info)
        
        if row.chrom1 == row.chrom2:
            max_ext = min(args.break_extend, int((row.start2 - row.stop1) / 3.0))
            r1 = (max(0, row.start1 - args.break_extend), row.stop1 + max_ext)
            r2 = (max(0, row.start2 - max_ext), row.stop2 + args.break_extend)
            if r1[1] > r2[0]:
                starts = [r1[0]]
                stops = [r2[1]]
                chroms = [row.chrom1]
            else:
                starts = [r1[0], r2[0]]
                stops = [r1[1], r2[1]]
                chroms = [row.chrom1, row.chrom2]
        else:
            r1 = (max(0, row.start1 - args.break_extend), row.stop1 + args.break_extend)
            r2 = (max(0, row.start2 - args.break_extend), row.stop2 + args.break_extend)
            starts = [r1[0], r2[0]]
            stops = [r1[1], r2[1]]
            chroms = [row.chrom1, row.chrom2]

        readpairs = tk_readpairs.get_readpairs2(in_bam, chroms, starts, stops,
            max_insert = max_insert, min_mapq = args.min_mapq)

        # Distal readpairs across the breakpoints
        dist_readpairs = filter(filter_fun(in_bam, row.chrom1, row.chrom2, r1, r2), readpairs)
            
        if len(dist_readpairs) > MAX_READPAIRS:
            sel = numpy.random.choice(len(dist_readpairs), MAX_READPAIRS)
        else:
            sel = np.arange(len(dist_readpairs))
        dist_readpairs = [dist_readpairs[ridx] for ridx in sel]

        res_arr = tk_readpairs.get_sv_lr(dist_readpairs, ins_logsf_fun, max_insert, chimera_rates)

        if len(res_arr) == 0:
            out_quals.append(row.qual)
            out_chroms1.append(row.chrom1)
            out_starts1.append(row.start1)
            out_stops1.append(row.stop1)
            out_chroms2.append(row.chrom2)
            out_starts2.append(row.start2)
            out_stops2.append(row.stop2)
            out_infos.append(row['info'])
        else:
            if args.best_only:
                res_arr = sorted(res_arr, key = lambda x: x[0], reverse = True)
                res_arr = [res_arr[0]]
                
            for (lr, num_split, num_pairs, sv_len, support_range, svt, support_readpairs) in res_arr:
                range1, range2 = support_range
                if num_split + num_pairs >= args.min_reads_to_call and lr >= args.min_lr_to_call and not range1 is None and not range2 is None:
                    out_quals.append(row.qual + args.rp_lr_multiplier * lr)
                    out_chroms1.append(row.chrom1)
                    out_starts1.append(range1[0])
                    out_stops1.append(range1[1])
                    out_chroms2.append(row.chrom2)
                    out_starts2.append(range2[0])
                    out_stops2.append(range2[1])
                    if svt != in_svt and in_svt != 'TRANS':
                        in_svt = 'UNK'
                else:
                    out_quals.append(row.qual)
                    out_chroms1.append(row.chrom1)
                    out_starts1.append(row.start1)
                    out_stops1.append(row.stop1)
                    out_chroms2.append(row.chrom2)
                    out_starts2.append(row.start2)
                    out_stops2.append(row.stop2)
                    
                out_infos.append(tk_sv_io.update_info(row['info'], ['NPAIRS', 'NSPLIT', 'RP_LR', 'RP_TYPE', 'TYPE'],
                                                      [num_pairs, num_split, lr, svt, in_svt]))

    in_bam.close()

    if args.best_only:
        out_names = [n for n in bedpe_df['name']]
    else:
        out_names = np.arange(len(bedpe_df))

    out_df = tk_sv_io.create_sv_df(out_chroms1, out_starts1, out_stops1,
                                   out_chroms2, out_starts2, out_stops2,
                                   out_names, out_quals, out_infos)
    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.sv_variants)


def join(args, outs, chunk_defs, chunk_outs):
    join_df = None
    for chunk in chunk_outs:
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(chunk.sv_variants)
        join_df = pd.concat([join_df, bedpe_df], ignore_index = True)

    if not args.best_only:
        join_df['name'] = np.arange(len(join_df))
        
    tk_sv_io.write_sv_df_to_bedpe(join_df, outs.sv_variants)


def filter_fun(in_bam, chrom1, chrom2, r1, r2):
    return lambda rp: (in_bam.getrname(rp.chrom1) == chrom1 and \
                       in_bam.getrname(rp.chrom2) == chrom2 and \
                       tk_readpairs.pos_overlaps(rp.read1.pos, r1) and \
                       tk_readpairs.pos_overlaps(rp.read2.pos, r2)) or \
                      (in_bam.getrname(rp.chrom1) == chrom2 and \
                       in_bam.getrname(rp.chrom2) == chrom1 and \
                       tk_readpairs.pos_overlaps(rp.read1.pos, r2) and \
                       tk_readpairs.pos_overlaps(rp.read2.pos, r1))
    
