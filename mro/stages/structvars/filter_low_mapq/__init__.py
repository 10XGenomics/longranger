#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
#

import os.path
import numpy as np
import tenkit.bam as tk_bam
import longranger.sv.io as tk_sv_io
import tenkit.pandas as pd
import longranger.sv.utils as tk_sv_utils

__MRO__ = """
stage FILTER_LOW_MAPQ(
    in  bam    possorted_bam,
    in  bedpe  sv_calls,
    in  float  max_frac_low_mapq,
    out bedpe  sv_calls,
    out bedpe  pileups,
    src py     "stages/structvars/filter_low_mapq",
) split using (
    in int start_idx,
    in int stop_idx,
)
"""

BREAK_EXT = 1000

def mapq_changed(read):
    try:
        om = read.opt('OM')
    except KeyError:
        return False
    if om == '':
        return False
    return int(om) != read.mapq


def get_frac_mapq_changed(in_bam, chrom1, start1, stop1, chrom2, start2, stop2, min_mapq=60):
    changed1 = []
    for read in in_bam.fetch(str(chrom1), start1, stop1):
        if not read.is_unmapped and not read.is_duplicate:
            changed1.append(mapq_changed(read) and read.mapq >= min_mapq)
                
    changed2 = []
    for read in in_bam.fetch(str(chrom2), start2, stop2):
        if not read.is_unmapped and not read.is_duplicate:
            changed2.append(mapq_changed(read) and read.mapq >= min_mapq)
                
    changed1 = np.array(changed1, dtype=np.bool)
    changed2 = np.array(changed2, dtype=np.bool)
    if len(changed1) == 0 or len(changed2) == 0:
        return 0.0

    # Compute the fraction of reads that changed mapq on each breakpoint.
    return np.max([np.mean(changed1), np.mean(changed2)])


def split(args):
    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_calls)

    nsvs = sv_df.shape[0]
    nbreaks_per_chunk = max(100, int(np.ceil(nsvs / 32.0))) # avoid overchunking
    nchunks = int(np.ceil(nsvs / float(nbreaks_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nbreaks_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nbreaks_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':6})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':6}]

    return {'chunks': chunk_defs, 'join': {'__mem_gb':6}}


def main(args, outs):
    pred_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_calls)
    pred_df = tk_sv_utils.get_dataframe_loc(pred_df, list(range(args.start_idx, args.stop_idx)))

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    frac_changed = np.zeros((len(pred_df), ), dtype=np.float)
    
    for i, (_, row) in enumerate(pred_df.iterrows()):
        frac_changed[i] = get_frac_mapq_changed(in_bam, row.chrom1, max(0, row.start1 - BREAK_EXT),row.stop1 + BREAK_EXT,
                                                row.chrom2, max(0, row.start2 - BREAK_EXT),row.stop2 + BREAK_EXT,
                                                min_mapq=60)

    pileups = pred_df[frac_changed > args.max_frac_low_mapq]
    pred_df = pred_df[frac_changed <= args.max_frac_low_mapq]

    tk_sv_io.write_sv_df_to_bedpe(pred_df, outs.sv_calls)
    tk_sv_io.write_sv_df_to_bedpe(pileups, outs.pileups)


def join(args, outs, chunk_defs, chunk_outs):
    out_calls = None
    out_pileups = None
    for c in chunk_outs:
        if not os.path.isfile(c.sv_calls):
            continue
        calls = tk_sv_io.read_sv_bedpe_to_df(c.sv_calls)
        pileups = tk_sv_io.read_sv_bedpe_to_df(c.pileups)
        out_calls = pd.concat([out_calls, calls], ignore_index=True)
        out_pileups = pd.concat([out_pileups, pileups], ignore_index=True)

    tk_sv_io.write_sv_df_to_bedpe(out_calls, outs.sv_calls)
    tk_sv_io.write_sv_df_to_bedpe(out_pileups, outs.pileups)
