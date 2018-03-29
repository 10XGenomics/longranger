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
import tenkit.hdf5 as tk_hdf5
import longranger.sv.utils as tk_sv_utils
from longranger.sv.constants import MAX_FRAG_SIZE

__MRO__ = """
stage FILTER_PILEUPS(
    in  bam    possorted_bam,
    in  bedpe  sv_calls,
    in  h5     hap_coverage,
    in  float  min_rel_depth,
    in  float  max_clipped_frac,
    out bedpe  sv_calls,
    out bedpe  pileups,
    src py     "stages/structvars/filter_pileups",
) split using (
    in int start_idx,
    in int stop_idx,
)
"""

# Look at coverage BREAK_EXT before and after the breakpoints to get a
# normalization factor.
BREAK_EXT = 1000
BIN_WIN = 1000
MAX_CLIPPED = 10
DEL_REL_COV = 0.8
DUP_REL_COV = 1.2


def has_too_many_clipped(in_bam, chrom, start, stop, min_mapq = 30, max_clipped_frac = 0.1):
    clip_starts = np.zeros((stop - start, ), dtype = np.int)

    for read in in_bam.fetch(str(chrom), int(start), int(stop)):
        if not read.is_unmapped and not read.is_duplicate and read.mapq >= min_mapq:
            align_start = read.pos - start
            align_end = read.aend - start
            cigar_tuples = read.cigar
            # Cigar starts with soft or hard clipping
            if cigar_tuples[0][0] in [4,5] and align_start >= 0:
                clip_starts[align_start] += 1
            # Cigar ends with soft or hard clipping
            if cigar_tuples[-1][0] in [4,5] and align_end < len(clip_starts):
                clip_starts[align_end] += 1

    return np.mean(clip_starts > max(MAX_CLIPPED, np.median(clip_starts))) > max_clipped_frac


def split(args):
    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_calls)

    nsvs = sv_df.shape[0]
    nbreaks_per_chunk = max(100, int(np.ceil(nsvs / 32.0))) # avoid overchunking
    nchunks = int(np.ceil(nsvs / float(nbreaks_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nbreaks_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nbreaks_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':12})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':12}]

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 12}}


def main(args, outs):
    pred_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_calls)
    pred_df = tk_sv_utils.get_dataframe_loc(pred_df, list(range(args.start_idx, args.stop_idx)))

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    cov_reader = tk_hdf5.DataFrameReader(args.hap_coverage)
    sel_cols = ['cov_q30_hap0', 'cov_q30_hap1', 'cov_q30_hap2']

    has_pileups = np.zeros((len(pred_df), ), dtype = np.bool)

    for i, (_, row) in enumerate(pred_df.iterrows()):
        has_clipped1 = has_too_many_clipped(in_bam, row.chrom1, max(0, row.start1 - BREAK_EXT),
                                            row.stop1 + BREAK_EXT,
                                            max_clipped_frac = args.max_clipped_frac)
        has_clipped2 = has_too_many_clipped(in_bam, row.chrom2, max(0, row.start2 - BREAK_EXT),
                                            row.stop2 + BREAK_EXT,
                                            max_clipped_frac = args.max_clipped_frac)
        has_clipped = has_clipped1 and has_clipped2

        if row.chrom1 != row.chrom2 or row.start2 - row.stop1 > MAX_FRAG_SIZE:
            has_pileups[i] = has_clipped
            continue

        cov = cov_reader.query((row.chrom1, max(0, row.start1 - BREAK_EXT),
                                row.stop2 + BREAK_EXT))
        cov['bin'] = np.array(cov['pos'] / BIN_WIN, dtype = np.int)
        if not 'coverage_deduped' in cov.columns:
            cov['coverage_deduped'] = cov[sel_cols].sum(axis = 1)
        cov_arr = np.array(cov.groupby('bin').mean()['coverage_deduped'])
        median_cov = np.median(cov_arr)

        # Rescue for deletions or duplications with breakpoints on the pileups
        sv_len = row.stop2 - row.start1
        side_cov = cov_reader.query((row.chrom1, max(0, row.start1 - BREAK_EXT - sv_len / 2),
                                     row.start1 - BREAK_EXT))
        side_cov = pd.concat([side_cov,
                              cov_reader.query((row.chrom2,
                                                row.stop2 + BREAK_EXT,
                                                row.stop2 + BREAK_EXT + sv_len / 2))],
                             ignore_index = True)
        if not 'coverage_deduped' in side_cov.columns:
            side_cov['coverage_deduped'] = side_cov[sel_cols].sum(axis = 1)

        # Ignore pileups, enough evidence for a large-scale copy number variant
        if np.median(cov.coverage_deduped) < DEL_REL_COV * np.median(side_cov.coverage_deduped):
            continue
        if np.median(cov.coverage_deduped) > DUP_REL_COV * np.median(side_cov.coverage_deduped):
            continue

        # Filter out the call if there are pileups very close to the breakpoints
        has_pileups[i] = len(cov_arr) > 4 and np.any(cov_arr[[0, 1, -2, -1]] > args.min_rel_depth * median_cov)
        has_pileups[i] = has_pileups[i] or has_clipped

    pileups = pred_df[has_pileups]
    pred_df = pred_df[np.logical_not(has_pileups)]

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
        out_calls = pd.concat([out_calls, calls], ignore_index = True)
        out_pileups = pd.concat([out_pileups, pileups], ignore_index = True)

    tk_sv_io.write_sv_df_to_bedpe(out_calls, outs.sv_calls)
    tk_sv_io.write_sv_df_to_bedpe(out_pileups, outs.pileups)
