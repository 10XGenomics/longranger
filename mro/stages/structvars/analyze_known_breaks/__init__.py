#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import os
import os.path
import json
import numpy as np
import tenkit.pandas as pd

import martian

import tenkit.bam as tk_bam
from longranger.sv.stats import *
from longranger.sv.utils import *
from longranger.sv.io import *
import tenkit.hdf5 as tk_hdf5
from longranger.sv.constants import MIN_FRAG_SIZE_WGS, MIN_READS_PER_FRAG_WGS, MIN_FRAG_SIZE_TARGET, MIN_READS_PER_FRAG_TARGET
from longranger.sv.constants import SV_FRAGMENT_LINK_DISTANCE_WGS, SV_FRAGMENT_LINK_DISTANCE_TARGET
from longranger.sv.constants import MAX_FRAG_SIZE, TARGET_COV_BIN

__MRO__ = """
stage ANALYZE_KNOWN_BREAKS(
    in  bam   input,
    in  bed   targets,
    in  int   nx,
    in  int   extend_win,
    in  bedpe gt_variants,
    in  h5    fragments,
    in  json  fragment_histogram,
    in  h5    barcodes,
    in  tsv   barcode_blacklist,
    in  json  coverage,
    out bedpe summary,
    src py    "stages/structvars/analyze_known_breaks",
) split using (
    in  int   start_idx,
    in  int   stop_idx,
)
"""

def split(args):
    if args.gt_variants is None:
        chunk_defs = [{'start_idx':0, 'stop_idx':0}]
        return {'chunks': chunk_defs}

    sv_df = read_sv_bedpe_to_df(args.gt_variants)
    check_sv_names(sv_df)
    nsvs = sv_df.shape[0]
    nbreaks_per_chunk = max(20, int(np.ceil(nsvs / 32.0))) # avoid overchunking
    nchunks = int(np.ceil(nsvs / float(nbreaks_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nbreaks_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nbreaks_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':8})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':8}]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    out_bedpe = None
    for c in chunk_outs:
        in_bedpe = read_sv_bedpe_to_df(c.summary)
        out_bedpe = pd.concat([out_bedpe, in_bedpe], ignore_index = True)
    write_sv_df_to_bedpe(out_bedpe, outs.summary)


def isfile(fn):
    return not fn is None and os.path.isfile(fn)

def main(args, outs):
    sv_df = read_sv_bedpe_to_df(args.gt_variants)
    sv_df = get_dataframe_loc(sv_df, list(range(args.start_idx, args.stop_idx)))

    if not isfile(args.fragments) or not isfile(args.fragment_histogram) or not isfile(args.barcodes) \
            or not isfile(args.barcode_blacklist) or not isfile(args.coverage):
        sv_df['qual'] = 0
        if 'info' in sv_df.columns:
            info_strs = [s for s in sv_df['info']]
        else:
            info_strs = ['.' for i in range(len(sv_df))]
        for i in range(len(info_strs)):
            info_strs[i] = update_info(info_strs[i], ['BCOV', 'NBCS1', 'NBCS2', 'NOOV'], [0, 0, 0, 0])
        sv_df['info'] = info_strs
        sv_df['strand1'] = '.'
        sv_df['strand2'] = '.'
        sv_df['filters'] = '.'
        write_sv_df_to_bedpe(sv_df, outs.summary)
        martian.log_info('One or more files needed for computing quality scores are missing.')
        return

    input_bam = tk_bam.create_bam_infile(args.input)
    genome_size = np.sum(np.array(input_bam.lengths))
    input_bam.close()

    frag_file = args.fragments
    frag_hist_file = args.fragment_histogram
    barcode_file = args.barcodes
    barcode_blacklist_file = args.barcode_blacklist

    if args.targets is None:
        target_coverage = None
        corr_factor = 1.0
        min_frag_size = MIN_FRAG_SIZE_WGS
        min_reads_per_frag = MIN_READS_PER_FRAG_WGS
        link_distance = SV_FRAGMENT_LINK_DISTANCE_WGS
    else:
        target_regions = bed_to_region_map(args.targets, merge = True)
        target_coverage = region_cum_coverage_map(target_regions, TARGET_COV_BIN)

        with open(args.coverage, 'r') as f:
            cov_sum = json.load(f)['target_info']
        if 'on_target_bases' in cov_sum:
            prob_off_target = 1 - cov_sum['on_target_bases'] / float(cov_sum['total_bases'])
        else:
            prob_off_target = 0.001

        corr_factor = off_target_amp_corr_factor(target_regions, prob_off_target, genome_size = genome_size)
        min_frag_size = MIN_FRAG_SIZE_TARGET
        min_reads_per_frag = MIN_READS_PER_FRAG_TARGET
        link_distance = SV_FRAGMENT_LINK_DISTANCE_TARGET

    frag_sizes, frag_counts, frag_prc, blacklist_barcodes, frag_filter_fun = get_frag_data(frag_hist_file, barcode_file,
        barcode_blacklist_file, nx = args.nx, min_frag_size = min_frag_size, min_reads_per_frag = min_reads_per_frag)

    min_sv_len = link_distance
    prob_store = {}

    required_cols = ['bc', 'bc_num_reads', 'bc_est_len', 'bc_mean_reads_per_fragment',
                     'chrom', 'start_pos', 'end_pos', 'num_reads', 'obs_len', 'est_len']
    fragment_reader = tk_hdf5.DataFrameReader(frag_file)

    out_df = None

    for (chrom1, chrom2), g in sv_df.groupby(['chrom1', 'chrom2']):
        query1 = (chrom1, max(0, np.min(g['start1']) - MAX_FRAG_SIZE), np.max(g['stop1']) + MAX_FRAG_SIZE)
        filt = frag_filter_fun(*query1)
        frags1 = filt(fragment_reader.query(query1, query_cols=required_cols, id_column = 'fragment_id'))

        query2 = (chrom2, max(0, np.min(g['start2']) - MAX_FRAG_SIZE), np.max(g['stop2']) + MAX_FRAG_SIZE)
        filt = frag_filter_fun(*query2)
        frags2 = filt(fragment_reader.query(query2, query_cols=required_cols, id_column = 'fragment_id'))

        lrs = np.zeros((len(g),), dtype = np.int)
        if 'info' in g.columns:
            info_strs = [s for s in g['info']]
        else:
            info_strs = ['.' for i in range(len(g))]

        for i, (s1, e1, s2, e2) in enumerate(zip(g['start1'], g['stop1'], g['start2'], g['stop2'])):
            if chrom1 == chrom2:
                middle = 0.5 * (e1 + s2)
                locus1 = (chrom1, max(0, s1 - args.extend_win), min(middle, e1 + args.extend_win))
                locus2 = (chrom2, max(middle, s2 - args.extend_win), e2 + args.extend_win)
            else:
                locus1 = (chrom1, max(0, s1 - args.extend_win), e1 + args.extend_win)
                locus2 = (chrom2, max(0, s2 - args.extend_win), e2 + args.extend_win)

            tmp_frags1 = overlap_frags(max(0, s1 - args.extend_win), e1 + args.extend_win, frags1)
            tmp_frags2 = overlap_frags(max(0, s2 - args.extend_win), e2 + args.extend_win, frags2)

            (lr, pb, nov, ndiscordant, new_start, new_stop, _) = get_lr_from_frag_df(frags1, frags2,
                locus1, locus2, prob_store, frag_sizes, frag_counts, blacklist_barcodes, target_coverage,
                corr_factor, genome_size = genome_size, min_dist = min_sv_len)
            lrs[i] = lr
            if new_start[0] is None or new_start[1] is None:
                new_start = (s1, e1)
            if new_stop[0] is None or new_stop[1] is None:
                new_stop = (s2, e2)
            info_strs[i] = update_info(info_strs[i], ['BCOV', 'NBCS1', 'NBCS2', 'NOOV', 'P_SV', 'START', 'STOP'],
                [nov, len(set(tmp_frags1.bc)), len(set(tmp_frags2.bc)), ndiscordant, pb,
                '-'.join([str(p) for p in new_start]), '-'.join([str(p) for p in new_stop])])

        g_cp = g.copy()
        g_cp['qual'] = np.maximum(lrs, 0)
        g_cp['info'] = info_strs
        g_cp['strand1'] = '+'
        g_cp['strand2'] = '+'
        g_cp['filters'] = '.'
        out_df = pd.concat([out_df, g_cp])
    write_sv_df_to_bedpe(out_df, outs.summary)
