#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import tenkit.pandas as pd
import longranger.sv.io as tk_sv_io
import longranger.sv.utils as tk_sv_utils
from longranger.sv.constants import MAX_FRAG_SIZE

__MRO__ = """
stage MERGE_SV_CALLS(
    in  bedpe  sv_variants1,
    in  bedpe  sv_variants2,
    in  int    min_overlap,
    in  int    max_dist,
    in  float  min_frac_overlap,
    in  bool   merge_dist,
    out bedpe  sv_variants,
    src py     "stages/structvars/merge_sv_calls",
)
"""

def merge_overlapping(sv_df, min_frac_over):
    last_start = None
    last_stop = None
    last_row_idx = 0
    last_row = None
    final_rows = []
    sv_types = set([])
    infos = [s for s in sv_df['info']]

    for i, (_, row) in enumerate(sv_df.iterrows()):
        if last_row is None:
            last_start = row.start1
            last_stop = row.stop2
            last_row = row
            last_row_idx = i
        else:
            ov = float(min(last_stop, row.stop2) - row.start1)
            frac_ov = min(ov / (last_stop - last_start), ov / (row.stop2 - row.start1))

            if frac_ov < min_frac_over:
                final_rows.append(last_row_idx)
                if len(sv_types) > 1:
                    infos[last_row_idx] = tk_sv_io.update_info(infos[last_row_idx],
                                                               ['TYPE'], ['UNK'])
                sv_types = set([])
                last_start = row.start1
                last_stop = row.stop2
                last_row = row
                last_row_idx = i
            else:
                if row.qual > last_row.qual:
                    last_row_idx = i
                    last_row = row
                    last_stop = row.stop2
        sv_types.add(tk_sv_io.get_sv_type(row.info))

    if not last_row is None:
        final_rows.append(last_row_idx)
        if len(sv_types) > 1:
            infos[last_row_idx] = tk_sv_io.update_info(infos[last_row_idx], ['TYPE'], ['UNK'])

    sv_df['info'] = infos
    return tk_sv_utils.get_dataframe_loc(sv_df, final_rows)


def main(args, outs):

    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_variants1)

    if not args.sv_variants2 is None:
        sv_df = pd.concat([sv_df, tk_sv_io.read_sv_bedpe_to_df(args.sv_variants2)], ignore_index = True)
    sv_df['name'] = np.arange(len(sv_df))

    if args.merge_dist:
        is_prox = np.ones((len(sv_df),), dtype=np.bool)
    else:
        is_prox = np.logical_and(sv_df.chrom1 == sv_df.chrom2, sv_df.start2 - sv_df.stop1 < MAX_FRAG_SIZE)
    prox_df = sv_df[is_prox]
    dist_df = sv_df[np.logical_not(is_prox)]

    prox_df.sort(['chrom1', 'chrom2'], inplace = True)

    res_df = None
    for _, tmp_df in prox_df.groupby(['chrom1', 'chrom2']):
        tmp_df.sort(['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2'], inplace = True)
        # cluster the loci in the group based on proximity
        groups = tk_sv_utils.get_break_groups(tmp_df, args.max_dist)

        if args.min_bc_overlap is None:
            # for each cluster, get the row with max qual
            # tmp_df.loc[g] gets the subset of tmp_df in the cluster.
            # then idxmax gets the max index
            out_df = tmp_df.loc[[tmp_df.loc[g]['qual'].idxmax() for g in groups]]
            out_df.sort(['start1', 'stop1', 'start2', 'stop2'], inplace = True)
        else:
            sel_idx = []
            for group in groups:
                sel_idx.extend(merge_bcs(tmp_df.loc[group], args.min_bc_overlap))
            out_df = tmp_df.loc[sel_idx]
            out_df = tk_sv_io.drop_info(out_df, ['BCS'])
            
        if not args.min_frac_overlap is None:
            out_df = merge_overlapping(out_df, args.min_frac_overlap)

        res_df = pd.concat([res_df, out_df], ignore_index = True)

    res_df = pd.concat([dist_df, res_df], ignore_index = True)
    tk_sv_io.write_sv_df_to_bedpe(res_df, outs.sv_variants)


def merge_bcs(call_df, min_bc_overlap):
    bc_clusters = {}
    for idx, row in call_df.iterrows():
        bcs = frozenset(tk_sv_io.extract_sv_info(row.info, ['BCS'])[0].split(','))
        match = None
        for cluster_id, (cluster, _, _) in bc_clusters.iteritems():
            print cluster, bcs
            if bcs_overlap(cluster, bcs, min_bc_overlap):
                match = cluster_id
                break
        if match is None:
            bc_clusters[len(bc_clusters)] = (bcs, idx, row)
        elif bc_clusters[match][2].qual < row.qual:
            bc_clusters[match] = (bcs, idx, row)         
    
    return [idx for (_, idx, _) in bc_clusters.values()]


def bcs_overlap(bcs1, bcs2, min_bc_overlap):
    if min_bc_overlap >= 1:
        return len(bcs1.intersection(bcs2)) >= min_bc_overlap
    return len(bcs1.intersection(bcs2)) / float(len(bcs1.union(bcs2))) >= min_bc_overlap
