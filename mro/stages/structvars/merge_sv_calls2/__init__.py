#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import tenkit.pandas as pd
import longranger.sv.io as tk_sv_io
import longranger.sv.utils as tk_sv_utils

__MRO__ = """
stage MERGE_SV_CALLS2(
    in  bedpe  sv_variants,
    in  bedpe  cnv_variants,
    in  int    max_dist,
    in  json   sv_summary,
    out bedpe  sv_variants,
    src py     "stages/structvars/merge_sv_calls2",
)
"""


def main(args, outs):

    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_variants)
    sv_df["info2"] = "SV"

    cnv_df = tk_sv_io.read_sv_bedpe_to_df(args.cnv_variants)
    cnv_df["info2"] = "CNV"

    sv_df = pd.concat([sv_df, cnv_df], ignore_index = True)
    sv_df['name'] = np.arange(len(sv_df))
    sv_df.sort(['chrom1', 'chrom2'], inplace = True)

    res_df = None
    for _, tmp_df in sv_df.groupby(['chrom1', 'chrom2']):
        tmp_df.sort(['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2'], inplace = True)
        # cluster the loci in the group based on proximity
        groups = tk_sv_utils.get_break_groups(tmp_df, args.max_dist)

        # for each cluster, get the row with max qual
        # tmp_df.loc[g] gets the subset of tmp_df in the cluster.
        # then idxmax gets the max index

        out_df = pd.DataFrame(columns=sv_df.columns)
        idx = 0
        for g in groups:
            row = tmp_df.loc[tmp_df.loc[g]['qual'].idxmax()]
            if (tmp_df.loc[g]['info2'] == 'SV').any():
                row = tmp_df.loc[(tmp_df.loc[g]['info2'] == 'SV').idxmax()]

            source = list(set(tmp_df.loc[g]['info2']))
            row['info'] += (";SOURCE=" + ",".join(source))
            out_df.loc[idx] = row
            idx += 1

        out_df.sort(['start1', 'stop1', 'start2', 'stop2'], inplace = True)
        res_df = pd.concat([res_df, out_df], ignore_index = True)

    tk_sv_io.write_sv_df_to_bedpe(res_df, outs.sv_variants)
