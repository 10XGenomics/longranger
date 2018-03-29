# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.pandas as pd
import longranger.sv.io as sv_io
import numpy as np
from numpy.random import choice

def has_inter_chromosomal(callset):
    return np.any(callset.chrom1 != callset.chrom2)

def get_priorities(callset):
    priorities = np.array([int(sv_io.extract_sv_info(info, ['SOURCE'])[0]) for info in callset['info']])
    callset['priority'] = priorities
    return callset

def group_overlapping(callset):
    """Group overlapping calls and return the groups.
    """
    groups = np.zeros((len(callset), ), dtype=np.int)
    last_chrom = None
    last_end = None
    group = 0

    for idx, (_, row) in enumerate(callset.iterrows()):
        if row is None or row.chrom1 != last_chrom or row.start1 > last_end:
            group += 1
            last_chrom = row.chrom1
            last_end = row.start2
        else:
            last_end = max(row.start2, last_end)
        groups[idx] = group

    callset['group'] = groups
    return callset.groupby('group')

def select_first():
    return lambda x:x.iloc[0]

def select_widest():
    # np.argmax returns the index in the dataframe NOT a positional index
    return lambda x: x.loc[np.argmax(x.start2 - x.start1)]

def select_highest_qual():
    return lambda x: x.loc[np.argmax(x.qual)]

def select_random():
    return lambda x: x.loc[choice(x.index, 1)]

def merge_overlapping(callsets, selection_fun=select_first()):
    """Merge overlapping calls and remove redundancies based on the specified function.

    - callsets: list of call dataframes
    - selection_fun: given a group of overlapping calls, this function decides
    which of these calls should be output. Default is to pick the first in each group.
    """

    # Make sure there are no inter-chromosomal calls
    assert np.all(not has_inter_chromosomal(callset) for callset in callsets)
    callset = pd.concat(callsets, ignore_index=True)
    # The selection functions might break if this is empty
    if callset.empty:
        return None
    callset.sort(['chrom1', 'start1', 'start2'], inplace=True)
    groups = group_overlapping(callset)

    # agg returns a multilevel dataframe. reset_index removes one level and
    # adds a group column, which we drop.
    return groups.agg(selection_fun).reset_index().drop('group', axis=1)

def main(args, outs):
    callsets = []

    if args.calls1 is not None:
        c = sv_io.read_sv_bedpe_to_df(args.calls1)
        callsets.append(c)

    if args.calls2 is not None:
        c = sv_io.read_sv_bedpe_to_df(args.calls2)
        callsets.append(c)

    if args.calls3 is not None:
        c = sv_io.read_sv_bedpe_to_df(args.calls3)
        callsets.append(c)

    # Select the highest qual
    merged = merge_overlapping(callsets, select_widest())
    sv_io.write_sv_df_to_bedpe(merged, outs.merged)
