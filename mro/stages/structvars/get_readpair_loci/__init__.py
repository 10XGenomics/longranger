#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Get candidate loci for short SVs by grouping together discordant reads

import cPickle
import os.path
import numpy as np
import tenkit.pandas as pd
import martian
from collections import defaultdict
import json
import tenkit.safe_json
import longranger.sv.io as tk_sv_io
import longranger.sv.utils as tk_sv_utils
import longranger.sv.readpairs as tk_readpairs
from longranger.sv.constants import MAX_INSERT_SIZE_PRC, MIN_SV_INSERT_SIZE_PRC
__MRO__ = """
stage GET_READPAIR_LOCI(
    in  bam    possorted_bam,
    in  string barcode_whitelist,
    in  json   basic_summary,
    in  int    min_mapq,
    in  json   insert_sizes,
    in  float  merge_range_factor,
    in  int    min_reads_to_call,
    in  int    min_lr_to_call,
    in  int    min_sv_len,
    in  int    max_sv_len,
    in  pickle ranges,
    out bedpe  sv_calls,
    out json   discordant_read_counts,
    src py    "stages/structvars/get_readpair_loci",
) split using (
    in  string chrom,
    in  int[]  starts,
    in  int[]  stops,
)
"""

# Maximum number of basepairs considered in a chunk.
MAX_BP_PER_CHUNK = 40000000

def split(args):
    with open(args.ranges, 'r') as f:
        ranges = cPickle.load(f)

    chunks = []
    for c, r in ranges.iteritems():
        starts = []
        stops = []
        tot_size = 0
        for (s, e) in r:
            if tot_size > MAX_BP_PER_CHUNK:
                chunks.append({'chrom':c, 'starts':starts, 'stops':stops, '__mem_gb':3})
                starts = []
                stops = []
                tot_size = 0
            starts.append(s)
            stops.append(e)
            tot_size += e - s
        chunks.append({'chrom':c, 'starts':starts, 'stops':stops, '__mem_gb':3})
    # Can't have an empty list of chunks
    if len(chunks) == 0:
        chunks.append({'chrom':None, 'starts':[], 'stops':[]})
    return {'chunks': chunks, 'join': {'__mem_gb': 3}}


def main(args, outs):
    if args.chrom is None or len(args.starts) == 0 or args.barcode_whitelist is None:
        tk_sv_io.write_sv_df_to_bedpe(None, outs.sv_calls)
        return

    max_insert, ins_logsf_fun = tk_sv_utils.get_insert_size_info(args.insert_sizes, MAX_INSERT_SIZE_PRC)
    if max_insert is None:
        martian.throw('No Q60 reads')

    # This is slightly bigger than the maximum "normal" insert
    min_call_insert, _ = tk_sv_utils.get_insert_size_info(args.insert_sizes, MIN_SV_INSERT_SIZE_PRC)
    min_sv_len = max(args.min_sv_len, min_call_insert)
    martian.log_info('Setting min_sv_len to {}'.format(min_sv_len))
    
    with open(args.basic_summary, 'r') as f:
        summary = json.load(f)

    chimera_rate_del = summary['far_chimera_rate']
    chimera_rate_inv = summary['far_chimera_rate'] + summary['same_dir_chimera_rate']
    chimera_rate_dup = summary['far_chimera_rate'] + summary['outward_dir_chimera_rate']
    chimera_rates = {tk_readpairs.DEL_STR:chimera_rate_del,
                     tk_readpairs.INV_STR:chimera_rate_inv,
                     tk_readpairs.TDUP_STR:chimera_rate_dup,
                     tk_readpairs.TRANS_STR:summary['far_chimera_rate']}

    df, read_counts, _ = tk_readpairs.get_discordant_loci(args.possorted_bam, chrom = str(args.chrom),
                                                          starts = args.starts, stops = args.stops,
                                                          min_mapq = args.min_mapq, min_insert = 0,
                                                          max_insert = max_insert,
                                                          max_merge_range = args.merge_range_factor * max_insert,
                                                          min_sv_len = min_sv_len, max_sv_len = args.max_sv_len,
                                                          ins_logsf_fun = ins_logsf_fun,
                                                          min_lr_to_call = args.min_lr_to_call,
                                                          min_reads_to_call = args.min_reads_to_call,
                                                          chimera_rate = chimera_rates, reads_as_qual = True)

    # Need to convert to dict because defaultdict doesn't get pickled properly
    read_counts['split'] = dict(read_counts['split'])
    read_counts['pair'] = dict(read_counts['pair'])
    tk_sv_io.write_sv_df_to_bedpe(df, outs.sv_calls)
    with open(outs.discordant_read_counts, 'w') as f:
        f.write(tenkit.safe_json.safe_jsonify(read_counts))


def join(args, outs, chunk_defs, chunk_outs):
    join_df = None
    read_counts = {}
    read_counts['split'] = defaultdict(int)
    read_counts['pair'] = defaultdict(int)

    for chunk in chunk_outs:
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(chunk.sv_calls)
        join_df = pd.concat([join_df, bedpe_df], ignore_index = True)

        if not os.path.isfile(chunk.discordant_read_counts):
            continue
        with open(chunk.discordant_read_counts, 'r') as f:
            counts = json.load(f)
        for t, c in counts['split'].iteritems():
            read_counts['split'][t] += c
        for t, c in counts['pair'].iteritems():
            read_counts['pair'][t] += c

    join_df['name'] = [str(i) for i in np.arange(len(join_df))]
    tk_sv_io.write_sv_df_to_bedpe(join_df, outs.sv_calls)

    read_counts['split'] = dict(read_counts['split'])
    read_counts['pair'] = dict(read_counts['pair'])

    with open(args.basic_summary, 'r') as f:
        num_reads = float(json.load(f)['num_reads']) / 2.0

    read_counts['frac_split'] = {}
    read_counts['frac_pair'] = {}
    for t, c in read_counts['split'].iteritems():
        read_counts['frac_split'][t] = c / num_reads
    for t, c in read_counts['pair'].iteritems():
        read_counts['frac_pair'][t] = c / num_reads

    with open(outs.discordant_read_counts, 'w') as f:
        f.write(tenkit.safe_json.safe_jsonify(read_counts))
