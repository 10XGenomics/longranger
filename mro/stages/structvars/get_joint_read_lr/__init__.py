#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import sys
import json
import numpy as np
import tenkit.pandas as pd
from itertools import groupby
import tenkit.bam as tk_bam
import longranger.sv.io as tk_sv_io
import longranger.sv.utils as tk_sv_utils
import longranger.sv.sv_call as sv_call
import longranger.sv.readpairs as tk_readpairs
import tenkit.tabix as tk_tabix
from longranger.sv.constants import MAX_DEL_READPAIRS

__MRO__ = """
stage GET_JOINT_READ_LR(
    in  bam     possorted_bam,
    in  json    basic_summary,
    in  bedpe   sv_calls,
    in  bedpe   sv_calls2,
    in  tsv.gz  fragment_phasing,
    in  json    insert_sizes,
    in  int     min_mapq,
    in  int     break_pad,
    in  float   min_lr,
    in  int     min_sv_len,
    in  int     em_iters,
    out bedpe   sv_calls,
    out bedpe   non_pass_sv_calls,
    src py     "stages/structvars/get_joint_read_lr",
) split using (
    in  int chunk_start,
    in  int chunk_end,
)
"""


MAX_DEL_BC_DEPTH = 1000

def read_bedpes(args):
    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.sv_calls)
    if not args.sv_calls2 is None:
        sv_df = pd.concat([sv_df, tk_sv_io.read_sv_bedpe_to_df(args.sv_calls2)],
                          ignore_index=True)
        sv_df['name'] = np.arange(len(sv_df))
    return sv_df


def split(args):
    sv_df = read_bedpes(args)
    nsvs = sv_df.shape[0]
    nsvs_per_chunk = max(100, int(np.ceil(nsvs / 200.0)))
    nchunks = int(np.ceil(nsvs / float(nsvs_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nsvs_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nsvs_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':9})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':16}]

    return {'chunks': chunk_defs}


def get_frag_coverage(frag_phasing, chrom, start, stop):
    bcs = set([])
    if not chrom in frag_phasing.contigs:
        return bcs

    for frag_line in frag_phasing.fetch(str(chrom), start, stop):
        frag = frag_line.strip().split('\t')
        bc = frag[6]
        bcs.add(bc)
    return bcs


def main(args, outs):
    sv_df = read_bedpes(args)
    sv_df = tk_sv_utils.get_dataframe_loc(sv_df, list(range(int(args.start_idx), int(args.stop_idx))))

    max_insert, ins_logsf_fun = tk_sv_utils.get_insert_size_info(args.insert_sizes)
    print >> sys.stderr, 'max insert', max_insert

    if max_insert is None:
        tk_sv_io.write_sv_df_to_bedpe(None, outs.sv_calls)
        tk_sv_io.write_sv_df_to_bedpe(sv_df, outs.non_pass_sv_calls)
        return

    with open(args.basic_summary, 'r') as f:
        summary = json.load(f)

    chimera_rate_del = summary['far_chimera_rate']
    chimera_rate_inv = summary['far_chimera_rate'] + summary['same_dir_chimera_rate']
    chimera_rate_dup = summary['far_chimera_rate'] + summary['outward_dir_chimera_rate']
    chimera_rates = {tk_readpairs.DEL_STR:chimera_rate_del,
                     tk_readpairs.INV_STR:chimera_rate_inv,
                     tk_readpairs.TDUP_STR:chimera_rate_dup,
                     tk_readpairs.TRANS_STR:summary['far_chimera_rate']}

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    frag_phasing = tk_tabix.create_tabix_infile(args.fragment_phasing)

    pass_calls = []
    non_pass_calls = []

    for i, (_, row) in enumerate(sv_df.iterrows()):
        sv_type = tk_sv_io.get_sv_type(row.info)

        middle = int(0.5 * (row.stop1 + row.start2))

        # Bail out on all non deletions
        if sv_type != tk_readpairs.DEL_STR:
            continue

        if row.chrom1 == row.chrom2:
            r1 = (max(0, row.start1 - args.break_pad), min(middle, row.stop1 + args.break_pad))
            r2 = (max(middle, row.start2 - args.break_pad), row.stop2 + args.break_pad)

            if row.start2 - row.stop1 > 4 * args.break_pad:
                starts = [r1[0], r2[0]]
                stops = [r1[1], r2[1]]
                chroms = [row.chrom1, row.chrom2]
            else:
                starts = [r1[0]]
                stops = [r2[1]]
                chroms = [row.chrom1]
        else:
            r1 = (max(0, row.start1 - args.break_pad), row.stop1 + args.break_pad)
            r2 = (max(0, row.start2 - args.break_pad), row.stop2 + args.break_pad)
            starts = [r1[0], r2[0]]
            stops = [r1[1], r2[1]]
            chroms = [row.chrom1, row.chrom2]

        bc_cov1 = len(get_frag_coverage(frag_phasing, row.chrom1, r1[0], r1[1]))
        bc_cov2 = len(get_frag_coverage(frag_phasing, row.chrom2, r2[0], r2[1]))
        if sv_type == tk_readpairs.DEL_STR and max(bc_cov1, bc_cov2) > MAX_DEL_BC_DEPTH:
            print >> sys.stderr, 'Too many barcodes in DEL candidate', row.chrom1, row.start1, row.stop2
            continue

        readpairs = tk_readpairs.get_readpairs(in_bam, chroms, starts, stops,
                                               max_insert=max_insert, min_mapq=args.min_mapq)

        normal_readpairs = [rp for rp in readpairs if rp.sv_type == tk_readpairs.NORMAL_STR]
        if len(normal_readpairs) > MAX_DEL_READPAIRS:
            sel = np.random.choice(len(normal_readpairs), MAX_DEL_READPAIRS)
        else:
            sel = np.arange(len(normal_readpairs))
        normal_readpairs = [normal_readpairs[ridx] for ridx in sel]

        # Distal readpairs across the breakpoints
        dist_readpairs = [rp for rp in readpairs if rp.sv_type == sv_type and ((tk_readpairs.pos_overlaps(rp.read1.pos, r1) and tk_readpairs.pos_overlaps(rp.read2.pos, r2)) or
                                                                               (tk_readpairs.pos_overlaps(rp.read1.pos, r2) and tk_readpairs.pos_overlaps(rp.read2.pos, r1)))]
        if len(dist_readpairs) > MAX_DEL_READPAIRS:
            sel = np.random.choice(len(dist_readpairs), MAX_DEL_READPAIRS)
        else:
            sel = np.arange(len(dist_readpairs))
        dist_readpairs = [dist_readpairs[ridx] for ridx in sel]

        dist_readpairs.extend(normal_readpairs)
        if sv_type == tk_readpairs.DEL_STR and len(starts) == 2:
            more_readpairs = tk_readpairs.get_readpairs(in_bam, [row.chrom1], [r1[1] + 1], [r2[0] - 1],
                                                        max_insert=max_insert,
                                                        min_mapq=args.min_mapq,
                                                        normal_only=True)
            if len(more_readpairs) > MAX_DEL_READPAIRS:
                sel = np.random.choice(len(more_readpairs), MAX_DEL_READPAIRS)
            else:
                sel = np.arange(len(more_readpairs))
            dist_readpairs.extend([more_readpairs[ridx] for ridx in sel if more_readpairs[ridx].sv_type == tk_readpairs.NORMAL_STR])

        readpairs = sorted(dist_readpairs, key = lambda x: x.barcode)
        read_groups = {}
        for bc, read_group_iter in groupby(dist_readpairs, lambda x: x.barcode):
            read_groups[bc] = list(read_group_iter)

        bc_set = set(read_groups.keys())
        bc_list = sorted(read_groups.keys())
        phase_set1 = tk_sv_utils.get_phase_set(frag_phasing, row.chrom1, r1[0], r1[1])
        phase_set2 = tk_sv_utils.get_phase_set(frag_phasing, row.chrom2, r2[0], r2[1])

        if len(bc_list) < 1:
            print >> sys.stderr, 'Not enough barcodes. Skipping'
            continue

        bc_phase_sets1 = tk_sv_utils.get_barcode_phase_probs(frag_phasing, row.chrom1, r1[0], r1[1], bc_set, in_ps = phase_set1)
        bc_phase_sets2 = tk_sv_utils.get_barcode_phase_probs(frag_phasing, row.chrom2, r2[0], r2[1], bc_set, in_ps = phase_set2)

        cand_breaks1 = np.arange(r1[0], r1[1] + 1, 5)
        cand_breaks2 = np.arange(r2[0], r2[1] + 1, 5)

        res = tk_readpairs.eval_sv_em(read_groups, cand_breaks1, cand_breaks2,
                                      sv_type, chimera_rates, phase_set1, phase_set2,
                                      bc_phase_sets1, bc_phase_sets2, max_insert, ins_logsf_fun,
                                      em_iters = args.em_iters)

        ((no_sv_max, sv_max, het_sv_max), max_locus, zygosity, max_hap, prior_hap_probs, hap_probs, support) = res

        lr = sv_max - no_sv_max if max_hap is None else het_sv_max - no_sv_max

        hap_probs1 = hap_probs[:, 0:2]
        hap_probs2 = hap_probs[:, 2:]

        new_call = sv_call.SvCall.from_em_results(row.chrom1, row.chrom2, phase_set1, phase_set2,
                                                  (no_sv_max, sv_max, het_sv_max), max_locus,
                                                  sv_call._SvType(sv_type, ('.', '.')), zygosity,
                                                  max_hap, support, (hap_probs1, hap_probs2, None))

        # the break interval is inclusive
        if lr >= args.min_lr and new_call.qual >= args.min_qv and new_call.break2[0] - new_call.break1[1] + 1 >= args.min_sv_len:
            pass_calls.append(new_call)
        else:
            # Leave breakpoints unchanged
            new_call.break1 = (row.start1, row.stop1)
            new_call.break2 = (row.start2, row.stop2)
            non_pass_calls.append(new_call)

    out_df = sv_call.SvCall.svs_to_dataframe(pass_calls)
    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.sv_calls)

    out_df = sv_call.SvCall.svs_to_dataframe(non_pass_calls)
    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.non_pass_sv_calls)
    in_bam.close()
    frag_phasing.close()


def join(args, outs, chunk_defs, chunk_outs):
    join_df = None
    non_pass_join_df = None
    for chunk in chunk_outs:
        df = tk_sv_io.read_sv_bedpe_to_df(chunk.sv_calls)
        non_pass_df = tk_sv_io.read_sv_bedpe_to_df(chunk.non_pass_sv_calls)
        join_df = pd.concat([join_df, df], ignore_index = True)
        non_pass_join_df = pd.concat([non_pass_join_df, non_pass_df], ignore_index = True)

    join_df['name'] = np.arange(len(join_df))
    tk_sv_io.write_sv_df_to_bedpe(join_df, outs.sv_calls)
    non_pass_join_df['name'] = np.arange(len(join_df), len(join_df) + len(non_pass_join_df))
    tk_sv_io.write_sv_df_to_bedpe(non_pass_join_df, outs.non_pass_sv_calls)
