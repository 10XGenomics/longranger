#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Finds pairs of genomic windows with significant barcode overlaps.
#

import sys
import cPickle
import numpy as np
import martian
from itertools import groupby, combinations

from longranger.sv.overlap_detector import OverlapDetector
import tenkit.bam as tk_bam
import tenkit.hdf5 as tk_hdf5
import longranger.sv.stats as tk_sv_stats
import longranger.sv.utils as tk_sv_utils
import longranger.sv.io as tk_sv_io
from longranger.sv.constants import SV_MIN_BARCODES, MIN_FRAG_SIZE, MAX_FRAG_SIZE, MIN_SV

__MRO__ = """
stage DETECT_OVERLAPS(
    in  bam      possorted_bam,
    in  json     fragment_histogram,
    in  h5       fragments,
    in  bed      targets,
    in  pickle[] bc_pos_mats,
    in  pickle   inv_bc_map,
    in  pickle   bc_counts,
    in  pickle   win_counts,
    in  pickle   loci,
    in  int      nx,
    in  int      min_overlap,
    in  int      min_call_dist,
    in  int      max_call_dist,
    in  float    max_logp,
    in  float    max_bcs_to_call,
    in  string   test, "eg. BINOM or EXP_COUNT/0.8 where the thing before / is the name of the test and the number after is a potential parameter"
    out pickle   overlap_loci,
    src py       "stages/structvars/detect_overlaps",
) split using (
    in  pickle[] bc_mat_list1,
    in  pickle[] bc_mat_list2,
    in  pickle   inv_bc_map,
    in  pickle   win_counts,
)
"""

STEP = 100
MIN_NUM_FRAGS = 100

def split(args):
    mem_request = 18

    bc_mat_filenames = args.bc_pos_mats

    with open(args.loci, 'rb') as f:
        loci = cPickle.load(f) # loci used for chunking bc count matrices

    # Compute the maximum number of barcodes in a window in order to make a call
    # based on the percentiles of these counts genome-wide.
    if args.max_bcs_to_call <= 1.0:
        # If less than one, interpret as a percentile
        with open(args.bc_counts, 'rb') as f:
            bc_counts = cPickle.load(f)
        bc_count_list = []
        for v in bc_counts.values():
            bc_count_list.extend(v)
        max_bcs = int(np.percentile(bc_count_list, args.max_bcs_to_call * 100.0))
    else:
        # Otherwise, interpret as a static cutoff
        max_bcs = args.max_bcs_to_call
    martian.log_info('Max BCs per window ' + str(max_bcs))

    nloci = len(loci)
    chunk_defs = []
    # First create one output chunk for each input chunk
    for i in range(nloci):
        chunk_defs.append({'bc_mat_list1': [bc_mat_filenames[i]], 'bc_mat_list2': [bc_mat_filenames[i]],
            'inv_bc_map': args.inv_bc_map, 'win_counts': args.win_counts, 'max_bcs':max_bcs, '__mem_gb':mem_request})

    if nloci > 1:
        # list of all pairs of loci. Each chunk of this stage will call SVs across
        # a list of pairs of loci.
        chunk_combs = list(combinations(range(nloci), 2))
        nchunks = min(max(100, 500 - nloci), len(chunk_combs))
        npairs_per_chunk = int(np.ceil(len(chunk_combs) / float(nchunks)))

        for i in range(nchunks):
            idx_list = list(range(i * npairs_per_chunk, min((i + 1) * npairs_per_chunk, len(chunk_combs))))
            if len(idx_list) == 0:
                break
            loci_list1 = [bc_mat_filenames[chunk_combs[c][0]] for c in idx_list]
            loci_list2 = [bc_mat_filenames[chunk_combs[c][1]] for c in idx_list]
            # Each chunk will call SVs across all pairs in zip(bc_mat_list1, bc_mat_list2).
            chunk_defs.append({'bc_mat_list1': loci_list1, 'bc_mat_list2': loci_list2,
                               'inv_bc_map': args.inv_bc_map, 'win_counts': args.win_counts, 'max_bcs':max_bcs, '__mem_gb':mem_request})
    return {'chunks': chunk_defs, 'join': {'__mem_gb':mem_request}}


def join(args, outs, chunk_defs, chunk_outs):
    all_overlap_loci = []
    for c in chunk_outs:
        with open(c.overlap_loci, 'rb') as f:
            all_overlap_loci.extend(cPickle.load(f))
    with open(outs.overlap_loci, 'wb') as f:
        cPickle.dump(all_overlap_loci, f)


def main(args, outs):

    input_bam = tk_bam.create_bam_infile(args.possorted_bam)

    genome_size = np.sum(input_bam.lengths)
    input_bam.close()

    if args.fragment_histogram is None:
        martian.log_info('Missing fragment histogram. No calls will be made')
        with open(outs.overlap_loci, 'wb') as f:
            cPickle.dump([], f)
        return

    frag_res = tk_sv_stats.read_frag_hist(args.fragment_histogram,
        MIN_FRAG_SIZE, adjust = 0)
    frag_sizes, frag_counts = frag_res

    with open(args.win_counts, 'rb') as f:
        read_counts = cPickle.load(f) # number of reads per bc

    overlap_loci = []

    is_nx_bc, bc_rank = tk_sv_utils.get_nx_bcs(read_counts, args.nx)
    if np.sum(is_nx_bc) < SV_MIN_BARCODES:
        martian.log_info('Too few barcodes for SV-calling. No calls will be made.')
        with open(outs.overlap_loci, 'wb') as f:
            cPickle.dump(overlap_loci, f)
        return
    hot_bcs = np.where(is_nx_bc)[0]

    # Group pairs by the first element, to avoid loading the same thing twice.
    grouped = groupby(zip(args.bc_mat_list1, args.bc_mat_list2), lambda x:x[0])

    test = eval('tk_sv_utils.' + args.test.split('/')[0])
    if len(args.test.split('/')) > 1:
        bc_ov_adjust = float(args.test.split('/')[1])
    else:
        bc_ov_adjust = 1.0

    if test == tk_sv_utils.EXP_COUNT2:
        fragments = tk_hdf5.read_data_frame_limited(args.fragments, query_cols = ['bc', 'obs_len'], max_rows=10000000)
        if len(fragments) < MIN_NUM_FRAGS:
            martian.log_info('Too few fragments ({}). No calls will be made.'.format(len(fragments)))
            with open(outs.overlap_loci, 'wb') as f:
                cPickle.dump(overlap_loci, f)
            return
        exp_bc_ov = OverlapDetector.precompute_bc_overlaps(fragments, test, step=STEP,
                                                           max_frag_size=MAX_FRAG_SIZE,
                                                           genome_size=genome_size)

    elif test == tk_sv_utils.EXP_COUNT3:
        fragments = tk_hdf5.read_data_frame_limited(args.fragments, query_cols = ['bc', 'obs_len', 'num_reads'], max_rows=10000000)
        fragments = fragments[fragments.num_reads > 1]

        fragments.drop(['num_reads'], axis=1, inplace=True)
        if len(fragments) < MIN_NUM_FRAGS:
            martian.log_info('Too few fragments ({}). No calls will be made.'.format(len(fragments)))
            with open(outs.overlap_loci, 'wb') as f:
                cPickle.dump(overlap_loci, f)
            return
        print >> sys.stderr, 'Precomputing expected barcode overlap distribution...'
        exp_bc_ov = OverlapDetector.precompute_bc_overlaps(fragments, test, step=STEP,
                                                           max_frag_size=MAX_FRAG_SIZE,
                                                           genome_size=genome_size)
    else:
        exp_bc_ov = None

    for bc_mat_filename1, bc_mat_filename_iter2 in grouped:
        with open(bc_mat_filename1, 'rb') as f:
            chrom1, starts1, stops1 = cPickle.load(f)
            bc_mat1 = cPickle.load(f)

        loci1 = (chrom1, starts1, stops1)
        if len(hot_bcs) > 0:
            bc_mat1 = bc_mat1.tocsr()[hot_bcs, :]
        bc_mat1 = bc_mat1.tolil()

        b2_list = [b[1] for b in bc_mat_filename_iter2]
        for bc_mat_filename2 in b2_list:
            with open(bc_mat_filename2, 'rb') as f:
                chrom2, starts2, stops2 = cPickle.load(f)
                bc_mat2 = cPickle.load(f)

            loci2 = (chrom2, starts2, stops2)
            if len(hot_bcs) > 0:
                bc_mat2 = bc_mat2.tocsr()[hot_bcs, :]
            bc_mat2 = bc_mat2.tolil()

            overlap_detector = OverlapDetector(bc_mat1, bc_mat2, loci1, loci2, min_reads = 1,
                max_bcs = args.max_bcs)

            if test == tk_sv_utils.BINOM:
                res = overlap_detector.get_overlaps(test, max_logp = args.max_logp, min_ov = args.min_overlap, min_dist = args.min_call_dist,
                                                    max_dist = args.max_call_dist, verbose = True)
            else:
                res = overlap_detector.get_overlaps(test, min_ov = args.min_overlap, exp_bc_ov=exp_bc_ov, step=STEP,
                                                    genome_size = genome_size, frag_sizes = frag_sizes, max_logp = args.max_logp,
                                                    frag_counts = frag_counts, min_dist = args.min_call_dist,
                                                    max_dist = args.max_call_dist, bc_ov_adjust=bc_ov_adjust)
            (out_loci1, out_loci2, _, _, _, _) = res

            if len(out_loci1[1]) == 0:
                continue

            for (start1, stop1, start2, stop2) in zip(out_loci1[1], out_loci1[2], out_loci2[1], out_loci2[2]):
                if chrom1 > chrom2 or (chrom1 == chrom2 and (start1, stop1) > (start2, stop2)):
                    overlap_loci.append((chrom2, start2, stop2, chrom1, start1, stop1))
                else:
                    overlap_loci.append((chrom1, start1, stop1, chrom2, start2, stop2))

    if test != tk_sv_utils.BINOM:
        prox_loci = [locus for locus in overlap_loci if
                     locus[0] == locus[3] and locus[4] - locus[2] < MAX_FRAG_SIZE]

        dist_loci = [locus for locus in overlap_loci if
                     locus[0] != locus[3] or locus[4] - locus[2] >= MAX_FRAG_SIZE]

        new_loci = merge_loci(dist_loci)
        new_loci.extend(merge_loci(prox_loci))
    else:
        new_loci = overlap_loci

    with open(outs.overlap_loci, 'wb') as f:
        cPickle.dump(new_loci, f)


def merge_loci(loci, start_idx = 0, merge = True):
    dist_df = tk_sv_io.create_sv_df([l[0] for l in loci],
                              [l[1] for l in loci],
                              [l[2] for l in loci],
                              [l[3] for l in loci],
                              [l[4] for l in loci],
                              [l[5] for l in loci],
                              np.arange(len(loci)),
                              [1 for l in loci])

    # Merge adjacent distal windows
    groups = tk_sv_utils.get_break_groups(dist_df, 1)

    new_loci = []
    for i, g in enumerate(groups):
        g_df = dist_df.loc[g]
        c1 = g_df.iloc[0].chrom1
        c2 = g_df.iloc[0].chrom2
        s1 = np.min(g_df.start1)
        e1 = np.max(g_df.stop1)
        s2 = np.min(g_df.start2)
        e2 = np.max(g_df.stop2)

        if c1 == c2 and s2 - e1 < MIN_SV:
            # Do not merge if this would result in a huge blob
            for gi in g:
                locus = (c1, dist_df.loc[gi].start1, dist_df.loc[gi].stop1,
                         c2, dist_df.loc[gi].start2, dist_df.loc[gi].stop2)
                new_loci.append(locus)
        else:
            locus = (c1, s1, e1, c2, s2, e2)
            new_loci.append(locus)

    return new_loci
