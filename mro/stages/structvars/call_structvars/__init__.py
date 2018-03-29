#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Call large scale structural variants.
#

import sys
import os
import os.path
import cPickle
import numpy as np
from numpy.random import choice, seed
import tenkit.pandas as pd
import bisect
import martian

import tenkit.hdf5 as tk_hdf5
import tenkit.tabix as tk_tabix
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
from longranger.sv.constants import MAX_FRAG_SIZE, MIN_FRAG_SIZE, MIN_SV, MAX_DIST_FROM_CAND

from itertools import product

import longranger.sv.utils as sv_utils
import longranger.sv.stats as tk_sv_stats
import longranger.sv.read_model as tk_sv_read_model
import longranger.sv.sv_call as tk_sv_call
import longranger.sv.io as tk_sv_io

__MRO__ = """
stage CALL_STRUCTVARS(
    in  bam      possorted_bam,
    in  pickle   overlap_loci,
    in  pickle   low_depth_loci,
    in  bool     recompute,
    in  int      sv_min_qv,
    in  int      min_mapq,
    in  int      min_bcs,
    in  h5       fragments,
    in  json     fragment_histogram,
    in  tsv.gz   fragment_phasing,
    in  tsv      barcode_blacklist,
    in  json     coverage,
    in  h5       coverage_details,
    in  bed      targets,
    in  float    p_ov_mol,
    in  int      max_cand_breaks,
    in  int      grid_len,
    out bedpe    sv_variants,
    src py       "stages/structvars/call_structvars",
) split using (
    in int start_idx,
    in int stop_idx,
)
"""

MAX_CHUNKS = 500

MAX_READ_LEN = 150

# Min and max reads in a molecule in order to consider for SV-calling
MIN_READS_PER_FRAG = 5
MAX_READS_PER_FRAG = 1000
# Ignore candidate loci if they are wider than this
MAX_REGION_LEN = 300000

# Extend distal breakpoints by this amount before searching for exact breakpoint
BREAK_EXT = 10000

MAX_READ_GROUPS = 1000 # Max barcodes to consider per locus
MAX_READS_TO_READ = 10000000 # Max num of reads to read per locus

# Fetch reads within this window around breakpoints in order to form "molecules"
FRAG_EXTEND = 100000


def prepare_loci(args):
    """Merge and sort input lists of candidate loci."""

    overlap_loci = []

    # Loci based on barcode overlaps. Type of SV is unknown.
    if not args.overlap_loci is None:
        with open(args.overlap_loci, 'rb') as f:
            loci = cPickle.load(f)
        overlap_loci.extend([(x[0], x[1], x[2], x[3], x[4], x[5], None) for x in loci])

    # Low depth loci. These will only be evaluated for deletions.
    if not args.low_depth_loci is None:
        del_calls = tk_sv_io.read_sv_bedpe_to_df(args.low_depth_loci)
        for _, row in del_calls.iterrows():
            overlap_loci.append((row.chrom1, row.start1, row.stop1,
                                 row.chrom2, row.start2, row.stop2, 'DEL'))

    # Loci based on read-pair support. These will only be evaluated for the
    # type of SV supported by the readpairs.
    if not args.rp_calls is None:
        rp_calls = tk_sv_io.read_sv_bedpe_to_df(args.rp_calls)
        for _, row in rp_calls.iterrows():
            sv_type = tk_sv_io.get_sv_type(row.info)
            if not sv_type in ['DEL', 'INV', 'DUP']:
                sv_type = None
            else:
                sv_type = [sv_type]
            overlap_loci.append((row.chrom1, row.start1, row.stop1,
                                 row.chrom2, row.start2, row.stop2, sv_type))

    # Sort by position and also get the sorted indices.
    sorted_overlap_loci = sorted(overlap_loci,
                                 key=lambda x: (x[0], x[1], x[2], x[3], x[4], x[5]))
    sorted_overlap_loci_idx = sorted(range(len(overlap_loci)),
                                     key=lambda x: (overlap_loci[x][0], overlap_loci[x][1],
                                                    overlap_loci[x][2], overlap_loci[x][3],
                                                    overlap_loci[x][4], overlap_loci[x][5]))

    # If there is a single source of candidate loci, coming from a BEDPE, then
    # keep track of the names in the BEDPE, so you can annotate SV-calls made
    # with the BEDPE line from which they came.
    if args.overlap_loci is None and args.low_depth_loci is None and not args.rp_calls is None:
        input_calls = tk_sv_io.read_sv_bedpe_to_df(args.rp_calls)
        input_names = list(input_calls['name'])
        input_names = [input_names[n] for n in sorted_overlap_loci_idx]
    else:
        input_names = None

    return sorted_overlap_loci, input_names


def split(args):
    overlap_loci, _ = prepare_loci(args)

    nloci = len(overlap_loci)
    nloci_per_chunk = max(2, int(np.ceil(nloci / MAX_CHUNKS))) # avoid overchunking
    nchunks = int(np.ceil(nloci / float(nloci_per_chunk)))
    chunk_defs = []

    mem_request = 11

    for i in range(nchunks):
        chunk_start = i * nloci_per_chunk
        chunk_end = min(nloci, (i + 1) * nloci_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':mem_request})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':5}]
    return {'chunks': chunk_defs, 'join': {'__mem_gb':mem_request}}


def join(args, outs, chunk_defs, chunk_outs):
    out_bedpe = None
    for c in chunk_outs:
        if not os.path.isfile(c.sv_variants):
            continue
        in_bedpe = tk_sv_io.read_sv_bedpe_to_df(c.sv_variants)
        if not in_bedpe is None:
            out_bedpe = pd.concat([out_bedpe, in_bedpe], ignore_index=True)

    if out_bedpe is None:
        col_names = ['chrom1', 'start1', 'stop1',
                     'chrom2', 'start2', 'stop2', 'name', 'qual',
                     'strand1', 'strand2', 'filters', 'info']
        out_bedpe = pd.DataFrame(columns=col_names)
    out_bedpe.names = np.arange(len(out_bedpe))

    out_bedpe = out_bedpe[out_bedpe.qual >= args.sv_min_qv]
    tk_sv_io.write_sv_df_to_bedpe(out_bedpe, outs.sv_variants)


def isfile(fn):
    return not fn is None and os.path.isfile(fn)


def get_reads(in_bam, chrom, start, stop, in_read_df=None,
              min_mapq=30, max_reads=500000, blacklist_barcodes=None):
    poses = []
    ends = []
    bcs = []

    if not in_read_df is None and len(in_read_df) > 0:
        ret_df = in_read_df.sort('pos')
        old_poses = np.array(ret_df['pos'])
        # Subtracting the read length is roughly right, ideally we should sort
        # by aend.
        # Loci are considered in an ordered fashion, so we should never fetch
        # reads "earlier" in the bam.
        start = max(old_poses[0], max(0, start - MAX_READ_LEN))
        if start >= old_poses[0] and start <= old_poses[-1]:
            start_idx = bisect.bisect_left(old_poses, start)
            if stop >= old_poses[0] and stop <= old_poses[-1]:
                stop_idx = min(len(ret_df), bisect.bisect(old_poses, stop))
            else:
                stop_idx = len(ret_df)
            # Remove all positions that are smaller than the input start
            ret_df = ret_df.iloc[start_idx:stop_idx]
            # Set the new start to the end of the input data frame.
            # Add an overlap of READ_LEN to capture reads that were right on
            # the boundary between the old and new data frame.
            start = max(0, old_poses[stop_idx - 1] - MAX_READ_LEN)
            stop = max(start, stop)
    else:
        ret_df = None

    for i, read in enumerate(in_bam.fetch(str(chrom), int(start), int(stop))):
        if i > max_reads:
            break
        bc = tk_io.get_read_barcode(read)
        if read.pos < start:
            continue
        if not blacklist_barcodes is None and bc in blacklist_barcodes:
            continue
        if not read.is_secondary and not read.is_duplicate and read.is_read1 and \
           not read.is_unmapped and read.mapq >= min_mapq and read.is_proper_pair and \
           not bc is None:
            poses.append(read.pos)
            ends.append(read.aend)
            bcs.append(tk_io.get_read_barcode(read))

    tmp_ret_df = pd.DataFrame({'chrom':chrom, 'pos':poses, 'aend':ends, 'bc':bcs})

    ret_df = pd.concat([ret_df, tmp_ret_df], ignore_index=True)
    ret_df.sort(['bc', 'pos'], inplace=True)
    return ret_df


def loci_close(locus1, locus2):
    """Return True if the two loci overlap after accounting for extentions."""
    (c1, s1, _, c2, _, e2) = locus1
    (c1b, s1b, _, c2b, _, _) = locus2
    if c1 == c1b and c2 == c2b and c1 == c2:
        return s1b >= s1 and s1b - FRAG_EXTEND <= e2 + FRAG_EXTEND
    else:
        return False


def main(args, outs):
    """SV calling on a subset of the input loci."""

    #### Prepare input files and parameters ####
    if not isfile(args.fragment_histogram) or not isfile(args.fragments) or \
       not isfile(args.fragment_phasing):
        martian.log_info('One or more files needed for SV-calling are missing. No calls will be made.')
        tk_sv_io.write_sv_df_to_bedpe(None, outs.sv_variants)
        return

    # Get candidate loci and subset to the loci for this chunk.
    overlap_loci, input_names = prepare_loci(args)

    overlap_loci = [overlap_loci[i] for i in range(int(args.start_idx), int(args.stop_idx))]
    if not input_names is None:
        input_names = [input_names[i] for i in range(int(args.start_idx), int(args.stop_idx))]

    # Get molecule size distribution.
    frag_res = tk_sv_stats.read_frag_hist(args.fragment_histogram, MIN_FRAG_SIZE)
    frag_sizes, frag_counts = frag_res

    # Get fragment phasing info. This will be used to get the barcode phasing.
    frag_phasing = tk_tabix.create_tabix_infile(args.fragment_phasing)

    # Estimate the Poisson reads rate alpha.
    fragment_df = tk_hdf5.read_data_frame_limited(args.fragments, query_cols=['obs_len', 'num_reads'], max_rows=20000000)
    fragment_df = fragment_df[fragment_df.num_reads > MIN_READS_PER_FRAG]
    alpha = np.median(np.array(fragment_df['num_reads']) / np.array(fragment_df['obs_len'], dtype=np.float))
    martian.log_info('Using alpha = {}'.format(alpha))

    sv_model = tk_sv_read_model.ReadModel(alpha, frag_sizes, frag_counts,
                                          p_ov_mol=args.p_ov_mol, step=1000)

    if not args.targets is None:
        msg = 'Read-based SV-calling from targeted data not supported.'
        martian.log_info(msg)
        return

    # Get a set of barcodes to remove from SV-calling.
    if not args.barcode_blacklist is None:
        tmp_df = pd.read_csv(args.barcode_blacklist, sep='\t', index_col=False)
        blacklist_barcodes = set(tmp_df.bc)
    else:
        blacklist_barcodes = set([])

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    #### End prepare input files and parameters ####

    old_locus = None
    old_reads = None

    res_arr = []
    group_ids = []

    # Iterate through all loci and evaluate them for the presence of SVs
    for locus_idx, (c1, s1, e1, c2, s2, e2, _) in enumerate(overlap_loci):

        print >> sys.stderr, 'Evaluating locus', c1, s1, e1, c2, s2, e2

        # Candidate locus too wide. Skip.
        if e1 - s1 > MAX_REGION_LEN or e2 - s2 > MAX_REGION_LEN:
            print >> sys.stderr, 'Locus too wide. Skipping.'
            continue

        # Candidate loci too close to each other. Skip.
        if c1 == c2 and s2 - e1 < 2 * BREAK_EXT:
            print >> sys.stderr, 'Breakpoints of too close. Skipping.'
            continue

        # Evaluate for proximal SVs (DEL, INV, DUP) if the distance between
        # breakpoints is < MAX_FRAG_SIZE. Otherwise, the event will be called
        # a translocation and we'll try to infer the signal orientation.
        if c1 == c2 and s2 - e1 < MAX_FRAG_SIZE:
            if not old_locus is None and loci_close(old_locus, (c1, s1, e1, c2, s2, e2)):
                in_read_df = old_reads
            else:
                in_read_df = None
            res, reads = call_proximal(sv_model, c1, s1, e1, s2, e2,
                                       in_bam, in_read_df, frag_phasing,
                                       blacklist_barcodes, args)

            old_locus = (c1, s1, e1, c2, s2, e2)
            old_reads = reads
        else:
            res = call_distal(sv_model, c1, max(0, s1 - BREAK_EXT), e1 + BREAK_EXT,
                              c2, max(0, s2 - BREAK_EXT), e2 + BREAK_EXT, in_bam,
                              frag_phasing, blacklist_barcodes, args)

        if not res is None:
            res_arr.extend(res)
            group_ids.extend(locus_idx * np.ones((len(res),), dtype=np.int))

    in_bam.close()

    out_df = tk_sv_call.SvCall.svs_to_dataframe(res_arr)

    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.sv_variants)


def get_read_group_endpoints(reads):
    read_poses = np.array(reads.pos)
    uniq_poses = [np.min(read_poses), np.max(read_poses)]
    p = np.argmax(np.diff(read_poses))
    if np.max(np.diff(read_poses)) > MIN_SV:
        uniq_poses.append(read_poses[p])
        uniq_poses.append(read_poses[p+1])
    uniq_poses = np.array(uniq_poses)
    return uniq_poses


def refine_breaks(read_groups, s1, e1, s2, e2,
                  grid_len=1000, max_cand_breaks=100):
    # Split each of the loci into small windows and count how many
    # molecules end or start in each region of the resulting grid.
    ticks1 = np.arange(s1, e1 + 1, grid_len)
    ticks2 = np.arange(s2, e2 + 1, grid_len)
    tick_counts = np.zeros((len(ticks1), len(ticks2)), dtype=np.int)

    for _, reads in read_groups:
        uniq_poses = get_read_group_endpoints(reads)

        #uniq_poses = np.ones((len(read_poses),), dtype = np.bool)
        #uniq_poses[1:] = np.diff(read_poses) > 0
        #read_poses = read_poses[uniq_poses]

        xposes = uniq_poses[np.logical_and(uniq_poses >= ticks1[0], uniq_poses < ticks1[-1])]
        yposes = uniq_poses[np.logical_and(uniq_poses >= ticks2[0], uniq_poses < ticks2[-1])]

        if len(xposes) >= 0 and len(yposes) >= 0:
            pairs = [(a, b) for (a, b) in product(xposes, yposes) if b - a > MIN_SV]
            xposes = np.array([(v[0] - ticks1[0]) / grid_len for v in pairs], dtype=np.int)
            yposes = np.array([(v[1] - ticks2[0]) / grid_len for v in pairs], dtype=np.int)
            tick_counts[xposes, yposes] += 1

    m = np.argmax(tick_counts)
    x, y = np.unravel_index(m, (len(ticks1), len(ticks2)))
    new_start1 = (x * grid_len + ticks1[0])
    new_start2 = (y * grid_len + ticks2[0])

    cand_breaks1 = []
    cand_breaks2 = []
    # Get more refined breaks within the specified window
    for _, reads in read_groups:
        uniq_poses = get_read_group_endpoints(reads)
        cand_breaks1.extend(uniq_poses[np.logical_and(uniq_poses >= new_start1,
                                                      uniq_poses <= new_start1 + grid_len)])
        cand_breaks2.extend(uniq_poses[np.logical_and(uniq_poses >= new_start2,
                                                      uniq_poses <= new_start2 + grid_len)])

    if len(cand_breaks1) > max_cand_breaks:
        cand_breaks1 = sorted(list(set(choice(cand_breaks1, max_cand_breaks, replace=False))))
    else:
        cand_breaks1 = sorted(list(set(cand_breaks1)))

    if len(cand_breaks2) > max_cand_breaks:
        cand_breaks2 = sorted(list(set(choice(cand_breaks2, max_cand_breaks, replace=False))))
    else:
        cand_breaks2 = sorted(list(set(cand_breaks2)))

    return cand_breaks1, cand_breaks2


def call_proximal(sv_model, c1, s1, e1, s2, e2,
                  in_bam, in_reads,
                  frag_phasing, blacklist_barcodes, args):
    """Evaluate a pair of loci (c1, s1, e1), (c1, s2, e2) for proximal SV calls.

    Return value: A tuple (results, reads).
    results is a list of SvCall objects.
    reads is a DataFrame with the fetched reads at the locus.
    """

    seed(1)

    # Get the phase set on each of the two loci
    ps1 = sv_utils.get_phase_set(frag_phasing, c1, max(0, s1 - BREAK_EXT), e1 + BREAK_EXT)
    ps2 = sv_utils.get_phase_set(frag_phasing, c1, max(0, s2 - BREAK_EXT), e2 + BREAK_EXT)

    # Fetch the reads on the entire locus (c1, s1, e2)
    # Add some padding to make sure you fetched all the reads from the
    # involved molecules.
    reads = get_reads(in_bam, c1, max(0, s1 - FRAG_EXTEND),
                      e2 + FRAG_EXTEND, in_read_df=in_reads,
                      min_mapq=args.min_mapq, max_reads=MAX_READS_TO_READ,
                      blacklist_barcodes=blacklist_barcodes)

    break_filter = filter_fun([-np.inf, np.inf], set(), set())
    read_groups = reads.groupby('bc').filter(break_filter).groupby('bc')

    if e1 - s1 > args.grid_len and e2 - s2 > args.grid_len:
        cand_breaks1, cand_breaks2 = refine_breaks(read_groups, max(0, s1 - BREAK_EXT), e1 + BREAK_EXT,
                                                   max(0, s2 - BREAK_EXT), e2 + BREAK_EXT,
                                                   grid_len=args.grid_len,
                                                   max_cand_breaks=args.max_cand_breaks)
    else:
        cand_breaks1, cand_breaks2 = refine_breaks(read_groups, s1 - args.grid_len, e1 + args.grid_len,
                                                   s2 - args.grid_len, e2 + args.grid_len,
                                                   grid_len=args.grid_len,
                                                   max_cand_breaks=args.max_cand_breaks)

    if len(cand_breaks1) == 0 or len(cand_breaks2) == 0:
        return (None, reads)

    # Restrict the read groups based on the new better estimate of the candidate
    # breakpoints.
    break_filter = filter_fun([cand_breaks1[0], cand_breaks2[-1]],
                              blacklist_barcodes, set([]))
    read_groups = read_groups.filter(break_filter).groupby('bc')

    loci_pairs = [(a, b) for a, b in product(cand_breaks1, cand_breaks2) if b - a > MIN_SV]
    if len(loci_pairs) == 0:
        return (None, reads)

    #read_groups = sorted(read_groups, key = lambda x: len(x[1]), reverse = True)
    #read_groups = read_groups[0:min(len(read_groups), MAX_READ_GROUPS)]

    read_group_dict = {}
    for bc_idx, (bc, read_group) in enumerate(read_groups):
        if bc_idx > MAX_READ_GROUPS:
            break
        read_group_dict[bc] = tk_sv_read_model.ReadInfo(read_group.chrom, read_group.pos, group_ids=None)

    bc_set = set(read_group_dict.keys())
    bc_phase_set_dict1 = sv_utils.get_barcode_phase_probs(frag_phasing, c1,
                                                          cand_breaks1[0], cand_breaks1[-1], bc_set, in_ps=ps1)
    bc_phase_set_dict2 = sv_utils.get_barcode_phase_probs(frag_phasing, c1,
                                                          cand_breaks2[0], cand_breaks2[-1], bc_set, in_ps=ps2)

    res = sv_model.em_it_away(loci_pairs, read_group_dict, ps1, ps2,
                              bc_phase_set_dict1, bc_phase_set_dict2,
                              em_iters=5)

    (max_lrs, max_locus, sv_type, zygosity, max_hap, support, posterior_probs) = res
    sv_call = tk_sv_call.SvCall.from_em_results(c1, c1, ps1, ps2, max_lrs,
                                                max_locus, sv_type, zygosity, max_hap,
                                                support, posterior_probs)

    return ([sv_call], reads)


def filter_fun(breaks, blacklist_barcodes, common_bcs):
    """Return a function for filtering groups of barcodes.

    The returned function will only keep reads that satisfy the following:
    - their barcode is not blacklisted
    - there are at least MIN_READS_PER_FRAG with the corresponding barcode and
    at most MAX_READS_PER_FRAG such reads
    - the reads with the corresponding barcode span more than MIN_FRAG_SIZE bases
    - the barcode is in the common_bcs and it does not fall entirely before
    breaks[0] or after breaks[1]"""

    return lambda x: not x.iloc[0].bc in blacklist_barcodes and \
        len(x) <= MAX_READS_PER_FRAG and len(x) >= MIN_READS_PER_FRAG and \
        np.max(x.pos) - np.min(x.pos) > MIN_FRAG_SIZE and \
        (x.iloc[0].bc in common_bcs or not(np.min(x.pos) > breaks[-1] or np.max(x.pos) < breaks[0]))


def filter_irrelevant_frags1(svt, cand_breaks1):
    """Return a filter function for molecules not participating in the SV"""
    if svt == 'TRANS_0e' or svt == 'TRANS_00':
        return lambda x: not(np.all(x.pos < cand_breaks1[0]) or \
                             np.all(x.pos > cand_breaks1[-1] + MAX_DIST_FROM_CAND))
    return lambda x: not(np.all(x.pos > cand_breaks1[-1]) or \
                         np.all(x.pos < cand_breaks1[0] - MAX_DIST_FROM_CAND))


def filter_irrelevant_frags2(svt, cand_breaks2):
    if svt == 'TRANS_e0' or svt == 'TRANS_00':
        return lambda x: not(np.all(x.pos < cand_breaks2[0]) or \
                             np.all(x.pos > cand_breaks2[-1] + MAX_DIST_FROM_CAND))
    return lambda x: not(np.all(x.pos > cand_breaks2[-1]) or \
                         np.all(x.pos < cand_breaks2[0] - MAX_DIST_FROM_CAND))


def call_distal(sv_model, c1, s1, e1, c2, s2, e2, in_bam,
                frag_phasing, blacklist_barcodes, args):
    """Evaluate a pair of loci (c1, s1, e1), (c2, s2, e2) for distal SV calls.

    Return value: a list of SvCall objects.
    """

    # Set the random seed so if we run this for the same locus multiple times we
    # always get the same answer.
    seed(1)

    # Fetch reads from each of the two loci
    reads1 = get_reads(in_bam, c1, max(0, s1 - FRAG_EXTEND), e1 + FRAG_EXTEND,
                       min_mapq=args.min_mapq, max_reads=MAX_READS_TO_READ,
                       blacklist_barcodes=blacklist_barcodes)
    reads1 = reads1.groupby('bc').filter(lambda x: len(x) <= MAX_READS_PER_FRAG and \
                                         len(x) >= MIN_READS_PER_FRAG / 2 and \
                                         np.max(x.pos) - np.min(x.pos) > MIN_FRAG_SIZE / 2)

    reads2 = get_reads(in_bam, c2, max(0, s2 - FRAG_EXTEND), e2 + FRAG_EXTEND,
                       min_mapq=args.min_mapq, max_reads=MAX_READS_TO_READ,
                       blacklist_barcodes=blacklist_barcodes)
    reads2 = reads2.groupby('bc').filter(lambda x: len(x) <= MAX_READS_PER_FRAG and \
                                         len(x) >= MIN_READS_PER_FRAG / 2 and \
                                         np.max(x.pos) - np.min(x.pos) > MIN_FRAG_SIZE / 2)

    common_bcs = (set(reads1.bc).intersection(set(reads2.bc)))
    if len(common_bcs) < args.min_bcs:
        return None

    common_reads1 = reads1[np.array([bc in common_bcs for bc in reads1.bc], dtype=np.bool)]
    common_reads2 = reads2[np.array([bc in common_bcs for bc in reads2.bc], dtype=np.bool)]

    # Split each of the loci into small windows and count how many
    # molecules end or start in each region of the resulting grid.
    ticks1 = np.arange(np.min(common_reads1.pos), np.max(common_reads1.pos), args.grid_len)
    ticks2 = np.arange(np.min(common_reads2.pos), np.max(common_reads2.pos), args.grid_len)
    tick_counts = np.zeros((len(ticks1), len(ticks2)), dtype=np.int)
    nticks = len(ticks1) * len(ticks2)
    # barcode x cell, 1 if there is barcode overlap in the cell for that barcode
    bc_tick_counts = np.zeros((len(common_bcs), nticks), dtype=np.bool)

    common_reads1 = common_reads1.groupby('bc')
    common_reads2 = common_reads2.groupby('bc')

    common_bcs = list(common_bcs)
    for i, bc in enumerate(common_bcs):
        r1 = common_reads1.get_group(bc)
        r2 = common_reads2.get_group(bc)
        x = (np.array([np.min(r1.pos), np.max(r1.pos)], dtype=np.int) - ticks1[0]) / args.grid_len
        y = (np.array([np.min(r2.pos), np.max(r2.pos)], dtype=np.int) - ticks2[0]) / args.grid_len
        x = x[np.logical_and(x >= 0, x < len(ticks1))]
        y = y[np.logical_and(y >= 0, y < len(ticks2))]
        vals = list(set(product(x, y)))
        a = np.array([v[0] for v in vals], dtype=np.int)
        b = np.array([v[1] for v in vals], dtype=np.int)
        tick_counts[a, b] += 1
        bc_tick_counts[i, np.ravel_multi_index((a, b), (len(ticks1), len(ticks2)))] = True

    rounds = 0
    final_res = []
    ps1 = sv_utils.get_phase_set(frag_phasing, c1, s1, e1)
    ps2 = sv_utils.get_phase_set(frag_phasing, c2, s2, e2)

    # Get the cell of the grid with the maximum barcode overlap, get the
    # barcodes/molecules involved, figure out the direction of the signal,
    # get a better estimate of the breakpoints, and compute the SV score at
    # these breakpoints. Then, subtract the contribution of the involved
    # barcodes from the barcode overlap grid and repeat.
    # In practice, we only do it once.
    while np.max(tick_counts) > 2 and len(common_reads1) > 2 and len(common_reads2) > 2 and rounds < 1:
        rounds += 1

        # Get region of maximum barcode overlap. m is flat index
        m = np.argmax(tick_counts)
        x, y = np.unravel_index(m, (len(ticks1), len(ticks2)))
        new_start1 = (x * args.grid_len + ticks1[0])
        new_start2 = (y * args.grid_len + ticks2[0])

        # Get the fragments overlapping the region of max barcode overlap
        common_reads1 = common_reads1.filter(lambda x: not(np.max(x.pos) < new_start1 - 5000 or
                                                           np.min(x.pos) > new_start1 + 5000))
        common_reads2 = common_reads2.filter(lambda x: not(np.max(x.pos) < new_start2 - 5000 or
                                                           np.min(x.pos) > new_start2 + 5000))
        # Get fragments with common barcodes that overlap the selected cell
        # This is a superset of np.where(bc_tick_counts[:, m])[0] because it doesn't require that
        # the fragment starts or ends within the cell.
        new_common_bcs = set(common_reads1.bc).intersection(set(common_reads2.bc))

        # Infer orientation based on ends of fragments
        sv_types, cand_breaks1, cand_breaks2 = get_orient(common_reads1.groupby('bc'),
                                                          common_reads2.groupby('bc'),
                                                          new_start1, new_start2,
                                                          new_common_bcs,
                                                          max_cand_breaks=args.max_cand_breaks)
        if len(cand_breaks1) == 0 or len(cand_breaks2) == 0 or len(sv_types) == 0:
            return None

        svt = sv_types[0]

        # Use the inferred orientation of the signal to remove barcodes that are
        # clearly irrelevant.
        tmp_read_groups1 = reads1.groupby('bc').filter(filter_irrelevant_frags1(svt, cand_breaks1))
        tmp_read_groups1['break'] = 1

        tmp_read_groups2 = reads2.groupby('bc').filter(filter_irrelevant_frags2(svt, cand_breaks2))
        tmp_read_groups2['break'] = 2

        read_groups = pd.concat([tmp_read_groups1, tmp_read_groups2], ignore_index=True)

        if len(read_groups) == 0:
            return None

        sel_bcs = set(read_groups.bc)
        if len(sel_bcs) > MAX_READ_GROUPS:
            sel_bcs = list(sel_bcs)
            sel_bcs = set([sel_bcs[i] for i in choice(len(sel_bcs), MAX_READ_GROUPS, replace=False)])

        read_groups = read_groups[[_b in sel_bcs for _b in read_groups.bc]]

        read_groups.sort('bc', inplace=True)
        bc_set = set(read_groups.bc)
        read_groups = read_groups.groupby('bc')

        bc_phase_set_dict1 = sv_utils.get_barcode_phase_probs(frag_phasing, c1, cand_breaks1[0], cand_breaks1[-1], bc_set, in_ps=ps1)
        bc_phase_set_dict2 = sv_utils.get_barcode_phase_probs(frag_phasing, c2, cand_breaks2[0], cand_breaks2[-1], bc_set, in_ps=ps2)

        loci_pairs = [(_a, _b) for _a, _b in product(cand_breaks1, cand_breaks2)]

        read_group_dict = {}
        for bc, read_group in read_groups:
            read_group_dict[bc] = tk_sv_read_model.ReadInfo(read_group.chrom, read_group.pos, group_ids=read_group['break'])

        res = sv_model.em_it_away(loci_pairs, read_group_dict, ps1, ps2,
                                  bc_phase_set_dict1, bc_phase_set_dict2,
                                  em_iters=5, proximal=False)

        (max_lrs, max_locus, sv_type, zygosity, max_hap, support, posterior_probs) = res
        sv_call = tk_sv_call.SvCall.from_em_results(c1, c2, ps1, ps2, max_lrs,
                                                    max_locus, sv_type, zygosity, max_hap,
                                                    support, posterior_probs)
        final_res.append(sv_call)

    return final_res


def frag_overlaps(frags, s, e):
    return np.logical_not(np.logical_or(np.min(frags.pos) > e, np.max(frags.pos) < s))
    #return np.logical_not(np.logical_or(frags.start_pos > e, frags.end_pos < s))


def get_orient(reads1, reads2, pos1, pos2, common_bcs, max_cand_breaks=100):
    """Estimate the orientation of a distal SV and get refined breakpoints.

    Inputs:
    - frags1/2: molecules on each of the two breakpoints
    - pos1/2: positions of max barcode overlap
    """

    counts = [0, 0, 0, 0]
    cand_breaks1 = []
    cand_breaks2 = []
    for bc in common_bcs:
        r1 = reads1.get_group(bc)
        r2 = reads2.get_group(bc)

        start1, stop1 = np.min(r1.pos), np.max(r1.pos)
        start2, stop2 = np.min(r2.pos), np.max(r2.pos)

        a = int(abs(pos1 - start1) > abs(stop1 - pos1))
        b = int(abs(pos2 - start2) > abs(stop2 - pos2))
        counts[a + 2 * b] += 1
        if a == 1:
            cand_breaks1.append(stop1)
        else:
            cand_breaks1.append(start1)
        if b == 1:
            cand_breaks2.append(stop2)
        else:
            cand_breaks2.append(start2)

    cand_breaks1 = np.array(sorted(cand_breaks1))
    cand_breaks2 = np.array(sorted(cand_breaks2))
    cand_breaks1 = cand_breaks1[np.logical_and(cand_breaks1 >= pos1 - 1000,
                                               cand_breaks1 <= pos1 + 1000)]
    cand_breaks2 = cand_breaks2[np.logical_and(cand_breaks2 >= pos2 - 1000,
                                               cand_breaks2 <= pos2 + 1000)]

    if len(cand_breaks1) > max_cand_breaks:
        cand_breaks1 = sorted(list(set(choice(cand_breaks1, max_cand_breaks, replace=False))))
    if len(cand_breaks2) > max_cand_breaks:
        cand_breaks2 = sorted(list(set(choice(cand_breaks2, max_cand_breaks, replace=False))))

    idx = np.argsort(counts)[::-1]
    idx = idx[0:1]

    types = ['TRANS_00', 'TRANS_e0', 'TRANS_0e', 'TRANS_ee']
    test_types = []
    for i in idx:
        if counts[i] > 0.25 * np.sum(counts):
            test_types.append(types[i])
    return test_types, cand_breaks1, cand_breaks2
