#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Call large scale structural variants based on molecule-level information.
# This should be used for targeted settings. For WGS, CALL_STRUCTVARS should
# be used.
#

import sys
import os
import os.path
import shutil
import cPickle
import numpy as np
import tenkit.pandas as pd
import json
from itertools import repeat

import tenkit.safe_json
import tenkit.tabix as tk_tabix
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.hdf5 as tk_hdf5
from longranger.sv.constants import (MIN_READS_PER_FRAG_TARGET,
                                     MIN_FRAG_SIZE_TARGET,
                                     SV_FRAGMENT_LINK_DISTANCE_TARGET,
                                     MAX_FRAG_SIZE, TARGET_COV_BIN)

import martian

import longranger.sv.stats as sv_stats
import longranger.sv.utils as sv_utils
import longranger.sv.io as sv_io
import longranger.sv.sv_call as sv_call

__MRO__ = """
stage CALL_STRUCTVARS_FRAG(
    in  bam    input              "possorted BAM (doesnt need to be phased)",
    in  int    min_mapq,
    in  pickle overlap_loci       "list of loci to consider from DETECT_OVERLAPS",
    in  int    sv_min_qv,
    in  int    min_call_dist      "min SV to call (max of this and frag_size_prc-implied will be used)",
    in  h5     fragments,
    in  json   fragment_histogram,
    in  tsv.gz fragment_phasing,
    in  tsv    barcode_blacklist,
    in  json   coverage,
    in  bed    targets,
    in  float  frag_size_prc      "percentile of fragment size distribution to use for computing min SV to call",
    in  int    grid_len           "window size for refining barcode overlap regions in the absence of targets",
    in  int    break_ext          "input candidate loci will be extended by +/-break_ext",
    in  float  p_ov_mol           "probability of barcode collisions",
    in  int    min_frag_size,
    in  int    min_reads_per_frag,
    out bedpe  sv_variants,
    out json   summary,
    src py     "stages/structvars/call_structvars_frag",
) split using (
    in  int    start_idx,
    in  int    stop_idx,
)

"""

FRAG_EXTEND = 100000

def get_reads(in_bam, chrom, start, stop, min_mapq=60):
    bcs = []
    poses = []

    for read in in_bam.fetch(str(chrom), start, stop):
        mapq = read.mapq
        if mapq < min_mapq or read.is_secondary or read.is_duplicate:
            continue
        bc = tk_io.get_read_barcode(read)
        if bc is None:
            continue

        bcs.append(bc)
        poses.append(read.pos)
    df = pd.DataFrame({'pos':poses, 'bc':bcs})
    df.sort('bc', inplace=True)
    return df


def get_frags_from_reads(in_bam, chrom1, start1, stop1, chrom2, start2, stop2,
                         min_mapq=60, min_sv_len=45000,
                         min_reads_per_frag=MIN_READS_PER_FRAG_TARGET,
                         min_frag_size=MIN_FRAG_SIZE_TARGET):
    """Reconstruct molecules around two loci."""

    if chrom1 == chrom2 and start2 - stop1 < MAX_FRAG_SIZE:
        # Hard case: The two loci are close enough that we could have molecules
        # spanning them.
        frag_starts = []
        frag_stops = []
        frag_reads = []
        frag_bcs = []

        reads = get_reads(in_bam, chrom1, max(0, start1 - FRAG_EXTEND),
                          stop2 + FRAG_EXTEND, min_mapq=min_mapq).groupby('bc')

        for bc, group in reads:
            poses = np.array(group.pos)
            # Split positions into groups separated by a gap > min_sv_len
            pos_diff = np.where(np.diff(poses) > min_sv_len)[0]
            new_starts = np.concatenate([np.array([0]), pos_diff + 1])
            new_stops = np.concatenate([pos_diff, np.array([len(poses) - 1])])
            frag_starts.extend(poses[new_starts])
            frag_stops.extend(poses[new_stops])
            frag_reads.extend(new_stops - new_starts + 1)
            frag_bcs.extend(repeat(bc, len(pos_diff) + 1))


        frags = pd.DataFrame({'bc':frag_bcs, 'start_pos':frag_starts, 'end_pos':frag_stops,
                              'num_reads':frag_reads})
        # Remove spanning fragments
        frags = frags[(frags.start_pos > stop1) | (frags.end_pos < start2)]
        frags1 = frags.copy()
        frags2 = frags.copy()
    else:
        reads1 = get_reads(in_bam, chrom1, max(0, start1 - FRAG_EXTEND),
                           stop1 + FRAG_EXTEND, min_mapq=min_mapq)
        reads2 = get_reads(in_bam, chrom2, max(0, start2 - FRAG_EXTEND),
                           stop2 + FRAG_EXTEND, min_mapq=min_mapq)

        frags1 = reads1.groupby('bc').agg(['min', 'max', 'count'])['pos'].reset_index()
        frags2 = reads2.groupby('bc').agg(['min', 'max', 'count'])['pos'].reset_index()

        frags1.columns = ['bc', 'start_pos', 'end_pos', 'num_reads']
        frags2.columns = ['bc', 'start_pos', 'end_pos', 'num_reads']

    frags1 = frags1[(frags1.end_pos > start1) & (frags1.start_pos < stop1) & \
                    (frags1.num_reads > min_reads_per_frag) & \
                    (frags1.end_pos - frags1.start_pos > min_frag_size)]
    frags2 = frags2[(frags2.end_pos > start2) & (frags2.start_pos < stop2) & \
                    (frags2.num_reads > min_reads_per_frag) & \
                    (frags2.end_pos - frags2.start_pos > min_frag_size)]
    frags1['chrom'] = chrom1
    frags2['chrom'] = chrom2

    return frags1, frags2


def split(args):

    with open(args.overlap_loci, 'rb') as f:
        overlap_loci = cPickle.load(f)

    nloci = len(overlap_loci)
    nloci_per_chunk = max(10, int(np.ceil(nloci / 500.0))) # avoid overchunking
    nchunks = int(np.ceil(nloci / float(nloci_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nloci_per_chunk
        chunk_end = min(nloci, (i + 1) * nloci_per_chunk)
        chunk_defs.append({'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':8})
    if len(chunk_defs) == 0:
        chunk_defs = [{'start_idx':0, 'stop_idx':0, '__mem_gb':8}]
    return {'chunks': chunk_defs, 'join': {'__mem_gb':16}}


def join(args, outs, chunk_defs, chunk_outs):
    out_bedpe = None
    for c in chunk_outs:
        if not os.path.isfile(c.sv_variants):
            continue
        in_bedpe = sv_io.read_sv_bedpe_to_df(c.sv_variants)
        if not in_bedpe is None:
            out_bedpe = pd.concat([out_bedpe, in_bedpe], ignore_index=True)

    if not out_bedpe is None:
        out_bedpe['name'] = np.arange(len(out_bedpe))
    sv_io.write_sv_df_to_bedpe(out_bedpe, outs.sv_variants)

    if chunk_outs[0] is not None and os.path.exists(chunk_outs[0].summary):
        shutil.copyfile(chunk_outs[0].summary, outs.summary)
    else:
        outs.summary = None


def isfile(fn):
    return not fn is None and os.path.isfile(fn)


def main(args, outs):
    if not isfile(args.fragment_histogram) \
        or not isfile(args.barcode_blacklist) or not isfile(args.coverage):
        martian.log_info('One or more files needed for SV-calling are missing. No calls will be made.')
        return

    in_bam = tk_bam.create_bam_infile(args.input)
    genome_size = np.sum(np.array(in_bam.lengths))

    frag_hist_file = args.fragment_histogram
    barcode_blacklist_file = args.barcode_blacklist

    if args.targets is None:
        martian.exit('You should use CALL_STRUCTVARS for WGS samples.')
    else:
        target_regions = sv_utils.bed_to_region_map(args.targets, merge=True)
        target_coverage = sv_utils.region_cum_coverage_map(target_regions, TARGET_COV_BIN)
        link_distance = SV_FRAGMENT_LINK_DISTANCE_TARGET
        with open(args.coverage, 'r') as f:
            cov_sum = json.load(f)['target_info']
        if 'on_target_bases' in cov_sum:
            prob_off_target = 1 - cov_sum['on_target_bases'] / float(cov_sum['total_bases'])
        else:
            prob_off_target = 0.001
        corr_factor = sv_stats.off_target_amp_corr_factor(target_regions, prob_off_target,
                                                          genome_size=genome_size)

    res = sv_stats.get_frag_data(frag_hist_file, barcode_blacklist_file,
                                 min_frag_size=0, frag_size_prc_cutoff=args.frag_size_prc)
    frag_sizes, frag_counts, frag_prc, blacklist_barcodes = res

    frag_phasing = tk_tabix.create_tabix_infile(args.fragment_phasing)

    min_sv_len = int(max(0.8 * args.min_call_dist, link_distance))
    if not frag_prc is None and not args.targets is None:
        min_sv_len = int(max(min_sv_len, frag_prc))
    martian.log_info('Calling SVs with min length: {:d}'.format(min_sv_len))

    fragment_df = tk_hdf5.read_data_frame(args.fragments, query_cols=['obs_len', 'num_reads'])
    fragment_df = fragment_df[fragment_df.num_reads > MIN_READS_PER_FRAG_TARGET]
    alpha = np.median(np.array(fragment_df['num_reads']) / np.array(fragment_df['obs_len'], dtype=np.float))
    martian.log_info('Using alpha = {}'.format(alpha))

    summary = {}
    summary['min_sv_len'] = min_sv_len
    with open(outs.summary, 'w') as out_fn:
        out_fn.write(tenkit.safe_json.safe_jsonify(summary, pretty=True))

    model = sv_stats.FragModel(frag_sizes, frag_counts, blacklist_barcodes,
                               target_coverage, cov_bin=TARGET_COV_BIN,
                               corr_factor=corr_factor, genome_size=genome_size,
                               target_regions=target_regions, alpha=alpha,
                               p_ov_mol=args.p_ov_mol)

    with open(args.overlap_loci, 'rb') as f:
        overlap_loci = cPickle.load(f)

    overlap_loci = [overlap_loci[i] for i in range(int(args.start_idx), int(args.stop_idx))]

    final_res = []

    for i, (c1, s1, e1, c2, s2, e2) in enumerate(overlap_loci):
        frags1, frags2 = get_frags_from_reads(in_bam, c1, max(0, s1 - args.break_ext), e1 + args.break_ext,
                                              c2, max(0, s2 - args.break_ext), e2 + args.break_ext,
                                              min_mapq=args.min_mapq,
                                              min_sv_len=SV_FRAGMENT_LINK_DISTANCE_TARGET,
                                              min_frag_size=args.min_frag_size,
                                              min_reads_per_frag=args.min_reads_per_frag)

        bc_set = set(frags1.bc).union(set(frags2.bc))

        ps1 = sv_utils.get_phase_set(frag_phasing, c1, s1, e1)
        ps2 = sv_utils.get_phase_set(frag_phasing, c2, s2, e2)
        bc_phase_set_dict1 = sv_utils.get_barcode_phase_probs(frag_phasing, c1, s1, e1, bc_set, in_ps=ps1)
        bc_phase_set_dict2 = sv_utils.get_barcode_phase_probs(frag_phasing, c2, s2, e2, bc_set, in_ps=ps2)

        print >> sys.stderr, 'Evaluating locus', c1, s1, e1, c2, s2, e2

        res = model.eval_sv(frags1, frags2, (c1, s1, e1), (c2, s2, e2),
                            min_dist=min_sv_len, ps1=ps1, ps2=ps2,
                            phase_set_dict1=bc_phase_set_dict1,
                            phase_set_dict2=bc_phase_set_dict2,
                            grid_len=args.grid_len)

        if not res is None and res.qual >= args.sv_min_qv:
            final_res.append(res)

    in_bam.close()

    out_df = sv_call.SvCall.svs_to_dataframe(final_res)
    sv_io.write_sv_df_to_bedpe(out_df, outs.sv_variants)
