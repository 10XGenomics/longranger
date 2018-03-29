#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Finds pairs of genomic windows with significant barcode overlaps.
#

import os.path
from itertools import groupby
import numpy as np
from scipy.stats import binom
import tenkit.pandas as pd
import tenkit.reference as tk_reference
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import longranger.sv.io as tk_sv_io
import tenkit.regions as tk_regions
from tenkit.constants import PARALLEL_LOCUS_SIZE
from tenkit.chunk_utils import generate_chrom_loci, pack_loci
from tenkit.coverage import get_hap_coverage

import martian

__MRO__ = """
stage GET_DEL_CANDIDATES(
    in  bam    possorted_bam,
    in  string sex,
    in  string reference_path,
    in  string barcode_whitelist,
    in  int    min_mapq,
    in  int    min_del_len,
    in  int    max_del_len,
    in  int    min_bad_region,
    in  float  transition_prob,
    in  float  het_read_prob,
    out bedpe  del_candidates,
    src py     "stages/structvars/get_del_candidates",
) split using (
    in  string[] loci,
)
"""

MIN_LOG_PROB = -10
DEL_PRIOR_PROB = 0.5
MIN_COV = 2
OK_COV = 4
PHASED_COV_WINDOW = 3001
MAX_COV = 200
MIN_GAP = 200

def split(args):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    if not args.reference_path is None and os.path.exists(tk_reference.get_primary_contigs(args.reference_path)):
        with open(tk_reference.get_primary_contigs(args.reference_path), 'r') as f:
            primary_contigs = set([line.strip() for line in f.readlines()])
    else:
        # Default is to include all contigs
        primary_contigs = set(in_bam.references)

    all_loci = []
    for (chrom_name, chrom_size) in zip(in_bam.references, in_bam.lengths):
        if chrom_name in primary_contigs and not (args.sex in ['f', 'female'] and chrom_name in ['Y', 'chrY']):
            all_loci.extend(generate_chrom_loci(None, chrom_name, chrom_size, PARALLEL_LOCUS_SIZE))
    in_bam.close()

    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci, '__mem_gb':8} for loci in locus_sets]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    out_df = None
    for chunk in chunk_outs:
        tmp_df = tk_sv_io.read_sv_bedpe_to_df(chunk.del_candidates)
        out_df = pd.concat([out_df, tmp_df], ignore_index=True)

    out_df['name'] = np.arange(len(out_df))
    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.del_candidates)


def get_candidate_del_loci(hap_cov, transition_prob=1e-2, het_read_prob=0.9):
    sel_cols = ['cov_q30_hap' + str(i) for i in range(3)]
    #hap_cov['total_cov'] = hap_cov[sel_cols].sum(axis=1)
    hap_cov['total_cov'] = hap_cov["cov_q30_hap0"] + hap_cov["cov_q30_hap1"]
    npos = len(hap_cov)

    ### Emission probabilities
    em_probs = np.ones((npos, 2)) * MIN_LOG_PROB
    # emission given no del
    em_probs[:, 0] = np.maximum(MIN_LOG_PROB, binom.logpmf(np.maximum(hap_cov.cov_q30_hap0, hap_cov.cov_q30_hap1),
                                                           np.array(hap_cov.total_cov), 0.65))
    em_probs[:, 1] = np.maximum(MIN_LOG_PROB, binom.logpmf(np.maximum(hap_cov.cov_q30_hap0, hap_cov.cov_q30_hap1),
                                                           np.array(hap_cov.total_cov), het_read_prob))
    ### Transition probabilities
    # [no-del -> no-del , del -> no-del,
    #  no-del -> del    , del -> del]
    trans_probs = np.array([[1 - transition_prob, transition_prob],
                            [transition_prob    , 1 - transition_prob]])
    trans_probs = np.log(trans_probs)

    ### Prior state probabilities (prob of starting on a del or no del)
    del_prior_prob = DEL_PRIOR_PROB
    priors = np.array([1 - del_prior_prob, del_prior_prob])
    priors = np.log(priors)

    max_probs = np.zeros((npos, 2))
    max_state = np.zeros((npos, 2))

    max_probs[0, :] = priors + em_probs[0, :]

    for i in range(1, npos):
        # max_probs[i-1, 0] is the probability of the most probable path
        # (i.e. hidden state sequence) that ends at position i-1 with a 0
        new_probs = max_probs[i - 1, :] + trans_probs
        max_probs[i, :] = em_probs[i, :] + np.max(new_probs, axis=1)
        max_state[i, :] = np.argmax(new_probs, axis=1)

    best_path = np.zeros((npos, ))
    best_path[-1] = np.argmax(max_probs[-1, :])
    for i in range(npos - 2, 0, -1):
        best_path[i] = max_state[i, best_path[i + 1]]

    return best_path


def group_bit_arr(arr, start):
    ret_regions = []
    pos = start
    for bit, group in groupby(arr):
        group_size = len(list(group))
        group_start = pos
        group_stop = group_start + group_size
        if bit:
            ret_regions.append((group_start, group_stop))
        pos += group_size
    return ret_regions


def main(args, outs):

    if args.barcode_whitelist is None:
        # write empty dataframe
        tk_sv_io.write_sv_df_to_bedpe(None, outs.del_candidates)
        martian.log_info('Data seem un-barcoded. No deletion candidates will be computed.')
        return

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    del_loci = []
    for (chrom, start, stop) in (tk_io.get_locus_info(l) for l in args.loci):
        cov_df = get_hap_coverage(in_bam, None, chrom, start, stop, cov_quals=[30])
        best_path = get_candidate_del_loci(cov_df, transition_prob=args.transition_prob, het_read_prob=args.het_read_prob)

        # Get regions with good coverage for a het del (not too high, not too low)
        bad_cov = np.logical_or(cov_df['total_cov'] < MIN_COV,
                                cov_df['total_cov'] > MAX_COV)
        bad_regions = tk_regions.Regions([ (s,e) for (s,e) in group_bit_arr(bad_cov, start) if e-s > args.min_bad_region])

        # Group the states of the HMM and exclude bad regions
        pos = start
        out_loci = []
        for bit, group in groupby(best_path):
            group_size = len(list(group))
            group_start = pos
            group_stop = group_start + group_size
            if bit and group_size >= args.min_del_len and group_size <= args.max_del_len and \
               not bad_regions.overlapping_regions(group_start, group_stop):
                out_loci.append((chrom, group_start, group_stop))
            pos += group_size

        # Get regions that look like hom dels
        hom_del_loci = group_bit_arr(cov_df['total_cov'] < MIN_COV, start)
        out_loci.extend([(chrom, s, e) for (s, e) in hom_del_loci])
        out_loci = sorted(out_loci)

        # Now merge deletion candidates that are separated by short non-dels
        if out_loci:
            new_out_loci = []
            last_locus = out_loci[0]
            for i, locus in enumerate(out_loci[1:]):
                if locus[1] - last_locus[2] > MIN_GAP:
                    new_out_loci.append(last_locus)
                    last_locus = locus
                else:
                    last_locus = (last_locus[0], min(locus[1], last_locus[1]), max(locus[2], last_locus[2]))
            new_out_loci.append(last_locus)

            del_loci.extend(new_out_loci)

    final_loci = [locus for locus in del_loci if locus[2] - locus[1] >= args.min_del_len and locus[2] - locus[1] <= args.max_del_len]
    info_strs = ['TYPE=DEL' for _ in final_loci]
    in_bam.close()

    chroms = [locus[0] for locus in final_loci]
    starts1 = np.array([locus[1] for locus in final_loci], dtype=np.int)
    starts2 = np.array([locus[2] for locus in final_loci], dtype=np.int)
    sv_df = tk_sv_io.create_sv_df(chroms, starts1, starts1 + 1,
                                  chroms, starts2, starts2 + 1,
                                  np.arange(len(chroms)), 1, info_strs = info_strs)
    tk_sv_io.write_sv_df_to_bedpe(sv_df, outs.del_candidates)
