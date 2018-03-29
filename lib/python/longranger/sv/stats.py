#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for computing fragment overlaps and calling large scale structural variants.
#
import numpy as np
from scipy.stats import poisson
from bisect import bisect
import json
import pandas as pd
import longranger.sv.sv_call as sv_call
import tenkit.stats as tk_stats
import tenkit.regions as tk_regions
from longranger.sv.constants import (MIN_LOG_PROB, MAX_FRAG_SIZE, TARGET_COV_BIN,
                                     SV_FRAG_SIZE_PRC, MAX_DIST_FROM_CAND)
from collections import defaultdict
from itertools import product


def merge_freq_bins(frag_sizes, frag_counts, merge_thresh):
    """Merges bins of a histogram with similar frequency.
    Args:
    - frag_sizes: Ends of bins.
    - frag_counts: counts in bins.
    - merge_thresh: maximum difference in the frequencies (not counts) of consecutive
    bins in order to merge.
    Return value:
    A tuple (new_frag_sizes, new_frag_counts)
    """
    frag_sizes = np.array(frag_sizes)
    frag_counts = np.array(frag_counts)
    if len(frag_sizes) > 0:
        new_frag_counts = [frag_counts[0]]
        new_frag_sizes = [frag_sizes[0]]
        for c, s in zip(frag_counts[1:], frag_sizes[1:]):
            if abs((c - new_frag_counts[-1]) / float(frag_counts.sum())) < merge_thresh:
                # difference in frequencies is low, so we can merge the bins
                new_frag_counts[-1] += c
                new_frag_sizes[-1] = s
            else:
                new_frag_counts.append(c)
                new_frag_sizes.append(s)
        return (np.array(new_frag_sizes), np.array(new_frag_counts))
    return (frag_sizes, frag_counts)


def get_frag_data(frag_hist_file, barcode_blacklist_file=None,
                  min_frag_size=1000, frag_size_prc_cutoff=SV_FRAG_SIZE_PRC):

    """Get the data structures required for fragment-based sv-calling.
    Args:
    - frag_hist_file: histogram json like the one produced by REPORT_SINGLE_PARTITION
    Keys are bins (fragment sizes), values are counts.
    - barcode_blacklist_file: tsv file containing coalescent barcodes.
    - min_frag_size: ignore all fragments smaller than this

    Return value:
    Tuple (frag_sizes, frag_counts, max_frag, blacklist_barcodes).
    - blacklist_barcodes: Set of blacklisted barcodes.
    - max_frag: the frag_size_prc_cutoff-th percentile of fragment sizes
    """

    with open(frag_hist_file, 'r') as f:
        frag_hist = json.load(f)

    frag_sizes = np.array([int(k) for k in frag_hist.keys()])
    frag_counts = np.array([frag_hist[str(k)] for k in frag_sizes])
    sort_idx = np.argsort(frag_sizes)
    frag_sizes = frag_sizes[sort_idx]
    frag_counts = frag_counts[sort_idx]
    sel = frag_sizes >= min_frag_size
    frag_sizes = frag_sizes[sel]
    frag_counts = frag_counts[sel]
    frag_size_prc = np.where(np.cumsum(frag_counts) > frag_size_prc_cutoff * np.sum(frag_counts))[0]
    if len(frag_size_prc) > 0:
        max_frag = frag_sizes[frag_size_prc[0]]
    else:
        max_frag = None

    if len(frag_sizes) > 0 and np.max(frag_sizes) < MAX_FRAG_SIZE:
        frag_sizes = np.concatenate((frag_sizes, np.array([MAX_FRAG_SIZE])))
        frag_counts = np.concatenate((frag_counts, np.array([int(0.01 * np.sum(frag_counts))])))

    frag_sizes, frag_counts = merge_freq_bins(frag_sizes, frag_counts, 1e-3)

    if not barcode_blacklist_file is None:
        tmp_df = pd.read_csv(barcode_blacklist_file, sep = '\t', index_col = False)
        blacklist_barcodes = set(tmp_df.bc)
    else:
        blacklist_barcodes = set([])

    return (frag_sizes, frag_counts, max_frag, blacklist_barcodes)


def off_target_amp_corr_factor(target_regions, prob_off_target, genome_size=3e9):
    """Correction factor for amp rate (alpha) in off target regions.
    Args:
    - target_regions: region map with target regions per chrom.
    - prob_off_target: fraction of reads off target
    """

    # Fraction of genome on target
    genome_target_ratio = np.sum([np.sum([a[1] - a[0] for a in r]) for r in target_regions.values()]) / float(genome_size)
    # Fraction of reads on target
    b = 1 - prob_off_target
    corr_factor = (genome_target_ratio / float(1 - genome_target_ratio)) * (float(1 - b) / b)
    return corr_factor


def get_max_bc_ov(target_regions, common_frags, grid_len=1000):
    """ Get region of maximum barcode overlap between the input molecules.

    If the molecules overlap targets, get the target with the maximum barcode
    overlap. Otherwise get the grid_len-sized window with the max barcode overlap.
    """
    bc_ov = defaultdict(int)
    min_start = np.min(common_frags.start_pos_left)
    max_stop = np.max(common_frags.end_pos_left)
    starts = np.arange(min_start, max_stop, grid_len)
    left_regions = tk_regions.Regions(list(zip(starts, starts + grid_len)))

    min_start = np.min(common_frags.start_pos_right)
    max_stop = np.max(common_frags.end_pos_right)
    starts = np.arange(min_start, max_stop, grid_len)
    right_regions = tk_regions.Regions(list(zip(starts, starts + grid_len)))

    for _, row in common_frags.iterrows():
        if row.chrom_left in target_regions:
            # Get all the target regions overlapped by the fragment
            left_targets = target_regions[row.chrom_left].overlapping_regions(row.start_pos_left, row.end_pos_left)
        else:
            left_targets = []
        left_targets.extend(left_regions.overlapping_regions(row.start_pos_left, row.end_pos_left))

        if row.chrom_right in target_regions:
            right_targets = target_regions[row.chrom_right].overlapping_regions(row.start_pos_right, row.end_pos_right)
        else:
            right_targets = []
        right_targets.extend(right_regions.overlapping_regions(row.start_pos_right, row.end_pos_right))

        for (left_start, left_end), (right_start, right_end) in product(left_targets, right_targets):
            bc_ov[(left_start, left_end, right_start, right_end)] += 1

    max_bc_ov = np.max(np.array(bc_ov.values()))
    for region, region_ov in sorted(bc_ov.iteritems()):
        if region_ov == max_bc_ov:
            return (region, region_ov)


def get_orient(left_pos, right_pos, common_frags):
    """Infer orientation of barcode ovelap signal."""

    counts = [0, 0, 0, 0]
    for _, row in common_frags.iterrows():
        a = abs(row.start_pos_left - left_pos) > abs(row.end_pos_left - left_pos)
        b = abs(row.start_pos_right - right_pos) > abs(row.end_pos_right - right_pos)
        counts[a + 2 * b] += 1

    orient_arr = [(sv_call.DOWNSTREAM, sv_call.DOWNSTREAM),
                  (sv_call.UPSTREAM, sv_call.DOWNSTREAM),
                  (sv_call.DOWNSTREAM, sv_call.UPSTREAM),
                  (sv_call.UPSTREAM, sv_call.UPSTREAM)]
    return orient_arr[np.argmax(counts)]


def safe_logaddexp(arr):
    if len(arr) > 1:
        return tk_stats.logaddexp(arr)
    if len(arr) == 1:
        return arr[0]
    return MIN_LOG_PROB


def emp_logpdf(arr, sizes, counts):
    idx = np.array([bisect(sizes, x) for x in arr])
    res = np.zeros((len(arr),))
    in_range = np.logical_and(arr >= sizes[0], arr < sizes[-1])
    res[np.logical_not(in_range)] = MIN_LOG_PROB
    res[in_range] = np.log(counts[idx[in_range]]) / float(np.sum(counts))
    return res


def frag_size_logpmf(min_frag_len, frag_sizes, frag_counts):
    """Returns a tuple (P, S) where:
    P is an array with log-probability of length being L for all L in 
    frag_sizes that are > min_frag_len
    S[i] is the length corresponding to P[i]"""
    if min_frag_len >= frag_sizes[-1]:
        return np.array([MIN_LOG_PROB]), np.array([min_frag_len])
    idx = bisect(frag_sizes, min_frag_len)
    return (np.log(frag_counts[idx:] / float(np.sum(frag_counts))), frag_sizes[idx:])


def safe_poisson_logpmf(n, alpha):
    res = np.zeros((len(n), ))
    if alpha == 0:
        res[n > 0] = MIN_LOG_PROB
    else:
        res = poisson.logpmf(n, alpha)
    return res


def safe_poisson_logpmf2(n, alphas):
    res = np.zeros((len(alphas), ))
    res[alphas > 0] = poisson.logpmf(n, alphas[alphas > 0])
    if n > 0:
        res[alphas == 0] = MIN_LOG_PROB
    return res


def read_generation_pmf(n, high_alpha, low_alpha, high_len, low_len):
    """Computes the log-probability of generating n reads from a region of
    length high_len + low_len, of which high_len have an amp rate of high_alpha
    and low_len have an amp rate of low_alpha."""

    assert(high_len >= 0)
    assert(low_len >= 0)
    assert(low_len + high_len > 0)
    if low_len == 0:
        return safe_poisson_logpmf(np.array([n]), high_alpha * high_len)[0]
    if high_len == 0:
        return safe_poisson_logpmf(np.array([n]), low_alpha * low_len)[0]
    # probability of getting 0..n reads from the high yield regions
    p_high = safe_poisson_logpmf(np.arange(n + 1), high_alpha * high_len)
    assert(np.all(np.logical_not(np.isnan(p_high))))
    # probability of getting n..0 reads from the low yield regions
    p_low = safe_poisson_logpmf(n - np.arange(n + 1), low_alpha * low_len)
    assert(np.all(np.logical_not(np.isnan(p_low))))
    # probs[i] is the probability of getting i reads from the high yield regions
    # and n - i reads from the low yield regions
    probs = p_high + p_low
    # sum over i
    return safe_logaddexp(probs)


def read_frag_hist(frag_hist_file, min_frag_size, adjust=0.01):
    with open(frag_hist_file, 'r') as f:
        frag_hist = json.load(f)

    frag_sizes = np.array([int(k) for k in frag_hist.keys()])
    frag_counts = np.array([frag_hist[str(k)] for k in frag_sizes])
    sort_idx = np.argsort(frag_sizes)
    frag_sizes = frag_sizes[sort_idx]
    frag_counts = frag_counts[sort_idx]
    sel = frag_sizes >= min_frag_size
    frag_sizes = frag_sizes[sel]
    frag_counts = frag_counts[sel]

    if adjust > 0 and len(frag_sizes) > 0 and np.max(frag_sizes) < MAX_FRAG_SIZE:
        frag_sizes = np.concatenate((frag_sizes, np.array([MAX_FRAG_SIZE])))
        frag_counts = np.concatenate((frag_counts, np.array([int(adjust * np.sum(frag_counts))])))

    # Handle edge case of no long fragments -- just supply dummy values
    if len(frag_sizes) == 0:
        frag_sizes = np.array([1000, 2000])
        frag_counts = np.array([100, 50])

    frag_sizes, frag_counts = merge_freq_bins(frag_sizes, frag_counts, 1e-3)
    return frag_sizes, frag_counts


def get_cov_bases(target_coverage, cov_bin, start, stop):
    """Number of bases of [start, stop) covered by the targets.
    Args:
    - target_coverage: an array with the number of bases covered in bins
    [0, cov_bin), [cov_bin + 1, 2 * cov_bin).
    - cov_bin: length of bins for reporting target_coverage
    - start, stop
    """
    first_bin = int(np.floor(start / float(cov_bin)))
    if first_bin >= len(target_coverage):
        return 0
    last_bin = min(len(target_coverage) - 1, int(np.floor((stop - 1) / float(cov_bin))))
    if first_bin > 0:
        prev_coverage = target_coverage[first_bin - 1]
    else:
        prev_coverage = 0
    res_coverage = target_coverage[last_bin] - prev_coverage

    # Assume that the coverage of each bin is evenly distributed across the bin and remove
    # fractional parts of the coverage of the first and last bins.
    first_bin_coverage = target_coverage[first_bin] - prev_coverage
    res_coverage -= (start - cov_bin * first_bin) * first_bin_coverage / float(cov_bin)
    last_bin_coverage = target_coverage[last_bin] - target_coverage[last_bin - 1]
    res_coverage -= max(0, cov_bin * (last_bin + 1) - stop) * last_bin_coverage / float(cov_bin)
    return min(int(res_coverage), stop - start)


class FragModel(object):
    """Class for evaluating SVs based on molecule-level information."""

    def __init__(self, frag_sizes, frag_counts, blacklist_barcodes,
                 target_coverage, cov_bin=1000, corr_factor=1.0,
                 genome_size=3e9, target_regions=None, alpha=1e-2, p_ov_mol=1e-10):
        self.frag_sizes = frag_sizes
        self.frag_counts = frag_counts
        self.blacklist_barcodes = blacklist_barcodes
        self.target_coverage = target_coverage
        self.corr_factor = corr_factor
        self.genome_size = genome_size
        self.target_regions = target_regions
        self.alpha = alpha
        self.p_ov_mol = np.log(p_ov_mol)
        self.cov_bin = cov_bin
        self.prob_store = {}


    def eval_sv(self, frags1, frags2, locus1, locus2, break_ext=1000,
                genome_size=3e9, min_dist=0, min_bc_support=2,
                ps1=None, ps2=None, phase_set_dict1=None, phase_set_dict2=None,
                grid_len=1000):

        if self.frag_sizes is None or len(self.frag_sizes) == 0 or len(self.frag_counts) == 0:
            return None

        chrom1, s1, e1 = locus1
        chrom2, s2, e2 = locus2

        def get_frag_lr(row):
            nr1 = row['num_reads_left']
            nr2 = row['num_reads_right']
            start_pos1 = row['start_pos_left']
            start_pos2 = row['start_pos_right']
            end_pos1 = row['end_pos_left']
            end_pos2 = row['end_pos_right']
            left_idx = (row['bc'], row['start_pos_left'], row['end_pos_left'])
            right_idx = (row['bc'], row['start_pos_right'], row['end_pos_right'])
            store1 = self.prob_store.get(left_idx, None)
            store2 = self.prob_store.get(right_idx, None)

            res = self.lr_target(nr1, nr2, chrom1, start_pos1, end_pos1,
                                 chrom2, start_pos2, end_pos2,
                                 min_dist=min_dist, store1=store1, store2=store2)

            out_lr, _, _, new_store1, new_store2, _, _, _, _ = res
            if left_idx not in self.prob_store and not new_store1 is None:
                self.prob_store[left_idx] = new_store1
            if right_idx not in self.prob_store and not new_store2 is None:
                self.prob_store[right_idx] = new_store2
            return out_lr

        # Get common fragments between the left and right loci
        merged_frags = pd.merge(frags1, frags2, on='bc', suffixes=['_left', '_right'], how='outer')

        common_bcs = set(frags1.bc).intersection(set(frags2.bc))
        if not common_bcs:
            return None

        #print >> sys.stderr, 'common bcs', len(common_bcs)

        # Get the target region with the maximum barcode overlap
        (start1, stop1, start2, stop2), _ = get_max_bc_ov(self.target_regions, merged_frags, grid_len=grid_len)

        # Infer the orientation of the signal
        orient = get_orient((start1 + stop1)/2, (start2 + stop2)/2, merged_frags)
        #print >> sys.stderr, 'Inferred location and orientation', start1, start2, orient

        # Whether the barcode supports the event. 1 means supporting, -1 means
        # opposing and 0 is no vote
        bc_support = np.zeros((len(merged_frags), ), dtype=np.int)
        hap_probs1 = np.ones((len(merged_frags), 2)) * 0.5
        hap_probs2 = np.ones((len(merged_frags), 2)) * 0.5
        bc_list = np.array(merged_frags.bc)

        for pidx, (_, row) in enumerate(merged_frags.iterrows()):
            if pd.isnull(row.start_pos_right):
                # Fragment spans the breakpoint
                if row.start_pos_left < start1 and row.end_pos_left > stop1:
                    bc_support[pidx] = -1
                if row.bc in phase_set_dict1:
                    hap_probs = phase_set_dict1[row.bc][2]
                    hap_probs1[pidx, 0] = hap_probs[0]
                    hap_probs1[pidx, 1] = hap_probs[1]
            elif pd.isnull(row.start_pos_left):
                if row.start_pos_right < start2 and row.end_pos_right > stop2:
                    bc_support[pidx] = -1
                if row.bc in phase_set_dict2:
                    hap_probs = phase_set_dict2[row.bc][2]
                    hap_probs2[pidx, 0] = hap_probs[0]
                    hap_probs2[pidx, 1] = hap_probs[1]
            else:
                if chrom1 != chrom2 or start2 - stop1 > MAX_FRAG_SIZE:
                    bc_support[pidx] = 1
                else:
                    if orient[0] == sv_call.UPSTREAM and row.end_pos_left - stop1 > MAX_DIST_FROM_CAND:
                        # signal extends in the wrong direction
                        bc_support[pidx] = -1
                    elif orient[0] == sv_call.DOWNSTREAM and start1 - row.start_pos_left > MAX_DIST_FROM_CAND:
                        bc_support[pidx] = -1
                    elif orient[1] == sv_call.UPSTREAM and row.end_pos_right - stop2 > MAX_DIST_FROM_CAND:
                        bc_support[pidx] = -1
                    elif orient[1] == sv_call.DOWNSTREAM and start2 - row.start_pos_right > MAX_DIST_FROM_CAND:
                        bc_support[pidx] = -1
                    else:
                        lr = get_frag_lr(row)
                        bc_support[pidx] = 1 if lr > 0 else (0 if lr == 0 else -1)

                if row.bc in phase_set_dict1:
                    hap_probs = phase_set_dict1[row.bc][2]
                    hap_probs1[pidx, 0] = hap_probs[0]
                    hap_probs1[pidx, 1] = hap_probs[1]
                if row.bc in phase_set_dict2:
                    hap_probs = phase_set_dict2[row.bc][2]
                    hap_probs2[pidx, 0] = hap_probs[0]
                    hap_probs2[pidx, 1] = hap_probs[1]

        # Support from each haplotype on the left breakpoint
        left_support = np.array([np.sum(hap_probs1[bc_support == 1, 0] > 0.99),
                                 np.sum(hap_probs1[bc_support == 1, 1] > 0.99)])
        right_support = np.array([np.sum(hap_probs2[bc_support == 1, 0] > 0.99),
                                  np.sum(hap_probs2[bc_support == 1, 1] > 0.99)])

        if np.sum(left_support >= min_bc_support) > 1 or np.sum(right_support >= min_bc_support) > 1:
            # There is support on both haplotypes. Call homozygous.
            haps = (sv_call.Haplotype.unknown, sv_call.Haplotype.unknown)
            zygosity = sv_call.Zygosity.hom
            support = np.sum(bc_support == 1)
            non_support = np.sum(bc_support == -1)
            correct_hap_frac = None
        else:
            zygosity = sv_call.Zygosity.het
            if ps1 == ps2:
                if np.max(left_support) >= min_bc_support:
                    # There is enough support on a single haplotype. Assign to that haplotype.
                    haps = (sv_call.Haplotype.hap1 if left_support[0] < left_support[1] else sv_call.Haplotype.hap0,
                            sv_call.Haplotype.hap1 if left_support[0] < left_support[1] else sv_call.Haplotype.hap0)
                    hap_probs = hap_probs1[:, int(haps[0].value)]
                else:
                    haps = (sv_call.Haplotype.unknown, sv_call.Haplotype.unknown)
                    hap_probs = None
                    correct_hap_frac = None
            else:
                if np.max(left_support) >= min_bc_support:
                    new_hap0 = sv_call.Haplotype.hap1 if left_support[0] < left_support[1] else sv_call.Haplotype.hap0
                else:
                    new_hap0 = sv_call.Haplotype.unknown

                if np.max(right_support) >= min_bc_support:
                    new_hap1 = sv_call.Haplotype.hap1 if right_support[0] < right_support[1] else sv_call.Haplotype.hap0
                else:
                    new_hap1 = sv_call.Haplotype.unknown

                haps = (new_hap0, new_hap1)

                if haps[0] == sv_call.Haplotype.unknown or haps[1] == sv_call.Haplotype.unknown:
                    hap_probs = None
                else:
                    hap_probs = np.minimum(hap_probs1[:, int(haps[0].value)], hap_probs2[:, int(haps[1].value)])

            if hap_probs is None:
                support = np.sum(bc_support == 1)
                non_support = np.sum(bc_support == -1)
                correct_hap_frac = None
            else:
                support = int(np.ceil(np.sum(hap_probs[bc_support == 1])))
                non_support = int(np.ceil(np.sum(hap_probs[bc_support == -1])))
                correct_hap_frac = tk_stats.robust_divide(support, np.sum(bc_support == 1))

        info = {}
        info['PS1'] = ps1
        info['PS2'] = ps2
        info['HAP_ALLELIC_FRAC'] = tk_stats.robust_divide(support, support + non_support)
        info['FRAC_HAP_SUPPORT'] = correct_hap_frac
        info['ALLELIC_FRAC'] = np.mean(bc_support == 1)

        info['BCS'] = ','.join(list(bc_list[bc_support == 1]))

        sv_type = sv_call._SvType('UNK', (orient[0], orient[1]))

        res = sv_call.SvCall(sv_type, break1=(start1, stop1),
                             break2=(start2, stop2),
                             chrom1=chrom1, chrom2=chrom2, qual=support,
                             zygosity=zygosity, hap1=haps[0], hap2=haps[1], info=info)
        #print >> sys.stderr, res

        return res


    def log_prob_frag_target(self, nr, start_pos, end_pos, target_coverage, mean_target_ratio):
        """ Compute P(nr reads in [start_pos, end_pos); alpha, alpha_corr_factor, mean_target_ratio)
        """
        assert(nr >= 2)
        obs_len = end_pos - start_pos
        prob_true_len, true_len = frag_size_logpmf(obs_len, self.frag_sizes, self.frag_counts)

        high_len = get_cov_bases(target_coverage, self.cov_bin, start_pos, end_pos)
        low_len = obs_len - high_len
        prob_reads1 = read_generation_pmf(nr - 2, self.alpha, self.corr_factor * self.alpha, high_len, low_len)
        prob_reads_off = safe_poisson_logpmf2(0, self.alpha * (true_len - obs_len) * mean_target_ratio)
        prob_reads_off += safe_poisson_logpmf2(0, self.corr_factor * self.alpha * (true_len - obs_len) * (1 - mean_target_ratio))

        prob_reads = prob_reads1 + 2 * np.log(self.alpha) + prob_reads_off + np.log(true_len - obs_len)
        p = np.maximum(MIN_LOG_PROB, safe_logaddexp(prob_reads + prob_true_len))
        return p, prob_reads1


    def mean_target_by_offset(self, start_pos, end_pos, target_coverage, len_bins):
        """Average fraction of on-target reads for a fragment with observed area [start_pos, end_pos)
        and has true length in the distribution given by len_bins.
        """
        # Get the possible bins of the true fragment length.
        obs_len = end_pos - start_pos
        idx = min(bisect(len_bins, obs_len), len(len_bins) - 1)
        len_bins = np.array(len_bins[idx:] - obs_len)

        def safe_div(a, b):
            if b == 0:
                return 0
            return a / float(b)

        average_cov = np.zeros((len(len_bins), ))
        average_cov_up = np.zeros((len(len_bins), ))
        average_cov_down = np.zeros((len(len_bins), ))
        for i, bl in enumerate(len_bins):
            cu = get_cov_bases(target_coverage, self.cov_bin, max(0, start_pos - bl), start_pos)
            cd = get_cov_bases(target_coverage, self.cov_bin, end_pos, end_pos + bl)
            average_cov_up[i] = safe_div(cu, start_pos - max(0, start_pos - bl))
            average_cov_down[i] = safe_div(cd, bl)
            average_cov[i] = safe_div(cu + cd, bl + start_pos - max(0, start_pos - bl))

        return average_cov, average_cov_up, average_cov_down, len_bins


    @staticmethod
    def create_prob_store_target(target_ratio, target_ratio_up, target_ratio_down,
                                 len_bins, on_prob, read_prob):
        return {'target_ratio':target_ratio, 'target_ratio_up':target_ratio_up,
                'target_ratio_down':target_ratio_down, 'len_bins':len_bins,
                'on_prob':on_prob, 'read_prob':read_prob}


    @staticmethod
    def extract_prob_store_target(store):
        return (store['target_ratio'], store['target_ratio_up'], store['target_ratio_down'],
                store['len_bins'], store['on_prob'], store['read_prob'])


    def prob_no_reads(self, tot_len, target_ratio):
        p = safe_poisson_logpmf2(0, self.alpha * tot_len * target_ratio)
        p += safe_poisson_logpmf2(0, self.corr_factor * self.alpha * tot_len * (1 - target_ratio))
        return p


    def estimate_extent(self, chrom, start_pos, end_pos, ext_bins, ext_prob=0.99):
        """Finds the maximum extension to the left and right of [start_pos, end_pos) such
        that the probability of observing no reads in the extension is > ext_prob.
        In other words, finds how far the true fragment length could have extended."""

        res = self.mean_target_by_offset(start_pos, end_pos, self.target_coverage.get(chrom, []),
                                         np.array(ext_bins) + (end_pos - start_pos))
        _, target_ratio_up, target_ratio_down, len_bins = res
        no_read_by_len = self.prob_no_reads(len_bins, target_ratio_down)
        pos = np.where(no_read_by_len > np.log(ext_prob))[0]
        if len(pos) > 0:
            ext_right = len_bins[pos[-1]]
            prob_ext_right = no_read_by_len[pos[-1]]
        else:
            ext_right = 0
            prob_ext_right = 0

        no_read_by_len = self.prob_no_reads(len_bins, target_ratio_up)
        pos = np.where(np.logical_and(len_bins < start_pos, no_read_by_len > np.log(ext_prob)))[0]
        if len(pos) > 0:
            ext_left = len_bins[pos[-1]]
            prob_ext_left = no_read_by_len[pos[-1]]
        else:
            ext_left = 0
            prob_ext_left = 0

        return (ext_left, ext_right, prob_ext_left, prob_ext_right)


    def merge_lr_target(self, nr1, nr2, chrom1, start_pos1, end_pos1,
                        chrom2, start_pos2, end_pos2, min_dist=0,
                        store1=None, store2=None):
        """Computes the log-probabilities that two fragments with the same barcode come from the same or different molecules.
        Args:
        - As in lr_target. It is assumed that the second fragment always has larger genomic coordinates than the first.
        """

        if store1 is None:
            # Compute average fraction of on-target reads for the different possible true fragment sizes
            # (averaged across offsets with respect to the original region)
            res = self.mean_target_by_offset(start_pos1, end_pos1, self.target_coverage.get(chrom1, []), self.frag_sizes)
            target_ratio1, target_ratio_up1, target_ratio_down1, len_bins1 = res

            p1, preads1 = self.log_prob_frag_target(nr1, start_pos1, end_pos1,
                                                    self.target_coverage.get(chrom1, []), target_ratio1)
            store1 = self.create_prob_store_target(target_ratio1, target_ratio_up1, target_ratio_down1, len_bins1, p1, preads1)
        else:
            (target_ratio1, target_ratio_up1, target_ratio_down1, len_bins1, p1, preads1) = self.extract_prob_store_target(store1)

        if store2 is None:
            res = self.mean_target_by_offset(start_pos2, end_pos2, self.target_coverage.get(chrom2, []), self.frag_sizes)
            target_ratio2, target_ratio_up2, target_ratio_down2, len_bins2 = res

            p2, preads2 = self.log_prob_frag_target(nr2, start_pos2, end_pos2,
                                                    self.target_coverage.get(chrom2, []), target_ratio2)
            store2 = self.create_prob_store_target(target_ratio2, target_ratio_up2, target_ratio_down2, len_bins2, p2, preads2)
        else:
            (target_ratio2, target_ratio_up2, target_ratio_down2, len_bins2, p2, preads2) = self.extract_prob_store_target(store2)

        p_two_frags = np.maximum(MIN_LOG_PROB, p1 + p2 + self.p_ov_mol)

        tot_obs = end_pos2 - start_pos1
        if chrom1 == chrom2 and tot_obs <= self.frag_sizes[-1]:
            high_len = get_cov_bases(self.target_coverage.get(chrom1, []), self.cov_bin, end_pos1, start_pos2)
            low_len = start_pos2 - end_pos1 - high_len
            prob_reads_middle = read_generation_pmf(0, self.alpha, self.corr_factor * self.alpha, high_len, low_len)

            prob_true_len, true_len = frag_size_logpmf(tot_obs, self.frag_sizes, self.frag_counts)

            res = self.mean_target_by_offset(start_pos1, end_pos2,
                                             self.target_coverage.get(chrom1, []), self.frag_sizes)
            target_ratio_middle, _, _, _ = res

            prob_reads_off = self.prob_no_reads(true_len - tot_obs, target_ratio_middle)
            prob_reads = preads1 + preads2 + 4 * np.log(self.alpha) + prob_reads_middle + prob_reads_off + np.log(true_len - tot_obs)

            p_one_frag = np.maximum(MIN_LOG_PROB, safe_logaddexp(prob_reads + prob_true_len) + np.log1p(-np.exp(self.p_ov_mol)))
        else:
            p_one_frag = MIN_LOG_PROB
            prob_reads_middle = MIN_LOG_PROB

        return (p_one_frag, p_two_frags, store1, store2, prob_reads_middle)


    def lr_target(self, nr1, nr2, chrom1, start_pos1, end_pos1,
                  chrom2, start_pos2, end_pos2, min_dist=0,
                  store1=None, store2=None, ext_prob=0.75):
        """Computes the log-LR of the two hypotheses that there was an SV vs there was no SV:
        logP(observed fragments | SV) - logP(observed fragments | no SV).
        - nr1, nr2: number of reads on each of the fragments.
        - (chrom1, start_pos1, end_pos1): Coordinates of the first fragment.
        - (chrom2, start_pos2, end_pos2): Coordinates of the second fragment. It is assumed that the second fragment always
          has larger genomic coordinates than the first.
        - alpha: reads/bp on targets
        - target_coverage: dictionary of binned target bases
        - cov_bin: bin size used for computing coverage.
        - alpha_corr_factor: correction factor for amp-rate. reads/bp off target are alpha * alpha_corr_factor
        - frag_sizes, frag_counts: histogram of fragment lengths
        - store1, store2: dictionaries
        - ext_prob: Determines how long the observed fragments will be extended in order to try to estimate
        the size of the true fragment under the assumption of an SV.

        Return values:
        (lr, store1, store2, ext_left1, ext_right1, ext_left2, ext_right2)
        - LR
        - updated store1, store2
        - extension of fragments to the left and right (ie. how far did the true fragment extended from
        the observed fragment)
        """

        assert(chrom1 != chrom2 or end_pos1 < start_pos2)

        res = self.merge_lr_target(nr1, nr2, chrom1, start_pos1, end_pos1,
                                   chrom2, start_pos2, end_pos2, min_dist=min_dist,
                                   store1=store1, store2=store2)
        p_one_frag, p_two_frags, store1, store2, _ = res
        preads1 = store1['read_prob']
        preads2 = store2['read_prob']

        if chrom1 == chrom2 and start_pos2 - end_pos1 <= min_dist:
            return (0.0, 0.0, 0.0, None, None, 0, 0, 0, 0)

        p_no_break = np.logaddexp(p_one_frag, p_two_frags)

        # Break
        obs_len = end_pos2 - start_pos2 + end_pos1 - start_pos1
        est_len = obs_len
        prob_reads_off = 0.0

        res = self.estimate_extent(chrom1, start_pos1, end_pos1,
                                   np.arange(1, 3 * TARGET_COV_BIN, TARGET_COV_BIN), ext_prob=ext_prob)
        ext_left1, ext_right1, _, prob_ext_right = res
        res = self.estimate_extent(chrom2, start_pos2, end_pos2,
                                   np.arange(1, 3 * TARGET_COV_BIN, TARGET_COV_BIN), ext_prob=ext_prob)
        ext_left2, ext_right2, prob_ext_left, _ = res
        est_len += ext_left2 + ext_right1
        prob_reads_off += prob_ext_left + prob_ext_right

        # This will ensure we don't get negative numbers in the computation of target_ratio_off2
        # below, while keeping the upstream and downstream vectors (off1 and off2) the same length.
        est_len = min(est_len, end_pos2)

        if est_len > self.frag_sizes[-1]:
            p_break_one_frag = MIN_LOG_PROB
        else:
            prob_true_len, true_len = frag_size_logpmf(est_len, self.frag_sizes, self.frag_counts)
            res = self.mean_target_by_offset(start_pos1, start_pos1 + est_len,
                                             self.target_coverage.get(chrom1, []), self.frag_sizes)
            target_ratio_off1, target_ratio_off_up1, _, _ = res
            res = self.mean_target_by_offset(end_pos2 - est_len, end_pos2,
                                             self.target_coverage.get(chrom2, []), self.frag_sizes)
            target_ratio_off2, _, target_ratio_off_down2, _ = res
            if chrom1 == chrom2:
                target_ratio = 0.5 * (target_ratio_off_up1 + target_ratio_off_down2)
            else:
                # Unclear where we should be looking to estimate the target fraction.
                target_ratio = 0.5 * (target_ratio_off1 + target_ratio_off2)

            prob_reads_off += self.prob_no_reads(true_len - est_len, target_ratio)
            prob_reads = preads1 + preads2 + 4 * np.log(self.alpha) + prob_reads_off + np.log(true_len - est_len)

            p_break_one_frag = np.maximum(MIN_LOG_PROB, safe_logaddexp(prob_reads + prob_true_len) + np.log1p(-np.exp(self.p_ov_mol)))
        p_break = np.logaddexp(p_break_one_frag, p_two_frags)

        lr = p_break - p_no_break

        return (lr, p_break, p_no_break, store1, store2, ext_left1, ext_right1, ext_left2, ext_right2)

