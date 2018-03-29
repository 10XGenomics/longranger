#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for SV calling
#

import numpy as np
import json
from bisect import bisect_left
import tenkit.pandas as pd
import tenkit.bio_io as tk_io
import tenkit.regions as tk_regions
from scipy.stats import binom
from collections import defaultdict

import longranger.sv.io as tk_sv_io
from longranger.sv.constants import MAX_FRAG_SIZE, MIN_LOG_PROB, MAX_PHASE_PROB

BINOM = 0
EXP_COUNT = 2
EXP_COUNT2 = 3
EXP_COUNT3 = 4
BINOM_EMP_BC_COUNT_BC_FREQ = 1
BINOM_LR_FRAG = 5

MAX_OUT_QUAL = 3000


def get_read_haps(readpairs):
    nhap0_r1 = np.sum(np.array([tk_io.get_read_haplotype(rp.read1) == 1 for rp in readpairs]))
    nhap1_r1 = np.sum(np.array([tk_io.get_read_haplotype(rp.read1) == 2 for rp in readpairs]))
    nhap0_r2 = np.sum(np.array([tk_io.get_read_haplotype(rp.read2) == 1 for rp in readpairs]))
    nhap1_r2 = np.sum(np.array([tk_io.get_read_haplotype(rp.read2) == 2 for rp in readpairs]))
    return (nhap0_r1, nhap1_r1, nhap0_r2, nhap1_r2)


def sort_and_merge(regions, extend = 0):
    """Sorts and merges list of tuples (chrom, start, stop).
    Args:
    - regions: list of tuples (chrom, start, stop).
    - extend: non-negative integer (if negative, will be set to 0). Extend regions to the left
    and right by that amount before sorting and merging.

    Return value:
    - list of tuples (chrom, start, stop).
    """
    extend = max(extend, 0)
    regions = [(r[0], max(0, r[1] - extend), r[2] + extend) for r in regions]
    regions = sorted(regions, key = lambda x:(x[0], x[1], x[2]))

    merged_regions = []
    prev_chrom = ''
    prev_start = -1
    prev_stop = -1
    for ridx, (chrom, start, stop) in enumerate(regions):
        if ridx > 0:
            if chrom != prev_chrom or start >= prev_stop:
                merged_regions.append((prev_chrom, prev_start, prev_stop))
                prev_chrom, prev_start, prev_stop = chrom, start, stop
            else:
                prev_stop = max(prev_stop, stop)
        else:
            prev_chrom, prev_start, prev_stop = chrom, start, stop
    if prev_chrom != '':
        merged_regions.append((prev_chrom, prev_start, prev_stop))

    return merged_regions


def log10_binom_pval(nov, bcs1, bcs2, nbcs):
    # logsf is log(1 - cdf)
    return binom.logsf(nov, bcs1, bcs2 / float(nbcs)) / np.log(10)


def log10_emp_pval(win_idx, bc_idx, n, freq):
    assert(len(win_idx) == len(bc_idx))
    if len(bc_idx) == 0:
        return np.zeros((1,))
    return np.bincount(win_idx, weights = binom.logsf(0, n, freq[bc_idx])) / np.log(10)


def pval_to_qual(log10_pval):
    return int(min(-10 * log10_pval, MAX_OUT_QUAL))


def get_nx_bcs(read_counts, n):
    n = min(n, 100.0) / 100.0
    nbcs = read_counts.size
    bc_rank = np.argsort(read_counts)[::-1]
    is_nx_bc = np.zeros((nbcs,), dtype = np.bool)
    if n > 0 and len(read_counts) > 0:
        read_cumsum = np.cumsum(read_counts[bc_rank])
        # Get index of first bc where sum of reads becomes > n% of reads
        nx_idx = np.where(read_cumsum >= n * read_counts.sum())[0][0]
        is_nx_bc[bc_rank[0:(nx_idx + 1)]] = True
    return (is_nx_bc, bc_rank)


def bc_map_from_inv(inv_bc_map, is_nx_bc):
    bc_map = {}
    for i, b in inv_bc_map.iteritems():
        if is_nx_bc[i]:
            bc_map[b] = i
    return bc_map


def  get_bcs_at_region(in_bam, chrom, start, end, min_mapq = 60, bc_map = None,
    other_chrom = None, other_start = None, other_end = None, read_to_bc = None):

    bc_list = {}

    for read in in_bam.fetch(str(chrom), int(start), int(end)):
        if read.mapq < min_mapq or read.is_duplicate or read.is_secondary:
            continue
        if read.pos < start or read.pos > end:
            continue

        if not other_start is None and not other_end is None:
            tag_names = [t[0] for t in read.tags]
            tag_vals = [t[1] for t in read.tags]
            if 'XA' in tag_names:
                idx = np.where(np.array(tag_names) == 'XA')[0]
                chrom = tag_vals[idx].split(',')[0]
                pos = int(tag_vals[idx].split(',')[1].strip('+-'))
                if chrom == other_chrom and pos >= other_start and pos <= other_end:
                    continue
        if read_to_bc is None:
            bc = tk_io.get_read_barcode(read)
        else:
            bc = read_to_bc[read.qname]

        if not(bc is None or bc == '') and (bc_map is None or bc in bc_map):
            if not bc in bc_list:
                bc_list[bc] = []
            bc_list[bc].append(read.qname)

    for bc in bc_list.keys():
        bc_list[bc] = list(set(bc_list[bc]))
    return bc_list


def loci_to_region_map(loci, merge = False):
    if loci is None:
        return None
    regions = {}
    for chrom, starts, stops in loci:
        if not chrom in regions:
            regions[chrom] = []
        regions[chrom].extend(zip(starts, stops))
    for chrom, region_list in regions.iteritems():
        regions[chrom] = tk_regions.Regions(regions = regions[chrom], merge = merge)
    return regions


def loci_to_named_region_map(loci, singletons = False):
    if loci is None:
        return None
    regions = defaultdict(list)
    for chrom, starts, stops, names in loci:
        if singletons:
            regions[chrom].append((starts, stops, names))
        else:
            regions[chrom].extend(zip(starts, stops, names))
    for chrom, region_list in regions.iteritems():
        regions[chrom] = tk_regions.NamedRegions(regions = regions[chrom])
    return regions


def bed_to_region_map(bed_file, merge = False, extend = 0):
    loci = []
    if bed_file is not None:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track'):
                    continue
                fields = line.strip().split('\t')
                loci.append((fields[0], [max(0, int(fields[1]) - extend)], [int(fields[2]) + extend]))
    return loci_to_region_map(loci, merge = merge)


def bedpe_to_region_map(bed_file, merge = False, extend = 0):
    loci = []
    if bed_file is not None:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                loci.append((fields[0], [max(0, int(fields[1]) - extend)], [int(fields[2]) + extend]))
                loci.append((fields[3], [max(0, int(fields[4]) - extend)], [int(fields[5]) + extend]))
    return loci_to_region_map(loci, merge = merge)


def bedpe_df_to_named_region_map(df):
    regions_map = defaultdict(list)
    for _, row in df.iterrows():
        assert row.chrom1 == row.chrom2
        regions_map[row.chrom1].append((row.start1, row.stop2, row['name']))

    for chrom, regions in regions_map.iteritems():
        regions_map[chrom] = tk_regions.NamedRegions(regions)
    return regions_map


def merge_region_maps(map1, map2):
    chroms = list(set(map1.keys()).union(set(map2.keys())))
    merged_map = {}
    for c in chroms:
        if c in map1:
            map1[c].merge(map2.get(c, tk_regions.Regions([])))
            merged_map[c] = map1[c]
        else:
            merged_map[c] = map2[c]
    return merged_map


def strictly_overlapping_regions(regions, chrom, start, stop):
    ov = regions[chrom].overlapping_regions(start, stop)
    return [(max(r[0], start), min(r[1], stop)) for r in ov]


def get_region_overlaps(regions1, regions2):
    """Finds overlaps between pairs of regions.
    regions1: Dictionary of NamedRegions objects as returned by loci_to_named_region_map.
    regions2: Similar.
    Return value:
    A dict from names in regions1 to sets of names of overlapping regions in regions2.
    """
    mapping = {}
    for chrom, reg in regions1.iteritems():
        if chrom in regions2:
            for r in reg:
                names = regions2[chrom].overlapping_region_names(r.start, r.end)
                mapping[r.name] = set(names)
    return mapping


def region_coverage(regions, bin_len):
    """Coverage of regions in bins.
    Bins a locus (eg. a chromosome) into bins and computes the number of bases in each
    bin that are covered by the given regions.
    Args:
    - regions: A Regions object with non-overlapping regions.
    - bin_len: length of bins
    Return value:
    - cov: An array with the number of bases covered in each bin. cov[i] is the coverage
    of bases [i * bin_len, (i + 1) * bin_len). The last bin is determined by the maximum
    basepair covered by the regions (so no empty bins past the last region).
    """
    max_bin_idx = int(np.ceil(float(np.max([r[1] for r in regions])) / bin_len))
    cov = np.zeros((max_bin_idx, ), dtype = np.int)

    for r in regions:
        covered_bins = np.arange(int(np.floor(r[0] / float(bin_len))),
            int(np.floor((r[1] - 1) / float(bin_len))) + 1)
        if len(covered_bins) > 1:
            covered_bases = np.ones((len(covered_bins), ), dtype = np.int) * int(bin_len)
            covered_bases[0] = (covered_bins[0] + 1) * bin_len - r[0]
            covered_bases[-1] = int(r[1] - covered_bins[1] * bin_len)
        else:
            covered_bases = np.array([r[1] - r[0]], dtype = np.int)
        cov[covered_bins] += covered_bases
    return cov


def region_cum_coverage_map(region_map, bin_len):
    cov_map = {}
    for c, regions in region_map.iteritems():
        cov_map[c] = np.cumsum(region_coverage(regions, bin_len))
    return cov_map


def overlap_breaks(pred_loci, true_loci, min_rel_overlap = 0.5):
    if not isinstance(pred_loci, pd.DataFrame):
        pred_df = tk_sv_io.read_sv_bedpe_to_df(pred_loci)
    else:
        pred_df = pred_loci

    if not isinstance(true_loci, pd.DataFrame):
        gt_df = tk_sv_io.read_sv_bedpe_to_df(true_loci)
    else:
        gt_df = true_loci

    gt_df.index = gt_df['name']
    gt_map = bedpe_df_to_named_region_map(gt_df)

    pred_to_matching_true = defaultdict(set)
    true_to_matching_pred = defaultdict(set)

    for _, row in pred_df.iterrows():
        if not row.chrom1 in gt_map:
            continue
        matches = gt_map[row.chrom1].overlapping_region_names(row.start1, row.stop2)
        this_len = row.stop2 - row.start1
        if len(matches) > 0:
            match_df = gt_df.loc[list(matches)]
            lengths = np.array(match_df.stop2 - match_df.start1, dtype=np.float)
            ov = np.minimum(match_df.stop2, row.stop2) - np.maximum(match_df.start1, row.start1)
            good_matches = np.logical_and(ov / float(this_len) > min_rel_overlap,
                                          ov / lengths > min_rel_overlap)

            if np.any(good_matches):
                true_names = match_df[good_matches]['name']
                for matching_true in true_names:
                    true_to_matching_pred[matching_true].add(row['name'])
                    pred_to_matching_true[row['name']].add(matching_true)

    return (pred_to_matching_true, true_to_matching_pred, {})


def compare_breaks(pred_loci, true_loci = None, max_dist = 100, window_loci = None):
    """
    pred_file: BEDPE file with sv calls or pandas DataFrame as returned by tk_sv_io.read_sv_bedpe_to_df
    true_file: BEDPE file with ground truth variants (or other set of variants against
        which pred_file will be compared)
    max_dist: maximum distance between a true and a predicted breakpoint in order to
        say that they overlap
    window_loci: list of tuples (chrom, starts, stops), where chrom is a chromosome name and
        starts/stops are lists/arrays of start and ending positions. If this is provided,
        true svs that completely fall within such a locus will be marked as "filtered" (i.e.
        not detectable). For example, these can be the windows used for detecting svs.
        An SV that lies completely within a single window cannot be detected.
    """

    if true_loci is None or pred_loci is None:
        return ({}, {}, set([]))

    ###### Read predicted breakpoints and extend them by max_dist
    pred_breaks1 = []
    pred_breaks2 = []
    if not isinstance(pred_loci, pd.DataFrame):
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(pred_loci)
    else:
        bedpe_df = pred_loci
    for n, row in bedpe_df.iterrows():
        break1 = (row.chrom1, max(0, row.start1 - max_dist), row.stop1 + max_dist, row['name'])
        break2 = (row.chrom2, max(0, row.start2 - max_dist), row.stop2 + max_dist, row['name'])
        if break1 > break2:
            break1, break2 = break2, break1
        pred_breaks1.append(break1)
        pred_breaks2.append(break2)

    pred_regions1 = loci_to_named_region_map(pred_breaks1, singletons = True)
    pred_regions2 = loci_to_named_region_map(pred_breaks2, singletons = True)

    ###### Read provided loci
    regions = loci_to_region_map(window_loci)

    ###### Read true svs
    filtered_svs = set([]) # set of true svs that are non-detectable
    if not isinstance(true_loci, pd.DataFrame):
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(true_loci)
    else:
        bedpe_df = true_loci
    true_breaks1 = []
    true_breaks2 = []
    for n, row in bedpe_df.iterrows():
        name = row['name']
        chrom1, start1, stop1 = row.chrom1, row.start1, row.stop1
        chrom2, start2, stop2 = row.chrom2, row.start2, row.stop2
        break1 = (chrom1, start1, stop1, name)
        break2 = (chrom2, start2, stop2, name)
        if break1 > break2:
            break1, break2 = break2, break1

        is_filtered = False
        if not regions is None and chrom1 == chrom2 and chrom1 in regions:
            # SV is filtered if both its breakpoints are on the same window
            ovs1 = regions[chrom1].overlapping_regions(start1, stop1)
            ovs2 = regions[chrom2].overlapping_regions(start2, stop2)
            if len(set(ovs1).intersection(set(ovs2))) > 0:
                is_filtered = True
                filtered_svs.add(name)
        if not is_filtered:
            true_breaks1.append(break1)
            true_breaks2.append(break2)

    true_regions1 = loci_to_named_region_map(true_breaks1, singletons = True)
    true_regions2 = loci_to_named_region_map(true_breaks2, singletons = True)

    ###### Get overlaps bewtween predicted and true breakpoints
    mapping_break1 = get_region_overlaps(pred_regions1, true_regions1)
    mapping_break2 = get_region_overlaps(pred_regions2, true_regions2)

    pred_to_matching_true = {}
    true_to_matching_pred = {}
    for pred_name, matched in mapping_break1.iteritems():
        if not pred_name in mapping_break2:
            # There was a match only for one of the breakpoints of this predicted sv.
            continue
        for true_name in matched:
            if true_name in mapping_break2[pred_name]:
                s1 = pred_to_matching_true.setdefault(pred_name, set([]))
                s1.add(true_name)
                s1 = true_to_matching_pred.setdefault(true_name, set([]))
                s1.add(pred_name)

    return (pred_to_matching_true, true_to_matching_pred, filtered_svs)


def cluster_loci(loci, win, max_range = np.inf):
    """
    Args:
    loci: list of tuples (chrom, start, stop, name)
    win: distance between loci to cluster together
    max_range: used to avoid overmerging chained loci. This is the maximum distance
    between a locus start and the beginning of the cluster (earliest start of loci in
    cluster).

    Return value:
    cluster_to_mem: dictionary from cluster index to a list of names of loci in the cluster.
    """

    win = max(win, 0)
    cluster = 0
    cluster_to_mem = {}
    cluster_to_mem[0] = []
    cluster_ranges = {}
    cluster_ranges[0] = {}
    mem_to_cluster = {}
    last_stop = 0
    first_start = np.inf
    loci = sorted(loci, key = lambda x: (x[0], x[1], x[2]))

    for idx, locus in enumerate(loci):
        name = locus[3]
        assert locus[2] >= locus[1], 'Invalid locus ' + str(locus)
        if idx > 0:
            prev_locus = loci[idx - 1]
            # Chromosome change or start of this locus is more than "win" aften end of previous CLUSTER.
            if locus[0] != prev_locus[0] or locus[1] - last_stop > win or locus[1] - first_start > max_range:
                cluster_ranges[cluster] = (prev_locus[0], first_start, last_stop)
                cluster += 1
                cluster_to_mem[cluster] = []
                last_stop = 0
                first_start = locus[1]
        else:
            first_start = locus[1]
        last_stop = max(locus[2], last_stop)
        cluster_to_mem[cluster].append(name)
        mem_to_cluster[name] = cluster
    if len(loci) > 0:
        cluster_ranges[cluster] = (loci[-1][0], first_start, last_stop)

    return cluster_to_mem, mem_to_cluster, cluster_ranges


def get_break_groups(bedpe_df, merge_win = 10000, max_range = np.inf):
    """A simplified version of merge_breaks"""

    if not isinstance(bedpe_df, pd.DataFrame):
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(bedpe_df)
    else:
        bedpe_df = pd.DataFrame(bedpe_df)

    breaks = []
    for i, (n, row) in enumerate(bedpe_df.iterrows()):
        breaks.append((row.chrom1, row.start1, row.stop1, (n, 1)))
        breaks.append((row.chrom2, row.start2, row.stop2, (n, 2)))
    _, mem_to_cluster, _ = cluster_loci(breaks, merge_win, max_range = max_range)

    cluster_pairs = defaultdict(list)
    for i, (n, row) in enumerate(bedpe_df.iterrows()):
        cluster_idx1 = mem_to_cluster[(n, 1)]
        cluster_idx2 = mem_to_cluster[(n, 2)]
        cluster_pairs[(cluster_idx1, cluster_idx2)].append(n)
    return cluster_pairs.values()


def get_break_groups_from_loci(loci, merge_win = 10000, max_range = np.inf):
        breaks = []
        for (c1, s1, e1, c2, s2, e2, n) in loci:
            breaks.append((c1, s1, e1, (n, 1)))
            breaks.append((c2, s2, e2, (n, 2)))
        _, mem_to_cluster, _ = cluster_loci(breaks, merge_win, max_range = max_range)

        cluster_pairs = defaultdict(list)
        for (c1, s1, e1, c2, s2, e2, n) in loci:
            cluster_idx1 = mem_to_cluster[(n, 1)]
            cluster_idx2 = mem_to_cluster[(n, 2)]
            cluster_pairs[(cluster_idx1, cluster_idx2)].append(n)
        return cluster_pairs.values()


def merge_breaks(bedpe_df, out_bedpe, merge_win = 10000, max_range = np.inf, max_nmates = np.inf,
                 cluster_qual_factor = 0.2):
    """Merges a set of SVs into a non-redundant set.
    Args:
    - bedpe_df: Either a bedpe file or a DataFrame like the one returned by
    tk_sv_io.read_sv_bedpe_to_df.
    - out_bedpe: Path to file where output will be written.
    - merge_win: Breakpoints will be merged if they are within this distance from
    each other. Two SVs will be merged if both their breakpoints can be merged.
    - max_range: See max_range field of cluster_loci.
    - max_nmates: Two extra info fields will be added to the output BEDPE, NMATES1,
    and NMATES2. NMATES1 is the number of mate breakpoints (after merging, so
    breakpoint clusters), of the first breakpoint of an SV.
    SVs whose breakpoints both exceed the max_nmates cutoff will not be included in the
    output.

    Return value:
    The output BEDPE.
    """
    if not isinstance(bedpe_df, pd.DataFrame):
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(bedpe_df)
    else:
        bedpe_df = pd.DataFrame(bedpe_df)
    breaks = []
    for i in range(bedpe_df.shape[0]):
        breaks.append((bedpe_df.iloc[i, 0], bedpe_df.iloc[i, 1], bedpe_df.iloc[i, 2], (bedpe_df.iloc[i, 6], 1)))
        breaks.append((bedpe_df.iloc[i, 3], bedpe_df.iloc[i, 4], bedpe_df.iloc[i, 5], (bedpe_df.iloc[i, 6], 2)))
    _, mem_to_cluster, _ = cluster_loci(breaks, merge_win, max_range = max_range)

    cluster_pairs = {}
    for i in range(bedpe_df.shape[0]):
        name = bedpe_df.iloc[i]['name']
        cluster_idx1 = mem_to_cluster[(name, 1)]
        cluster_idx2 = mem_to_cluster[(name, 2)]
        if not (cluster_idx1, cluster_idx2) in cluster_pairs:
            cluster_pairs[(cluster_idx1, cluster_idx2)] = [i]
        else:
            old_pair = cluster_pairs[(cluster_idx1, cluster_idx2)][0]
            # Make sure the old and the new pair have breaks on the same chromosomes
            assert(bedpe_df.iloc[old_pair, 0] == bedpe_df.iloc[i, 0])
            assert(bedpe_df.iloc[old_pair, 3] == bedpe_df.iloc[i, 3])
            cluster_pairs[(cluster_idx1, cluster_idx2)].append(i)

    new_cluster_pairs = {}
    cluster_dist_ratio = {}
    for p, pos_list in cluster_pairs.iteritems():
        pos_arr = np.array(pos_list)
        tmp_df = get_dataframe_loc(bedpe_df, pos_arr)
        quals = np.array(tmp_df.qual)
        best_call = pos_arr[np.argmax(quals)]
        close_calls = np.where(quals >= cluster_qual_factor * np.max(quals))[0]
        close_df = get_dataframe_loc(tmp_df, close_calls)

        same_chrom = bedpe_df.iloc[best_call]['chrom2'] == bedpe_df.iloc[best_call]['chrom1']
        min_break_dist = np.min(close_df.start2) - np.max(close_df.stop1)
        max_break_dist = bedpe_df.iloc[best_call]['start2'] - bedpe_df.iloc[best_call]['stop1']

        new_cluster_pairs[p] = best_call
        if not same_chrom or max_break_dist > MAX_FRAG_SIZE:
            cluster_dist_ratio[p] = '.'
        elif min_break_dist <= 0:
            cluster_dist_ratio[p] = float('NaN')
        else:
            cluster_dist_ratio[p] = float(max_break_dist) / min_break_dist

    cluster_pairs = new_cluster_pairs

    def clusters_close(i, j):
        chrom1, start1, stop1 = bedpe_df.iloc[i, 0], bedpe_df.iloc[i, 1], bedpe_df.iloc[i, 2]
        chrom2, start2, stop2 = bedpe_df.iloc[i, 3], bedpe_df.iloc[i, 4], bedpe_df.iloc[i, 5]
        next_chrom1, next_start1, next_stop1 = bedpe_df.iloc[j, 0], bedpe_df.iloc[j, 1], bedpe_df.iloc[j, 2]
        next_chrom2, next_start2, next_stop2 = bedpe_df.iloc[j, 3], bedpe_df.iloc[j, 4], bedpe_df.iloc[j, 5]
        dist1 = max(next_start1 - stop1, start1 - next_stop1)
        dist2 = max(next_start2 - stop2, start2 - next_stop2)
        return (chrom1 == next_chrom1 and chrom2 == next_chrom2 and dist1 <= merge_win and dist2 <= merge_win)

    # The "chain-breaking" in cluster_loci might still leave some redundancy.
    # In particular, we might leave some almost touching clusters that were
    # separated only because of chain-breaking. Do a second round of clustering
    # where you go through consecutive pairs of cluster and merge them if they're merge-able.
    new_cluster_pairs = {}
    for (cluster1, cluster2) in sorted(cluster_pairs.keys()):
        if cluster_pairs[(cluster1, cluster2)] == -1:
            continue
        # Consider all neighboring clusters after this cluster.
        # Notice that the cluster indices are sorted by genomic coordinates.
        neigh_clusters = [(cluster1, cluster2 + 1), (cluster1 + 1, cluster2 - 1),
                          (cluster1 + 1, cluster2), (cluster1 + 1, cluster2 + 1)]
        idx = cluster_pairs[(cluster1, cluster2)]
        # Best cluster among neighboring clusters
        max_cluster = ((cluster1, cluster2), idx)
        for next_cluster1, next_cluster2 in neigh_clusters:
            if not (next_cluster1, next_cluster2) in cluster_pairs:
                continue
            if cluster_pairs[(next_cluster1, next_cluster2)] == -1:
                continue
            next_idx = cluster_pairs[(next_cluster1, next_cluster2)]
            if clusters_close(idx, next_idx):
                cluster_pairs[(next_cluster1, next_cluster2)] = -1
                if bedpe_df.iloc[idx]['qual'] < bedpe_df.iloc[next_idx]['qual']:
                    max_cluster = ((next_cluster1, next_cluster2), next_idx)
        new_cluster_pairs[max_cluster[0]] = max_cluster[1]

    cluster_pairs = new_cluster_pairs

    # Now compute the number of mate breakpoints for each cluster
    num_mates = {}
    for (cluster1, cluster2) in cluster_pairs.keys():
        if not cluster1 in num_mates:
            num_mates[cluster1] = 0
        if not cluster2 in num_mates:
            num_mates[cluster2] = 0
        num_mates[cluster1] += 1
        if cluster2 != cluster1:
            num_mates[cluster2] += 1

    sel_loc = []
    new_info_strs = []
    for (cluster1, cluster2) in sorted(cluster_pairs.keys()):
        sv_loc = cluster_pairs[(cluster1, cluster2)]
        if num_mates[cluster1] > max_nmates and num_mates[cluster2] > max_nmates:
            continue
        sel_loc.append(sv_loc)
        new_info_strs.append(tk_sv_io.update_info(bedpe_df.iloc[sv_loc]['info'], ['NMATES1', 'NMATES2', 'RESOLUTION'],
                                         [num_mates[cluster1], num_mates[cluster2], cluster_dist_ratio[(cluster1, cluster2)]]))
    if len(sel_loc) > 0:
        bedpe_df = bedpe_df.iloc[sel_loc]
        bedpe_df['info'] = new_info_strs
    else:
        bedpe_df = pd.DataFrame(columns = bedpe_df.columns)
    if not out_bedpe is None:
        tk_sv_io.write_sv_df_to_bedpe(bedpe_df, out_bedpe)

    return bedpe_df


def merge_multiple_breaks(in_bedpes, out_bedpe, merge_win = 10000, max_range = np.inf):
    assert(len(in_bedpes) > 0)
    in_bedpe_df = None
    for bi, bedpe in enumerate(in_bedpes):
        bedpe_df = tk_sv_io.read_sv_bedpe_to_df(bedpe)
        assert(bedpe_df.shape[1] > 11)
        bedpe_df = bedpe_df.iloc[:, 0:12]
        # Make sure that all names from all files are unique
        bedpe_df['name'] = [str(n) + '_' + str(bi) for n in bedpe_df['name']]
        in_bedpe_df = pd.concat([in_bedpe_df, bedpe_df], ignore_index = True)

    return merge_breaks(in_bedpe_df, out_bedpe, merge_win = merge_win, max_range = max_range)


def compare_multiple_breaks(in_bedpes, sample_names, out_bedpe, merge_win = 0, max_range = np.inf):
    """Compares multiple BEDPE files.
    Args:
    - in_bedpes: A list of BEDPE files to compare.
    - sample_names: A list of the same size with unique names for the input samples.
    - out_bedpe: Where union BEDPE will be written.

    Return value:
    A DataFrame with the union of calls and information about which
    calls are present in which input files. This DataFrame will have one entry per call in the
    union and will include (among other things) columns <sample>_qual, <sample>_filtered,
    <sample>_correct, and <sample>_dist for each of the input BEDPEs.
    """

    assert(len(sample_names) == len(in_bedpes))

    # Merge all the input files. This will get rid of redundant entries.
    # The quality in the output will be the maximum quality across all files.
    merged_df = merge_multiple_breaks(in_bedpes, out_bedpe, merge_win = merge_win, max_range = max_range)
    num_merged = len(merged_df)

    # Map the name of each entry in the union to its index in the DataFrame.
    name_to_ind = {}
    for i, n in enumerate(merged_df['name']):
        name_to_ind[n] = i

    new_filters = [set([]) for i in range(num_merged)]
    new_matches = [set([]) for i in range(num_merged)]

    # For each of the input BEDPEs find which of the entries in the union it
    # overlaps. This is somewhat duplicated work, but it's simpler this way.
    for sample, bedpe in zip(sample_names, in_bedpes):
        in_df = tk_sv_io.read_sv_bedpe_to_df(bedpe)
        name_to_ind2 = {}
        for i, n in enumerate(in_df['name']):
            name_to_ind2[n] = i

        matched_qual = np.zeros((num_merged, ), dtype = np.int)
        is_correct = np.zeros((num_merged, ), dtype = np.bool)
        is_filtered = np.zeros((num_merged, ), dtype = np.bool)
        tmp_dist = np.zeros((num_merged, ), dtype = np.int)
        matched_names = ['' for i in range(num_merged)]

        # merged_to_this will be a dictionary from a name in the union to a set
        # of names in the input bedpe
        merged_to_this, _, _ = compare_breaks(merged_df, bedpe, max_dist = merge_win)
        for name1, name2_set in merged_to_this.iteritems():
            ind1 = name_to_ind[name1]
            matched_names[ind1] = ';'.join([str(s) for s in name2_set])
            for name2 in name2_set:
                ind2 = name_to_ind2[name2]
                matched_qual[ind1] = max(matched_qual[ind1], in_df.iloc[ind2]['qual'])
                match = tk_sv_io.extract_sv_info(in_df.iloc[ind2]['info'], ['MATCHES'])[0]
                is_match_correct = (match != '.' and match != '' and not match is None)
                if is_match_correct:
                    new_matches[ind1].add(match)
                    # Never set back to False if it was set to true.
                    is_correct[ind1] = True
                is_filtered[ind1] = in_df.iloc[ind2]['filters'] != '.'
                if in_df.iloc[ind2]['filters'] != '.':
                    new_filters[ind1] = new_filters[ind1].union(set(in_df.iloc[ind2]['filters'].split(';')))
                if in_df.iloc[ind2]['chrom1'] != in_df.iloc[ind2]['chrom2']:
                    tmp_dist[ind1] = -1
                else:
                    tmp_dist[ind1] = in_df.iloc[ind2]['start2'] - in_df.iloc[ind2]['stop1']

        merged_df[str(sample) + '_matches'] = matched_names
        merged_df[str(sample) + '_qual'] = matched_qual
        merged_df[str(sample) + '_correct'] = is_correct
        merged_df[str(sample) + '_filtered'] = is_filtered
        merged_df[str(sample) + '_dist'] = tmp_dist

    info_strs = ['.' for i in range(num_merged)]
    filter_strs = ['.' for i in range(num_merged)]
    for i in range(num_merged):
        match_str = ','.join(new_matches[i]) if len(new_matches[i]) > 0 else '.'
        info_strs[i] = tk_sv_io.update_info('.', ['MATCHES'], [match_str])
        filter_strs[i] = ';'.join(new_filters[i]) if len(new_filters[i]) > 0 else '.'

    merged_df['qual'] = np.array(np.max(merged_df[[str(s) + '_qual' for s in sample_names]], axis = 1), dtype = np.int)
    merged_df['filters'] = filter_strs
    merged_df['info'] = info_strs
    merged_df.sort(['qual', 'chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2'],
        ascending = [0, 1, 1, 1, 1, 1, 1], inplace = True)

    return merged_df


def get_dataframe_loc(df, loc):
    if len(loc) > 0:
        return df.iloc[loc]
    else:
        return pd.DataFrame(columns = df.columns)


def get_insert_size_info(insert_size_json, prc=0.9):
    with open(insert_size_json, 'r') as f:
        insert_size_dict = json.load(f)

    if not '60' in insert_size_dict or len(insert_size_dict['60']) == 0:
        return None, None

    insert_size_dict = insert_size_dict['60']
    insert_sizes = np.array([int(k) for k in insert_size_dict.keys()])
    insert_counts = np.array([insert_size_dict[str(k)] for k in insert_sizes])
    sort_idx = np.argsort(insert_sizes)
    insert_sizes = insert_sizes[sort_idx]
    insert_cum_counts = np.cumsum(insert_counts[sort_idx])

    # cum_arr[i] will be the empirical probability of observing
    # insert sizes >= i.
    # Using a fixed size array can save a ton of time downstream
    # since getting the empirical distribution is an operation
    # done very frequently in the short read sv-caller.
    cum_arr = np.zeros((insert_sizes[-1], ))
    for i in range(len(cum_arr)):
        if i <= insert_sizes[0]:
            bidx = 0
        else:
            bidx = bisect_left(insert_sizes, i)
        cum_arr[i] = np.log(1.0 - insert_cum_counts[bidx] / float(insert_cum_counts[-1]))

    assert len(cum_arr) > 0

    def emp_logsf(x, cum):
        x = np.array(x, dtype=np.int)
        res = np.ones((len(x),)) * MIN_LOG_PROB
        res[x < len(cum) - 1] = cum[x[x < len(cum) - 1]]
        return res

    pidx = np.where(insert_cum_counts > prc * insert_cum_counts[-1])[0][0]
    max_insert = insert_sizes[pidx]
    return max_insert, lambda x: emp_logsf(x, cum_arr)


def get_phase_set(frag_phasing, chrom, start, stop):
    """Get the closest phase set to a specified region.
    All overlapping phase sets are assumed to be at distance 0 from the region.
    frag_phasing is a handle for a tabix file.

    Why use the tabix frag_phasing file to get the overlapping
    phase sets instead of the phase_block h5?
    If the point (chrom, start, stop) overlaps a phase set,
    then this should not make a difference (we get the overlapping
    phase set anyway).
    BUT consider the case where X = (chrom, start, stop) is between phase
    sets. There might exist fragments that span X and a nearby phase set.
    Fetching on frag_phasing will get all the fragments that overlap X.
    Then we can use this information to assign a phase set to X, even
    if X was outside a phase block.
    """
    ps_dist = np.inf
    final_ps = None
    if not chrom in frag_phasing.contigs:
        return None

    for frag_line in frag_phasing.fetch(str(chrom), start, stop):
        frag = frag_line.strip().split('\t')
        ps, ps_start, ps_stop = frag[3], int(frag[4]), int(frag[5])
        if not ps is None:
            if start > ps_stop:
                new_ps_dist = start - ps_stop
            elif ps_start > stop:
                new_ps_dist = ps_start - stop
            else:
                new_ps_dist = 0
            if new_ps_dist < ps_dist:
                final_ps = ps
                ps_dist = new_ps_dist
            # If you find an overlapping phase set, you're done.
            if new_ps_dist == 0:
                break
    return final_ps


def get_barcode_phase_probs(frag_phasing, chrom, start, stop, bcs, in_ps=None):
    """Returns the phasing probabilities at a given locus for a set of barcodes.
    bcs is a set or dict.
    Return value:
    dict from barcodes to tuples (phase_set, prob_hap0, prob_hap1).
    """

    if not chrom in frag_phasing.contigs:
        return {}

    out_bcs = {}
    # Find a fragment that overlaps the position and has the same barcode.
    for frag_line in frag_phasing.fetch(str(chrom), start, stop):
        frag = frag_line.strip().split('\t')
        ps, ps_start, ps_stop = frag[3], int(frag[4]), int(frag[5])
        bc, p_hap1, p_hap2 = frag[6], float(frag[7]), float(frag[8])
        p_hap1 = max(1 - MAX_PHASE_PROB, min(MAX_PHASE_PROB, p_hap1))
        p_hap2 = max(1 - MAX_PHASE_PROB, min(MAX_PHASE_PROB, p_hap2))
        if bc in bcs and (in_ps is None or in_ps == ps):
            new_dist = max(0, max(ps_start - stop, start - ps_stop))
            if not bc in out_bcs or out_bcs[bc][1] > new_dist:
                out_bcs[bc] = (ps, new_dist, (p_hap1, p_hap2))
    return out_bcs
