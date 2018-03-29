#!/usr/bin/env python
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

import sys
import numpy as np
import tenkit.pandas as pd
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.regions as tk_regions
from longranger.sv.utils import *

import martian


def merge_bc_counts(tot_bc_counts, bc_counts):
    for b, read_list in bc_counts.iteritems():
        if not b in tot_bc_counts:
            tot_bc_counts[b] = set([])
        tot_bc_counts[b] = tot_bc_counts[b].union(read_list)
    return tot_bc_counts


class BreakpointAnalyzer():
    """Class for finding barcode and read overlaps around given SVs"""

    def __init__(self, in_bam_filename, bc_freq, bc_map, read_freq = None,
                 regions_file = None, extend = 0):
        """
        in_bam_filename: BAM file
        bc_freq: numpy array with barcode frequencies
        bc_map: dict from barcodes to indices (in bc_freq and read_freq)
        read_freq: numpy array with read frequencies
        regions_file: BED file with target regions (or None). If provided, only
           reads from these regions will be counted (after extending these regions
            by "extend")
        extend: How much to extend target regions. Should be a non-negative integer.
        """
        self.bam = in_bam_filename
        self.bc_freq = bc_freq
        self.bc_map = bc_map
        self.read_freq = read_freq
        if regions_file is None:
            self.regions = None
        else:
            self.regions = {}
            extend = max(extend, 0)
            bed_iterator = tk_io.get_bed_iterator(regions_file)
            for chrom, start, stop in bed_iterator:
                start = max(0, start - extend)
                stop = stop + extend
                if not chrom in self.regions:
                    self.regions[chrom] = []
                self.regions[chrom].append((start, stop))

            for chrom, region_list in self.regions.iteritems():
                self.regions[chrom] = tk_regions.Regions(self.regions[chrom])


    def get_closest_regions(self, chrom, start, stop, ext_start, ext_stop):
        """Get a list of target regions that are within (ext_start, ext_stop),
        sorted by their distance from (start, stop)"""


        scan_regions_tmp = self.regions[chrom].overlapping_regions(ext_start, ext_stop)
        scan_regions = []
        for (scan_start, scan_stop) in scan_regions_tmp:
            # Distance will be 0 if regions are overlapping.
            dist = max(0, max(scan_start - stop, start - scan_stop))
            scan_regions.append((scan_start, scan_stop, dist))
        if len(scan_regions) == 0:
            scan_regions = [(0, 0)]
        else:
            scan_regions = [(r[0], r[1]) for r in sorted(scan_regions, key = lambda x: (x[2], x[0], x[1]))]
        return scan_regions


    def get_bcs_around_break(self, in_bam, chrom, start, stop, win_left, win_right,
        min_reads = 100, min_mapq = 60):

        ext_start = max(0, int(start - win_left))
        ext_stop = int(stop + win_right)

        if self.regions is None:
            scan_regions = [(ext_start, ext_stop)]
        else:
            if not chrom in self.regions:
                scan_regions = [(0, 0)]
            else:
                scan_regions = self.get_closest_regions(chrom, start, stop, ext_start, ext_stop)

        tot_reads = 0
        tot_bc_counts = {}
        out_start = scan_regions[0][0]
        out_stop = scan_regions[0][1]
        for (scan_start, scan_stop) in scan_regions:
            bc_counts = get_bcs_at_region(in_bam, chrom, scan_start, scan_stop, min_mapq = min_mapq, bc_map = self.bc_map)
            tot_bc_counts = merge_bc_counts(tot_bc_counts, bc_counts)
            tot_reads = tot_reads + np.sum([len(b) for b in bc_counts.values()])
            out_start = min(out_start, scan_start)
            out_stop = max(out_stop, scan_stop)
            if tot_reads > min_reads:
                break
        return (tot_bc_counts, out_start, out_stop)


    def get_stats_at_breaks(self, breakpoints, win_left, win_right,
        min_reads = 100, min_mapq = 60, method = BINOM_EMP_BC_COUNT_BC_FREQ,
        outward_only = False):

        columns = ['chrom1', 'start1', 'stop1', 'extStart1', 'extStop1',
            'chrom2', 'start2', 'stop2', 'extStart2', 'extStop2',
            'bcOv', 'nbcs1', 'nbcs2', 'readOv', 'nreads1', 'nreads2', 'binomQual', 'qual',
            'bcs', 'bcFreqs']

        stat_df = pd.DataFrame(columns = columns, index = np.arange(len(breakpoints)))

        in_bam = tk_bam.create_bam_infile(self.bam)

        nbcs = len(self.bc_map)

        for bidx, breakpoint in enumerate(breakpoints):
            chrom1, start1, stop1, chrom2, start2, stop2 = breakpoint[0:6]
            if (chrom1, start1, stop1) > (chrom2, start2, stop2):
                chrom1, start1, stop1, chrom2, start2, stop2 = chrom2, start2, stop2, chrom1, start1, stop1

            bc1_counts, ext_start1, ext_stop1 = self.get_bcs_around_break(in_bam, chrom1, start1, stop1,
                win_left, (0 if outward_only else win_right), min_reads = min_reads, min_mapq = min_mapq)
            bc2_counts, ext_start2, ext_stop2  = self.get_bcs_around_break(in_bam, chrom2, start2, stop2,
                (0 if outward_only else win_left), win_right, min_reads = min_reads, min_mapq = min_mapq)

            bc_ov = set(bc1_counts.keys()).intersection(set(bc2_counts.keys()))

            nbcs1, nbcs2 = len(bc1_counts), len(bc2_counts)
            read_list1 = []
            for b, reads in bc1_counts.iteritems():
                read_list1.extend(reads)
            read_list1 = set(read_list1)
            read_list2 = []
            for b, reads in bc2_counts.iteritems():
                read_list2.extend(reads)
            read_list2 = set(read_list2)
            nreads1 = len(read_list1)
            nreads2 = len(read_list2)
            read_ov = len(read_list1.intersection(read_list2))

            if len(bc_ov) > 0:
                bc_idx = np.array([self.bc_map[b] for b in bc_ov]).flatten()
                win_idx = np.zeros(bc_idx.shape, dtype = np.int)
                if method == BINOM_EMP_BC_COUNT_BC_FREQ or method == 5:
                    pov = log10_emp_pval(win_idx, bc_idx, max(nbcs1, nbcs2), self.bc_freq)
                elif method == BINOM:
                    pov = log10_binom_pval(len(bc_ov), nbcs1, nbcs2, nbcs)
                else:
                    martian.throw('Unsupported method for quality computation.')
                qual = pval_to_qual(pov)
                binom_qual = pval_to_qual(log10_binom_pval(len(bc_ov), nbcs1, nbcs2, nbcs))
            else:
                qual = 0
                binom_qual = 0
            bcs = ','.join(list(bc_ov))
            bc_freqs = ','.join(['{:.2f}'.format(-np.log10(self.bc_freq[self.bc_map[b]])) for b in bc_ov])

            stat_df.loc[bidx] = [chrom1, int(start1), int(stop1), ext_start1, ext_stop1,
                chrom2, int(start2), int(stop2), ext_start2, ext_stop2,
                len(bc_ov), nbcs1, nbcs2, read_ov, nreads1, nreads2, binom_qual, qual, bcs, bc_freqs]

        in_bam.close()

        return stat_df
