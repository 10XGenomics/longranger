#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import tenkit.bam as tk_bam
import longranger.test as lr_test
from longranger.sv.breakpoint_analyzer import *

bam_filename = lr_test.in_path('test_count_reads.bam')

class TestBreakpointAnalyzer(lr_test.UnitTestBase):
    def setUp(self):
        bc_map = {}
        bc_map['4-GCAGTTAGAGAAAT'] = 0
        bc_map['1-GCTCCTGTATGGCG'] = 1
        bc_map['1-GATGAAGTACTGAA'] = 2
        bc_map['8-ACTTTCGTTAATCT'] = 3
        bc_map['8-GGGTAGTCAGTAAG'] = 4
        bc_map['2-TCCCGTTCCTGGAT'] = 5
        bc_map['6-CGTCAATCTTGGCA'] = 6
        bc_map['3-TGTCGAGTCCGCTG'] = 7
        bc_map['6-GCGAAGTCCCTAAG'] = 8
        bc_map['3-TCAGTGGTCCAATC'] = 9
        bc_map['2-CAGAAAGTCTTGCA'] = 10
        bc_map['6-ATGCGTAGTTTCTA'] = 11
        bc_map['6-ACTCAGCAGACATA'] = 12
        bc_map['1-GGGACATCTCCACC'] = 13
        bc_map['2-TCCTTATCCTGGAT'] = 14
        bc_map['3-GTCGTAAGTGACAT'] = 15
        bc_map['2-TCCCGTTCCTGGAT'] = 16
        bc_map['1-CATTCTCATCGTCA'] = 17
        bc_map['8-GTTCTTTCTTCGAG'] = 18
        bc_map['5-CGTCAAGTTAGACA'] = 19

        bc_freq = np.ones((len(bc_map), )) * 0.01
        read_freq = np.ones((len(bc_map), )) * 0.01
        self.targets = lr_test.in_path('test_breakpoint_analyzer_targets.bed')
        self.target_analyzer = BreakpointAnalyzer(bam_filename, bc_freq, bc_map,
            read_freq = read_freq, regions_file = self.targets, extend = 0)
        self.target_analyzer_100 = BreakpointAnalyzer(bam_filename, bc_freq, bc_map,
            read_freq = read_freq, regions_file = self.targets, extend = 100)
        self.analyzer = BreakpointAnalyzer(bam_filename, bc_freq, bc_map,
            read_freq = read_freq, regions_file = None, extend = 0)
        self.analyzer_100 = BreakpointAnalyzer(bam_filename, bc_freq, bc_map,
            read_freq = read_freq, regions_file = None, extend = 100)


    def test_get_closest_regions(self):
        regions = self.target_analyzer.get_closest_regions('chr19', 287513, 287516, 287300, 287528)
        self.assertEqual(regions, [(287513, 287516), (287518, 287520),
            (287503, 287507), (287523, 287526), (287400, 287407)])

        regions = self.target_analyzer.get_closest_regions('chr19', 100, 200, 400, 500)
        self.assertEqual(regions, [(0, 0)])

        regions = self.target_analyzer_100.get_closest_regions('chr19', 287513, 287516, 287300, 287528)
        self.assertEqual(regions, [(287300, 287708)])


    def test_get_bcs_around_break(self):
        # This is a little tricky because sams/bams are 1-based but the coordinates pysam uses
        # are 0-based. So the coordinates below would correspond to
        # samtools view bam_filename chr19:(bstart+1):(bstop+1).
        bstart = 287605
        bstop = 287608
        bam = tk_bam.create_bam_infile(bam_filename)
        bcs, start, stop = self.analyzer.get_bcs_around_break(bam, 'chr19', bstart, bstop, 0, 0,
            min_reads = 100, min_mapq = 60)
        self.assertEqual(start, bstart)
        self.assertEqual(stop, bstop)
        self.assertEqual(len(bcs), 3)

        bcs, start, stop = self.analyzer.get_bcs_around_break(bam, 'chr19', bstart, bstop, 10, 0,
            min_reads = 100, min_mapq = 60)
        self.assertEqual(start, bstart - 10)
        self.assertEqual(stop, bstop)
        self.assertEqual(len(bcs), 5)

        bcs, start, stop = self.target_analyzer.get_bcs_around_break(bam, 'chr19', bstart, bstop, 0, 0,
            min_reads = 100, min_mapq = 60)
        self.assertEqual(start, 287606)
        self.assertEqual(stop, 287608)
        self.assertEqual(len(bcs), 1)

        bcs, start, stop = self.target_analyzer.get_bcs_around_break(bam, 'chr19', bstart, bstop, 10, 10,
            min_reads = 100, min_mapq = 60)
        self.assertEqual(start, 287606)
        self.assertEqual(stop, 287608)
        self.assertEqual(len(bcs), 1)

        bcs, start, stop = self.target_analyzer.get_bcs_around_break(bam, 'chr19', 287513, 287516, 100, 100,
            min_reads = 100, min_mapq = 60)
        # Overlapping regions are: [(287513, 287516), (287518, 287520),
        #    (287503, 287507), (287523, 287526), (287606, 287608)]
        self.assertEqual(start, 287503)
        self.assertEqual(stop, 287608)
        self.assertEqual(len(bcs), 3)

        bcs, start, stop = self.target_analyzer.get_bcs_around_break(bam, 'chr19', 287513, 287516, 100, 100,
            min_reads = 0, min_mapq = 60)
        self.assertEqual(start, 287513)
        self.assertEqual(stop, 287516)
        self.assertEqual(len(bcs), 1)
        bam.close()


    def test_get_stats_at_breaks(self):
        breakpoints = [('chr19', 287605, 287608, 'chr19', 287605, 287608),
                       ('chr19', 287605, 287608, 'chr19', 287512, 287516),]
        stat_df = self.analyzer.get_stats_at_breaks(breakpoints, 0, 0)
        self.assertEqual(stat_df.iloc[0]['bcOv'], 3)
        self.assertEqual(stat_df.iloc[1]['bcOv'], 1)
