#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import longranger.test as lr_test
import tenkit.regions as tk_regions
from longranger.sv.constants import MIN_LOG_PROB
from numpy.random import random_integers

import longranger.sv.stats as tk_sv_stats
from longranger.sv.utils import region_cum_coverage_map

class TestSVStats(lr_test.UnitTestBase):
    def setUp(self):
        pass


    def test_get_cov_bases(self):
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 10)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 0, 5), 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 10, 20), 0)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 10, 25), 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 25, 27), 2)
        # The coverage in [32, 34) is int((5/10) * 2)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 32, 34), 1)
        # int((5/10) * 8 + (5/10) * 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 32, 44), 6)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 32, 100), 9)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 0, 0), 0)
        # int((5/10) * 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 45, 50), 2)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 10, 56, 80), 0)

        cov_map = region_cum_coverage_map(regions, 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 5, 0, 5), 5)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 5, 10, 20), 0)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 5, 32, 34), 0)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 5, 32, 45), 10)
        self.assertEqual(tk_sv_stats.get_cov_bases(cov_map['chr1'], 5, 45, 50), 0)


    def test_mean_target_by_offset(self):
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 5)['chr1']
        c, cu, cd, bins = tk_sv_stats.mean_target_by_offset(20, 30, cov_map, 5, np.array([10]))
        self.assertEqual(list(c), [0])
        self.assertEqual(list(cu), [0])
        self.assertEqual(list(cd), [0])
        self.assertEqual(list(bins), [0])
        c, cu, cd, bins = tk_sv_stats.mean_target_by_offset(20, 30, cov_map, 5, np.array([11, 14, 20, 30, 60]))
        self.assertEqual(list(c), [0, 0, 0.25, 0.5, 2/7.0])
        self.assertEqual(list(cu), [0, 0, 0, 0.5, 0.5])
        self.assertEqual(list(cd), [0, 0, 0.5, 0.5, 0.2])
        self.assertEqual(list(bins), [1, 4, 10, 20, 50])

        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 10)['chr1']
        c, cu, cd, bins = tk_sv_stats.mean_target_by_offset(20, 30, cov_map, 10, np.array([10]))
        self.assertEqual(list(c), [0])
        self.assertEqual(list(cu), [0])
        self.assertEqual(list(cd), [0])
        self.assertEqual(list(bins), [0])
        c, cu, cd, bins = tk_sv_stats.mean_target_by_offset(20, 30, cov_map, 10, np.array([11, 15, 20]))
        self.assertEqual(list(c), [0, 0.2, 0.25])
        self.assertEqual(list(cu), [0, 0, 0])
        self.assertEqual(list(cd), [0, 0.4, 0.5])
        self.assertEqual(list(bins), [1, 5, 10])
        c, cu, cd, bins = tk_sv_stats.mean_target_by_offset(40, 55, cov_map, 10, np.array([20, 100]))
        self.assertEqual(list(c), [0.2, 25/(40.0 + 85.0)])
        self.assertEqual(list(cu), [0.4, 25/40.0])
        self.assertEqual(list(cd), [0, 0])
        self.assertEqual(list(bins), [5, 85])


    def test_estimate_extent(self):
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 5)
        self.assertEqual(tk_sv_stats.estimate_extent('chr1', 0, 10, cov_map, 5, [], 0.0001, 0.5), (0, 0, 0, 0))
        self.assertEqual(tk_sv_stats.estimate_extent('chr1', 0, 10, cov_map, 5, [10], 0.001, 0.1), (0, 10, 0, -0.001))
        res = tk_sv_stats.estimate_extent('chr1', 20, 30, cov_map, 5, [10], 0.001, 0.1)
        np.testing.assert_array_almost_equal(res, [10, 10, -0.001, np.log(np.exp(-0.001 * 5) * np.exp(-0.0001 * 5))])


    def merge_freq_bins(self):
        sizes = [0, 1, 2, 3, 4]
        counts = [1, 1, 1, 1, 1]
        new_sizes, new_counts = tk_sv_stats.merge_freq_bins(sizes, counts, 1)
        self.assertEqual(list(new_sizes), [4])
        self.assertEqual(list(new_counts), [5])
        new_sizes, new_counts = tk_sv_stats.merge_freq_bins(sizes, counts, 0.2)
        self.assertEqual(list(new_sizes), [1, 3, 4])
        self.assertEqual(list(new_counts), [2, 2, 1])
        new_sizes, new_counts = tk_sv_stats.merge_freq_bins(sizes, counts, 0)
        self.assertEqual(list(new_sizes), list(sizes))
        self.assertEqual(list(new_counts), list(counts))

        sizes = [0, 1, 2, 3, 4]
        counts = [1, 0, 1, 1, 1]
        new_sizes, new_counts = tk_sv_stats.merge_freq_bins(sizes, counts, 0.3)
        self.assertEqual(list(new_sizes), [3, 4])
        self.assertEqual(list(new_counts), [3, 1])

        sizes = np.arange(100)
        counts = random_integers(0, 100, 100)
        new_sizes, new_counts = tk_sv_stats.merge_freq_bins(sizes, counts, 0.1)
        self.assertEqual(np.sum(new_counts), np.sum(counts))


    def log_prob_frag(self):
        nr = 5
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 10)['chr1']
        alpha = 0.01
        prob_off_target = 1.0
        frag_sizes = np.array([10, 20, 30])
        frag_counts = np.array([30, 20, 10])
        mean_target_ratio = np.array([1.0, 1.0])

        # Test that the target and wgs versions produce the same results when
        # the correction factor for off-target regions is 1 (so the amp rate
        # is the same as on-target).
        p1, pr1 = tk_sv_stats.log_prob_frag_target(nr, cov_map, 10, 0, 10, alpha,
            prob_off_target, mean_target_ratio, frag_sizes, frag_counts)
        p2, pr2 = tk_sv_stats.log_prob_frag(nr, 0, 10, alpha, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        np.testing.assert_almost_equal(pr1, pr2, decimal = 8)

        # Test that the probability of generating reads off-target when the
        # correction factor is 0 is MIN_LOG_PROB.
        frag_sizes = np.array([10, 20])
        frag_counts = np.array([10, 20])
        p1, pr1 = tk_sv_stats.log_prob_frag_target(nr, cov_map, 10, 10, 20, alpha, 0,
            mean_target_ratio, frag_sizes, frag_counts)
        self.assertEqual(p1, MIN_LOG_PROB)
        self.assertEqual(pr1, MIN_LOG_PROB)

        # Test that target probs are smaller when corr_factor < 1 and
        # larger when corr_factor > 1.
        p1, pr1 = tk_sv_stats.log_prob_frag_target(nr, cov_map, 10, 10, 20, alpha,
            0.5, 0.5 * mean_target_ratio, frag_sizes, frag_counts)
        p2, pr2 = tk_sv_stats.log_prob_frag(nr, 10, 20, alpha, frag_sizes, frag_counts)
        assert(p1 < p2)
        assert(pr1 < pr2)

        p1, pr1 = tk_sv_stats.log_prob_frag_target(nr, cov_map, 10, 10, 20, alpha,
                2, 0.5 * mean_target_ratio, frag_sizes, frag_counts)
        p2, pr2 = tk_sv_stats.log_prob_frag(nr, 10, 20, alpha, frag_sizes, frag_counts)
        assert(p1 > p2)
        assert(pr1 > pr2)


    def read_generation_pmf(self):
        n = 10
        p1 = tk_sv_stats.read_generation_pmf(n, 0.2, 0.1, 10, 5)
        p2 = tk_sv_stats.read_generation_pmf(n, 0.1, 0.2, 5, 10)
        self.assertEqual(p1, p2)

        # Probability of generating reads from the low is 0
        # so it doesn't matter how long the low regions are.
        n = 10
        p1 = tk_sv_stats.read_generation_pmf(n, 0.2, 0, 10, 5)
        p2 = tk_sv_stats.read_generation_pmf(n, 0.2, 0, 10, 100)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        n = 0
        p1 = tk_sv_stats.read_generation_pmf(n, 0.2, 0, 10, 5)
        p2 = tk_sv_stats.read_generation_pmf(n, 0.2, 0, 10, 100)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)

        p1 = tk_sv_stats.read_generation_pmf(0, 0.2, 0.1, 10, 5)
        np.testing.assert_almost_equal(p1, np.log(np.exp(-2) * np.exp(-0.5)))


    def off_target_amp_corr_factor(self):
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        self.assertEqual(tk_sv_stats.off_target_amp_corr_factor(regions, 0.5, 60), 1.0)
        np.testing.assert_almost_equal(tk_sv_stats.off_target_amp_corr_factor(regions, 0.1, 3e6), 0, decimal = 5)


    def test_merge_lr(self):
        nr = 5
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 10)
        alpha = 0.01
        frag_sizes = np.array([10, 20, 30])
        frag_counts = np.array([30, 20, 10])

        # Test that the target and wgs versions produce the same results when
        # the correction factor for off-target regions is 1 (so the amp rate
        # is the same as on-target).
        p1, pr1, _, _, po1 = tk_sv_stats.merge_lr_target(nr, nr, 'chr1', 0, 2, 'chr1', 8, 10,
            alpha, cov_map, 10, 1.0, frag_sizes, frag_counts)
        p2, pr2, _, _, po2 = tk_sv_stats.merge_lr_wgs(nr, nr, 'chr1', 0, 2, 'chr1', 8, 10,
            alpha, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        np.testing.assert_almost_equal(pr1, pr2, decimal = 8)
        np.testing.assert_almost_equal(po1, po2, decimal = 8)

        p1, pr1, _, _, po1 = tk_sv_stats.merge_lr_target(nr, nr, 'chr1', 0, 2, 'chr1', 45, 48,
            alpha, cov_map, 10, 1.0, frag_sizes, frag_counts)
        p2, pr2, _, _, po2 = tk_sv_stats.merge_lr_wgs(nr, nr, 'chr1', 0, 2, 'chr1', 45, 48,
            alpha, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        np.testing.assert_almost_equal(pr1, pr2, decimal = 8)
        np.testing.assert_almost_equal(po1, po2, decimal = 8)
        self.assertEqual(p1, MIN_LOG_PROB)

        # Longer fragments: p_one_frag is larger
        frag_sizes = np.array([10, 20, 30])
        frag_counts = np.array([10, 50, 60])
        p1, pr1, _, _, po1 = tk_sv_stats.merge_lr_target(5, 5, 'chr1', 0, 10, 'chr1', 15, 20,
            2, cov_map, 10, 1.0, frag_sizes, frag_counts)
        p2, pr2, _, _, po2 = tk_sv_stats.merge_lr_wgs(5, 5, 'chr1', 0, 10, 'chr1', 15, 20,
            2, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        np.testing.assert_almost_equal(pr1, pr2, decimal = 8)
        np.testing.assert_almost_equal(po1, po2, decimal = 8)
        assert(p1 > pr1)

        # Short fragments: p_one_frag is smaller
        frag_sizes = np.array([10, 20, 30])
        frag_counts = np.array([10, 5, 1])
        p1, pr1, _, _, po1 = tk_sv_stats.merge_lr_target(5, 5, 'chr1', 0, 10, 'chr1', 15, 20,
            2, cov_map, 10, 1.0, frag_sizes, frag_counts)
        p2, pr2, _, _, po2 = tk_sv_stats.merge_lr_wgs(5, 5, 'chr1', 0, 10, 'chr1', 15, 20,
            2, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(p1, p2, decimal = 8)
        np.testing.assert_almost_equal(pr1, pr2, decimal = 8)
        np.testing.assert_almost_equal(po1, po2, decimal = 8)
        assert(p1 < pr1)


    def test_lr(self):
        nr = 5
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10000], [50000, 60000], [75000, 85000]])
        cov_map = region_cum_coverage_map(regions, 1000)

        # Test that the code is force-linking fragments that are too close to each other
        # if they are on the same chromosome.
        alpha = 0.01
        nr = 4
        frag_sizes = np.array([20001, 40002, 80002, 100000])
        frag_counts = np.array([300, 100, 100, 100])
        # fragments should be merged
        lr1, _, _, s1, s2, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr1', 22000, 32000,
                                                  alpha, cov_map, 1.0, frag_sizes, frag_counts,
                                                  min_dist = 10000)
        self.assertEqual(lr1, 0.0)
        self.assertEqual(s1, None)
        self.assertEqual(s2, None)
        # fragments should not be merged because they are on different chromosomes
        lr1, _, _, s1, s2, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr2', 22000, 32000,
            alpha, cov_map, 1.0, frag_sizes, frag_counts)
        assert(lr1 > 0.0)
        assert(s1 != None)
        lr2, _, _, s1, s2, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr2', 0, 20000, 'chr1', 22000, 32000,
            alpha, cov_map, 1.0, frag_sizes, frag_counts)
        self.assertEqual(lr1, lr2)

        # The only possible explanation is that these are two independent molecules in
        # the same GEM (too long to be one molecule with SV or without).
        alpha = 0.001
        frag_sizes = np.array([20001])
        frag_counts = np.array([300])
        lr1, _, _, _, _, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, cov_map, 1.0, frag_sizes, frag_counts)
        lr2, _, _, _, _, _ = tk_sv_stats.lr_wgs(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(lr1, lr2, decimal = 8)
        np.testing.assert_almost_equal(lr1, 0.0, decimal = 8)
        # ... and irrespective of prior.
        lr1, _, _, _, _, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, cov_map, 1.0, frag_sizes, frag_counts, p_prior = np.log(0.01))
        lr2, _, _, _, _, _ = tk_sv_stats.lr_wgs(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, frag_sizes, frag_counts, p_prior = np.log(0.01))
        np.testing.assert_almost_equal(lr1, lr2, decimal = 8)
        np.testing.assert_almost_equal(lr1, 0.0, decimal = 8)

        # Given high amp-rate, it is unlikely that this was a single fragment.
        alpha = 0.01
        nr = 200
        frag_sizes = np.array([20001, 40002, 60000])
        frag_counts = np.array([300, 100, 100])
        lr1, _, _, _, _, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, cov_map, 1.0, frag_sizes, frag_counts)
        lr2, _, _, _, _, _ = tk_sv_stats.lr_wgs(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, frag_sizes, frag_counts)
        np.testing.assert_almost_equal(lr1, lr2, decimal = 8)
        np.testing.assert_almost_equal(lr1, 0.0, decimal = 8)

         # Absence of reads should be more likely in the targeted case.
        alpha = 0.01
        nr = 4
        frag_sizes = np.array([20001, 80002, 100000])
        frag_counts = np.array([300, 100, 100])
        lr1, _, _, _, _, _, _, _, _ = tk_sv_stats.lr_target(nr, nr, 'chr1', 0, 20000, 'chr1', 72000, 82000,
            alpha, cov_map, 0.01, frag_sizes, frag_counts)
        lr2, _, _, _, _, _ = tk_sv_stats.lr_wgs(nr, nr, 'chr1', 0, 20000, 'chr1', 52000, 72000,
            alpha, frag_sizes, frag_counts)
        assert(lr1 > lr2)
