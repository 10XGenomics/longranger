#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import longranger.test as lr_test
import numpy as np
import scipy.sparse as sp
from longranger.sv.overlap_detector import OverlapDetector, filter_mat

class TestOverlapDetector(lr_test.UnitTestBase):
    def setUp(self):
        pass

    def test_filter_mat(self):
        mat = [[0, 1, 2, 3], [10, 1, 1, 1], [0, 0, 1, 5]]
        bc_mat = sp.lil_matrix(np.array(mat))

        out_mat, sel_wins, read_counts = filter_mat(bc_mat)
        correct_mat = [[0, 1, 1, 1], [1, 1, 1, 1], [0, 0, 1, 1]]
        assert(np.all(out_mat.todense() == correct_mat))
        assert(list(sel_wins) == list(np.arange(4)))
        assert(list(read_counts) == [10, 2, 4, 9])

        out_mat, sel_wins, read_counts = filter_mat(bc_mat, max_bcs = 2)
        correct_mat = [[0, 1], [1, 1], [0, 0]]
        assert(np.all(out_mat.todense() == correct_mat))
        assert(list(sel_wins) == [0, 1])
        assert(list(read_counts) == [10, 2])

        out_mat, sel_wins, read_counts = filter_mat(bc_mat, min_reads = 2)
        correct_mat = [[0, 0, 1, 1], [1, 0, 0, 0], [0, 0, 0, 1]]
        assert(np.all(out_mat.todense() == correct_mat))
        assert(list(sel_wins) == [0, 1, 2, 3])
        assert(list(read_counts) == [10, 2, 4, 9])

        out_mat, sel_wins, read_counts = filter_mat(bc_mat, min_reads = 0, max_bcs = 0)
        correct_mat = [[], [], []]
        assert(np.all(out_mat.todense() == correct_mat))
        assert(list(sel_wins) == [])
        assert(list(read_counts) == [])


    def test_empty_bc_mat(self):
        """Ensure reasonable behavior when one of the matrices is empty,
        or becomes empty after filtering."""
        mat1 = sp.lil_matrix(np.array([[0, 1, 2, 3], [10, 1, 1, 1], [0, 0, 1, 5]]))
        mat2 = sp.lil_matrix((3, 4))
        loci1 = (('chr1', np.arange(100, 500, 100), np.arange(200, 600, 100)))
        loci2 = (('chr2', np.arange(100, 500, 100), np.arange(200, 600, 100)))
        overlap_detector = OverlapDetector(mat1, mat2, loci1, loci2)
        (out_loci1, out_loci2, p_vals, overlaps, nbcs1, nbcs2) = overlap_detector.get_overlaps(0, 0, verbose = False)
        assert(len(p_vals) == 0)
        assert(len(overlaps) == 0)
        assert(len(nbcs1) == 0)
        assert(len(nbcs2) == 0)
        assert(len(out_loci1[1]) == 0)
        assert(len(out_loci1[2]) == 0)
        assert(len(out_loci2[1]) == 0)
        assert(len(out_loci2[2]) == 0)

        overlap_detector = OverlapDetector(mat1, mat2, loci1, loci2, min_reads = 0, max_bcs = 0)
        (out_loci1, out_loci2, p_vals, overlaps, nbcs1, nbcs2) = overlap_detector.get_overlaps(0, 0, verbose = False)
        assert(len(p_vals) == 0)
        assert(len(overlaps) == 0)
        assert(len(nbcs1) == 0)
        assert(len(nbcs2) == 0)
        assert(len(out_loci1[1]) == 0)
        assert(len(out_loci1[2]) == 0)
        assert(len(out_loci2[1]) == 0)
        assert(len(out_loci2[2]) == 0)


    def test_zero_dist(self):
        """Ensure that all pairs are returned when dist is set to 0."""
        mat1 = sp.lil_matrix(np.ones((3, 4)))
        loci1 = (('chr1', np.arange(100, 500, 100), np.arange(200, 600, 100)))
        overlap_detector = OverlapDetector(mat1, mat1, loci1, loci1)
        (out_loci1, out_loci2, _, _, _, _) = overlap_detector.get_overlaps(0, 0, min_dist = 0, max_dist = 0, verbose = False)
        assert(list(out_loci1[1]) == [100, 200, 300])
        assert(list(out_loci1[2]) == [200, 300, 400])
        assert(list(out_loci2[1]) == [200, 300, 400])
        assert(list(out_loci2[2]) == [300, 400, 500])

        (out_loci1, out_loci2, _, _, _, _) = overlap_detector.get_overlaps(0, 0, min_dist = 0, max_dist = 1000, verbose = False)
        assert(list(out_loci1[1]) == [100, 100, 100, 200, 200, 300])
        assert(list(out_loci1[2]) == [200, 200, 200, 300, 300, 400])
        assert(list(out_loci2[1]) == [200, 300, 400, 300, 400, 400])
        assert(list(out_loci2[2]) == [300, 400, 500, 400, 500, 500])


    def test_min_overlap(self):
        """Test that increasing min_overlap does not mess up locus indices."""
        mat1 = sp.lil_matrix(np.array([[5, 0, 0, 5], [0, 1, 1, 1], [3, 2, 3, 2]]))
        loci1 = (('chr1', np.arange(100, 500, 100), np.arange(200, 600, 100)))
        overlap_detector = OverlapDetector(mat1, mat1, loci1, loci1)
        res = overlap_detector.get_overlaps(0, 0, min_ov = 2, min_dist = 0,
                                            max_dist = 1000, verbose = False)
        (out_loci1, out_loci2, _, _, _, _) = res
        assert(list(out_loci1[1]) == [100, 200, 200, 300])
        assert(list(out_loci1[2]) == [200, 300, 300, 400])
        assert(list(out_loci2[1]) == [400, 300, 400, 400])
        assert(list(out_loci2[2]) == [500, 400, 500, 500])


    def test_passing_indices_around(self):
        mat1 = np.zeros((100, 10))
        mat1[0:6, 0] = 1
        mat1[0:3, 2] = 1
        mat1[3:6, 4] = 1
        mat1[0:3, 6] = 1
        mat1[3:6, 8] = 1
        # Binomial logp-value of windows 0-2 and 0-6: -4
        # Binomial logp-value of windows 2-6 and 4-8: -inf
        # Binomial logp-value of windows 0-4 and 0-8: -4

        mat1 = sp.lil_matrix(mat1)
        loci1 = (('chr1', np.arange(100, 1100, 100), np.arange(200, 1200, 100)))

        overlap_detector = OverlapDetector(mat1, mat1, loci1, loci1)
        res = overlap_detector.get_overlaps(method = 0, max_logp = 0, min_ov = 1,
                                            min_dist = 0, max_dist = 1000, verbose = False)
        print res
        (out_loci1, out_loci2, p_vals, _, _, _) = res
        assert(list(out_loci1[1]) == [100, 100, 100, 100,  300, 500])
        assert(list(out_loci1[2]) == [200, 200, 200, 200,  400, 600])
        assert(list(out_loci2[1]) == [300, 500, 700, 900,  700, 900])
        assert(list(out_loci2[2]) == [400, 600, 800, 1000, 800, 1000])
        self.assertApproxEqual(p_vals[0], -4, precision = 1e-6)
        self.assertApproxEqual(p_vals[1], -4, precision = 1e-6)
