#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Test the demultiplexer

import numpy as np
import tenkit.test as tk_test
from .. import *


class TestPrepareSvCallingRanges(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_get_good_ranges(self):
        ranges = get_good_ranges([], 100)
        self.assertEqual(ranges, [])
        
        ranges = get_good_ranges(np.arange(0, 100), 10)
        self.assertEqual(ranges, [(0, 100)])
        
        ranges = get_good_ranges(np.arange(0, 100), 101)
        self.assertEqual(ranges, [])
        
        ranges = get_good_ranges(np.array([0, 1, 4, 5]), 3)
        self.assertEqual(ranges, [(0, 6)])
        
        ranges = get_good_ranges(np.array([0, 1, 4, 5]), 1)
        self.assertEqual(ranges, [(0, 2), (4, 6)])
        
        ranges = get_good_ranges(np.array([1, 2, 3, 8, 9, 10]), 3)
        self.assertEqual(ranges, [(1, 4), (8, 11)])

        ranges = get_good_ranges(np.array([1, 2, 8, 9, 10]), 3)
        self.assertEqual(ranges, [(8, 11)])
        
        ranges = get_good_ranges(np.array([1, 2, 3, 10]), 3)
        self.assertEqual(ranges, [(1, 4)])

    
    def test_get_good_ranges_old(self):
        # No bad regions
        ranges = get_good_ranges_old(np.array([]), 1000)
        self.assertEqual(ranges, [(0, 1000)])
        
        # No bad regions, min region smaller than chr
        ranges = get_good_ranges_old(np.array([]), 100, 1000)
        self.assertEqual(ranges, [])

        # No good regions
        ranges = get_good_ranges_old(np.arange(0, 1000), 1000)
        self.assertEqual(ranges, [])
        
        ranges = get_good_ranges_old(np.arange(0, 1000), 100)
        self.assertEqual(ranges, [])

        # Bad regions extending beyond end of chromosome
        bad_poses = np.array([0, 1, 2, 10, 11, 12, 13, 20, 21, 100, 101, 102])
        ranges = get_good_ranges_old(bad_poses, 100, 1)
        self.assertEqual(ranges, [(3, 10), (14, 20), (22, 100)])
        
        bad_poses = np.array([100, 101, 102])
        ranges = get_good_ranges_old(bad_poses, 100, 1)
        self.assertEqual(ranges, [(0, 100)])
        
        # Bad regions starting at 0
        bad_poses = np.array([0, 1, 2, 10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 1000, 1)
        self.assertEqual(ranges, [(3, 10), (14, 20), (22, 1000)])
        
        bad_poses = np.array([0, 1, 2, 10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 1000, 4)
        self.assertEqual(ranges, [(0, 10), (14, 1000)])
        
        # Testing minimum size of bad region
        bad_poses = np.array([0, 1, 2, 10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 1000, 5)
        self.assertEqual(ranges, [(0, 1000)])
        
        bad_poses = np.array([0, 1, 2, 10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 2, 5)
        self.assertEqual(ranges, [])
        
        bad_poses = np.array([10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 1000, 5)
        self.assertEqual(ranges, [(0, 1000)])
        
        ranges = get_good_ranges_old(bad_poses, 1000, 3)
        self.assertEqual(ranges, [(0, 10), (14, 1000)])
        
        # Other tests
        bad_poses = np.array([10, 11, 12, 13, 20, 21])
        ranges = get_good_ranges_old(bad_poses, 1000, 1)
        self.assertEqual(ranges, [(0, 10), (14, 20), (22, 1000)])

