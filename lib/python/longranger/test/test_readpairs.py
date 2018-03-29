#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import os
import os.path
import numpy as np
import tenkit.bam as tk_bam
import longranger.test as lr_test
from scipy.stats import norm
from longranger.test import TEST_FILE_IN_DIR
import longranger.sv.readpairs as tk_readpairs

def infile(path):
    return os.path.join(TEST_FILE_IN_DIR, 'readpair_support', path)

NO_PAIR_FILE = infile('test_no_pairs.bam')
SPLIT_FILE = infile('test_split.bam')
PAIR_FILE = infile('test_pairs.bam')

class TestReadPair(lr_test.UnitTestBase):
    def setUp(self):
        pass

    def test_get_ranges(self):
        mat = np.zeros((10, 10))
        mat[0:3, 0:3] = 1
        mat[5:10, 0:3] = 1
        mat[5:6, 0:5] = 1
        mat[3:4, 7:9] = 1
        a, b = mat.nonzero()
        self.assertEqual(tk_readpairs.get_ranges(a, b, max_dist = 1), [((0,3), (0,3)), ((3,4), (7,9)), ((5,10), (0,5))])
        self.assertEqual(tk_readpairs.get_ranges(a, b, max_dist = 10), [((0,10), (0,9))])

        mat = np.diag(np.ones((10,)))
        a, b = mat.nonzero()
        self.assertEqual(tk_readpairs.get_ranges(a, b, 1), [((0,10),(0,10))])

        mat[5,5] = 0
        a, b = mat.nonzero()
        self.assertEqual(tk_readpairs.get_ranges(a, b, 1), [((0,5), (0,5)), ((6,10), (6,10))])


    def test_overlaps(self):
        assert(tk_readpairs.overlaps((0, 5), (0, 10), 0))
        assert(tk_readpairs.overlaps((0, 5), (6, 10), 1))
        assert(tk_readpairs.overlaps((5, 10), (0, 6), 0))
        assert(tk_readpairs.overlaps((5, 10), (6, 7), 0))
        assert(tk_readpairs.overlaps((5, 10), (15, 20), 10))
        assert(not tk_readpairs.overlaps((5, 10), (15, 20), 1))


    def test_no_readpairs(self):
        in_bam = tk_bam.create_bam_infile(NO_PAIR_FILE)
        for idx, read in enumerate(in_bam):
            if idx % 2 == 0:
                read1 = read
            else:
                read2 = read
                assert(not tk_readpairs.is_split(read1, read2) and not tk_readpairs.is_pair(read1, read2))
        in_bam.close()


    def test_split(self):
        max_ins = 386
        in_bam = tk_bam.create_bam_infile(SPLIT_FILE)
        for idx, read in enumerate(in_bam):
            if idx % 2 == 0:
                read1 = read
            else:
                read2 = read
                rp = tk_readpairs.ReadPair(read1, read2)
                if idx == 1: # HWI-D00684:18:HBEJCADXX:1:1214:3039:46742
                    self.assertEqual(tk_readpairs.map_qstart(read1), 40)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 0)
                    assert(not tk_readpairs.is_split_consecutive(read1, read2))
                    assert(not rp.is_first_half_first())
                    self.assertEqual(rp.sv_type, tk_readpairs.INV_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (189703752 - tk_readpairs.SPLIT_LEN, 189703752 + tk_readpairs.SPLIT_LEN + 1))
                    self.assertEqual(range2, (189790033 - tk_readpairs.SPLIT_LEN, 189790033 + tk_readpairs.SPLIT_LEN + 1))
                    probs = rp.get_break_lr([range1[0] - 3 * tk_readpairs.SPLIT_LEN, range1[0] + 1, range1[0]],
                        [range2[0] - 3 * tk_readpairs.SPLIT_LEN, range2[0], range2[0]], 0, lambda x: 1.0)
                    np.testing.assert_almost_equal(probs[0], 0.0)
                elif idx == 11: # HISEQ-002:213:HCGHVADXX:2:2213:5344:43278
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 39)
                    assert(tk_readpairs.is_split_consecutive(read1, read2))
                    assert(rp.is_first_half_first())
                    self.assertEqual(rp.sv_type, tk_readpairs.DEL_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (162512094 + 40 - tk_readpairs.SPLIT_LEN,
                                              162512094 + 40 + tk_readpairs.SPLIT_LEN + 1))
                    self.assertEqual(range2, (162626334 - tk_readpairs.SPLIT_LEN,
                                              162626334 + tk_readpairs.SPLIT_LEN + 1))
                elif idx == 13: # HISEQ-002:213:HCGHVADXX:2:2213:5344:43278 (both reads are split)
                    self.assertEqual(tk_readpairs.map_qstart(read1), 48)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 0)
                    assert(tk_readpairs.is_split_consecutive(read1, read2))
                    assert(not rp.is_first_half_first())
                    self.assertEqual(rp.sv_type, tk_readpairs.DEL_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (162512094 + 40 - tk_readpairs.SPLIT_LEN,
                                              162512094 + 40 + tk_readpairs.SPLIT_LEN + 1))
                    self.assertEqual(range2, (162626334 - tk_readpairs.SPLIT_LEN,
                                              162626334 + tk_readpairs.SPLIT_LEN + 1))
                elif idx == 21: # HISEQ-002:203:HBERLADXX:2:1211:3679:57658
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 39)
                    assert(not tk_readpairs.is_split_consecutive(read1, read2))
                    assert(rp.is_first_half_first())
                    self.assertEqual(rp.sv_type, tk_readpairs.INV_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (39226167 - tk_readpairs.SPLIT_LEN,
                                              39226167 + tk_readpairs.SPLIT_LEN + 1))
                    self.assertEqual(range2, (39392916 + 49 - tk_readpairs.SPLIT_LEN,
                                              39392916 + 49 + tk_readpairs.SPLIT_LEN + 1))
                elif idx == 23: # HISEQ-002:203:HBERLADXX:2:2213:8161:45667
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 31)
                    assert(tk_readpairs.is_split_consecutive(read1, read2))
                    assert(rp.is_first_half_first())
                    self.assertEqual(rp.sv_type, tk_readpairs.DEL_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (39225159 + 33 - tk_readpairs.SPLIT_LEN,
                                              39225159 + 33 + tk_readpairs.SPLIT_LEN + 1))
                    self.assertEqual(range2, (39399227 - tk_readpairs.SPLIT_LEN,
                                              39399227 + tk_readpairs.SPLIT_LEN + 1))
        in_bam.close()


    def test_pairs(self):
        mean_ins = 142
        std_ins = 81.5431220178
        max_ins = 386
        in_bam = tk_bam.create_bam_infile(PAIR_FILE)
        for idx, read in enumerate(in_bam):
            if idx % 2 == 0:
                read1 = read
            else:
                read2 = read
                rp = tk_readpairs.ReadPair(read1, read2)
                if idx == 1: # HISEQ-002:203:HBERLADXX:2:1109:13682:22490
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 0)
                    self.assertEqual(rp.sv_type, tk_readpairs.DEL_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    # Maximum distance between reads is max_ins - 2 * 88
                    self.assertEqual(range1, (189704422 + 88, 189704422 + max_ins - 88 + 1))
                    self.assertEqual(range2, (189783402 - max_ins + 2 * 88, 189783402 + 1))
                if idx == 13: # HISEQ-002:203:HBERLADXX:2:2201:3706:93604
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 0)
                    self.assertEqual(rp.sv_type, tk_readpairs.INV_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    self.assertEqual(range1, (162503343 - max_ins + 2 * 88, 162503343 + 1))
                    self.assertEqual(range2, (162627709 - max_ins + 2 * 88, 162627709 + 1))
                    self.assertEqual(tk_readpairs.get_pair_break_dist(read1, 162503243), 188)
                    self.assertEqual(tk_readpairs.get_pair_break_dist(read2, 162627609), 154)
                if idx == 19: # HISEQ-002:203:HBERLADXX:1:2214:1366:24373
                    self.assertEqual(tk_readpairs.map_qstart(read1), 0)
                    self.assertEqual(tk_readpairs.map_qstart(read2), 0)
                    self.assertEqual(rp.sv_type, tk_readpairs.DEL_STR)
                    range1, range2 = rp.get_break_ranges(max_ins)
                    # +8 because this read has 8bp of deletions with respect to the reference
                    # so the last aligned position on the reference is start + 88 + 8
                    self.assertEqual(range1, (78967094 + 88 + 8, 78967094 + max_ins - 88 + 1 + 8))
                    self.assertEqual(range2, (79036484 - max_ins + 2 * 88, 79036484 + 1))
                    self.assertEqual(tk_readpairs.get_pair_break_dist(read1, 78967294), 200)
                    self.assertEqual(tk_readpairs.get_pair_break_dist(read2, 79036284), 288)

                    # The ranges are single points. range2[1] is outside the valid range.
                    probs = rp.get_break_lr(np.arange(range1[0], range1[1]),
                                            np.ones((range1[1] - range1[0],)) * range2[1],
                                            0, lambda x: 0.0)
                    assert(np.all(probs == 0))
                    # probabilities returned should be decreasing because we "used up" all the
                    # insert size on the side of read2
                    probs = rp.get_break_lr(np.arange(range1[0], range1[1]),
                                            np.ones((range1[1] - range1[0],)) * range2[1],
                                            max_ins, lambda x: norm.logsf(x, mean_ins, std_ins))
                    self.assertEqual(list(probs), list(sorted(probs)[::-1]))
                    # probabilities returned should be increasing because the mean insert size is
                    # pretty large.
                    probs = rp.get_break_lr(np.arange(range1[0], range1[1]),
                                            np.ones((range1[1] - range1[0],)) * range2[0], 10000,
                                            lambda x: norm.logsf(x, 10000, std_ins))
                    self.assertEqual(list(probs), list(sorted(probs)))
        in_bam.close()
