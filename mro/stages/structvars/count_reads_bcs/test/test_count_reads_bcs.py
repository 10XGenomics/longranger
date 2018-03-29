#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import scipy.sparse
import math
import tenkit.test as tk_test
import tenkit.bio_io as tk_io
from .. import *
import martian

TEST_BAM = tk_test.in_path('test_count_reads.bam')
TEST_TARGETS = tk_test.in_path('test_targets_for_counting.txt')

martian.test_initialize(tk_test.out_path(""))

class TestCountReadsBcs(tk_test.UnitTestBase):
    def setUp(test):
        pass


    def test_get_non_overlapping_wins(self):
        starts = np.arange(0, 12, 2)
        stops = np.arange(2, 14, 2)
        sel = get_non_overlapping_wins(starts, stops)
        assert(np.all(sel == np.arange(6)))

        starts = np.arange(6)
        stops = np.arange(2, 8)
        sel = get_non_overlapping_wins(starts, stops)
        assert(list(sel) == [0, 2, 4])
        
        starts = np.array([0, 1])
        stops = np.array([5, 3])
        sel = get_non_overlapping_wins(starts, stops)


    def test_merge_bc_mat(self):
        test_mat = np.array([[0, 1, 2, 3, 2, 5], [3, 5, 4, 6, 2, 1]])
        bc_mat = sp.lil_matrix(test_mat)
        starts = np.arange(0, 12, 2)
        stops = np.arange(1, 13, 2)
        
        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 0, 0)
        assert(np.all(merged_mat.todense() == bc_mat.todense()))
        assert(np.all(new_starts == starts))
        assert(np.all(new_stops == stops))

        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, -100, 0)
        assert(np.all(merged_mat.todense() == bc_mat.todense()))
        assert(np.all(new_starts == starts))
        assert(np.all(new_stops == stops))
        
        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 4, 0)
        assert(np.all(merged_mat.todense() == bc_mat.todense()))
        assert(np.all(new_starts == starts))
        assert(np.all(new_stops == stops))

        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 4, -10)
        assert(np.all(merged_mat.todense() == bc_mat.todense()))
        assert(np.all(new_starts == starts))
        assert(np.all(new_stops == stops))

        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 4, 2)
        correct_mat = np.array([[1, 2, 3, 7], [8, 4, 6, 3]])
        correct_starts = np.array([0, 4, 6, 8])
        correct_stops = np.array([3, 5, 7, 11])
        assert(np.all(correct_mat == merged_mat.todense()))
        assert(len(new_starts) == correct_mat.shape[1])
        assert(len(new_stops) == correct_mat.shape[1])
        assert(np.all(new_starts == correct_starts))
        assert(np.all(new_stops == correct_stops))
        
        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 50, 2)
        correct_mat = np.array([[1, 5, 7], [8, 10, 3]])
        correct_starts = np.array([0, 4, 8])
        correct_stops = np.array([3, 7, 11])
        assert(np.all(correct_mat == merged_mat.todense()))
        assert(len(new_starts) == correct_mat.shape[1])
        assert(len(new_stops) == correct_mat.shape[1])
        assert(np.all(new_starts == correct_starts))
        assert(np.all(new_stops == correct_stops))
        
        merged_mat, new_starts, new_stops = merge_bc_mat(bc_mat, starts, stops, 50, 20)
        correct_mat = np.array([[13], [21]])
        correct_starts = np.array([0])
        correct_stops = np.array([11])
        assert(np.all(correct_mat == merged_mat.todense()))
        assert(len(new_starts) == correct_mat.shape[1])
        assert(len(new_stops) == correct_mat.shape[1])
        assert(np.all(new_starts == correct_starts))
        assert(np.all(new_stops == correct_stops))
        

    def get_local_bcs(self, bam_filename, chrom, starts, stops, blacklist = set([])):
        bc_map = {}
        inv_bc_map = {}
        input_bam = tk_bam.create_bam_infile(bam_filename)
        for start, stop in zip(starts, stops):
            for read in input_bam.fetch(str(chrom), start, stop):
                if read.is_duplicate or read.is_secondary:
                    continue
                bc = tk_io.get_read_barcode(read)
                if not(bc is None or bc == '') and not bc in bc_map and not bc in blacklist:
                    inv_bc_map[len(bc_map)] = bc
                    bc_map[bc] = len(bc_map)
        return (bc_map, inv_bc_map)

    
    def test_create_bc_matrix(self):
        starts = [287500, 287600]
        stops = [287600, 287700]
        bc_map, inv_bc_map = self.get_local_bcs(TEST_BAM, 'chr19', starts, stops)
        
        correct_bc_mat = [[2, 1, 1, 2, 1, 1, 1, 1, 1]]
        correct_bc_mat = np.array(correct_bc_mat).T
        
        bc_mat, new_starts, new_stops = create_bc_matrix_step(TEST_BAM, 'chr19', 
            [287500], [287700], 1000, 1000, bc_map, 
            min_mapq = 60, read1_only = False, no_split = False)
        assert(np.all(bc_mat.todense() == correct_bc_mat))
        assert(list(new_starts) == [287500])
        assert(list(new_stops) == [287700])
        
        bc_mat2, new_starts, new_stops = create_bc_matrix_step(TEST_BAM, 'chr19', 
            starts, stops, 1000, 1000, bc_map, 
            min_mapq = 60, read1_only = False, no_split = False)
        assert(list(new_starts) == starts)
        assert(list(new_stops) == stops)
        assert(np.all(bc_mat2.todense().sum(axis = 1) == bc_mat.todense()))

        bc_mat2, new_starts, new_stops = create_bc_matrix_step(TEST_BAM, 'chr19', 
            starts, stops, 10, 10, bc_map, 
            min_mapq = 60, read1_only = False, no_split = False)
        assert(list(new_starts) == list(np.arange(287500, 287700, 10)))
        assert(list(new_stops) == list(np.arange(287510, 287710, 10)))
        assert(np.all(bc_mat2.todense().sum(axis = 1) == bc_mat.todense()))
        
        bc_mat2, new_starts, new_stops = create_bc_matrix_step(TEST_BAM, 'chr19', 
            starts, stops, 10, 2, bc_map, 
            min_mapq = 60, read1_only = False, no_split = False)
        correct_starts = np.arange(287500, 287700, 2)
        correct_stops = np.concatenate((np.minimum(np.arange(287500, 287600, 2) + 10, 287600),
            np.minimum(np.arange(287600, 287700, 2) + 10, 287700)))
        assert(list(new_starts) == list(correct_starts))
        assert(list(new_stops) == list(correct_stops))
        sel_win = get_non_overlapping_wins(new_starts, new_stops)
        assert(np.all(bc_mat2[:, sel_win].todense().sum(axis = 1) == bc_mat.todense()))

