#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import scipy.sparse
import math
import tenkit.test as tk_test
from .. import *
import martian

TEST_SPLIT_BAM = tk_test.in_path('test_split.bam')

martian.test_initialize(tk_test.out_path(""))

class TestGetReadpairLoci(tk_test.UnitTestBase):
    def setUp(test):
        pass

    def test_get_discordant_loci(self):
        # Reads 3,4,7,11 are neither secondary nor read1 so they should never be considered.
        # Reads 5, 8, 9, 10 are split
        # Reads 1, 2, 6 are rp and read1

        loci = get_discordant_loci(TEST_SPLIT_BAM, min_insert = 0, max_insert = 300, min_sv_len = 300)
        # Only reads 5,6,8,9,10 are included (1-based read indices)
        self.assertEqual(len(loci), 10)
        self.assertEqual(loci[0], ('chr20', 60173, 60773, (0, 1)))
        self.assertEqual(loci[1], ('chr20', 60609, 61209, (0, 2)))
        self.assertEqual(loci[2], ('chr20', 60329, 60929, (1, 1)))
        self.assertEqual(loci[3], ('chr20', 60794, 61394, (1, 2)))
        
        # Same as above plus read 1.
        loci = get_discordant_loci(TEST_SPLIT_BAM, min_insert = 100, max_insert = 300, min_sv_len = 300)
        self.assertEqual(len(loci), 12)
        
        # Only reads 8,9,10
        loci = get_discordant_loci(TEST_SPLIT_BAM, min_insert = 0, max_insert = 600, min_sv_len = 600)
        self.assertEqual(len(loci), 6)
        # Note that 'chr20' < 'chr5'
        self.assertEqual(loci[0], ('chr20', 60555, 61755, (0, 1)))
        self.assertEqual(loci[1], ('chr5', 132882251, 132883451, (0, 2)))
        
        loci = get_discordant_loci(TEST_SPLIT_BAM, min_insert = 0, max_insert = 300, 
            min_sv_len = 300, min_mapq = 60)
        # Only reads 5,6,8,9,10 are included (1-based read indices)
        self.assertEqual(len(loci), 4)
