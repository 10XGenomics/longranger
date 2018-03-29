#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for testing analyze_run.py
#

from pyfasta import Fasta
import tenkit.bam as tk_bam
import tenkit.test as tk_test
from tenkit.regions import Regions
from .. import *

## Input files
bam_in_file = tk_test.in_path('test_analyze_bias.bam')
fasta_dir = tk_test.in_path('fasta/')

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        self.bam_in = tk_bam.create_bam_infile(bam_in_file)

    def test_get_depth_info(self):
        ref_fasta = Fasta(fasta_dir + 'test/chr0.fa')
        chr0 = ref_fasta['chr0']
        confident_regions = Regions([(0,10000000)])

        reads = list(self.bam_in)
        r = get_depth_info(reads, "chr0", 0, len(chr0), None, confident_regions)
        (depth_df, summary_depth_info, confident_depth_info, target_info, target_cov) = r

        reads_dd = filter(lambda x: not x.is_duplicate, reads)
        r_dd = get_depth_info(reads_dd, "chr0", 0, len(chr0), None, confident_regions)
        (dd_depth_df, summary_depth_info_deduped, confident_depth_info, target_info, target_cov) = r_dd

        self.assertEqual(summary_depth_info, {0: 10, 1: 10, 2: 10, 3: 10})
        self.assertEqual(summary_depth_info_deduped, {0: 10, 1: 20, 2: 10})
        self.assertEqual(target_info, {})


        r = get_depth_info(reads, "chr0", 0, len(chr0), Regions([(5, 15)]), confident_regions)
        (target_depth_df, summary_depth_info, confident_depth_info, target_info, target_cov) = r

        self.assertEqual(summary_depth_info, {2: 5, 3: 5})

        self.assertEqual(len(target_depth_df), 10)
        self.assertEqual(len(target_cov), 1)
        self.assertEqual(target_cov['mean'][0], 2.5)
        self.assertEqual(sum(target_depth_df.coverage), target_info['on_target_bases'])


        r_dd = get_depth_info(reads_dd, "chr0", 0, len(chr0), Regions([(5, 15)]), confident_regions)
        (target_depth_df, summary_depth_info_deduped, confident_depth_info, target_info, target_cov) = r_dd

        self.assertEqual(summary_depth_info_deduped, {1: 5, 2: 5})



if __name__ == '__main__':
    tk_test.run_tests()
