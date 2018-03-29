#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import os
import os.path
import longranger.test as lr_test
import tenkit.bio_io as tk_io
from longranger.sv.phase_utils import *

TEST_FILE_DIR = lr_test.in_path('sv_phasing')

class TestSvPhaseUtils(lr_test.UnitTestBase):
    def setUp(self):
        pass

    def test_select_best_hap(self):
        test_sv_phasing_file = os.path.join(TEST_FILE_DIR, 'test_sv_phasing.tsv')
        sv_phasing_df = select_best_hap(test_sv_phasing_file)
        self.assertEqual(sv_phasing_df.loc[1010, 1].called_hap, '1')
        self.assertEqual(sv_phasing_df.loc[1010, 2].called_hap, '1')
        self.assertEqual(sv_phasing_df.loc[1195, 2].called_hap, '0')
        self.assertEqual(len(sv_phasing_df.loc[1195]), 1)
        self.assertEqual(not 92 in sv_phasing_df.index.levels[0], True)
        self.assertEqual(sv_phasing_df.loc[2034, 1].called_hap, '0')
        self.assertEqual(sv_phasing_df.loc[2034, 2].called_hap, '1')
        sv_phasing_df = select_best_hap(test_sv_phasing_file, True)
        self.assertEqual(sv_phasing_df.loc[92, 2].called_hap, '.')


    def test_match_to_gt_haps(self):
        gt_vcf_in = tk_io.VariantFileReader(os.path.join(TEST_FILE_DIR, 'phasing_gt.vcf.gz'))
        res = match_to_gt_haps(gt_vcf_in, gt_vcf_in, 'chr1', 0, 31)
        self.assertEqual(res[0], 0)
        res = match_to_gt_haps(gt_vcf_in, gt_vcf_in, 'chr1', 500, 531)
        assert(res is None)
        vcf_in = tk_io.VariantFileReader(os.path.join(TEST_FILE_DIR, 'phasing.vcf.gz'))
        res = match_to_gt_haps(gt_vcf_in, vcf_in, 'chr1', 0, 31)
        self.assertEqual(res[0], 0)
        res = match_to_gt_haps(gt_vcf_in, vcf_in, 'chr1', 31, 61)
        self.assertEqual(res[0], 1)
        res = match_to_gt_haps(gt_vcf_in, vcf_in, 'chr1', 61, 120)
        assert(res is None)
