#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
import os.path
import tenkit.test as tk_test
from .. import *


IN_BAM= tk_test.in_path('test_attach_bc_vars.bam')
OUT_VCF = tk_test.out_path('test_call_variants.vcf')

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_call_variants(self):
        test_locus = "chr1:10000..20000"
        args = martian.Record({
                'input': IN_BAM,
                'locus': test_locus,
                'reference_path': 'hg19',
                'targets_file': None,
                'restrict_locus': None,
                'coverage': None,
                'max_coverage': None,
                'variant_mode':'freebayes'})
        outs = martian.Record({'default': OUT_VCF})

        main(args, outs)
        self.assertTrue(os.path.exists(OUT_VCF))
