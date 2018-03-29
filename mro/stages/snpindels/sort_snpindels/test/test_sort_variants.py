#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Simple test of variant sorting
#
import tenkit.test as tk_test
import tenkit.bio_io as tk_io
from .. import *

var_filename = tk_test.in_path('test_canonicalizer.vcf')
sort_var_filename = tk_test.out_path('default.vcf.gz')

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_sort_variants(self):
        args = martian.Record({"input": var_filename})
        outs = martian.Record({"default": sort_var_filename})
        main_sort_variants(args, outs)

        vfr = tk_io.VariantFileReader(sort_var_filename)
        var_pos = [ r.POS for r in vfr.record_getter() ]

        sort_var_pos = list(var_pos)
        sort_var_pos.sort()

        self.assertEqual(var_pos, sort_var_pos)

if __name__ == '__main__':
    tk_test.run_tests()
