#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Test infer fragments
import pysam
import tenkit.test as tk_test
import tenkit.seq as tk_seq
from .. import *

import martian


class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_infer_fragment(self):

        f = [('a', 1, 10), ('a', 100, 110), ('a', 109, 150), ('a', 609, 700), ('a', 1300, 1400)]
        frags = infer_fragments(f, 500)

        self.assertEqual(frags, [('a', 1, 700, 4), ('a', 1300, 1400, 1)])
