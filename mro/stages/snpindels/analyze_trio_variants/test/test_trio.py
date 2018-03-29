#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Test trio matching coe
#

import tenkit.test as tk_test
from .. import trio

class TestTrios(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_mark_mins(self):
        vals = [(1,1), (0,0), (0,1), None]
        at_min = trio.mark_mins(vals, lambda x: x[0])
        self.assertEqual(at_min, [False, True, True, False])

        vals = [None, None]
        at_min = trio.mark_mins(vals, lambda x: x[0])
        self.assertEqual(at_min, [False, False])


    def test_parallel_iter1(self):
        i1 = [(('a', 0), 1), (('a', 1), 7), (('a', 10), 2)]
        i2 = [(('a', 0), 1), (('a', 1), 7), (('a', 10), 2)]
        i3 = [(('a', 0), 1), (('a', 1), 7), (('a', 10), 2), (('a', 12), 4)]

        v = list(trio.parallel_iter([i1,i2,i3], lambda x: x[0]))

        self.assertEqual(len(v), 4)
        self.assertEqual(v[-1], [None, None, (('a', 12), 4)])


    def test_parallel_iter2(self):
        i1 = [               (('a', 1), 7), (('a', 10), 2)]
        i2 = [(('a', 0), 1),                (('a', 10), 2)]
        i3 = [               (('a', 1), 7),                (('a', 12), 4)]

        v = list(trio.parallel_iter([i1,i2,i3], lambda x: x[0]))
        print v

        self.assertEqual(len(v), 4)
        self.assertEqual(v[0],  [None, (('a', 0), 1), None])
        self.assertEqual(v[-1], [None, None, (('a', 12), 4)])


    def test_compat(self):
        self.assertFalse(trio.compat(('T', 'T'), False, 0, ('TG', 'TG'), False, 0))
        self.assertFalse(trio.compat(('T', 'T'), True,  0, ('TG', 'TG'), False, 0))


    def test_hmm(self):
        # define some variants:
        # (genotype, phased, phase_set)
        aa0 = (['A','A'], True,   1)
        tt0 = (['T','T'], True,   1)
        at0 = (['A','T'], True,   1)
        at1 = (['A','T'], False,  1)
        at2 = (['A','T'], True,   2)

        # Test 1: switch penalty works properly 
        # (i.e. maintain phase blocks, but don't worry about unphased variants)
        (states, errors, scores) = trio.infer_trio_pattern([
            [aa0,at0,at0],
            [aa0,at0,at0],
            [at0,aa0,at0], # flipped, but same phase set
            [at1,aa0,at1], # unphased - no problen
            [aa0,at2,at2], # new phase set
        ])

        self.assertEqual(states, [(0,1,0), (0,1,0), (0,1,0), (1,1,1), (1,1,0)])
        self.assertEqual(errors, [set(), set(), set(['p2_error']), set(), set()])

        # Test 2: detect mendelian violations
        (states, errors, scores) = trio.infer_trio_pattern([
            [aa0,aa0,at0],
            [aa0,aa0,tt0],
        ])

        self.assertIn('mendelian_violation', errors[0])
        self.assertIn('mendelian_violation', errors[1])
