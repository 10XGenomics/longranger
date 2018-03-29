#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Base class for tests.  Allows shared testing tools to be built
# as well as global specification of locations.
#

import unittest
import os
import os.path
import subprocess

CODE_PATH = os.path.dirname(os.path.abspath(__file__)) + '/'
TEST_FILE_IN_DIR  = CODE_PATH + 'test_files/inputs/'
TEST_FILE_OUT_DIR = CODE_PATH + 'test_files/outputs/'


def in_path(f):
    ''' Construct a test input file path from a relative path to a test file '''
    return os.path.join(TEST_FILE_IN_DIR, f)

def out_path(f):
    ''' Construct a test output fule path from a relative path to an output file '''
    return os.path.join(TEST_FILE_OUT_DIR, f)

class UnitTestBase(unittest.TestCase):
    def setUp(self):
        self.clear_directory()

    def clear_directory(self):
        if os.path.exists(TEST_FILE_OUT_DIR):
            subprocess.call(["rm", "-Rf", TEST_FILE_OUT_DIR])
        if not os.path.exists(TEST_FILE_OUT_DIR):
            os.mkdir(TEST_FILE_OUT_DIR)

    def tearDown(self):
        self.clear_directory()

    def assertApproxEqual(self, v1, v2, precision=1e-6):
        if v1 == v2:
            return

        err = float(abs(v1 - v2)) / v1
        self.assertLess(err, precision)

    def assertEqualDataFrames(self, df1, df2):
        self.assertEqual(df1.shape, df2.shape)
        self.assertTrue((df1.sort(axis=1) == df2.sort(axis=1)).all().all())
