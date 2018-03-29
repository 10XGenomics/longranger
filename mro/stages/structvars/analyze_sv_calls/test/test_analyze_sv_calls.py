#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Test the demultiplexer

import subprocess
import os
import os.path
import numpy as np
import json
import tenkit.pandas as pd
import tenkit.test as tk_test
from tenkit.constants import TEST_FILE_IN_DIR, TEST_FILE_OUT_DIR
from .. import *

import martian

def infile(path):
    return os.path.join(TEST_FILE_IN_DIR, 'analyze_sv_calls', path)

def outfile(path):
    return os.path.join(TEST_FILE_OUT_DIR, path)

mrs_job_name = 'test_analyze_sv_calls'
base_out_dir = outfile('')
job_dir = os.path.join(base_out_dir, mrs_job_name,  'ANALYZE_SV_CALLS', 'fork0', 'files')

class TestAnalyzeSvCalls(tk_test.UnitTestBase):

    def setUp(self):
        #self.clear_directory()
        martian.test_initialize(outfile(''))

        self.args = {
            'max_frac_black':0.1,
            'min_dist_from_black':10000,
            'seg_dup_min_dist':10000,
            'min_sv_len':1,
            'qv':50,
            'variants':infile('test_analyze_sv_calls2.bedpe'),
            'gt_variants':infile('test_analyze_sv_calls_gt2.bedpe'),
            'sv_blacklist_regions':infile('test_empty_blacklist.bed'),
            'seg_dups':infile('test_empty_seg_dups.bedpe'),
        }


    def write_mro(self):
        tpl = """
            @include "_structvar_caller_stages.mro"
            call ANALYZE_SV_CALLS(
                detect_dists = [10000, 50000],
                reference_path = null,
                sv_support_fragments  = null,
                fragments             = null,
                coverage              = null,
                target_dists = [100, 10000],
                max_frac_black = %(max_frac_black)s,
                min_read_support = 0,
                min_dist_from_black = %(min_dist_from_black)s,
                seg_dup_min_dist = %(seg_dup_min_dist)s,
                min_sv_len = %(min_sv_len)s,
                targets = %(targets)s,
                variants = "%(variants)s",
                gt_variants = "%(gt_variants)s",
                sv_blacklist_regions = "%(sv_blacklist_regions)s",
                seg_dups = "%(seg_dups)s",
                min_call_qv_wgs = %(qv)s,
                min_call_qv_target = %(qv)s,
                max_nmates            = null,
                max_break_res         = null,
                min_rel_depth         = null,
            )
        """

        fn = os.path.join(base_out_dir, 'test.mro')
        with open(fn, 'w') as f:
            f.write(tpl % self.args)
        return fn


    def run_stage_wgs(self):
        self.args['targets'] = 'null'
        mro_file = self.write_mro()
        proc = subprocess.Popen(['mrs', mro_file, mrs_job_name], cwd = base_out_dir, stdout = subprocess.PIPE)
        out, err = proc.communicate()

        if proc.returncode != 0:
            print out
            raise Exception("mrs failed on during test")

    
    def run_stage_targets(self):
        self.args['targets'] = '"{}"'.format(infile('test_targets2.bed'))
        mro_file = self.write_mro()
        proc = subprocess.Popen(['mrs', mro_file, mrs_job_name], cwd = base_out_dir, stdout = subprocess.PIPE)
        out, err = proc.communicate()

        if proc.returncode != 0:
            print out
            raise Exception("mrs failed on during test")


    def test_wgs(self):
        self.clear_directory()
        self.run_stage_wgs()

        with open(os.path.join(job_dir, "summary.json")) as summary_file:
            summary = json.load(summary_file)

        calls = pd.read_csv(os.path.join(job_dir, 'call_tsv.tsv'), sep = '\t', header = 0, index_col = False)

        # 34 calls, 3 have two matches each
        self.assertEqual(len(calls), 37)

        for _, row in calls.iterrows():
            if np.isnan(row.match_dist):
                assert(np.isnan(row.match))

        self.assertEqual(summary['Q50_num_calls'], 34)
        self.assertEqual(summary['Q50_calledINV_num_calls'], 4)
        self.assertEqual(summary['Q50_calledDEL_num_calls'], 18)
        
        np.testing.assert_almost_equal(summary['Q50_frac_read_support'], 11/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_mean_break_res'], 1235.294118, decimal = 6)
                
        self.assertEqual(summary['tier1_DEL_num_gt'], 4)
        self.assertEqual(summary['tier2_DEL_num_gt'], 6)
        self.assertEqual(summary['tier1_INV_num_gt'], 3)
        
        np.testing.assert_almost_equal(summary['Q50_tier1_DEL_sensitivity'], 3/4.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_DEL_sensitivity'], 5/6.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier1_INV_sensitivity'], 2/3.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier1_NA_sensitivity'], 14/16.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_NA_sensitivity'], 17/19.0, decimal = 10)

        # PPV using the full ground truth
        np.testing.assert_almost_equal(summary['Q50_tier2_NA_PPV'], 19/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_calledDEL_PPV'], 4/18.0, decimal = 10)
        
        # PPV using only tier1 deletions as ground truth
        np.testing.assert_almost_equal(summary['Q50_tier1_DEL_PPV'], 3/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier1_DEL_type_accuracy'], 1, decimal = 10)

        np.testing.assert_almost_equal(summary['Q50_tier2_INV_PPV'], 3/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_INV_type_accuracy'], 1/3.0, decimal = 10)
        
        np.testing.assert_almost_equal(summary['Q50_tier2_TRANS_PPV'], 12/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_TRANS_type_accuracy'], 1, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_TRANS_orient_accuracy'], 5/11.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier2_TRANS_orient_accuracy_missing'], 5/12.0, decimal = 10)        


    def test_targets(self):
        self.clear_directory()
        self.run_stage_targets()

        with open(os.path.join(job_dir, "summary.json")) as summary_file:
            summary = json.load(summary_file)

        calls = pd.read_csv(os.path.join(job_dir, 'call_tsv.tsv'), sep = '\t', header = 0, index_col = False)
        self.assertEqual(len(calls), 37)

        for _, row in calls.iterrows():
            if np.isnan(row.match_dist):
                assert(np.isnan(row.match))

        self.assertEqual(summary['Q50_num_calls'], 34)
        np.testing.assert_almost_equal(summary['Q50_frac_read_support'], 11/34.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_mean_break_res'], 1235.294118, decimal = 6)

        self.assertEqual(summary['tier1_DEL_num_gt_feasible100'], 1)
        self.assertEqual(summary['tier2_DEL_num_gt_feasible100'], 1)
        self.assertEqual(summary['tier1_INV_num_gt_feasible100'], 2)
        self.assertEqual(summary['tier1_DEL_num_gt_feasible10000'], 1)
        self.assertEqual(summary['tier2_DEL_num_gt_feasible10000'], 1)
        self.assertEqual(summary['tier1_INV_num_gt_feasible10000'], 3)
        self.assertEqual(summary['tier2_TRANS_num_gt_feasible10000'], 1)

        np.testing.assert_almost_equal(summary['Q50_tier1_INV_sensitivity_feasible100'], 0.5, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier1_INV_sensitivity_feasible10000'], 2/3.0, decimal = 10)
        np.testing.assert_almost_equal(summary['Q50_tier1_TRANS_sensitivity_feasible10000'], 1, decimal = 10)
