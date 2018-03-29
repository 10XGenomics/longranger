#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Test the demultiplexer

import subprocess
import os
import os.path
import tenkit.pandas as p
import numpy as np
import json
import tenkit.test as tk_test
import tenkit.hdf5
from tenkit.constants import TEST_FILE_IN_DIR, TEST_FILE_OUT_DIR
from .. import *

import martian

def infile(path):
    return os.path.join(TEST_FILE_IN_DIR, "analyze_snpindel_calls", path)

def outfile(path):
    return os.path.join(TEST_FILE_OUT_DIR, path)

mrs_job_name = "analyze_snpindel_calls"
base_out_dir = outfile("")
job_dir = os.path.join(base_out_dir, mrs_job_name,  "ANALYZE_SNPINDEL_CALLS", "fork0", "files")

class TestFunctions(tk_test.UnitTestBase):

    def setUp(self):
        self.clear_directory()
        martian.test_initialize(outfile(""))

        self.args = {
                'reference_path': 'hg19',
                'bam_file': infile("hg19_small.bam"),
                'ground_truth': infile("gt_sort.vcf.gz"),
                'input': infile("call_sort.vcf.gz"),
                'genes_file': infile("small_refseq_genes.tsv"),
                'long_switch_penalty_multiple': 5 }

    def write_mro(self, args, locus):
        tpl = """
            @include "_snpindel_caller_stages.mro"
            call ANALYZE_SNPINDEL_CALLS(
                reference_path = "%(reference_path)s",
                bam_file = "%(bam_file)s",
                coverage = null,
                genes_file = "%(genes_file)s",
                ground_truth = "%(ground_truth)s",
                input = "%(input)s",
                restrict_locus = "%(locus)s",
                targets_file = null,
                long_switch_penalty_multiple = %(long_switch_penalty_multiple)d,
                vc_conf_regions = null,
                fragment_phasing = null,
            )
        """
        args['locus'] = locus
        fn = os.path.join(base_out_dir, "test.mro")
        with open(fn, 'w') as f:
            f.write(tpl % args)

        return fn


    def run_stage(self, args, locus):
        mro_file = self.write_mro(args, locus)
        proc = subprocess.Popen(['mrs', mro_file, mrs_job_name], cwd=base_out_dir, stdout=subprocess.PIPE)
        out, err = proc.communicate()

        if proc.returncode != 0:
            print out
            raise Exception("mrs failed on during test")

    def test_gene_finder(self):
        gf = GeneFinder(infile("small_refseq_genes.tsv"))

        simple_regs = gf.simplify_region_set([(1,10, 'a'), (12,20,'a'), (13, 50, 'a'), (14, 60, 'a')])
        self.assertEqual(simple_regs, [(1, 10, 'a'), (12, 60, 'a')])


    def test_small(self):
        self.clear_directory()
        # Run a test on a tiny locus and compare with hand calculated metrics
        self.run_stage(self.args, "chr1:1000000..1800000")

        with open(os.path.join(job_dir, "summary.json")) as summary_file:
            summary = json.load(summary_file)

        # Check basic phasing metrics (verified with independent code)
        self.assertEqual(summary['N50_phase_block'], 65312)
        self.assertEqual(summary['longest_phase_block'], 66139)
        np.testing.assert_almost_equal(summary['mean_phase_block'], 2379.26373626)
        np.testing.assert_almost_equal(summary['fract_genes_phased'], 0.2727272727272727)
        np.testing.assert_almost_equal(summary['fract_genes_completely_phased'], 0.09090909090909091)
        np.testing.assert_almost_equal(summary['prob_snp_phased_in_gene'], 0.29032214854999666, decimal = 5)

        np.testing.assert_almost_equal(summary['short_switch_error'], 2.0/27, decimal = 10)
        np.testing.assert_almost_equal(summary['long_switch_error'], 0.0, decimal = 10)


    def test_basic(self):
        self.clear_directory()
        self.run_stage(self.args, "chr1:0..10000000")

        # Load all the output files
        variants = tenkit.hdf5.read_data_frame(os.path.join(job_dir, "variants.h5"))
        genes = p.read_csv(os.path.join(job_dir, "gene_stats.csv"))

        with open(os.path.join(job_dir, "summary.json")) as summary_file:
            summary = json.load(summary_file)

        print variants.head()

        # Check basic phasing metrics (verified with independent code)
        self.assertEqual(summary['N50_phase_block'], 362897)
        self.assertEqual(summary['longest_phase_block'], 1588263)
        np.testing.assert_almost_equal(summary['mean_phase_block'], 7067.99817851)

        # Hand calculated snp-weighted probabilities
        np.testing.assert_almost_equal(summary['fract_genes_phased'], 0.44155844155844154, decimal = 10)
        np.testing.assert_almost_equal(summary['fract_genes_completely_phased'], 0.1038961038961039, decimal = 10)
        self.assertEqual(summary['prob_snp_phased_in_gene'], 0.36608192963189135)
        self.assertEqual(summary['prob_snp_correct_in_gene'], 0.950884655460024)

        # Hand test some variant metrics
        gene_name = "KIAA1751"
        gene = genes[genes.gene == gene_name]

        self.assertEqual(gene.pair_phased_rate.values[0], float(79*78)/(143*142))
        np.testing.assert_almost_equal(gene.pair_correct_rate.values[0], float(77*76/2)/(78*77/2))

        self.assertEqual(gene.start.values[0], 1884751)
        self.assertEqual(gene.end.values[0], 1935276)


        # Hand test some variant metrics
        gene_name = "AJAP1"
        gene = genes[genes.gene == gene_name]

        np.testing.assert_almost_equal(gene.pair_phased_rate.values[0], float(149*148)/(184*183))
        np.testing.assert_almost_equal(gene.pair_correct_rate.values[0], float((143*142)/2 + 1) / (145*144/2))

        self.assertEqual(gene.start.values[0], 4715104)
        self.assertEqual(gene.end.values[0], 4843851)

        sample_variant_df = p.DataFrame({'in_obs':[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ],
                                        'in_gt':  [ 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 ],
                                   'variant_type':['S','S','D','I','S','C','D','I','D','S','I','I','D','S','C','D'],
                                    'FILTER':     [ 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 ],
                                'variant_length': [ 0 , 0 , -1, 1 , 0 , 1 , -1, 1 , -6, 0 , 6 , 7 , 5 , 0 , 6 , 1 ]})

        stats = {}

        filter_states = ["unfiltered", "filtered"]
        var_types = ["any", "snp", "insertion", "deletion", "complex", "insertion_lt5bp", "deletion_lt5bp"]

        for filter_condition in filter_states:
            for var_type in var_types:
                compute_snpindel_stats(stats, filter_condition, var_type, None, sample_variant_df, True)

        print stats
        #not covering all, but covering each dimension at least once
        self.assertEqual(stats['tps_filtered'], 2)
        self.assertEqual(stats['tps_unfiltered'], 4)
        self.assertEqual(stats['tps_filtered_snp'], 1)
        self.assertEqual(stats['tps_filtered_deletion'], 1)
        self.assertEqual(stats['tps_filtered_complex'], 0)
        self.assertEqual(stats['sensitivity_unfiltered'], 0.5)
        self.assertEqual(stats['ppv_unfiltered'], 0.5)
        self.assertEqual(stats['sensitivity_filtered'], 0.25)
        self.assertEqual(stats['ppv_filtered'], 0.5)
        self.assertEqual(stats['fps_filtered'], 2)
        self.assertEqual(stats['fns_unfiltered_insertion_lt5bp'], 1)
        self.assertEqual(stats['tps_filtered_deletion_lt5bp'], 1)
        self.assertEqual(stats['ppv_unfiltered_snp'], 0.5)
        self.assertEqual(stats['sensitivity_unfiltered_complex'], 0.5)
        self.assertEqual(stats['fns_unfiltered'], 4)
        self.assertEqual(stats['fns_filtered_deletion_lt5bp'], 1)

        # Test that some blocks of short switches are caught correctly
        ss = variants[variants.pos.isin([1529457,1529511])]
        self.assertTrue(ss.short_switch.all())

        # Used to be a long switch - now should be nothinh
        non_long = variants[variants.pos == 1529950]
        self.assertTrue((non_long.short_switch == False).all())
        self.assertTrue((non_long.long_switch == False).all())

        switch_fixes = variants[np.logical_and(variants.pos >= 3173228, variants.pos <= 3181098)]
        self.assertTrue((switch_fixes.long_switch == False).all())
        self.assertEqual(switch_fixes.short_switch.sum(), 5)
