#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
from .. import *
import tenkit.test as tk_test
import tenkit.seq as tk_seq
import tenkit.hdf5

bam_in_file = tk_test.in_path('test_map_rate.bam')

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        self.bam_in = tk_bam.create_bam_infile(bam_in_file)

    def _test_summary_metrics(self):

        read_info_out = tk_test.out_path("read_info.h5")

        insert_size_dists, nearest_targ_dists, summary_metrics, bc_table, mapq_counts, insert_size_hist = \
                compute_basic_stats(self.bam_in, {}, read_h5_out=read_info_out)

        summary = summary_metrics

        self.assertEqual(summary["mapped_bases"], 6)
        self.assertEqual(summary["mean_dup_rate"], 1.0)
        self.assertEqual(summary["num_reads"], 6)
        self.assertEqual(summary["total_bases"], 12)

        p = tenkit.hdf5.read_data_frame(read_info_out)
        self.assertEqual(p.shape[0], 3)

    def test_barcode_counts(self):
        bam_bc_file = tk_test.in_path("attach_bcs/attach_bcs_output.bam")
        read_info_out = tk_test.out_path("read_info.h5")
        barcode_whitelist = tk_seq.load_barcode_whitelist("737K-april-2014")
        bam_in = tk_bam.create_bam_infile(bam_bc_file)
        r = compute_basic_stats(bam_in,
                {}, 2000,
                bam_in.references,
                barcode_whitelist=barcode_whitelist,
                read_h5_out=read_info_out)
        # insert_size_dists, nearest_targ_dists, summary_metrics, bc_table, mapq_counts, insert_size_hist = r
        misc_sm, bc_sms = r

        # Look at the barcode results -- there should be a raw bc count for each read pair
        # n_raw_bcs = bc_table["count"].sum()
        n_reads = len([ x for x in tk_bam.create_bam_infile(bam_bc_file) ])

        # self.assertEqual(n_raw_bcs, n_reads / 2)

        # Load the per-cluster table -- there should be a row for each read pair
        read_info = tenkit.hdf5.read_data_frame(read_info_out)

        self.assertEqual(read_info.shape[0], n_reads / 2)

    def test_targets(self):
        bam_bc_file = tk_test.in_path("namesort_test.bam")
        read_info_out = tk_test.out_path("read_info.h5")
        barcode_whitelist = tk_seq.load_barcode_whitelist("737K-april-2014")

        targets_filename = tk_test.in_path('agilent_kinome_targs.bed')
        targets_file = open(targets_filename, 'r')
        target_regions = tk_io.get_target_regions(targets_file)

        bam_in = tk_bam.create_bam_infile(bam_bc_file)
        r = compute_basic_stats(bam_in,
                target_regions, 1000,
                bam_in.references,
                barcode_whitelist=barcode_whitelist,
                read_h5_out=read_info_out)
        # insert_size_dists, nearest_targ_dists, summary_metrics, bc_table, mapq_counts, insert_size_hist = r
        misc_sm, bc_sms = r

        nearest_targ_dists = bc_sms.get('nearest_targ_dists')
        maxTargetDist = max(nearest_targ_dists.get_summarizer(60).dict.keys())
        minTargetDist = min(nearest_targ_dists.get_summarizer(60).dict.keys())

        self.assertEqual(minTargetDist, 130)
        self.assertEqual(maxTargetDist, 10000)

    def test_target_finding(self):
        # Check a few targets by hand

        targets_filename = tk_test.in_path('agilent_kinome_targs.bed')
        targets_file = open(targets_filename, 'r')
        target_regions = tk_io.get_target_regions(targets_file)

        chr1_regions = target_regions['chr1']
        chr1_list = chr1_regions.get_region_list()
        test_reg = chr1_list[0]

        dist = get_read_regions_dist(0, test_reg.start-10, chr1_regions)
        self.assertEqual(dist, 10)

        dist = get_read_regions_dist(0, test_reg.start+10, chr1_regions)
        self.assertEqual(dist, 0)

        dist = get_read_regions_dist(test_reg.end+2, test_reg.end+3, chr1_regions)
        self.assertEqual(dist, 2)

        dist = get_read_regions_dist(test_reg.end-2, test_reg.end+3, chr1_regions)
        self.assertEqual(dist, 0)



