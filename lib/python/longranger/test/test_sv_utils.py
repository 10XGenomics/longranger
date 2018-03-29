#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
import longranger.test as lr_test
import tenkit.regions as tk_regions
from longranger.sv.utils import *
import longranger.sv.io

class TestSvUtils(lr_test.UnitTestBase):
    def setUp(self):
        pass


    def test_sort_and_merge(self):
        in_filename = lr_test.in_path('test_merge_regions.bed')
        regions = []
        with open(in_filename, 'r') as f:
            for line in f:
                fields = line.strip().split()
                regions.append((fields[0], int(fields[1]), int(fields[2])))

        out_regions = sort_and_merge(regions, 1000)
        # File created using BedTool's slopBed and mergeBed.
        out_filename = lr_test.in_path('test_merge_regions_d1000.bed')
        with open(out_filename, 'r') as f:
            for idx, line in enumerate(f):
                fields = line.strip().split()
                assert((fields[0], int(fields[1]), int(fields[2])) == out_regions[idx])

        out_regions = sort_and_merge(regions, 0)
        out_filename = lr_test.in_path('test_merge_regions_d0.bed')
        with open(out_filename, 'r') as f:
            for idx, line in enumerate(f):
                fields = line.strip().split()
                assert((fields[0], int(fields[1]), int(fields[2])) == out_regions[idx])


    def test_bc_map_from_inv(self):
        inv_bc_map = {0:'a', 1:'b', 3:'c'}
        bc_map = bc_map_from_inv(inv_bc_map, np.ones((4,), dtype = np.bool))
        assert(bc_map['a'] == 0)
        assert(bc_map['b'] == 1)
        assert(bc_map['c'] == 3)
        assert(len(bc_map) == 3)
        bc_map = bc_map_from_inv(inv_bc_map, np.array([0, 0, 0, 1], dtype = np.bool).flatten())
        assert(len(bc_map) == 1)
        assert(bc_map['c'] == 3)


    def test_get_nx_bc(self):
        read_counts = np.array([1,2,3,4,5,6,7,8,9,10])
        res, rank = get_nx_bcs(read_counts, 0)
        assert(np.all(res == False))
        assert(list(rank) == list((read_counts - 1)[::-1]))
        res, rank = get_nx_bcs(read_counts, -10)
        assert(np.all(res == False))
        res, rank = get_nx_bcs(read_counts, 0.0001)
        assert(list(res) == [False, False, False, False, False, False, False, False, False, True])
        res, rank = get_nx_bcs(read_counts, 100)
        assert(np.all(res))
        res, rank = get_nx_bcs(read_counts, 110)
        assert(np.all(res))
        res, rank = get_nx_bcs(read_counts, 50)
        assert(list(res) == [False, False, False, False, False, False, True, True, True, True])


    def test_loci_to_region_map(self):
        loci = [('chr1', [1000], [2000]), ('chr1', [2100], [2200]),
            ('chr2', [100, 500], [200, 700]),
            ('chr3', [1000, 900], [2000, 3000]), ('chr3', [2500], [2600])]
        regions = loci_to_region_map(loci)
        assert(regions['chr1'].starts == [1000, 2100])
        assert(regions['chr1'].ends == [2000, 2200])
        assert(regions['chr2'].starts == [100, 500])
        assert(regions['chr2'].ends == [200, 700])
        assert(regions['chr3'].starts == [900, 1000, 2500])
        assert(regions['chr3'].ends == [3000, 2000, 2600])


    def test_cluster_loci(self):
        loci = [('chr1', 1000, 2000, 1), ('chr2', 1000, 2000, 2), ('chr1', 2100, 2200, 3),
            ('chr2', 800, 950, 4), ('chr3', 1000, 2000, 5), ('chr3', 900, 3000, 6),
            ('chr3', 2500, 2600, 7)]

        clusters, _ = cluster_loci(loci, 0)
        assert(clusters[0] == [1])
        assert(clusters[1] == [3])
        assert(clusters[2] == [4])
        assert(clusters[3] == [2])
        assert(clusters[4] == [6, 5, 7])

        clustersB, _ = cluster_loci(loci, -10)
        assert(clustersB == clusters)

        clusters, _ = cluster_loci(loci, 50)
        assert(clusters[0] == [1])
        assert(clusters[1] == [3])
        assert(clusters[2] == [4, 2])
        assert(clusters[3] == [6, 5, 7])

        clusters, _ = cluster_loci(loci, 500)
        assert(clusters[0] == [1, 3])
        assert(clusters[1] == [4, 2])
        assert(clusters[2] == [6, 5, 7])

        clusters, _ = cluster_loci(loci, 500, 1000)
        assert(clusters[0] == [1])
        assert(clusters[1] == [3])
        assert(clusters[2] == [4, 2])
        assert(clusters[3] == [6, 5])
        assert(clusters[4] == [7])


    # Cases that should be tested
    # True bbbb-----------bbbb
    # Pred b--b
    # Pred   bbbb----------bbbb
    # Pred       bbbb-----------bbbb
    def test_compare_breaks(self):
        true_bed_filename = lr_test.in_path('test_compare_breaks_true.bedpe')
        pred_bed_filename = lr_test.in_path('test_compare_breaks_pred.bedpe')
        pred_to_match, true_to_match, is_sv_filtered = compare_breaks(pred_bed_filename)
        assert(len(pred_to_match) == 0)
        assert(len(true_to_match) == 0)
        assert(len(is_sv_filtered) == 0)

        pred_to_match, true_to_match, is_sv_filtered  = compare_breaks(pred_bed_filename, true_bed_filename, max_dist = 0)
        assert(len(is_sv_filtered) == 0)
        assert(len(true_to_match) == 1)
        assert(true_to_match['true_1'] == set(['pred_3']))
        assert(len(pred_to_match) == 1)
        assert(pred_to_match['pred_3'] == set(['true_1']))

        pred_to_match, true_to_match, is_sv_filtered  = compare_breaks(pred_bed_filename, true_bed_filename,
                                                                              max_dist = 500)
        assert(len(is_sv_filtered) == 0)
        assert(len(true_to_match) == 2)
        assert(true_to_match['true_1'] == set(['pred_3', 'pred_4']))
        assert(true_to_match['true_2'] == set(['pred_6']))
        assert(len(pred_to_match) == 3)
        assert(pred_to_match['pred_6'] == set(['true_2']))
        assert(pred_to_match['pred_3'] == set(['true_1']))
        assert(pred_to_match['pred_4'] == set(['true_1']))

        pred_to_match, true_to_match, is_sv_filtered = compare_breaks(pred_bed_filename, true_bed_filename,
            max_dist = 500, window_loci = [('chr4', [1000], [3000])])
        assert(len(is_sv_filtered) == 1)
        assert(len(true_to_match) == 1)
        assert(true_to_match['true_1'] == set(['pred_3', 'pred_4']))
        assert(len(pred_to_match) == 2)
        assert(pred_to_match['pred_3'] == set(['true_1']))
        assert(pred_to_match['pred_4'] == set(['true_1']))


    def test_merge_breaks(self):
        true_bed_filename = lr_test.in_path('test_merge_breaks.bedpe')
        out_bed_filename = lr_test.out_path('test_merge_breaks_out.bedpe')
        merge_breaks(true_bed_filename, out_bed_filename, merge_win = 50)
        self.assert_files_same(out_bed_filename, lr_test.in_path('test_merge_breaks_win50.bedpe'), False)
        merge_breaks(true_bed_filename, out_bed_filename, merge_win = 500)
        self.assert_files_same(out_bed_filename, lr_test.in_path('test_merge_breaks_win500.bedpe'), False)
        merge_breaks(true_bed_filename, out_bed_filename, merge_win = 500, max_nmates = 0)
        with open(out_bed_filename, 'r') as f:
            lines = [line for line in f.readlines() if not line.startswith('#')]
        assert(len(lines) == 0)

        res_df = merge_breaks(true_bed_filename, out_bed_filename, merge_win = 500, max_range = 1000)
        self.assertEqual(set(res_df['name']), set(['F', 'B', 'D']))
        res_df = merge_breaks(true_bed_filename, out_bed_filename, merge_win = 50, max_range = 1000)
        self.assertEqual(set(res_df['name']), set(['F', 'B', 'D', 'A', 'C']))

        res_df = merge_breaks(lr_test.in_path('test_merge_breaks_chain.bedpe'), out_bed_filename,
            merge_win = 1000, max_range = 1000)
        self.assertEqual(list(res_df['name']), ['B', 'D'])
        self.assertEqual(list(res_df['info']), ['NMATES1=2;NMATES2=1', 'NMATES1=2;NMATES2=1'])


    def compare_dfs_without_names(self, file1, file2):
        df1 = longranger.sv.io.read_sv_bedpe_to_df(file1).drop('name', 1)
        df2 = longranger.sv.io.read_sv_bedpe_to_df(file2).drop('name', 1)
        assert(np.all(df1 == df2))


    def test_merge_multiple_breaks(self):
        true_bed_filename = lr_test.in_path('test_merge_breaks.bedpe')
        out_bed_filename = lr_test.out_path('test_merge_breaks_out.bedpe')
        merge_multiple_breaks([true_bed_filename, true_bed_filename], out_bed_filename, merge_win = 50)
        # Remove the names, because these might not match
        self.compare_dfs_without_names(out_bed_filename, lr_test.in_path('test_merge_breaks_win50.bedpe'))

        true_bed_filename1 = lr_test.in_path('test_merge_breaks1.bedpe')
        true_bed_filename2 = lr_test.in_path('test_merge_breaks2.bedpe')
        merge_multiple_breaks([true_bed_filename1, true_bed_filename2], out_bed_filename, merge_win = 50)
        self.compare_dfs_without_names(out_bed_filename, lr_test.in_path('test_merge_breaks_win50.bedpe'))
        merge_multiple_breaks([true_bed_filename1, true_bed_filename2], out_bed_filename, merge_win = 500)
        self.compare_dfs_without_names(out_bed_filename, lr_test.in_path('test_merge_breaks_win500.bedpe'))


    def test_compare_multiple_breaks(self):
        filenames = ['test_merge_breaks.bedpe', 'test_merge_breaks1.bedpe',
                     'test_merge_breaks2.bedpe']
        in_bedpes = [lr_test.in_path(s) for s in filenames]
        merged_df = compare_multiple_breaks(in_bedpes, [0, 1, 2],
            lr_test.out_path('test_compare_breaks_out.bedpe'))
        merged_df = merged_df.sort(['qual', 'chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2'],
            ascending = [0, 1, 1, 1, 1, 1, 1])
        assert(np.all(merged_df['0_filtered'] == False))
        assert(np.all(merged_df['1_filtered'] == False))
        assert(np.all(merged_df['2_filtered'] == False))
        assert(np.all(merged_df['0_correct'] == False))
        assert(np.all(merged_df['1_correct'] == False))
        assert(np.all(merged_df['2_correct'] == False))
        self.assertEqual(list(merged_df['0_qual']), [50, 20, 20, 10, 10])
        self.assertEqual(list(merged_df['1_qual']), [0, 20, 0, 10, 10])
        self.assertEqual(list(merged_df['2_qual']), [50, 0, 20, 5, 0])
        self.assertEqual(list(merged_df['0_dist']), [-1, -1, 100, -1, -1])
        self.assertEqual(list(merged_df['1_dist']), [0, -1, 0, -1, -1])
        self.assertEqual(list(merged_df['2_dist']), [-1, 0, 100, -1, 0])


    def test_merge_region_maps(self):
        map1 = {}
        map1['chr1'] = tk_regions.Regions([(0, 10), (30, 40)])
        map1['chr2'] = tk_regions.Regions([(0, 20)])
        map2 = {}
        map2['chr1'] = tk_regions.Regions([(0, 15), (50, 60)])
        merged_map = merge_region_maps(map1, map2)
        self.assertEqual(list(merged_map['chr1'].starts), [0, 30, 50])
        self.assertEqual(list(merged_map['chr1'].ends), [15, 40, 60])
        self.assertEqual(list(merged_map['chr2'].starts), [0])
        self.assertEqual(list(merged_map['chr2'].ends), [20])


    def test_region_coverage(self):
        regions = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        self.assertEqual(list(region_coverage(regions, 5)), [5, 5, 0, 0, 5, 5, 0, 5, 5])
        self.assertEqual(list(region_coverage(regions, 10)), [10, 0, 10, 5, 5])
        regions = tk_regions.Regions([[0, 1], [1, 2], [2, 3]])
        self.assertEqual(list(region_coverage(regions, 10)), [3])
        self.assertEqual(list(region_coverage(regions, 1)), [1, 1, 1])
        regions = tk_regions.Regions([[0, 0], [1, 1], [2, 2]])
        self.assertEqual(list(region_coverage(regions, 1)), [0, 0])
        self.assertEqual(list(region_coverage(regions, 10)), [0])
        regions = tk_regions.Regions([[0, 1], [100, 110]])
        self.assertEqual(list(region_coverage(regions, 200)), [11])


    def test_region_cum_coverage_map(self):
        regions = {}
        regions['chr1'] = tk_regions.Regions([[0, 10], [20, 30], [35, 45]])
        cov_map = region_cum_coverage_map(regions, 10)
        self.assertEqual(list(cov_map['chr1']), [10, 10, 20, 25, 30])
        cov_map = region_cum_coverage_map(regions, 5)
        self.assertEqual(list(cov_map['chr1']), [5, 10, 10, 10, 15, 20, 20, 25, 30])
        regions['chr1'] = tk_regions.Regions([[0, 0], [1, 1], [2, 2]])
        cov_map = region_cum_coverage_map(regions, 1)
        self.assertEqual(list(cov_map['chr1']), [0, 0])


    def assert_files_same(self, file1, file2, match_order = True):
        with open(file1, 'r') as f1:
            lines1 = [line for line in f1.readlines() if not line.startswith('#')]
        with open(file2, 'r') as f2:
            lines2 = [line for line in f2.readlines() if not line.startswith('#')]
        if match_order:
            assert(lines1 == lines2)
        else:
            assert(set(lines1) == set(lines2))


if __name__ == '__main__':
    lr_test.run_tests()
