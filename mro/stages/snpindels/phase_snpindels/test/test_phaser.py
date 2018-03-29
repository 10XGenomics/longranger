#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for testing phaser.py
#

import scipy.sparse
import math
import tenkit.test as tk_test
import tenkit.bio_io as tk_io
from tenkit.bio_io import VariantFileReader, VariantFileWriter
from ..phaser import *

CONTIG_FILE = tk_test.in_path('phasing/small_fragments.h5')
SNP_INPUT_VCF = tk_test.in_path('test_phasing_snps_sorted.vcf.gz')
INDEL_INPUT_VCF = tk_test.in_path('test_phasing_indels_sorted.vcf.gz')
OUTPUT_VCF = tk_test.out_path('test_phasing_output.vcf')
OUTPUT_TSV = tk_test.out_path("test_bc_hap_out.tsv")

import martian
martian.test_initialize(tk_test.out_path(""))

class TestPhaserBig(tk_test.UnitTestBase):
    def setUp(self):
        contig_file = tk_test.in_path("phasing/big_fragments.h5")
        snp_vfr = VariantFileReader(tk_test.in_path("phasing/default.vcf.gz"))
        self.p = Phaser(snp_vfr, contig_file, 'chr21', 10000000, 10500000)

    def test_calc_hap1_hap2_log_prob(self):
        # Generate a random assignment
        n_snps = self.p.hap1_matrix.shape[1]
        print n_snps
        base_assign1 = [i%2 for i in xrange(n_snps)]
        base_assign2 = [(i+1)%2 for i in xrange(n_snps)]

        for test_pos in xrange(n_snps):

            if test_pos % 5 == 0:
                base_assign = base_assign2
            else:
                base_assign = base_assign1

            # Calculate the log-prob for each config, the make sure calc_hap1_hap2 gives the same thing
            hap1_cols = [i for (i,v) in enumerate(base_assign) if v == 0]
            hap2_cols = [i for (i,v) in enumerate(base_assign) if v == 1]
            hap1_log_probs = self.p.hap1_matrix[:,hap1_cols].sum(axis=1).A1 + self.p.hap2_matrix[:,hap2_cols].sum(axis=1).A1
            hap2_log_probs = self.p.hap1_matrix[:,hap2_cols].sum(axis=1).A1 + self.p.hap2_matrix[:,hap1_cols].sum(axis=1).A1
            mix_log_probs = self.p.mix_matrix.sum(axis=1).A1
            base_prob = self.p.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_log_probs)

            r = self.p.calc_hap1_hap2_logprob(base_assign, test_pos)

    def test_hap_block(self):
        p = self.p
        block_size = 11
        n_snps = self.p.hap1_matrix.shape[1]

        for test_pos in xrange(n_snps - block_size):
            (assign, score) = p.get_best_hap_block(p.hap1_matrix, p.hap2_matrix, p.mix_matrix, test_pos, test_pos+block_size)
            print (assign, score)
            (assign_beam, score_beam) = p.get_best_hap_block_beam(p.hap1_matrix, p.hap2_matrix, p.mix_matrix, test_pos, test_pos+block_size, beam=10)
            print (assign_beam, score_beam)

            np.testing.assert_almost_equal(score, score_beam)


class TestPhaser(tk_test.UnitTestBase):
    def setUp(self):
        fragment_phasing = tk_test.out_path("phasing/fragment_phasing.tsv")
        snp_vfr = VariantFileReader(SNP_INPUT_VCF)
        self.mix_prior = 0.0001
        self.p = Phaser(snp_vfr, CONTIG_FILE, 'chr1', 0, 10, bc_mix_prob=self.mix_prior, min_junction_hap_conf=0.99, min_var_hap_conf=0.99, hap_block_size=5)

    def test_get_bc_contig_info(self):
        (sub_bc_indexes, bc_regions, starts, stops, bcs, _) = self.p.get_bc_contig_info('chr1', 0, 10, CONTIG_FILE)
        self.assertEqual(set(bc_regions.keys()), set(['AAAAAAAAAAAA', 'CCCCCCCCCCCC', 'GGGGGGGGGGGG']))
        self.assertEqual(set(sub_bc_indexes.keys()), set(['AAAAAAAAAAAA_0_4', 'CCCCCCCCCCCC_0_4', 'GGGGGGGGGGGG_0_4', 'GGGGGGGGGGGG_4_10']))

    def test_get_specific_contig(self):
        (sub_bc_indexes, bc_regions, starts, stops, bcs, _) = self.p.get_bc_contig_info('chr1', 0, 10, CONTIG_FILE)
        self.assertEqual(self.p.get_specific_contig(bc_regions, 'AAAAAAAAAAAA', 2), (0,4))
        self.assertEqual(self.p.get_specific_contig(bc_regions, 'GGGGGGGGGGGG', 2), (0,4))
        self.assertEqual(self.p.get_specific_contig(bc_regions, 'GGGGGGGGGGGG', 6), (4,10))

    def test_get_bc_variant_matrices(self):
        (sub_bc_indexes, bc_regions, starts, stops, bcs, _) = self.p.get_bc_contig_info('chr1', 0, 10, CONTIG_FILE)

        prob_var_wrong = math.pow(10.0, -255.0/10.0)
        hap1_prob = (1.0 - math.pow(10.0, -3.9))*(1.0 - math.pow(10.0, -3.9))
        hap2_prob = math.pow(10.0, -3.9)*math.pow(10.0, -3.9)
        mix_prob = 0.5*0.5
        non_error_prob = self.mix_prior*mix_prob + (1.0 - self.mix_prior)*(0.5*hap1_prob + 0.5*hap2_prob)
        self.assertApproxEqual(self.p.hap1_matrix[sub_bc_indexes['AAAAAAAAAAAA_0_4'],self.p.pos_indexes[2]], (1 - prob_var_wrong)*math.log(hap1_prob/non_error_prob))
        self.assertApproxEqual(self.p.hap2_matrix[sub_bc_indexes['AAAAAAAAAAAA_0_4'],self.p.pos_indexes[2]], (1 - prob_var_wrong)*math.log(hap2_prob/non_error_prob))
        self.assertApproxEqual(self.p.mix_matrix[sub_bc_indexes['AAAAAAAAAAAA_0_4'],self.p.pos_indexes[2]], (1 - prob_var_wrong)*math.log(mix_prob/non_error_prob))

        prob_var_wrong = math.pow(10.0, -255.0/10.0)
        hap1_prob = (1.0 - math.pow(10.0, -3.9))
        hap2_prob = math.pow(10.0, -3.9)
        mix_prob = 0.5
        non_error_prob = self.mix_prior*mix_prob + (1.0 - self.mix_prior)*(0.5*hap1_prob + 0.5*hap2_prob)
        self.assertApproxEqual(self.p.hap1_matrix[sub_bc_indexes['GGGGGGGGGGGG_0_4'],self.p.pos_indexes[3]], (1 - prob_var_wrong)*math.log(hap1_prob/non_error_prob))
        self.assertApproxEqual(self.p.hap2_matrix[sub_bc_indexes['GGGGGGGGGGGG_0_4'],self.p.pos_indexes[3]], (1 - prob_var_wrong)*math.log(hap2_prob/non_error_prob))
        self.assertApproxEqual(self.p.mix_matrix[sub_bc_indexes['GGGGGGGGGGGG_0_4'],self.p.pos_indexes[3]],  (1 - prob_var_wrong)*math.log(mix_prob/non_error_prob))

    def test_call_haps(self):
        out_vcf = open(OUTPUT_VCF, 'w')
        vfw = VariantFileWriter(out_vcf, template_file=open(SNP_INPUT_VCF, 'r'))
        out_bc_haps = open(OUTPUT_TSV, 'w')
        self.p.call_haps(vfw, out_bc_haps)
        out_vcf.close()
        out_bc_haps.close()
        vfr = VariantFileReader(OUTPUT_VCF)
        hap_calls = {}
        for record in vfr.record_getter():
            chrom = tk_io.get_record_chrom(record)
            pos = tk_io.get_record_pos(record) - 1
            genotype, phased = tk_io.get_record_genotype_phased(record)
            hap_calls[(chrom, pos)] = genotype
            self.assertTrue(phased)

        print hap_calls
        self.assertTrue((hap_calls[('chr1', 2)] == [1,2] and hap_calls[('chr1', 3)] == [1,0]) or hap_calls[('chr1', 2)] == [2,1] and hap_calls[('chr1', 3)] == [0,1])

    def test_output_haps(self):
        pass

    def test_get_all_best_hap_blocks(self):
        m1 = scipy.sparse.csc_matrix([[1,1,1],[0,0,0], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,0,0],[1,1,1], [0,0,0]])
        m3 = scipy.sparse.csc_matrix([[1,1,1],[1,1,1], [1,1,1]])

        p = self.p
        p.hap1_matrix = m1
        p.hap2_matrix = m2
        p.mix_matrix = m3
        p.hap_block_size = 2
        p.block_buffer_size = 0
        all_best_hap_assigns = self.p.get_all_best_hap_blocks()
        self.assertTrue(all_best_hap_assigns == [(0,0), (0,)] or all_best_hap_assigns == [(1,1), (1,)])

        m1 = scipy.sparse.csc_matrix([[1,0,0],[0,5,1], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,1,0],[2,0,1], [1,1,1]])
        m3 = scipy.sparse.csc_matrix([[3,1,0],[0,0,1], [2,1,1]])

        p.hap1_matrix = m1
        p.hap2_matrix = m2
        p.mix_matrix = m3
        p.hap_block_size = 2
        all_best_hap_assigns = self.p.get_all_best_hap_blocks()
        print all_best_hap_assigns
        assign0 = all_best_hap_assigns[0]
        self.assertTrue(assign0 == (0,1) or assign0 == (1,0))

    def test_get_best_hap_block(self):
        m1 = scipy.sparse.csc_matrix([[1,1,1],[0,0,0], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,0,0],[1,1,1], [0,0,0]])
        m3 = scipy.sparse.csc_matrix([[1,1,1],[1,1,1], [1,1,1]])

        (best_hap_assigns, sc) = self.p.get_best_hap_block(m1, m2, m3, 0, 3)
        self.assertTrue(best_hap_assigns==(0,0,0) or best_hap_assigns==(1,1,1))

        m1 = scipy.sparse.csc_matrix([[1,0,1],[0,2,0], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,1,0],[1,0,1], [0,0,0]])
        m3 = scipy.sparse.csc_matrix([[1,1,1],[1,1,1], [1,1,1]])
        (best_hap_assigns, sc) = self.p.get_best_hap_block(m1, m2, m3, 0, 3)
        self.assertTrue(best_hap_assigns==(0,1,0) or best_hap_assigns==(1,0,1))

    def test_find_bc_inds(self):
        m1 = scipy.sparse.csc_matrix([[0,1,0],[0,0,1], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,1,0],[0,0,1], [1,1,1]])
        m3 = scipy.sparse.csc_matrix([[0,1,0],[0,0,1], [1,1,1]])

        bc_inds = self.p.find_bc_inds(m1, m2, m3, 0, 3)
        self.assertEqual(type(bc_inds), type(np.array([])))
        self.assertEqual(set(bc_inds), set([0,1,2]))

        bc_inds = self.p.find_bc_inds(m1, m2, m3, 0, 2)
        self.assertEqual(set(bc_inds), set([0,2]))

    def test_create_initial_log_probs(self):
        pass

    def test_update_log_probs(self):
        m1 = scipy.sparse.csc_matrix([[1,1,0],[0,0,1], [1,1,1]])
        m2 = scipy.sparse.csc_matrix([[0,1,0],[0,2,1], [1,1,1]])
        m3 = scipy.sparse.csc_matrix([[3,1,0],[0,0,1], [2,1,1]])

        bc_inds = self.p.find_bc_inds(m1, m2, m3, 0, 3)
        log_probs = self.p.create_initial_log_probs(3)
        log_probs = self.p.update_log_probs(m1[bc_inds,], m2[bc_inds,], m3[bc_inds,], 0, log_probs)
        self.assertEqual(set(log_probs.keys()), set([(0,), (1,)]))
        (lp1, lp2, lp3) = log_probs[(0,)]
        self.assertEqual(lp1[0,0], 1)
        self.assertEqual(lp1[1,0], 0)
        self.assertEqual(lp1[2,0], 1)
        self.assertEqual(lp2[0,0], 0)
        self.assertEqual(lp2[1,0], 0)
        self.assertEqual(lp2[2,0], 1)
        self.assertEqual(lp3[0,0], 3)
        self.assertEqual(lp3[1,0], 0)
        self.assertEqual(lp3[2,0], 2)

        (lp1, lp2, lp3) = log_probs[(1,)]
        self.assertEqual(lp2[0,0], 1)
        self.assertEqual(lp2[1,0], 0)
        self.assertEqual(lp2[2,0], 1)
        self.assertEqual(lp1[0,0], 0)
        self.assertEqual(lp1[1,0], 0)
        self.assertEqual(lp1[2,0], 1)
        self.assertEqual(lp3[0,0], 3)
        self.assertEqual(lp3[1,0], 0)
        self.assertEqual(lp3[2,0], 2)

        log_probs = self.p.update_log_probs(m1[bc_inds,], m2[bc_inds,], m3[bc_inds,], 1, log_probs)
        self.assertEqual(set(log_probs.keys()), set([(0,0), (0,1), (1,0), (1,1)]))
        (lp1, lp2, lp3) = log_probs[(0,1)]
        self.assertEqual(lp1[0,0], 2)
        self.assertEqual(lp1[1,0], 2)
        self.assertEqual(lp1[2,0], 2)

    def test_calc_total_log_prob(self):
        hap1_log_probs = np.array([[-2],[-3]])
        hap2_log_probs = np.array([[-4],[-2]])
        mix_log_probs = np.array([[-1],[-1]])

        corr_prob = np.sum(np.log(self.mix_prior*np.exp(mix_log_probs) + (1.0 - self.mix_prior)*(0.5*np.exp(hap1_log_probs) + 0.5*np.exp(hap2_log_probs))))
        prob = self.p.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_log_probs)
        self.assertEqual(round(prob,5), round(corr_prob, 5))

    def test_calc_bc_block_matrices(self):
        p = self.p
        (num_bcs, num_poses) = p.hap1_matrix.shape
        (hap_block1, hap_block2, mix_blocks) = p.calc_bc_block_matrices([[0]*num_poses])
        print hap_block1
        print p.hap1_matrix.sum(axis=1)
        np.testing.assert_allclose(hap_block1[:,0].todense(), p.hap1_matrix.sum(axis=1))


        (hap_block1, hap_block2, mix_blocks) = p.calc_bc_block_matrices([[0]])
        np.testing.assert_allclose(hap_block1[:,0].todense(), p.hap1_matrix[:,0].todense())

    def stitch_blocks_together(self):
        pass

    def test_reassign_vars(self):
        pass

    def test_get_phased_vars(self):
        pass

    def test_calc_var_hap_conf(self):
        pass

    def test_calc_hap1_hap2_logprob(self):
        pass

    def test_get_phase_block_ids(self):
        pass

    def test_calc_junction_conf(self):
        pass

    def test_get_sub_bc(self):
        pass
