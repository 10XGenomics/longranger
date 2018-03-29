#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Class for detecting barcode overlaps across genomic windows.
#
import sys
import numpy as np
from numpy.random import choice
import scipy.sparse as sp
from longranger.sv.utils import BINOM, EXP_COUNT, EXP_COUNT2, EXP_COUNT3, log10_binom_pval
import tenkit.stats as tk_stats
from longranger.sv.constants import MAX_FRAG_SIZE

MAX_SAMPLE_BCS = 100000

def filter_mat(in_bc_mat, max_bcs = np.inf, min_reads = 1):
    read_counts = np.array(in_bc_mat.sum(axis = 0)).flatten()
    if min_reads > 0:
        bc_mat = sp.lil_matrix(in_bc_mat.tolil() >= float(min_reads), dtype = in_bc_mat.dtype)
    else:
        bc_mat = in_bc_mat
    bc_counts = np.array(bc_mat.sum(axis = 0)).flatten()

    sel_wins = np.where(bc_counts <= max_bcs)[0]
    bc_mat = bc_mat.tocsc()
    if len(sel_wins) > 0:
        bc_mat = bc_mat[:, sel_wins].tolil()
    else:
        bc_mat = sp.lil_matrix((bc_mat.shape[0], 0), dtype = bc_mat.dtype)
    read_counts = read_counts[sel_wins]

    return (bc_mat, sel_wins, read_counts)


class OverlapDetector(object):
    def __init__(self, bc_mat1, bc_mat2, loci1, loci2, min_reads = 1,
        max_bcs = np.inf, nx = 90):
        """
        bc_mat1/2: matrices of barcodes x windows with read counts
        loci1/2: tuples (chrom, starts, stops) where chrom is a chromosome name and starts and stops
        are lists or numpy arrays with starting and ending positions of each window.
        min_reads: minimum number of reads supporting a barcode in a window for the barcode to be
        considered as appearing in the window
        max_bcs: maximum number of barcodes in a window to attempt to make a call (too high numbers
        can cause false positives).
        """
        self.nbcs = bc_mat1.shape[0]
        assert self.nbcs == bc_mat2.shape[0], 'Number of barcodes must be equal in the two input matrices.'
        assert len(loci1) == 3, 'loci1 must be a tuple (chrom, starts, stops)'
        assert len(loci2) == 3, 'loci2 must be a tuple (chrom, starts, stops)'
        self.bc_mat1, self.sel_wins1, self.read_counts1 = filter_mat(bc_mat1, max_bcs, min_reads)
        self.bc_mat2, self.sel_wins2, self.read_counts2 = filter_mat(bc_mat2, max_bcs, min_reads)
        self.chrom1, starts1, stops1 = loci1
        self.chrom2, starts2, stops2 = loci2
        assert len(starts1) == bc_mat1.shape[1], 'Number of windows in bc_mat1 is not equal to length of starts1.'
        assert len(starts2) == bc_mat2.shape[1], 'Number of windows in bc_mat2 is not equal to length of starts2.'
        assert len(stops1) == bc_mat1.shape[1], 'Number of windows in bc_mat1 is not equal to length of stops1.'
        assert len(stops2) == bc_mat2.shape[1], 'Number of windows in bc_mat2 is not equal to length of stops2.'
        self.starts1 = np.array(starts1, dtype = np.int).flatten()[self.sel_wins1]
        self.starts2 = np.array(starts2, dtype = np.int).flatten()[self.sel_wins2]
        self.stops1 = np.array(stops1, dtype = np.int).flatten()[self.sel_wins1]
        self.stops2 = np.array(stops2, dtype = np.int).flatten()[self.sel_wins2]
        self.bc_counts1 = np.array(self.bc_mat1.sum(axis = 0)).flatten() # Num BCs in each window
        self.bc_counts2 = np.array(self.bc_mat2.sum(axis = 0)).flatten()


    @staticmethod
    def precompute_bc_overlaps(frags, method, step=1000, max_frag_size=MAX_FRAG_SIZE, genome_size=3e9):

        nbins = int(np.ceil(max_frag_size / float(step))) + 1
        exp_bc_ov = np.zeros((nbins,))
        total_bcs = len(set(frags.bc))

        if method == EXP_COUNT2:
            frags['bin'] = np.minimum(np.array(frags.obs_len / float(step), dtype=np.int), nbins - 1)
            frags.sort('bin', inplace=True, ascending=False)

            bcs_so_far = set()
            # Start from the largest group and keep track of all unique barcodes you've seen so far
            for size_bin, group in sorted(frags.groupby('bin'), reverse=True):
                bcs_so_far = bcs_so_far.union(set(group.bc))
                # Fraction of barcodes that have at least one fragment in that bin or in a
                # larger bin.
                exp_bc_ov[size_bin] = tk_stats.robust_divide(len(bcs_so_far), total_bcs)

            # Fill in any remaining small bins
            exp_bc_ov[0:size_bin] = 1.0

        elif method == EXP_COUNT3:
            # Can't sort categorical data
            frags['bc'] = frags['bc'].astype('str')
            all_bcs = set(frags.bc)
            if len(all_bcs) > MAX_SAMPLE_BCS:
                sel_bcs = choice(list(all_bcs), MAX_SAMPLE_BCS, replace=False)
                frags = frags[frags.bc.isin(sel_bcs)]
                
            frags['bin'] = np.minimum(np.array(frags['obs_len'] / float(step), dtype=np.int), nbins - 1)
            frags.sort('bc', inplace=True)
            
            exp_bc_ov = np.zeros((nbins,))

            # We want to compute the following, averaged across barcodes:
            # P(overlap at X + d | barcode present at locus X) =
            # sum_(m:L(m)>d) (P(molecule at X is m)P(overlap at X + d | m is present at X)
            # where the sum is over molecules from barcode b with length > d.
            # The first probability is proportional to the length of m L(m).
            # The second is (L(m) - d)/L(m).
            # So after simplifying we get:
            # sum_(m:L(m)>d)[L(m) - d)] / sum_m(L(m))
            # However, simply ignoring the d factor seems to work better in practice.
            for idx, (bc, bc_frags) in enumerate(frags.groupby('bc')):
                #counts[b] is the number of molecules in bin b
                counts = np.bincount(bc_frags['bin'], minlength=nbins)
                #total_len[b] is the total length of molecules in bin b
                total_len = np.bincount(bc_frags['bin'], minlength=nbins, weights=bc_frags['obs_len'])
                counts = total_len - (counts * float(step))
                exp_bc_ov += np.cumsum(counts[::-1]) / float(np.sum(total_len))

                #counts = np.cumsum(counts[::-1]) * np.arange(len(exp_bc_ov) - 1, -1, -1) * float(step)
                #exp_bc_ov += (np.cumsum(total_len[::-1]) - counts) / float(np.sum(total_len))
                
            exp_bc_ov /= idx
            exp_bc_ov = exp_bc_ov[::-1]

        else:
            raise Exception('Invalid method for computing overlaps')
                
        return exp_bc_ov

    
    def get_overlaps(self, method, max_logp=0, min_ov=1, min_dist=0, step=0,
                     max_dist=None, genome_size=0, exp_bc_ov=None,
                     frag_sizes=None, frag_counts=None, verbose=True, bc_ov_adjust=1.0):

        assert max_dist is None or min_dist <= max_dist, 'max_dist must be None or >= min_dist'

        # ov_mat[i, j] is the number of overlapping barcodes between window i
        # in bc_mat1 and window j in bc_mat2
        ov_mat = self.bc_mat1.T * self.bc_mat2

        if verbose:
            print >> sys.stderr, 'Looking for barcode overlaps...'

        if method == BINOM:
            res = self.get_overlaps_binom(ov_mat, max_logp, min_overlap=min_ov,
                                          min_dist=min_dist, max_dist=max_dist)
            (idx1, idx2, p_vals, overlaps, nbcs1, nbcs2) = res
        elif method == EXP_COUNT:
            assert(genome_size > 0)
            assert(not frag_sizes is None)
            assert(not frag_counts is None)
            res = self.get_overlaps_adaptive(ov_mat, genome_size, frag_sizes, frag_counts,
                                             min_overlap=min_ov, max_logp=max_logp,
                                             min_dist=min_dist, max_dist=max_dist,
                                             method=method)
            (idx1, idx2, overlaps, nbcs1, nbcs2) = res
            p_vals = overlaps

        elif method == EXP_COUNT2 or method == EXP_COUNT3:
            assert not exp_bc_ov is None
            res = self.get_overlaps_adaptive(ov_mat, genome_size, None, None,
                                             step=step, exp_bc_ov=exp_bc_ov,
                                             min_overlap=min_ov, max_logp=max_logp,
                                             min_dist=min_dist, max_dist=max_dist,
                                             method=method, bc_ov_adjust=bc_ov_adjust)
            (idx1, idx2, overlaps, nbcs1, nbcs2) = res
            p_vals = overlaps
        else:
            raise Exception('Unexpected name of test in get_overlaps')
            
        loci1 = (self.chrom1, self.starts1[idx1], self.stops1[idx1])
        loci2 = (self.chrom2, self.starts2[idx2], self.stops2[idx2])

        if verbose:
            print >> sys.stderr, '#Overlapping regions', len(idx1)

        return (loci1, loci2, p_vals, overlaps, nbcs1, nbcs2)


    @staticmethod
    def expected_overlap(dists, genome_size, frag_sizes, frag_counts):
        exp_counts = np.zeros((len(dists),))
        for s, c in zip(frag_sizes, frag_counts):
            exp_counts += c * np.maximum(0, s - dists)
        return exp_counts / float(genome_size)


    def get_overlaps_adaptive(self, ov_mat, genome_size, frag_sizes, frag_counts, step=0, exp_bc_ov=None,
                              max_logp=0, min_overlap=0, min_dist=0, max_dist=np.inf,
                              method=EXP_COUNT, bc_ov_adjust=1.0):
        """
        ov_mat: sparse matrix window x window. ov_mat[i,j] is the number of
        BCs overlapping between the i-th window of the first chunk and the j-th
        window of the second chunk.
        """

        nwin = ov_mat.shape[0]
        loci1 = []
        loci2 = []
        overlaps = []
        nbcs1 = []
        nbcs2 = []

        # Average number of barcodes covering a locus
        # bg = self.expected_overlap([0], genome_size, frag_sizes, frag_counts)[0]

        # Splitting prevents memory problems.
        for i in range(nwin):
            if self.bc_counts1[i] < 1:
                continue
            # overlaps between the i-th window in bc_mat1 and all windows in bc_mat2
            ov_arr = ov_mat[i, :]

            # c is the indices of windows of chunk 2 that overlap with the i-th
            # window of chunk 1. v is the number of corresponding overlapping counts
            r, c, v = sp.find(ov_arr)
            if self.chrom1 == self.chrom2:
                start_dists = np.maximum(0, self.starts2[c] - self.stops1[i])
            else:
                start_dists = np.ones((len(c),)) * MAX_FRAG_SIZE
            if method == EXP_COUNT:
                exp_ov = self.expected_overlap(start_dists, genome_size, frag_sizes, frag_counts)
                # exp_ov *= np.maximum(1.0, np.maximum(self.bc_counts1[i], self.bc_counts2[c]) / float(bg))
            elif method == EXP_COUNT2 or method == EXP_COUNT3:
                start_dists = np.minimum(MAX_FRAG_SIZE, start_dists)
                exp_ov = exp_bc_ov[np.array(start_dists / step, dtype=np.int)]
                # Adjust the expected count. bc_ov_adjust can be used to relax our search (if < 1.0).
                exp_ov *= (bc_ov_adjust * np.minimum(self.bc_counts1[i], self.bc_counts2[c]))
            else:
                raise Exception('Unexpected name of test in get_overlaps_adaptive')
            
            good = v >= np.maximum(min_overlap, exp_ov)
            binom_p_arr = log10_binom_pval(v, self.bc_counts1[i], self.bc_counts2[c], self.nbcs)
            good = np.logical_and(good, binom_p_arr <= max_logp)

            if self.chrom1 == self.chrom2:
                good = np.logical_and(good,
                                      np.logical_and(start_dists >= min_dist,
                                                     start_dists <= max_dist))

            loci1.extend([i for count in range(sum(good))])
            loci2.extend(c[good].tolist())
            overlaps.extend(v[good].tolist())
            nbcs1.extend([self.bc_counts1[i] for count in range(sum(good))])
            nbcs2.extend(self.bc_counts2[c[good]].tolist())

        return (loci1, loci2, overlaps, nbcs1, nbcs2)


    def get_overlaps_binom(self, ov_mat, max_logp, min_overlap=0,
                           min_dist=0, max_dist=None):
        """
        ov_mat: sparse matrix window x window. ov_mat[i,j] is the number of
        BCs overlapping between the i-th window of the first chunk and the j-th
        window of the second chunk.
        bc_counts: flat (1D) arrays with the total number of BCs in each of the
        windows of chunks 1 and 2.
        nbcs: total number of barcodes
        max_p: maximum p-value for reporting overlaps.

        Return value:
        A tuple (locus_1_poses, locus_2_poses, p_vals) with the indices of
        overlapping windows and the corresponding p-values.
        """

        nwin = ov_mat.shape[0]
        loci1 = []
        loci2 = []
        p_vals = []
        overlaps = []
        nbcs1 = []
        nbcs2 = []

        # Splitting prevents memory problems.
        for i in range(nwin):
            if self.bc_counts1[i] < 1:
                continue
            ov_arr = ov_mat[i, :]
            # c is the indices of windows of chunk 2 that overlap with the i-th
            # window of chunk 1. v is the number of corresponding overlapping counts
            r, c, v = sp.find(ov_arr)
            good_dist = v >= min_overlap
            if not min_dist is None and not max_dist is None and max_dist >= 0 and self.chrom1 == self.chrom2:
                start_dist = self.starts2[c] - self.stops1[i]
                good_dist = np.logical_and(good_dist,
                                           np.logical_and(start_dist >= min_dist, start_dist <= max_dist))
            c = c[good_dist]
            v = v[good_dist]
            binom_p_arr = log10_binom_pval(v, self.bc_counts1[i], self.bc_counts2[c], self.nbcs)
            good = binom_p_arr <= max_logp
            p_vals.extend(binom_p_arr[good].tolist())
            loci1.extend([i for count in range(sum(good))])
            loci2.extend(c[good].tolist())
            overlaps.extend(v[good].tolist())
            nbcs1.extend([self.bc_counts1[i] for count in range(sum(good))])
            nbcs2.extend(self.bc_counts2[c[good]].tolist())

        self.bc_mat1 = self.bc_mat1.tolil()
        self.bc_mat2 = self.bc_mat2.tolil()

        return (loci1, loci2, p_vals, overlaps, nbcs1, nbcs2)
