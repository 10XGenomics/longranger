#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for computing fragment overlaps and calling large scale structural variants.
#
import sys
import numpy as np
from bisect import bisect
from itertools import product
import longranger.sv.stats as tk_sv_stats
from longranger.sv.constants import MIN_LOG_PROB, MAX_MOL_PROB, MAX_PHASE_PROB, MAX_DIST_FROM_CAND
import longranger.sv.sv_call as sv_call

class ReadInfo(object):
    def __init__(self, chroms, poses, group_ids=None):
        assert len(chroms) == len(poses)
        assert group_ids is None or len(poses) == len(group_ids)
        assert len(set(chroms)) <= 2
        self.is_grouped = not group_ids is None
        if not self.is_grouped:
            group_ids = np.ones((len(poses),), dtype=np.int)
        sorted_poses = sorted(zip(group_ids, chroms, poses))
        self.group_ids = np.array([gi for gi, chrom, pos in sorted_poses], dtype=np.int)
        self.poses = np.array([pos for gi, chrom, pos in sorted_poses], dtype=np.int)
        self.chroms = np.array([chrom for gi, chrom, pos in sorted_poses])


    def __len__(self):
        return len(self.poses)


    def split_in_two(self):
        if self.is_grouped:
            idx1 = self.group_ids == 1
        else:
            assert len(self.poses) >= 4

            uniq_chroms = list(set(list(self.chroms)))
            if len(uniq_chroms) == 2:
                idx1 = self.chroms == uniq_chroms[0]
            else:
                # Find the biggest gap
                max_gap_idx = np.argmax(np.diff(self.poses))
                idx1 = np.arange(0, len(self.poses)) <= max_gap_idx

        idx2 = np.logical_not(idx1)
        reads1 = ReadInfo(self.chroms[idx1], self.poses[idx1], self.group_ids[idx1])
        reads2 = ReadInfo(self.chroms[idx2], self.poses[idx2], self.group_ids[idx2])

        return (reads1, reads2)


class ReadModel(object):
    def __init__(self, alpha, frag_sizes, frag_counts, p_ov_mol=1e-9, step=1000):
        self.alpha = alpha
        self.p_ov_mol = p_ov_mol
        self.frag_sizes = frag_sizes
        self.frag_counts = frag_counts
        self.step = step
        self.size_log_probs = self.compute_log_len(alpha, frag_sizes, frag_counts, step)
        self.prior_hap_probs1 = None
        self.prior_hap_probs2 = None
        self.prior_mol_probs = None


    @staticmethod
    def compute_log_len(alpha, frag_sizes, frag_counts, step=1000):
        size_bins = np.arange(frag_sizes[0], max(frag_sizes[0] + step, frag_sizes[-1]), step)
        log_probs = np.zeros((len(size_bins),))

        for i, size in enumerate(size_bins):
            prob_true_len, true_len = tk_sv_stats.frag_size_logpmf(size, frag_sizes, frag_counts)
            #norm_factor = tk_sv_stats.safe_logaddexp(prob_true_len)
            log_probs[i] = tk_sv_stats.safe_logaddexp(prob_true_len - alpha * true_len) #- norm_factor
        return log_probs


    def mol_prob(self, nreads, obs_len):
        """Log-probability of observing n reads spanning an observed length of obs_len.

        Inputs:
        - nreads: number of observed reads
        - obs_len: distance between the first and last read
        """

        bin_idx = int(obs_len / float(self.step))
        size_prob = self.size_log_probs[bin_idx] if bin_idx < len(self.size_log_probs) else MIN_LOG_PROB
        return max(MIN_LOG_PROB, nreads * self.alpha + size_prob)


    def mol_prob_arr(self, nreads, obs_len):
        output = np.ones((len(obs_len), )) * MIN_LOG_PROB
        bin_idx = np.array(obs_len / float(self.step), dtype=np.int)
        size_probs = self.size_log_probs[bin_idx[bin_idx < len(self.size_log_probs)]]
        output[bin_idx < len(self.size_log_probs)] = np.maximum(MIN_LOG_PROB, nreads * self.alpha + size_probs)
        return output


    def get_ref_prob(self, reads1, reads2):
        """Log-probability under the reference model"""

        if reads1 is None or len(reads1) == 0:
            p_mol1 = float('NaN')
            p_mol2 = self.mol_prob(len(reads2), np.max(reads2.poses) - np.min(reads2.poses))
            p_joint = p_mol2
        elif reads2 is None or len(reads2) == 0:
            p_mol1 = self.mol_prob(len(reads1), np.max(reads1.poses) - np.min(reads1.poses))
            p_mol2 = float('NaN')
            p_joint = p_mol1
        else:
            p_mol1 = self.mol_prob(len(reads1), np.max(reads1.poses) - np.min(reads1.poses))
            p_mol2 = self.mol_prob(len(reads2), np.max(reads2.poses) - np.min(reads2.poses))

            chrom1 = reads1.chroms[0]
            chrom2 = reads2.chroms[0]
            if chrom1 == chrom2:
                p_joint = self.mol_prob(len(reads1) + len(reads2),
                                        np.max(reads2.poses) - np.min(reads1.poses))
            else:
                p_joint = MIN_LOG_PROB

        return (p_joint, p_mol1, p_mol2)


    def get_dup_prob_single_mol(self, reads, dup_start, dup_end, ref_prob,
                                left_idx, right_idx):

        start, end = np.min(reads.poses), np.max(reads.poses)
        dup_len = dup_end - dup_start

        if len(reads) <= 2:
            return ref_prob

        if end <= dup_start or start >= dup_end:
            return ref_prob

        if start <= dup_start and end >= dup_end:
            # molecule length is increased by the length of the dup
            return self.mol_prob(len(reads), end - start + dup_len)

        if start >= dup_start and end <= dup_end:
            if start - dup_start > MAX_DIST_FROM_CAND or dup_end - end > MAX_DIST_FROM_CAND:
                return ref_prob

            # all reads contained inside the duplication
            gap_pos = np.argmax(np.diff(reads.poses))
            new_len = reads.poses[gap_pos] - dup_start + dup_end - reads.poses[gap_pos + 1]
            return max(ref_prob, self.mol_prob(len(reads), new_len))

        return ref_prob


    def get_del_prob_single_mol(self, reads, del_start, del_end, ref_prob,
                                left_idx, right_idx):

        start, end = np.min(reads.poses), np.max(reads.poses)
        del_len = del_end - del_start

        if end <= del_start or start >= del_end:
            return ref_prob

        if left_idx != right_idx:
            return MIN_LOG_PROB

        return self.mol_prob(len(reads), end - start - del_len)


    def get_inv_prob_single_mol(self, reads, inv_start, inv_end, ref_prob,
                                left_idx, right_idx):
        start, end = np.min(reads.poses), np.max(reads.poses)

        if end <= inv_start or start >= inv_end:
            return ref_prob

        if start <= inv_start and end >= inv_end:
            return ref_prob

        if start >= inv_start and end <= inv_end:
            return ref_prob

        if inv_end > end:
            assert inv_start > start
            # x1 ... xi X    xj ... xn Y
            # whre X = inv_start, Y = inv_end, xj = read_poses[left_idx]
            # The length of the inverted read sequence is
            # xn - x1 - (xn - X) + (Y - xj) = X + Y - x1 - xj
            #if read_poses[left_idx] - inv_start < MAX_DIST_FROM_CAND:
            #    return ref_prob
            return self.mol_prob(len(reads), inv_start + inv_end - start - reads.poses[left_idx])

        # X x1 ... xi  Y xj ... xn
        # The length of the inverted read sequence is
        # xn - x1 - (Y - x1) + xi - X = xn + xi - X - Y
        assert inv_end < end
        #if inv_end - read_poses[right_idx - 1] < MAX_DIST_FROM_CAND:
        #    return ref_prob
        return self.mol_prob(len(reads), end + reads.poses[right_idx - 1] - inv_start - inv_end)


    def get_inv_dup_out_prob_single_mol(self, reads, inv_start, inv_end, ref_prob,
                                        left_idx, right_idx):
        start, end = np.min(reads.poses), np.max(reads.poses)

        if end <= inv_start or start >= inv_end:
            return ref_prob

        # This is different from an inversion and more like what
        # happens in a duplication.
        if start <= inv_start and end >= inv_end:
            return self.mol_prob(len(reads), end - start + inv_end - inv_start)

        if start >= inv_start and end <= inv_end:
            return ref_prob

        if start <= inv_end <= end:
            return ref_prob

        # Similar to inversion
        assert inv_start > start
        return self.mol_prob(len(reads), inv_start + inv_end - start - reads.poses[left_idx])


    def get_inv_dup_in_prob_single_mol(self, reads, inv_start, inv_end, ref_prob,
                                        left_idx, right_idx):
        start, end = np.min(reads.poses), np.max(reads.poses)

        if end <= inv_start or start >= inv_end:
            return ref_prob

        # This is different from an inversion and more like what
        # happens in a duplication.
        if start <= inv_start and end >= inv_end:
            return self.mol_prob(len(reads), end - start + inv_end - inv_start)

        if start >= inv_start and end <= inv_end:
            return ref_prob

        if start <= inv_start <= end:
            return ref_prob

        # Similar to inversion
        assert inv_end < end
        return self.mol_prob(len(reads), end + reads.poses[right_idx - 1] - inv_start - inv_end)


    def get_trans_prob(self, reads1, reads2, loci,
                       left_orient, right_orient,
                       ref_prob1, ref_prob2, ref_prob_joint):

        nloci = len(loci)
        left_breaks = [a for (a, b) in loci]
        right_breaks = [b for (a, b) in loci]

        left_probs = np.ones((nloci, )) * ref_prob1
        right_probs = np.ones((nloci, )) * ref_prob2
        joint_probs = np.ones((nloci, )) * ref_prob_joint

        if reads1 is None or len(reads1) == 0:
            right_start, right_end = np.min(reads2.poses), np.max(reads2.poses)
            bad_loci = np.logical_and(right_start < right_breaks,
                                      right_breaks < right_end)
            right_probs[bad_loci] = MIN_LOG_PROB
            joint_probs[bad_loci] = MIN_LOG_PROB

        elif reads2 is None or len(reads2) == 0:
            left_start, left_end = np.min(reads1.poses), np.max(reads1.poses)
            bad_loci = np.logical_and(left_start < left_breaks,
                                      left_breaks < left_end)
            left_probs[bad_loci] = MIN_LOG_PROB
            joint_probs[bad_loci] = MIN_LOG_PROB
        else:
            left_start, left_end = np.min(reads1.poses), np.max(reads1.poses)
            right_start, right_end = np.min(reads2.poses), np.max(reads2.poses)

            left_crosses = np.logical_and(right_start < right_breaks,
                                          right_breaks < right_end)
            right_crosses = np.logical_and(left_start < left_breaks,
                                           left_breaks < left_end)

            left_probs[left_crosses] = MIN_LOG_PROB
            right_probs[right_crosses] = MIN_LOG_PROB
            joint_probs[np.logical_or(left_crosses, right_crosses)] = MIN_LOG_PROB

            if left_orient == sv_call.UPSTREAM:
                bad_left = left_end > left_breaks
            else:
                bad_left = left_start < left_breaks
            joint_probs[bad_left] = MIN_LOG_PROB

            if right_orient == sv_call.UPSTREAM:
                bad_right = right_end > right_breaks
            else:
                bad_right = right_start < right_breaks
            joint_probs[bad_right] = MIN_LOG_PROB

            good_joint = np.logical_and(bad_left == False,
                                        bad_right == False)
            left_mol_lens = np.maximum(left_breaks - left_start, left_end - left_breaks)
            right_mol_lens = np.maximum(right_breaks - right_start, right_end - right_breaks)
            total_lens = left_mol_lens[good_joint] + right_mol_lens[good_joint]
            joint_probs[good_joint] = self.mol_prob_arr(len(reads1) + len(reads2), total_lens)

        return (left_probs, right_probs, joint_probs)


    def get_sv_prob(self, sv_type, loci, reads, reads1, reads2,
                    ref_prob1, ref_prob2, ref_prob_joint,
                    read_idx=None):

        read_poses_joint = reads.poses

        if sv_type in sv_call.DISTAL_SV_TYPES:
            left_orient = sv_type.orientation[0]
            right_orient = sv_type.orientation[1]
            all_probs = self.get_trans_prob(reads1, reads2, loci, left_orient, right_orient,
                                            ref_prob1, ref_prob2, ref_prob_joint)
            sv_probs1, sv_probs2, sv_probs_joint = all_probs
        else:
            nloci = len(loci)
            nreads1 = len(reads1)

            sv_probs1 = np.zeros((nloci, ))
            sv_probs2 = np.zeros((nloci, ))
            sv_probs_joint = np.zeros((nloci, ))

            sv_fun_dict = {sv_call.DEL:self.get_del_prob_single_mol,
                           sv_call.INV:self.get_inv_prob_single_mol,
                           sv_call.INV_DUP_OUT:self.get_inv_dup_out_prob_single_mol,
                           sv_call.INV_DUP_IN:self.get_inv_dup_in_prob_single_mol,
                           sv_call.DUP:self.get_dup_prob_single_mol}

            stored_vals = {}
            if read_idx is None:
                read_idx = {}

            for i, (sv_start, sv_end) in enumerate(loci):
                if sv_start in read_idx:
                    left_idx = read_idx[sv_start]
                else:
                    left_idx = bisect(read_poses_joint, sv_start)
                    read_idx[sv_start] = left_idx

                adj_left_idx = left_idx - nreads1

                if sv_end in read_idx:
                    right_idx = read_idx[sv_end]
                else:
                    right_idx = bisect(read_poses_joint, sv_end)
                    read_idx[sv_end] = right_idx

                adj_right_idx = right_idx - nreads1

                if (left_idx, right_idx) in stored_vals:
                    p1, p2, p0 = stored_vals[(left_idx, right_idx)]
                    sv_probs1[i] = p1
                    sv_probs2[i] = p2
                    sv_probs_joint[i] = p0
                else:
                    sv_fun = sv_fun_dict[sv_type]
                    sv_probs1[i] = sv_fun(reads1, sv_start, sv_end, ref_prob1,
                                          left_idx, min(right_idx, nreads1))
                    sv_probs2[i] = sv_fun(reads2, sv_start, sv_end, ref_prob2,
                                          adj_left_idx, adj_right_idx)
                    sv_probs_joint[i] = sv_fun(reads, sv_start, sv_end, ref_prob_joint,
                                               left_idx, right_idx)
                    stored_vals[(left_idx, right_idx)] = (sv_probs1[i], sv_probs2[i], sv_probs_joint[i])
        return (sv_probs_joint, sv_probs1, sv_probs2)


    def em_it_away(self, loci, read_groups, phase_set1, phase_set2,
                   bc_phase_set_dict1=None, bc_phase_set_dict2=None,
                   em_iters=5, proximal=True):

        nbcs = len(read_groups)
        nloci = len(loci)

        sv_types = sv_call.PROXIMAL_SV_TYPES if proximal else sv_call.DISTAL_SV_TYPES
        nsv_types = len(sv_types)

        # Get prior haplotype probabilities
        probs = self.init_hap_probs(read_groups, bc_phase_set_dict1, bc_phase_set_dict2)
        self.prior_hap_probs1, self.prior_hap_probs2 = probs

        # Get prior probabilities of reads coming from a single molecule
        self.prior_mol_probs = self.init_mol_probs(read_groups, self.p_ov_mol,
                                                   self.frag_sizes, self.frag_counts)

        # Log-likelihoods under the reference model, assuming all data come from a
        # single molecule.
        no_sv_probs = np.zeros((nbcs,))
        # Log-likelihoods for the two (most likely) splits of the reads into two
        # molecules
        no_sv_probs_mol1 = np.zeros((nbcs,))
        no_sv_probs_mol2 = np.zeros((nbcs,))

        sv_probs = np.zeros((nbcs, nloci, nsv_types))
        sv_probs_mol1 = np.zeros((nbcs, nloci, nsv_types))
        sv_probs_mol2 = np.zeros((nbcs, nloci, nsv_types))

        for i, (_, reads) in enumerate(sorted(read_groups.iteritems())):
            reads1, reads2 = reads.split_in_two()
            probs = self.get_ref_prob(reads1, reads2)
            no_sv_probs[i] = probs[0]
            no_sv_probs_mol1[i] = probs[1]
            no_sv_probs_mol2[i] = probs[2]

            read_idx = {}
            for sv_idx, svt in enumerate(sv_types):
                probs = self.get_sv_prob(svt, loci, reads, reads1, reads2,
                                         no_sv_probs_mol1[i], no_sv_probs_mol2[i],
                                         no_sv_probs[i], read_idx)
                sv_probs[i, :, sv_idx] = probs[0]
                sv_probs_mol1[i, :, sv_idx] = probs[1]
                sv_probs_mol2[i, :, sv_idx] = probs[2]

        # het_sv_probs[b, n, svt] is the probability of observing the data
        # for barcode b under the assumption that there is an sv of type svt at locus n
        # (or the n-th pair of loci), on haplotypes h1 and h2
        het_sv_probs = np.zeros((nbcs, nloci, nsv_types, 2, 2))
        if phase_set1 == phase_set2:
            het_sv_probs[:, :, :, 0, 1] = -np.inf
            het_sv_probs[:, :, :, 1, 0] = -np.inf

        posterior_hap_prob_dict = {}

        # Need to optimize separately by SV type otherwise we can get stuck in the
        # wrong local optimum.
        for sv_idx, svt in enumerate(sv_types):

            print >> sys.stderr, 'Evaluating sv type', svt

            for hap1, hap2 in product([0, 1], [0, 1]):
                if phase_set1 == phase_set2 and hap1 != hap2:
                    continue

                hap_probs1 = np.array(self.prior_hap_probs1)
                hap_probs2 = np.array(self.prior_hap_probs2)
                mol_probs = np.array(self.prior_mol_probs)

                for em_iter in range(em_iters):
                    print >> sys.stderr, 'EM iteration', em_iter

                    no_left_reads = np.zeros((len(read_groups),), dtype=np.bool)
                    no_right_reads = np.zeros((len(read_groups),), dtype=np.bool)

                    for i, (_, reads) in enumerate(sorted(read_groups.iteritems())):
                        reads1, reads2 = reads.split_in_two()

                        if reads1 is None or len(reads1) == 0:
                            no_left_reads[i] = True
                            hp_mol1 = hap_probs2[i, hap2]
                            probs_mol1 = np.logaddexp(np.log(hp_mol1) + sv_probs[i, :, sv_idx],
                                                      np.log(1 - hp_mol1) + no_sv_probs[i])
                            het_sv_probs[i, :, sv_idx, hap1, hap2] = probs_mol1
                        elif reads2 is None or len(reads2) == 0:
                            no_right_reads[i] = True
                            hp_mol1 = hap_probs1[i, hap1]
                            probs_mol1 = np.logaddexp(np.log(hp_mol1) + sv_probs[i, :, sv_idx],
                                                      np.log(1 - hp_mol1) + no_sv_probs[i])
                            het_sv_probs[i, :, sv_idx, hap1, hap2] = probs_mol1
                        else:
                            hp_mol1 = min(hap_probs1[i, hap1], hap_probs2[i, hap2])
                            if phase_set1 == phase_set2:
                                hp_mol2 = hp_mol1
                            else:
                                hp_mol2 = hap_probs1[i, hap1] * hap_probs2[i, hap2]
                            probs_mol1 = np.logaddexp(np.log(hp_mol1) + np.log(mol_probs[i]) + sv_probs[i, :, sv_idx],
                                                      np.log(1 - hp_mol1) + np.log(mol_probs[i]) + no_sv_probs[i])
                            probs_mol2 = np.logaddexp(np.log(hp_mol2) + np.log(1 - mol_probs[i]) + sv_probs_mol1[i, :, sv_idx] + sv_probs_mol2[i, :, sv_idx],
                                                      np.log(1 - hp_mol2) + np.log(1 - mol_probs[i]) + no_sv_probs_mol1[i] + no_sv_probs_mol2[i])
                            het_sv_probs[i, :, sv_idx, hap1, hap2] = np.logaddexp(probs_mol1, probs_mol2)

                    res = self.posterior_hap_probs(no_sv_probs, no_sv_probs_mol1, no_sv_probs_mol2,
                                                   sv_probs, sv_probs_mol1, sv_probs_mol2,
                                                   het_sv_probs, sv_idx, hap1, hap2,
                                                   no_left_reads, no_right_reads)
                    new_hap_probs1, new_hap_probs2, new_mol_probs = res
                    max_change = max(np.max(np.abs(new_hap_probs1 - hap_probs1)),
                                     np.max(np.abs(new_hap_probs2 - hap_probs2)))
                    print >> sys.stderr, 'Max change in phase probabilities', max_change
                    if max_change < 1e-4:
                        break
                    hap_probs1 = new_hap_probs1
                    hap_probs2 = new_hap_probs2
                    mol_probs = new_mol_probs

                posterior_hap_prob_dict[(svt, hap1, hap2)] = (hap_probs1, hap_probs2, mol_probs)

        res = self.get_best_model(no_sv_probs, sv_probs, het_sv_probs, sv_types=sv_types,
                                  no_sv_probs_mol1=no_sv_probs_mol1, no_sv_probs_mol2=no_sv_probs_mol2,
                                  sv_probs_mol1=sv_probs_mol1, sv_probs_mol2=sv_probs_mol2,
                                  posterior_hap_prob_dict=posterior_hap_prob_dict,
                                  no_left_reads=no_left_reads, no_right_reads=no_right_reads)
        het_sv_max_type_idx, (het_sv_max, het_sv_max_idx, het_sv_max_hap1, het_sv_max_hap2) = res

        best_combo = (sv_types[het_sv_max_type_idx], het_sv_max_hap1, het_sv_max_hap2)
        hap_probs1 = posterior_hap_prob_dict[best_combo][0]
        hap_probs2 = posterior_hap_prob_dict[best_combo][1]
        mol_probs = posterior_hap_prob_dict[best_combo][2]

        for i, (_, reads) in enumerate(sorted(read_groups.iteritems())):
            if not no_left_reads[i] and not no_right_reads[i]:
                no_sv_probs[i] = np.logaddexp(np.log(mol_probs[i]) + no_sv_probs[i],
                                              np.log(1 - mol_probs[i]) + no_sv_probs_mol1[i] + no_sv_probs_mol2[i])
                for sv_idx in range(nsv_types):
                    sv_probs[i, :, sv_idx] = np.logaddexp(np.log(mol_probs[i]) + sv_probs[i, :, sv_idx],
                                                          np.log(1 - mol_probs[i]) + sv_probs_mol1[i, :, sv_idx] + sv_probs_mol2[i, :, sv_idx])

        # Find optimal model
        # Sum over barcodes
        total_no_sv_probs = np.sum(no_sv_probs, axis=0)
        no_sv_max = np.max(total_no_sv_probs)

        total_sv_probs = np.sum(sv_probs, axis=0)
        sv_max_idx = np.argmax(total_sv_probs)
        sv_max_idx, sv_max_type_idx = np.unravel_index(sv_max_idx, total_sv_probs.shape)
        sv_max = total_sv_probs[sv_max_idx, sv_max_type_idx]

        if sv_max >= het_sv_max:
            # Homozygous event
            support = sv_probs[:, sv_max_idx, sv_max_type_idx] - no_sv_probs
            zygosity = sv_call.Zygosity.hom
            max_hap = (sv_call.Haplotype.unknown, sv_call.Haplotype.unknown)
            max_locus = loci[sv_max_idx]
            max_svt = sv_types[sv_max_type_idx]
        else:
            support = sv_probs[:, het_sv_max_idx, het_sv_max_type_idx] - no_sv_probs
            zygosity = sv_call.Zygosity.het
            max_hap = (sv_call.Haplotype.hap1 if het_sv_max_hap1 else sv_call.Haplotype.hap0,
                       sv_call.Haplotype.hap1 if het_sv_max_hap2 else sv_call.Haplotype.hap0)
            max_locus = loci[het_sv_max_idx]
            max_svt = sv_types[het_sv_max_type_idx]

        return ((no_sv_max, sv_max, het_sv_max), max_locus, max_svt,
                zygosity, max_hap, support,
                (hap_probs1, hap_probs2, mol_probs))

    @staticmethod
    def get_best_model(no_sv_probs, sv_probs, het_sv_probs, sv_types = None,
                       no_sv_probs_mol1 = None, no_sv_probs_mol2 = None,
                       sv_probs_mol1 = None, sv_probs_mol2 = None,
                       posterior_hap_prob_dict = None, no_left_reads = None,
                       no_right_reads = None):
        # sum over barcodes
        total_het_sv_probs = np.sum(het_sv_probs, axis=0)
        nsv_types = len(sv_types)
        nloci = total_het_sv_probs.shape[0]

        # total barcode support for each sv type
        support_bcs = np.zeros((nsv_types, ), dtype=np.int)
        # best combo of locus and haplotypes for each sv type
        best_combos = []

        for sv_idx in range(nsv_types):
            # get maximum likelihood position and haplotypes for this sv type
            het_sv_max_idx = np.argmax(total_het_sv_probs[:, sv_idx, :, :])
            res = np.unravel_index(het_sv_max_idx, (nloci, 2, 2))
            het_sv_max_idx, het_sv_max_hap1, het_sv_max_hap2 = res
            het_sv_max = np.max(total_het_sv_probs[:, sv_idx, :, :])
            best_combos.append((het_sv_max, het_sv_max_idx, het_sv_max_hap1, het_sv_max_hap2))

            best_combo = (sv_types[sv_idx], het_sv_max_hap1, het_sv_max_hap2)
            hap_probs1 = posterior_hap_prob_dict[best_combo][0]
            hap_probs2 = posterior_hap_prob_dict[best_combo][1]
            mol_probs = posterior_hap_prob_dict[best_combo][2]

            # make a copy of these arrays
            no_sv_probs_cp = np.array(no_sv_probs)
            sv_probs_cp = np.array(sv_probs[:, het_sv_max_idx, sv_idx])

            good_loci = np.logical_and(no_left_reads == False, no_right_reads == False)
            no_sv_probs_cp[good_loci] = np.logaddexp(np.log(mol_probs[good_loci]) + no_sv_probs[good_loci],
                                                     np.log(1 - mol_probs[good_loci]) + no_sv_probs_mol1[good_loci] + no_sv_probs_mol2[good_loci])
            sv_probs_cp[good_loci] = np.logaddexp(np.log(mol_probs[good_loci]) + sv_probs_cp[good_loci],
                                                  np.log(1 - mol_probs[good_loci]) + sv_probs_mol1[good_loci, het_sv_max_idx, sv_idx] + \
                                                  sv_probs_mol2[good_loci, het_sv_max_idx, sv_idx])

            hap_probs_sv = np.minimum(hap_probs1[:, het_sv_max_hap1],
                                      hap_probs2[:, het_sv_max_hap2])

            # get the total support at this position
            support = sv_probs_cp - no_sv_probs_cp
            support_bcs[sv_idx] = np.sum(hap_probs_sv[support > 0])

        # pick the model that has the maximum support at its maximum likelihood position
        sel_sv_idx = np.argmax(support_bcs)
        return (sel_sv_idx, best_combos[sel_sv_idx])

    @staticmethod
    def init_hap_probs(read_groups, bc_phase_set_dict1, bc_phase_set_dict2):
        nbcs = len(read_groups)
        # hap_probs1[bc, 0] and hap_probs[bc, 1] are going to be the
        # probabilities of barcode bc being on haplotypes 0 or 1 respectively
        # around the first breakpoint
        hap_probs1 = np.zeros((nbcs, 2))
        hap_probs2 = np.zeros((nbcs, 2))

        #bcs_on_hap0 = np.sum([v[2][0] > v[2][1] for v in bc_phase_set_dict1.values()])
        default_hap0_prob = 0.5 #bcs_on_hap0 / float(len(bc_phase_set_dict))

        for i, (bc, _) in enumerate(sorted(read_groups.iteritems())):
            if bc in bc_phase_set_dict1:
                hap_probs1[i, 0] = max(1 - MAX_PHASE_PROB,
                                       min(MAX_PHASE_PROB, bc_phase_set_dict1[bc][2][0]))
                hap_probs1[i, 1] = max(1 - MAX_PHASE_PROB, min(MAX_PHASE_PROB, bc_phase_set_dict1[bc][2][1]))
            else:
                hap_probs1[i, 0] = default_hap0_prob
                hap_probs1[i, 1] = 1 - default_hap0_prob

            if bc in bc_phase_set_dict2:
                hap_probs2[i, 0] = max(1 - MAX_PHASE_PROB,
                                       min(MAX_PHASE_PROB, bc_phase_set_dict2[bc][2][0]))
                hap_probs2[i, 1] = max(1 - MAX_PHASE_PROB,
                                       min(MAX_PHASE_PROB, bc_phase_set_dict2[bc][2][1]))
            else:
                hap_probs2[i, 0] = default_hap0_prob
                hap_probs2[i, 1] = 1 - default_hap0_prob

        return (hap_probs1, hap_probs2)


    @staticmethod
    def init_mol_probs(read_groups, p_ov_mol, frag_sizes, frag_counts):
        nbcs = len(read_groups)
        mol_probs = np.zeros((nbcs, ))

        for i, (_, reads) in enumerate(sorted(read_groups.iteritems())):
            # probability that data come from single molecule =
            # prob mol length > observed length *
            # prob there was no molecule overlap (i.e. no two
            # molecules from the same locus with the same barcode)
            if len(set(reads.chroms)) > 1:
                prob_len = MIN_LOG_PROB
            else:
                obs_len = np.max(reads.poses) - np.min(reads.poses)
                prob_true_len, _ = tk_sv_stats.frag_size_logpmf(obs_len, frag_sizes, frag_counts)
                prob_len = tk_sv_stats.safe_logaddexp(prob_true_len)
            mol_probs[i] = prob_len + np.log(1 - p_ov_mol)
        return np.maximum(1 - MAX_MOL_PROB, np.minimum(MAX_MOL_PROB, np.exp(mol_probs)))


    def posterior_hap_probs(self, no_sv_probs, no_sv_probs_mol1, no_sv_probs_mol2,
                            sv_probs, sv_probs_mol1, sv_probs_mol2,
                            het_sv_probs, svt_idx, best_hap_idx1, best_hap_idx2,
                            no_left_reads=False, no_right_reads=False):
        # Find optimal breakpoints and phasing
        # Sum across barcodes
        bc_probs = np.sum(het_sv_probs[:, :, svt_idx, best_hap_idx1, best_hap_idx2], axis=0)
        best_locus_idx = np.argmax(bc_probs)

        new_hap_probs1 = np.array(self.prior_hap_probs1)
        new_hap_probs2 = np.array(self.prior_hap_probs2)
        new_mol_probs = np.array(self.prior_mol_probs)

        for i in range(het_sv_probs.shape[0]):
            # Mixed hap, don't touch.
            if 1.0 - self.prior_hap_probs1[i, 0] - self.prior_hap_probs1[i, 1] > 0.1 or \
                1.0 - self.prior_hap_probs2[i, 0] - self.prior_hap_probs2[i, 1] > 0.1:
                continue
            if no_left_reads[i] or no_right_reads[i]:
                hp_mol1 = self.prior_hap_probs2[i, best_hap_idx2] if no_left_reads[i] else self.prior_hap_probs1[i, best_hap_idx1]
                p_on_hap_mol1 = sv_probs[i, best_locus_idx, svt_idx] + np.log(hp_mol1)
                p_off_hap_mol1 = no_sv_probs[i] + np.log(1 - hp_mol1)
                norm_factor = tk_sv_stats.safe_logaddexp([p_on_hap_mol1, p_off_hap_mol1])
                sel_hap = best_hap_idx2 if no_left_reads[i] else best_hap_idx1
                sel_probs = new_hap_probs2 if no_left_reads[i] else new_hap_probs1
                sel_probs[i, sel_hap] = np.exp(p_on_hap_mol1 - norm_factor)
                sel_probs[i, 1 - sel_hap] = 1 - sel_probs[i, sel_hap]
            else:
                hp_mol1 = min(self.prior_hap_probs1[i, best_hap_idx1],
                              self.prior_hap_probs2[i, best_hap_idx2])
                hp_mol2 = self.prior_hap_probs1[i, best_hap_idx1] * self.prior_hap_probs2[i, best_hap_idx2]
                p_on_hap_mol1 = sv_probs[i, best_locus_idx, svt_idx] + np.log(hp_mol1) + np.log(self.prior_mol_probs[i])
                p_on_hap_mol2 = sv_probs_mol1[i, best_locus_idx, svt_idx] + \
                                sv_probs_mol2[i, best_locus_idx, svt_idx] + \
                                np.log(hp_mol2) + np.log(1 - self.prior_mol_probs[i])
                p_off_hap_mol1 = no_sv_probs[i] + np.log(1 - hp_mol1) + np.log(self.prior_mol_probs[i])
                p_off_hap_mol2 = no_sv_probs_mol1[i] + no_sv_probs_mol2[i] + \
                                 np.log(1 - hp_mol2) + np.log(1 - self.prior_mol_probs[i])

                norm_factor = tk_sv_stats.safe_logaddexp([p_on_hap_mol1, p_on_hap_mol2,
                                                          p_off_hap_mol1, p_off_hap_mol2])
                new_hap_probs1[i, best_hap_idx1] = np.exp(np.logaddexp(p_on_hap_mol1, p_on_hap_mol2) - norm_factor)
                new_hap_probs1[i, 1 - best_hap_idx1] = 1 - new_hap_probs1[i, best_hap_idx1]
                new_hap_probs2[i, best_hap_idx2] = np.exp(np.logaddexp(p_on_hap_mol1, p_on_hap_mol2) - norm_factor)
                new_hap_probs2[i, 1 - best_hap_idx2] = 1 - new_hap_probs2[i, best_hap_idx2]

            new_mol_probs[i] = np.exp(np.logaddexp(p_on_hap_mol1, p_off_hap_mol1) - norm_factor)

        return new_hap_probs1, new_hap_probs2, new_mol_probs
