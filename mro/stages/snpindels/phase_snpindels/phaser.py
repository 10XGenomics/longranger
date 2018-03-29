#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import operator
import scipy.sparse
import scipy.misc
import numpy as np
import math
import tenkit.hdf5
import tenkit.bio_io as tk_io
from tenkit.regions import Regions
from tenkit.constants import ILLUMINA_QUAL_OFFSET
import tenkit.stats as tk_stats
import martian
# If we find a local swap during reassignment optmization, go back this many variants to see if
# any other nearby positions are affected
REASSIGN_REWIND = 50

# Don't let the multiple rewinds get out of control
# only go back this far before giving up
REASSIGN_OVERRIDE_COUNT = 200

# Width of beam for local block beam search
DEFAULT_BEAM_WIDTH = 430

# If a variant is not phasing consistent it might be homozygous.
ALLELE_FREQ_HOMOZYGOUS_MIN = 0.7

class Phaser(object):
    """Constructs and outputs inferred haplotypes.
    """
    def __init__(self, vfr, fragments, chrom, start, stop, bc_mix_prob=0.0005, min_junction_hap_conf=0.999, min_var_hap_conf=0.999, block_buffer_size=2, hap_block_size=10, max_reassign_rounds=2, vc_mode = "call"):
        self.vfr = vfr
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.beam_width = DEFAULT_BEAM_WIDTH
        # caching for use with reassign_vars
        self.cache_overall_assign = None
        self.cache_hap_cols = None

        self.vc_mode = vc_mode

        self.bc_mix_prob = bc_mix_prob
        self.bc_not_mix_prob = 1.0-bc_mix_prob
        self.block_buffer_size = block_buffer_size
        self.max_junction_hap_error = 1.0 - min_junction_hap_conf
        self.max_var_hap_error = 1.0 - min_var_hap_conf
        self.mix_log_prior = np.log(bc_mix_prob)
        self.not_mix_log_prior = np.log(self.bc_not_mix_prob)
        self.hap_log_prior = np.log(0.5) + self.not_mix_log_prior
        self.hap_block_size = hap_block_size
        self.max_reassign_rounds = max_reassign_rounds
        self.snp_mixed_log_prior = math.log(0.001)
        self.snp_phased_log_prior = math.log(0.999)
        self.indel_mixed_log_prior = math.log(0.1)
        self.indel_phased_log_prior = math.log(0.9)

        # Create bc_indexes and bc_regions
        (sub_bc_indexes, bc_regions, self.contig_starts, self.contig_stops, self.contig_bcs, self.num_reads, self.molecule_ids) = self.get_bc_contig_info(chrom, start, stop, fragments)

        # First just record the positions seen
        poses_seen = set()
        poses_vartype_seen = set()
        for record in tk_io.get_variant_iterator_pos(vfr, None, tk_io.create_locus_info(self.chrom, self.start, self.stop)):
            if self.try_to_phase_variant(record):
                var_pos = tk_io.get_record_pos(record) - 1
                alts = tk_io.get_record_alt_alleles(record)
                ref = tk_io.get_record_ref(record)
                snp = True
                for alt in alts:
                    if tk_io.get_var_type(ref, alt) != "S":
                        snp = False
                if not(var_pos in poses_seen): # note, if you have multiple vcf entries with the same position, they will get the same phasing, known issue
                    poses_seen.add(var_pos)
                    poses_vartype_seen.add((var_pos, snp))

        self.pos_indexes = {}
        self.snp_positions = []
        self.pos_is_snp = []
        for (pos_index, (pos, is_snp)) in enumerate(sorted(list(poses_vartype_seen), key=lambda tup: tup[0])):
            self.pos_indexes[pos] = pos_index
            self.pos_is_snp.append(is_snp)
            self.snp_positions.append(pos)
        self.snp_positions = np.array(self.snp_positions, dtype=np.int32)

        # Alllocate the sparse matrices
        num_bcs = len(sub_bc_indexes.keys())
        num_poses = len(self.pos_indexes.keys())
        self.num_bcs = num_bcs
        self.num_poses = num_poses

        martian.log_info("Setting up phasing matrices. Variants: %d, Barcodes: %d" % (num_poses, num_bcs))
        if num_poses > 0:
            self.hap1_matrix, self.hap2_matrix, self.mix_matrix = self.get_bc_variant_matrices(num_bcs, num_poses, bc_regions, sub_bc_indexes)
        else:
            self.hap1_matrix = scipy.sparse.csc_matrix((num_bcs, num_poses), dtype=np.float32)
            self.hap2_matrix = scipy.sparse.csc_matrix((num_bcs, num_poses), dtype=np.float32)
            self.mix_matrix = scipy.sparse.csc_matrix((num_bcs, num_poses), dtype=np.float32)

    def try_to_phase_variant(self, record):
        ''' Decide whether or not to try and phase a variant.
            Currently this just ignores homozygous variants '''
        if tk_io.get_record_passes_filters(record):
            (genotype, phased) = tk_io.get_record_genotype_phased(record)
            if len(genotype) == 1 or genotype[0] == genotype[1]:
                return False
            return True
        else:
            return False



    def get_bc_contig_info(self, chrom, start, stop, fragments):
        """ Gets the information on called contigs by barcode
        """
        bc_contigs = {}
        bc_regions = {}
        sub_bc_indexes = {}
        barcode_index = 0
        bcs = []
        starts = []
        stops = []
        reads = []
        mol_ids = []

        if fragments is not None:
            query = [(chrom, start, stop)]
            chunk_frags = tenkit.hdf5.read_data_frame_indexed(fragments, query, query_cols=['chrom', 'start_pos', 'end_pos', 'bc', 'num_reads', 'molecule_id'])

            for (k,r) in chunk_frags.iterrows():
                start_stops = bc_contigs.setdefault(r.bc, [])
                start_stops.append((r.start_pos, r.end_pos))

                sub_bc = self.get_sub_bc(r.bc, r.start_pos, r.end_pos)
                sub_bc_indexes[sub_bc] = barcode_index
                barcode_index += 1
                bcs.append(r.bc)
                starts.append(r.start_pos)
                stops.append(r.end_pos)
                reads.append(r.num_reads)
                mol_ids.append(r.molecule_id)

            for (bc, start_stops) in bc_contigs.iteritems():
                bc_regions[bc] = Regions(regions=start_stops)

        return (sub_bc_indexes, bc_regions, np.array(starts, dtype=int), np.array(stops, dtype=int), bcs, np.array(reads,dtype=int), np.array(mol_ids, dtype=int))

    def get_specific_contig(self, bc_regions, bc, var_pos):
        """ Gets the specific contig the variant is contained on
        """
        if not bc_regions.has_key(bc):
            return None

        return bc_regions[bc].get_region_containing_point(var_pos)


    def get_bc_variant_matrices(self, num_bcs, num_poses, bc_regions, sub_bc_indexes):
        """ Returns a sparse matrix (BC x Variant Position) containing the probability
        that the barcode could have been drawn from the given haplotype at that position
        hap = 0 or 1
        if mixed = True, ignores hap and calculates the probability that the bc was drawn from a mix
        of the two haplotypes at that position
        """
        hap1_matrix = scipy.sparse.lil_matrix((num_bcs, num_poses), dtype=np.float32)
        hap2_matrix = scipy.sparse.lil_matrix((num_bcs, num_poses), dtype=np.float32)
        mix_matrix = scipy.sparse.lil_matrix((num_bcs, num_poses), dtype=np.float32)

        for record in tk_io.get_variant_iterator_pos(self.vfr, None, tk_io.create_locus_info(self.chrom, self.start, self.stop)):

            # Skip homozygous variants
            if not self.try_to_phase_variant(record):
                continue

            barcodes = tk_io.get_record_barcodes(record)
            qual = tk_io.get_record_qual(record)
            var_pos = tk_io.get_record_pos(record) - 1
            (genotype, phased) = tk_io.get_record_genotype_phased(record)

            # Only work with variants we have chosen to phase:
            # some filtering happens in __init__
            if not self.pos_indexes.has_key(var_pos):
                continue

            def get_bc(pos):
                if len(barcodes) > pos:
                    return barcodes[pos]
                else:
                    return []

            bcs_allele1 = get_bc(genotype[0])
            bcs_allele2 = get_bc(genotype[1])

            # Handle missing qual -- set at Q40
            if qual is None:
                qual = 40.0

            # Cap quality at a reasonable value
            qual = min(50.0, max(1.0, qual))
            prob_var_wrong = math.pow(10.0, -qual/10.0)

            bc_quals1 = {}
            for bc in bcs_allele1:
                (bc_seq, bc_quals) = bc.split('_', 1)
                bc_quals = bc_quals.split('_')

                bc_quals1[bc_seq] = [math.log(math.pow(10.0, -max(1, (int(bc_qual) - ILLUMINA_QUAL_OFFSET))/10.0)) for bc_qual in bc_quals]

            bc_quals2 = {}
            for bc in bcs_allele2:
                (bc_seq, bc_quals) = bc.split('_', 1)
                bc_quals = bc_quals.split('_')

                bc_quals2[bc_seq] = [math.log(math.pow(10.0, -max(1, (int(bc_qual) - ILLUMINA_QUAL_OFFSET))/10.0)) for bc_qual in bc_quals]
            # Double-check below
            seen_bcs = set()
            for bc, quals1 in bc_quals1.iteritems():
                ctg = self.get_specific_contig(bc_regions, bc, var_pos)
                if ctg is None:
                    continue
                else:
                    (start, stop) = ctg

                sub_bc = self.get_sub_bc(bc, start, stop)
                seen_bcs.add(bc)
                quals2 = bc_quals2.get(bc, [])
                total_length = len(quals1) + len(quals2)
                log_hap_prob_1 = sum([math.log(1.0 - math.exp(x)) for x in quals1]) + sum(quals2)
                log_hap_prob_2 = sum([math.log(1.0 - math.exp(x)) for x in quals2]) + sum(quals1)
                log_hap_prob_mix = sum([math.log(0.5) for j in xrange(total_length)])

                #log_mix_prob = math.log(self.bc_mix_prob)
                #log_not_mix_prob = math.log(0.5*(1.0 -self.bc_mix_prob))
                #log_sum_probs = scipy.misc.logsumexp([log_mix_prob + log_hap_prob_mix, log_not_mix_prob + log_hap_prob_1, log_not_mix_prob + log_hap_prob_2])

                log_prob_wrong = math.log(prob_var_wrong)
                log_prob_notwrong = math.log(1.0 - prob_var_wrong)

                hap1_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_1)
                hap2_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_2)
                mix_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_mix)

            for bc, quals2 in bc_quals2.iteritems():
                if bc in seen_bcs:
                    continue
                ctg = self.get_specific_contig(bc_regions, bc, var_pos)
                if ctg is None:
                    continue
                else:
                    (start, stop) = ctg

                sub_bc = self.get_sub_bc(bc, start, stop)
                quals1 = bc_quals1.get(bc, [])
                total_length = len(quals1) + len(quals2)
                log_hap_prob_1 = sum([math.log(1.0 - math.exp(x)) for x in quals1]) + sum(quals2)
                log_hap_prob_2 = sum([math.log(1.0 - math.exp(x)) for x in quals2]) + sum(quals1)
                log_hap_prob_mix = sum([math.log(0.5) for j in xrange(total_length)])
                #sum_probs = self.bc_mix_prob*hap_prob_mix + (1.0 - self.bc_mix_prob)*(0.5*hap_prob_1 + 0.5*hap_prob_2)
                log_prob_wrong = math.log(prob_var_wrong)
                log_prob_notwrong = math.log(1.0 - prob_var_wrong)

                hap1_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_1)
                hap2_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_2)
                mix_matrix[sub_bc_indexes[sub_bc], self.pos_indexes[var_pos]] = np.logaddexp(log_prob_wrong, log_prob_notwrong + log_hap_prob_mix)

        return (hap1_matrix.tocsc(), hap2_matrix.tocsc(), mix_matrix.tocsc())


    def _call_haps(self):
        """ Calls haplotypes and outputs the variant file
        """

        # Short circuit 0 SNP case.
        if self.num_poses == 0:
            return ([], {}, {}, set(), {})

        all_best_hap_blocks = self.get_all_best_hap_blocks()
        (hap1_block_matrix, hap2_block_matrix, mix_block_matrix) = self.calc_bc_block_matrices(all_best_hap_blocks)
        overall_hap_assign = self.stitch_blocks_together(hap1_block_matrix, hap2_block_matrix, mix_block_matrix, all_best_hap_blocks)

        changes_made = True
        for j in xrange(self.max_reassign_rounds):
            martian.log_info("Doing reassign round %s" % str(j))
            (overall_hap_assign, phase_perrs, changes_made) = self.reassign_vars(overall_hap_assign)
            if not(changes_made):
                break

        (phase_block_ids, junction_qvs) = self.get_phase_block_ids(self.hap1_matrix, self.hap2_matrix, self.mix_matrix, overall_hap_assign)
        phased_vars = self.get_phased_vars(phase_perrs, phase_block_ids)
        return (overall_hap_assign, phase_block_ids, junction_qvs, phased_vars, phase_perrs)

    def call_haps(self, vfw, fragment_phasing):
        (overall_hap_assign, phase_block_ids, junction_qvs, phased_vars, phase_perrs) = self._call_haps()
        phase_set_labels = self.write_fragment_phasing(overall_hap_assign, phase_block_ids, phased_vars, fragment_phasing)
        self.output_haps(self.vfr, vfw, overall_hap_assign, phase_block_ids, phase_set_labels, junction_qvs, phased_vars, phase_perrs)

    def write_fragment_phasing(self, overall_hap_assign, phase_block_ids, phased_vars, fragment_phasing):
        # keep track of phase set labels
        phase_set_labels = {}

        overall_hap_assign = np.array(overall_hap_assign)

        # Make groups of variants based on phase_block ids
        phase_block_vars = {}
        phase_block_bounds = {}

        for (var_idx, phase_block) in phase_block_ids.iteritems():
            var_list = phase_block_vars.setdefault(phase_block, [])
            var_list.append(var_idx)

        # keep track of the first / last variants in each block
        for (phase_block_id, var_list) in phase_block_vars.items():
            start_pos_idx = min(var_list)
            end_pos_idx = max(var_list)
            phase_block_bounds[phase_block_id] = (start_pos_idx, end_pos_idx)

        phase_block_ids_sorted = sorted(phase_block_vars.keys())

        num_blocks = len(phase_block_ids_sorted)

        for i in range(num_blocks):
            phase_block_id = phase_block_ids_sorted[i]
            var_list = phase_block_vars[phase_block_id]

            # determine set of fragments that cover one of these variants
            # also find the closest SNPs on either side, to make sure we don't overstep
            # if we're extending the phase blocks
            (start_pos_idx, end_pos_idx) = phase_block_bounds[phase_block_id]
            prev_end_idx = phase_block_bounds[phase_block_ids_sorted[i-1]][1] if i > 0 else None
            next_start_idx = phase_block_bounds[phase_block_ids_sorted[i+1]][0] if i < num_blocks-1 else None

            ps_start = self.snp_positions[start_pos_idx]
            ps_end = self.snp_positions[end_pos_idx]

            # now extend blocks if possible
            max_extend = 50000 # TODO move to constants

            if prev_end_idx is not None:
                dist_to_prev_snp = ps_start - self.snp_positions[prev_end_idx]
                start_extend = min(max_extend, dist_to_prev_snp / 2)
            else:
                start_extend = max_extend

            if next_start_idx is not None:
                dist_to_next_snp = self.snp_positions[next_start_idx] - ps_end
                end_extend = min(max_extend, dist_to_next_snp / 2)
            else:
                end_extend = max_extend

            # don't go beyond the chunk boundaries
            # and subtract one from the end to avoid overlaps
            ps_start = max(ps_start - start_extend, self.start)
            ps_end = min(ps_end + end_extend, self.stop) - 1

            #phase_block_name = self.get_phase_block_name(phase_block_id)
            phase_block_name = ps_start + 1
            phase_set_labels[phase_block_id] = phase_block_name

            assigned_vars = []
            for idx in range(start_pos_idx, end_pos_idx+1):
                assigned_vars.append((idx, overall_hap_assign[idx]))

            bc_inds = self.find_bc_inds(self.hap1_matrix, self.hap2_matrix, self.mix_matrix, start_pos_idx, end_pos_idx+1)

            if len(bc_inds) == 0 or len(assigned_vars) == 0:
                continue

            #restrict_assignments = overall_hap_assign[phased_cols]
            #restrict_hap1_matrix = np.array(self.hap1_matrix[bc_inds,:][:, phased_cols].self.log_space_todense())
            #restrict_hap2_matrix = np.array(self.hap2_matrix[bc_inds,:][:, phased_cols].self.log_space_todense())
            #restrict_mix_matrix = np.array(self.mix_matrix[bc_inds,:][:, phased_cols].self.log_space_todense())
            p = self.calc_fragment_hap_posterior(assigned_vars, bc_inds, self.hap1_matrix, self.hap2_matrix, self.mix_matrix)

            sorted_bc_inds = sorted([(self.contig_starts[ind], ind, chunk_ind) for (chunk_ind, ind) in enumerate(bc_inds)])

            for (f_start, bc_ind, chunk_idx) in sorted_bc_inds:

                bc = self.contig_bcs[bc_ind]
                frag_start = self.contig_starts[bc_ind]
                frag_end = self.contig_stops[bc_ind]
                h1 = p[chunk_idx, 0]
                h2 = p[chunk_idx, 1]
                hm = p[chunk_idx, 2]
                frag_reads = self.num_reads[bc_ind]
                mol_id = self.molecule_ids[bc_ind]
                # Write out phasing details:
                # ['chrom', 'frag_start', 'frag_end', 'phase_set', 'ps_start', 'ps_end', 'bc', 'h0', 'h1', 'hmix', 'reads']
                line = [self.chrom, frag_start, frag_end, phase_block_name, ps_start, ps_end, bc, h1, h2, hm, frag_reads, mol_id]
                line_str = "\t".join([str(x) for x in line])
                fragment_phasing.write(line_str + "\n")

        return phase_set_labels

    def output_haps(self, vfr, vfw, overall_hap_assign, phase_block_ids, phase_set_labels, junction_qvs, phased_vars, phase_perrs):
        """ Outputs the phased variants
        """
        records_to_write = []

        for record in tk_io.get_variant_iterator_pos(vfr, None, tk_io.create_locus_info(self.chrom, self.start, self.stop)):

            var_pos = tk_io.get_record_pos(record) - 1
            barcodes = tk_io.get_record_barcodes(record)

            # If we didn't attempt to phase the variant, deal with it here.
            if not self.pos_indexes.has_key(var_pos):
                # Set reasonable phasing marks on this variant
                if not tk_io.get_record_has_filter('UNSUPPORTED_GENOTYPE',record):
                    input_genotype, input_phased = tk_io.get_record_genotype_phased(record)
                    # Set homozygotes to phase
                    if len(input_genotype) > 1 and input_genotype[0] == input_genotype[1]:
                        tk_io.set_record_genotype_phased(record, input_genotype, True)
                    elif len(input_genotype) > 1:
                        # heterozygotes that were not phased should be marked as unphased, overwriting the incoming
                        # phasing status
                        tk_io.set_record_genotype_phased(record, input_genotype, False)

                records_to_write.append(record)
                continue
            else:
                input_genotype, input_phased = tk_io.get_record_genotype_phased(record)

            pos_index = self.pos_indexes[var_pos]
            hap_assign = overall_hap_assign[pos_index]
            phase_block_id = phase_block_ids[pos_index]

            if pos_index in phased_vars and hap_assign != -1:
                phased = True
            else:
                phased = False
                # If we're not phased, used a fixed hap_assign, to
                # eliminate noise when diffing the VCF
                if hap_assign != -1:
                    hap_assign = 0

            phase_qual = int(-10.0*math.log10(max(phase_perrs[pos_index], 10.0**(-25.5))))
            junction_qual = junction_qvs[pos_index]

            if hap_assign == 0 or hap_assign == -1:
                phased_genotype = [input_genotype[0], input_genotype[1]]
            else:
                phased_genotype = [input_genotype[1], input_genotype[0]]

            tk_io.set_record_genotype_phased(record, phased_genotype, phased)
            #tk_io.set_record_phase_set(record, self.get_phase_block_name(phase_block_id))
            tk_io.set_record_phase_set(record, phase_set_labels[phase_block_id])
            tk_io.set_record_phase_qual(record, phase_qual)
            tk_io.set_record_junction_qual(record, junction_qual)
            if (hap_assign == -1 or phased == False)  and ((self.vc_mode == 'call') or (self.vc_mode == "precalled_plus" and "TENX" in record.INFO)):
                genotype, homozygous = self.is_homozygous_candidate(record)
                posthpb = tk_io.get_record_post_homopolymer_bases(record)[0].upper()
                posthpc = tk_io.get_record_post_homopolymer_counts(record)
                altbase = tk_io.get_record_alt_alleles(record)[0][-1].upper() #last base of first alt
                if homozygous:
                    record.INFO['AC'] = 2
                    tk_io.set_record_genotype_phased(record, [genotype, genotype], True)
                elif hap_assign == -1:
                    tk_io.set_record_filters(record, ["10X_PHASING_INCONSISTENT"])
                elif tk_io.get_var_type(tk_io.get_record_ref(record),tk_io.get_record_alt_alleles(record)[0]) == "I" and posthpc >= 4 and posthpb == altbase:
                    tk_io.set_record_filters(record, ["10X_HOMOPOLYMER_UNPHASED_INSERTION"])
                else:
                    pass

            if barcodes is not None:
                tk_io.set_record_barcodes(record, barcodes)

            records_to_write.append(record)

        def get_next_het(i, direction):
            ''' Get the next het variant closest to i in given direction '''
            while True:
                if i >= len(records_to_write) or i < 0:
                    return None
                variant_i = records_to_write[i]
                if not self.pos_indexes.has_key(tk_io.get_record_pos(variant_i) - 1):
                    i = i + direction
                    continue
                else:
                    return variant_i

        # Run through the records and put a phase block onto the variants that were not phased
        # Use the enclosing phase block for variants within a block
        # Use a new phase block for a variant between two phase blocks
        for idx, record in enumerate(records_to_write):
            pos = tk_io.get_record_pos(record) - 1
            if not self.pos_indexes.has_key(pos):
                phased_after = get_next_het(idx, 1)
                phased_before = get_next_het(idx, -1)

                # default case: give this record its own phase set
                ps = pos + 1 # switch back to 1-indexing for PS

                # if record is between two het variants with the same phase set, use that phase set
                if phased_after is not None and phased_before is not None:
                    block_after = tk_io.get_record_phase_set(phased_after)
                    block_before = tk_io.get_record_phase_set(phased_before)
                    if block_after == block_before:
                        ps = block_after

                tk_io.set_record_phase_set(record, ps)

        for record in records_to_write:
            vfw.write_record(record)

    def get_phase_block_name(self, phase_block_id):
        pos_based_phase_block_id = self.snp_positions[phase_block_id] + 1 # switch to 1-indexing
        return pos_based_phase_block_id

    def get_all_best_hap_blocks(self):
        """ Returns the best haplotype assignments over all blocks
        Potential improvement:
        Calculate on larger (overlapping) blocks and only retain the inference on the non-overlapping portions
        (so that edge variants are used to infer haplotypes of centers of blocks, but not the edges)
        """
        hap1_matrix = self.hap1_matrix
        hap2_matrix = self.hap2_matrix
        mix_matrix = self.mix_matrix
        block_size = self.hap_block_size

        (garb, tot_pos_inds) = hap1_matrix.shape
        all_best_hap_blocks = []
        for start_pos_index in xrange(0, tot_pos_inds, block_size - 2*self.block_buffer_size):
            if start_pos_index == 0:
                keep_start_index = 0
            else:
                keep_start_index = self.block_buffer_size

            stop_pos_index = min(start_pos_index + block_size, tot_pos_inds)
            keep_stop_index = block_size - self.block_buffer_size

            if stop_pos_index == tot_pos_inds:
                keep_stop_index = block_size
            else:
                keep_stop_index = block_size - self.block_buffer_size

            (hap_assigns, score) = self.get_best_hap_block_beam(hap1_matrix, hap2_matrix, mix_matrix, start_pos_index, stop_pos_index, self.beam_width)
            all_best_hap_blocks.append(hap_assigns[keep_start_index:keep_stop_index])

            # Should be handling in for loop, but indexes are hard...
            if start_pos_index + keep_stop_index >= tot_pos_inds:
                break

        return all_best_hap_blocks

    def is_homozygous_candidate(self, record):
        alt_allele_counts = tk_io.get_record_alt_allele_counts(record)
        ref_allele_counts = tk_io.get_record_ref_allele_count(record)
        sample = record.samples[0]
        if (len(sample.gt_alleles) > 1) and (ref_allele_counts is not None) and (alt_allele_counts is not None):
            genotype_1 = int(sample.gt_alleles[0])
            genotype_2 = int(sample.gt_alleles[1])
            assert(genotype_1 != genotype_2)
            ref = tk_io.get_record_ref(record)
            alt_alleles = tk_io.get_record_alt_alleles(record)
            alleles = [ref] + alt_alleles
            allele_counts = [ref_allele_counts] + alt_allele_counts
            allele1 = alleles[genotype_1]
            allele2 = alleles[genotype_2]
            allele_count1 = allele_counts[genotype_1]
            allele_count2 = allele_counts[genotype_2]
            allele1_perc = tk_stats.robust_divide(float(allele_count1), float(allele_count2) + float(allele_count1))
            allele2_perc = tk_stats.robust_divide(float(allele_count2), float(allele_count1) + float(allele_count2))
            if allele1_perc >= ALLELE_FREQ_HOMOZYGOUS_MIN:
                if allele1 != ref:
                    return (genotype_1, True)
                else:
                    return (None, False)
            elif allele2_perc >= ALLELE_FREQ_HOMOZYGOUS_MIN:
                if allele2 != ref:
                    return (genotype_2, True)
                else:
                    return (None, False)
        return (None, False)

    def get_best_hap_block_beam(self, hap1_matrix, hap2_matrix, mix_matrix, start_pos_index, end_pos_index, beam=DEFAULT_BEAM_WIDTH):
        bc_inds = self.find_bc_inds(hap1_matrix, hap2_matrix, mix_matrix, start_pos_index, end_pos_index)

        if len(bc_inds) == 0:
            return ((0,)*(end_pos_index - start_pos_index), 0.0)

        restrict_hap1_matrix = np.array(hap1_matrix[bc_inds, start_pos_index:end_pos_index].todense())
        restrict_hap2_matrix = np.array(hap2_matrix[bc_inds, start_pos_index:end_pos_index].todense())
        restrict_mix_matrix = np.array(mix_matrix[bc_inds, start_pos_index:end_pos_index].todense())
        restrict_is_snp = self.pos_is_snp[start_pos_index:end_pos_index]
        # select an ordering for the SNPs.  the most confidently phased SNPs come first
        snp_confidence = np.abs(restrict_hap1_matrix - restrict_hap2_matrix).sum(axis=0)
        snp_conf_index = sorted([(v,idx) for (idx,v) in enumerate(snp_confidence)])
        snp_order = range(len(snp_conf_index))#[idx for (v,idx) in snp_conf_index]
        num_bcs = len(bc_inds)
        num_poses = end_pos_index - start_pos_index
        log_probs = self.create_initial_log_probs(num_bcs)

        for pos_index in snp_order:
            log_probs = self.update_log_probs(restrict_hap1_matrix, restrict_hap2_matrix, restrict_mix_matrix, pos_index, log_probs, restrict_is_snp, False)
            # Don't bother scoring if we don't need to eliminate any options
            #if len(log_probs) < beam:
            #    continue

            scored_hap_assigns = []
            for item in log_probs.iteritems():
                (hap_assigns, (hap1_log_probs, hap2_log_probs, mix_log_probs, running_prior)) = item
                total_log_prob = self.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_log_probs, running_prior)
                scored_hap_assigns.append((total_log_prob, item))

            # Let only the beam highest scoring options through to the next round
            scored_hap_assigns.sort(reverse=True)
            best_score = scored_hap_assigns[0][0]
            index_cutoff = 0
            while index_cutoff < min(beam, len(scored_hap_assigns)) and scored_hap_assigns[0][0] - scored_hap_assigns[index_cutoff][0] < 40:
                index_cutoff += 1
            scored_hap_assigns = scored_hap_assigns[:index_cutoff]
            log_probs = {assign:probs for (score, (assign, probs)) in scored_hap_assigns}

        # Pick winning assignment and put back into original order
        (best_score, (hap_assigns, (hap1_log_probs, hap2_log_probs, mix_log_probs, running_prior))) = scored_hap_assigns[0]
        reorder_hap_assigns = [0] * num_poses
        for call_idx, orig_idx in enumerate(snp_order):
            reorder_hap_assigns[orig_idx] = hap_assigns[call_idx]
        return (tuple(reorder_hap_assigns), best_score)

    def find_bc_inds(self, hap1_matrix, hap2_matrix, mix_matrix, start_pos_index, end_pos_index):
        """ Returns the barcode indices that computation needs to be performed on for the relevant
        positions.
        """
        # The non-zero values should be the same for hap1_matrix, hap2_matrix, mix_matrix
        bcs, garb = hap1_matrix[:, start_pos_index:end_pos_index].nonzero()
        return np.unique(bcs)

    def create_initial_log_probs(self, num_bcs):
        """ Creates the initial data structure to hold the bc log probabilities for a given hap assignment
        Stored as:
        key: tuple of haplotype assignments (0 or 1)
        values: 3-tuple: (hap1_values, hap2_values, mix_values)
        """
        return {():(np.zeros((1, num_bcs)), np.zeros((1,num_bcs)), np.zeros((1,num_bcs)), 0)}


    def update_log_probs(self, hap1_matrix, hap2_matrix, mix_matrix, pos_index, log_probs, is_snp, block_mode):
        """ Updates the new log probabilities
        """
        new_log_probs = {}

        # The very first haplotype assignment is degenerate so just set it to zero
        # this halves the search space
        #hap_alternatives = 2
        if block_mode:
            hap_alternatives = [0,1]
        elif log_probs.has_key(()):
            hap_alternatives = [0]
        else:
            hap_alternatives = [0,1,-1] # -1 denotes that the variant is being assigned not a heterozygous variant
        # -1 last because in the search it should first by tried to be phased
        for new_hap_assign in hap_alternatives:
            for old_assigns, (old_hap1_log_probs, old_hap2_log_probs, old_mix_log_probs, old_running_prior) in log_probs.iteritems():
                new_assigns = old_assigns + (new_hap_assign,)
                new_running_prior = old_running_prior
                if is_snp is None:
                    new_running_prior = 0
                elif new_assigns >= 0:
                    if is_snp[pos_index]:
                        new_running_prior += self.snp_phased_log_prior
                    else:
                        new_running_prior += self.indel_phased_log_prior
                else:
                    if is_snp[pos_index]:
                        new_running_prior += self.snp_mixed_log_prior
                    else:
                        new_running_prior += self.indel_mixed_log_prior

                if new_hap_assign == 0:
                    new_hap1_log_probs = old_hap1_log_probs + hap1_matrix[:, pos_index]
                    new_hap2_log_probs = old_hap2_log_probs + hap2_matrix[:, pos_index]
                elif new_hap_assign == 1:
                    new_hap1_log_probs = old_hap1_log_probs + hap2_matrix[:, pos_index]
                    new_hap2_log_probs = old_hap2_log_probs + hap1_matrix[:, pos_index]
                else:
                    new_hap1_log_probs = old_hap1_log_probs + mix_matrix[:, pos_index]
                    new_hap2_log_probs = old_hap2_log_probs + mix_matrix[:, pos_index]

                new_mix_log_probs = old_mix_log_probs + mix_matrix[:, pos_index]

                new_log_probs[new_assigns] = (new_hap1_log_probs, new_hap2_log_probs, new_mix_log_probs, new_running_prior)

        return new_log_probs

    def calc_total_log_prob(self, hap1_log_probs, hap2_log_probs, mix_log_probs, running_prior):
        """ Given the bc log prob vectors for haplotype assignment, calculates the overall probability
        """
        if len(hap1_log_probs) == 0:
            return 0.0
        log_probs = np.logaddexp(self.hap_log_prior + np.logaddexp(hap1_log_probs,  hap2_log_probs), self.mix_log_prior + mix_log_probs)
        log_probs = np.sum(log_probs)
        log_probs += running_prior
        return log_probs

    def calc_fragment_hap_posterior(self, assigned_vars, bc_inds, allele1_log_probs, allele2_log_probs, mix_log_probs):
        """ Calculate the posterior distribution of the haplotype state of each fragment """

        bc_mix_prob = self.bc_mix_prob

        p1 = np.zeros(len(bc_inds), dtype=np.float64)
        p2 = np.zeros(len(bc_inds), dtype=np.float64)
        p_mix = np.zeros(len(bc_inds), dtype=np.float64)

        l1m = np.log(0.5 * (1.0 - bc_mix_prob))
        l_mix = np.log(bc_mix_prob)

        for (col_idx, assignment) in assigned_vars:
            if assignment == 0:
                p1 += allele1_log_probs[bc_inds,col_idx].toarray().squeeze()
                p2 += allele2_log_probs[bc_inds,col_idx].toarray().squeeze()
            elif assignment == 1:
                p1 += allele2_log_probs[bc_inds,col_idx].toarray().squeeze()
                p2 += allele1_log_probs[bc_inds,col_idx].toarray().squeeze()
            else:
                p1 += mix_log_probs[bc_inds, col_idx].toarray().squeeze()
                p2 += mix_log_probs[bc_inds, col_idx].toarray().squeeze()

            p_mix += mix_log_probs[bc_inds, col_idx].toarray().squeeze()
        p1 += l1m
        p2 += l1m
        p_mix += l_mix

        norm = scipy.misc.logsumexp([p1, p2, p_mix], axis=0)
        dists = np.exp(np.array([p1 - norm, p2 - norm, p_mix - norm])).transpose()
        return dists


    def calc_bc_block_matrices(self, block_assigns):
        """ Calculates the log probability matrices for blocks
        """
        hap1_matrix = self.hap1_matrix
        hap2_matrix = self.hap2_matrix
        mix_matrix = self.mix_matrix

        (num_bcs, num_poses) = hap1_matrix.shape
        num_blocks = len(block_assigns)
        hap1_block_matrix = scipy.sparse.lil_matrix((num_bcs, num_blocks))
        hap2_block_matrix = scipy.sparse.lil_matrix((num_bcs, num_blocks))
        mix_block_matrix = scipy.sparse.lil_matrix((num_bcs, num_blocks))
        pos_offset = 0
        for block_index in xrange(num_blocks):
            for assign_index in xrange(len(block_assigns[block_index])):
                pos_index = pos_offset + assign_index
                bc_inds = self.find_bc_inds(hap1_matrix, hap2_matrix, mix_matrix, pos_index, pos_index + 1)
                if len(bc_inds) == 0:
                    continue
                if block_assigns[block_index][assign_index] == 0:
                    hap1_block_matrix[bc_inds, block_index] += hap1_matrix[bc_inds, pos_index]
                    hap2_block_matrix[bc_inds, block_index] += hap2_matrix[bc_inds, pos_index]
                elif block_assigns[block_index][assign_index] == 1:
                    hap1_block_matrix[bc_inds, block_index] += hap2_matrix[bc_inds, pos_index]
                    hap2_block_matrix[bc_inds, block_index] += hap1_matrix[bc_inds, pos_index]
                else:
                    hap1_block_matrix[bc_inds, block_index] += mix_matrix[bc_inds, pos_index]
                    hap2_block_matrix[bc_inds, block_index] += mix_matrix[bc_inds, pos_index]
                mix_block_matrix[bc_inds, block_index] += mix_matrix[bc_inds, pos_index]

            pos_offset += len(block_assigns[block_index])

        return (hap1_block_matrix.tocsc(), hap2_block_matrix.tocsc(), mix_block_matrix.tocsc())


    def invert(self, assignment):
        if assignment >= 0:
            return 1-assignment
        else:
            return assignment

    def stitch_blocks_together(self, hap1_matrix, hap2_matrix, mix_matrix, block_assigns):
        """ Greedily stitches blocks together
        """
        modified_block_assigns = [block_assigns[0]]
        last_hap_assign = 0
        for block_index in xrange(len(block_assigns) - 1):
            bc_inds = self.find_bc_inds(hap1_matrix, hap2_matrix, mix_matrix, block_index, block_index + 1)
            num_bcs = len(bc_inds)

            if num_bcs == 0:
                modified_block_assigns.append(block_assigns[block_index + 1])
                continue

            num_poses = 2

            restrict_hap1_matrix = np.array(hap1_matrix[bc_inds, block_index:block_index + 2].todense())
            restrict_hap2_matrix = np.array(hap2_matrix[bc_inds, block_index:block_index + 2].todense())
            restrict_mix_matrix = np.array(mix_matrix[bc_inds, block_index:block_index + 2].todense())

            log_probs = self.create_initial_log_probs(num_bcs)
            for pos_index in xrange(num_poses):
                log_probs = self.update_log_probs(restrict_hap1_matrix, restrict_hap2_matrix, restrict_mix_matrix, pos_index, log_probs, None, True)

            hap1_log_probs, hap2_log_probs, mix_log_probs, prior = log_probs[(last_hap_assign, 0)]
            total_log_prob_0 = self.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_log_probs, 0)
            hap1_log_probs, hap2_log_probs, mix_log_probs, prior = log_probs[(last_hap_assign, 1)]
            total_log_prob_1 = self.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_log_probs, 0) # last argument is assignments used for priors of 0, 1, -1 but since they will be the same number of -1's in both of these it does not matter

            if total_log_prob_0 > total_log_prob_1:
                modified_block_assigns.append(block_assigns[block_index + 1])
                last_hap_assign = 0
            else:
                modified_block_assigns.append(tuple([self.invert(x) for x in block_assigns[block_index + 1]]))
                last_hap_assign = 1

        return reduce(operator.concat, modified_block_assigns)


    def reassign_vars(self, overall_assign):
        """ Greedily reassigns haplotype of variants if switching haplotype improves log
        probability.
        """

        (garb, tot_pos_inds) = self.hap1_matrix.shape
        change_made = False
        new_overall_assign = list(overall_assign)
        phase_perrs = {}
        pos_index = 0
        furthest_pos = 0
        furthest_pos_counter = 0
        trial_counter = 0

        while pos_index < tot_pos_inds:

            if pos_index > furthest_pos:
                furthest_pos = pos_index
                furthest_pos_counter = trial_counter

            (hap1_log_prob, hap2_log_prob, mix_log_prob) = self.calc_hap1_hap2_logprob(new_overall_assign, pos_index)
            max_log_prob = max(hap1_log_prob, hap2_log_prob, mix_log_prob)
            hap1_prob = math.exp(hap1_log_prob - max_log_prob)
            hap2_prob = math.exp(hap2_log_prob - max_log_prob)
            mix_prob = math.exp(mix_log_prob - max_log_prob)

            if hap1_log_prob > hap2_log_prob and hap1_log_prob > mix_log_prob:
                phase_perr = (hap2_prob + mix_prob)/(hap1_prob + hap2_prob + mix_prob)
                phase_perrs[pos_index] = phase_perr
                new_overall_assign[pos_index] = 0
                if overall_assign[pos_index] != 0 and phase_perr < 0.2:
                    change_made = True
                    pos_index = max(0, pos_index-REASSIGN_REWIND)
            elif hap2_log_prob > mix_log_prob:
                phase_perr = (hap1_prob + mix_prob)/(hap1_prob + hap2_prob + mix_prob)
                phase_perrs[pos_index] = phase_perr
                new_overall_assign[pos_index] = 1
                if overall_assign[pos_index] != 1 and phase_perr < 0.2:
                    change_made = True
                    pos_index = max(0, pos_index-REASSIGN_REWIND)
            else:
                phase_perr = (hap1_prob + hap2_prob)/(hap1_prob + hap2_prob + mix_prob)
                phase_perrs[pos_index] = phase_perr
                new_overall_assign[pos_index] = -1
                if overall_assign[pos_index] != -1 and phase_perr < 0.2:
                    change_made = True
                    pos_index = max(0, pos_index-REASSIGN_REWIND)

            pos_index += 1
            trial_counter += 1

            if trial_counter - furthest_pos_counter > REASSIGN_OVERRIDE_COUNT:
                pos_index = furthest_pos + 1

        return (new_overall_assign, phase_perrs, change_made)

    def get_phased_vars(self, phase_perrs, phase_block_ids):
        """ Gets the set of pos_indexes for variants that are phased beyond the given confidence, and are in phase blocks with >1 SNP
        """
        phased_vars = set()
        num_snps = len(phase_perrs)

        for pos_index, phase_perr in phase_perrs.iteritems():
            phase_block = phase_block_ids[pos_index]

            # Phasing confidence is high
            if phase_perr < self.max_var_hap_error:
                # Phase block size > 1
                if (pos_index > 0 and phase_block == phase_block_ids[pos_index-1]) or \
                        (pos_index < num_snps - 1 and phase_block == phase_block_ids[pos_index+1]):
                    phased_vars.add(pos_index)

        return phased_vars


    def get_hap_cols(self, overall_assign):
        if overall_assign == self.cache_overall_assign:
            return self.cache_hap_cols
        else:
            # Compute haplotype columns
            hap1_poses = [i for (i,v) in enumerate(overall_assign) if v == 0]
            hap2_poses = [i for (i,v) in enumerate(overall_assign) if v == 1]
            mix_poses = [i for (i,v) in enumerate(overall_assign) if v == -1]
            if len(hap1_poses) == 0 and len(mix_poses) == 0:
                hap1_col = self.hap2_matrix[:,hap2_poses].sum(axis=1)
                hap2_col = self.hap1_matrix[:,hap2_poses].sum(axis=1)
            elif len(hap2_poses) == 0 and len(mix_poses) == 0:
                hap1_col = self.hap1_matrix[:,hap1_poses].sum(axis=1)
                hap2_col = self.hap2_matrix[:,hap1_poses].sum(axis=1)
            elif len(hap1_poses) == 0:
                hap1_col = self.hap2_matrix[:,hap2_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
                hap2_col = self.hap1_matrix[:,hap2_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
            elif len(hap2_poses) == 0:
                hap1_col = self.hap1_matrix[:,hap1_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
                hap2_col = self.hap2_matrix[:,hap1_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
            elif len(mix_poses) == 0:
                hap1_col = self.hap1_matrix[:,hap1_poses].sum(axis=1) + self.hap2_matrix[:,hap2_poses].sum(axis=1)
                hap2_col = self.hap2_matrix[:,hap1_poses].sum(axis=1) + self.hap1_matrix[:,hap2_poses].sum(axis=1)
            else:
                hap1_col = self.hap1_matrix[:,hap1_poses].sum(axis=1) + self.hap2_matrix[:,hap2_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
                hap2_col = self.hap2_matrix[:,hap1_poses].sum(axis=1) + self.hap1_matrix[:,hap2_poses].sum(axis=1) + self.mix_matrix[:,mix_poses].sum(axis=1)
            mix_col = self.mix_matrix.sum(axis=1)

            hap1_col = np.array(hap1_col).flatten()
            hap2_col = np.array(hap2_col).flatten()
            mix_col =  np.array(mix_col).flatten()

            self.cache_overall_assign = [a for a in overall_assign]
            self.cache_hap_cols = (hap1_col, hap2_col, mix_col)
            return self.cache_hap_cols


    def calc_hap1_hap2_logprob(self, overall_assign, pos_index):
        """ Calculates the probability that the variant was on haplotype 1 or haplotype 2 given the assignment at all other
        locations.
        """
        bc_inds = self.find_bc_inds(self.hap1_matrix, self.hap2_matrix, self.mix_matrix, pos_index, pos_index + 1)

        if len(bc_inds) == 0:
            return (0.0, 0.0, 0.0)

        hap1_col, hap2_col, mix_col = self.get_hap_cols(overall_assign)
        hap1_col = hap1_col[bc_inds]
        hap2_col = hap2_col[bc_inds]
        mix_col  = mix_col[bc_inds]



        # subtract out the values from this position
        if overall_assign[pos_index] == 0:
            pos_hap1_col = np.matrix(self.hap1_matrix[:, pos_index].todense()).A1
            pos_hap2_col = np.matrix(self.hap2_matrix[:, pos_index].todense()).A1
        else:
            pos_hap1_col = np.matrix(self.hap2_matrix[:, pos_index].todense()).A1
            pos_hap2_col = np.matrix(self.hap1_matrix[:, pos_index].todense()).A1
        pos_mix_col = np.matrix(self.mix_matrix[:,pos_index].todense()).A1

        pos_hap1_col = pos_hap1_col[bc_inds]
        pos_hap2_col = pos_hap2_col[bc_inds]
        pos_mix_col = pos_mix_col[bc_inds]

        if overall_assign[pos_index] >= 0:
            hap1_hap2_col = hap1_col - pos_hap1_col + pos_hap2_col
            hap2_hap1_col = hap2_col - pos_hap2_col + pos_hap1_col
            hap1_mix_col = hap1_col - pos_hap1_col + pos_mix_col
            hap2_mix_col = hap2_col - pos_hap2_col + pos_mix_col
            hap1_hap1_col = hap1_col
            hap2_hap2_col = hap2_col
            if self.pos_is_snp[pos_index]:
                old_prior = self.snp_phased_log_prior
            else:
                old_prior = self.indel_phased_log_prior
        else:
            hap1_hap2_col = hap1_col - pos_mix_col + pos_hap2_col
            hap2_hap1_col = hap2_col - pos_mix_col + pos_hap1_col
            hap1_mix_col = hap1_col
            hap2_mix_col = hap2_col
            hap1_hap1_col = hap1_col - pos_mix_col + pos_hap1_col
            hap2_hap2_col = hap2_col - pos_mix_col + pos_hap2_col
            if self.pos_is_snp[pos_index]:
                old_prior = self.snp_mixed_log_prior
            else:
                old_prior = self.indel_mixed_log_prior

        keep_total_log_prob = self.calc_total_log_prob(hap1_col, hap2_col, mix_col, old_prior)
        phase_prior = 0
        mix_prior = 0
        if self.pos_is_snp[pos_index]:
            phase_prior = self.snp_phased_log_prior
            mix_prior = self.snp_mixed_log_prior
        else:
            phase_prior = self.indel_phased_log_prior
            mix_prior = self.indel_mixed_log_prior
        swap_total_log_prob = self.calc_total_log_prob(hap1_hap2_col, hap2_hap1_col, mix_col, phase_prior)
        if overall_assign[pos_index] == -1:
            swap2_total_log_prob = self.calc_total_log_prob(hap1_hap1_col, hap2_hap2_col, mix_col, phase_prior)
        else:
            swap2_total_log_prob = self.calc_total_log_prob(hap1_mix_col, hap2_mix_col, mix_col,  mix_prior)

        if overall_assign[pos_index] == 0:
            return (keep_total_log_prob, swap_total_log_prob, swap2_total_log_prob)
        elif overall_assign[pos_index] == 1:
            return (swap_total_log_prob, keep_total_log_prob, swap2_total_log_prob)
        else:
            return (swap_total_log_prob, swap2_total_log_prob, keep_total_log_prob)

    def get_phase_block_ids(self, hap1_matrix, hap2_matrix, mix_matrix, overall_assign):
        """ Breaks the phasing into phase blocks confidently phased above the given confidence threshold
        """
        phase_block_ids = {}
        junction_qvs = {}
        (tot_bc_inds, tot_pos_inds) = hap1_matrix.shape

        left_hap1_log_probs = np.zeros((tot_bc_inds, 1), dtype=np.float)
        left_hap2_log_probs = np.zeros((tot_bc_inds, 1), dtype=np.float)

        right_hap1_log_probs = np.zeros((tot_bc_inds, 1), dtype=np.float)
        right_hap2_log_probs = np.zeros((tot_bc_inds, 1), dtype=np.float)

        mix_log_probs = np.zeros((tot_bc_inds, 1), dtype=np.float32)

        # First populate the right side log probs
        for pos_index in xrange(0, tot_pos_inds):
            hap1_col = hap1_matrix[:,pos_index]
            hap2_col = hap2_matrix[:,pos_index]

            if overall_assign[pos_index] == 0:
                right_hap1_log_probs[hap1_col.nonzero()] += hap1_col.data
                right_hap2_log_probs[hap2_col.nonzero()] += hap2_col.data
            else:
                right_hap1_log_probs[hap2_col.nonzero()] += hap2_col.data
                right_hap2_log_probs[hap1_col.nonzero()] += hap1_col.data
            mix_col = mix_matrix[:,pos_index]
            mix_log_probs[mix_col.nonzero()] += mix_col.data

        block_id = 0
        phase_block_ids[0] = block_id
        junction_qvs[0] = 255

        for pos_index in xrange(1, tot_pos_inds):
            hap1_col = hap1_matrix[:,pos_index - 1]
            hap2_col = hap2_matrix[:,pos_index - 1]
            if overall_assign[pos_index-1] == 0:
                left_hap1_log_probs[hap1_col.nonzero()] += hap1_col.data
                left_hap2_log_probs[hap2_col.nonzero()] += hap2_col.data

                right_hap1_log_probs[hap1_col.nonzero()] -= hap1_col.data
                right_hap2_log_probs[hap2_col.nonzero()] -= hap2_col.data
            else:
                left_hap1_log_probs[hap2_col.nonzero()] += hap2_col.data
                left_hap2_log_probs[hap1_col.nonzero()] += hap1_col.data

                right_hap1_log_probs[hap2_col.nonzero()] -= hap2_col.data
                right_hap2_log_probs[hap1_col.nonzero()] -= hap1_col.data

            # Select out only the row of the score matrix that can help link this column
            snp_pos = self.snp_positions[pos_index]
            vrows = np.logical_and(self.contig_starts <= snp_pos, self.contig_stops > snp_pos)

            l_hap1_v = left_hap1_log_probs[vrows]
            l_hap2_v = left_hap2_log_probs[vrows]
            r_hap1_v = right_hap1_log_probs[vrows]
            r_hap2_v = right_hap2_log_probs[vrows]
            mix_v = mix_log_probs[vrows]

            hap1_log_probs = l_hap1_v + r_hap1_v
            hap2_log_probs = l_hap2_v + r_hap2_v

            switch_hap1_log_probs = l_hap1_v + r_hap2_v
            switch_hap2_log_probs = l_hap2_v + r_hap1_v

            total_log_prob = self.calc_total_log_prob(hap1_log_probs, hap2_log_probs, mix_v, 0)
            switch_total_log_prob = self.calc_total_log_prob(switch_hap1_log_probs, switch_hap2_log_probs, mix_v, 0)

            max_logprob = max(total_log_prob, switch_total_log_prob)
            total_prob = math.exp(total_log_prob - max_logprob)
            switch_total_log_prob = switch_total_log_prob - max_logprob

            junction_error_log_prob = switch_total_log_prob - np.log(total_prob + math.exp(switch_total_log_prob))

            if junction_error_log_prob > np.log(self.max_junction_hap_error):
                block_id = pos_index
            phase_block_ids[pos_index] = block_id
            junction_qvs[pos_index] = int(min(255.0, round(-10.0 * junction_error_log_prob / math.log(10) )))

        return (phase_block_ids, junction_qvs)

    def get_sub_bc(self, bc, start, stop):
        return '_'.join([bc, str(start), str(stop)])
