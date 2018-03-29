#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Trio tools
#

import itertools
import numpy as np
import tenkit.bio_io as tk_io

# For each variant, there are three dimensions in the state space:
# -> whether parent 1 transmitted hap 0 or 1
# -> whether parent 2 transmitted hap 0 or 1
# -> whether order of child haps is (P1|P2) or (P2|P1)
# Thus, 8 possible states are searched over.
state_idx = list(itertools.product((0,1),(0,1),(0,1)))
states = range(len(state_idx))

class TrioRow:
    def __init__(self, trio_phasing_row):
        self.alleles =    [v[0] for v in trio_phasing_row]
        self.phased =     [v[1] for v in trio_phasing_row]
        self.phase_sets = [v[2] for v in trio_phasing_row]

# Event masks
P1_ERROR = 0x1
P2_ERROR = 0x2

P1_UNPHASED = 0x10
P2_UNPHASED = 0x20
CH_UNPHASED = 0x40

P1_LONG_SWITCH = 0x100
P2_LONG_SWITCH = 0x200
CH_LONG_SWITCH = 0x400

MENDELIAN_VIOLATION = 0x1000 # a Mendelian violation in the trio at this site

error_bits = [
    ('p1_error', P1_ERROR),
    ('p2_error', P2_ERROR),

    ('p1_unphased', P1_UNPHASED),
    ('p2_unphased', P2_UNPHASED),
    ('ch_unphased', CH_UNPHASED),

    ('p1_long_switch', P1_LONG_SWITCH),
    ('p2_long_switch', P2_LONG_SWITCH),
    ('ch_long_switch', CH_LONG_SWITCH),

    ('mendelian_violation', MENDELIAN_VIOLATION),
]

def print_error_tags(v):
    str = ""
    for (name, bit) in error_bits:
        if v & bit:
            str += name + " "

    return str

def get_error_tags(v):
    tags = set()
    for (name, bit) in error_bits:
        if v & bit:
            tags.add(name)

    return tags

MATCH_ERROR = 1 # penalty for a single mismatch between parent / child haps
SWITCH_ERROR = 25 # penalty for switching transmitted allele within a single phase block

def is_mendelian_violation(p1_alleles, p2_alleles, ch_alleles):
    p1_shared = set(p1_alleles) & set(ch_alleles)
    p2_shared = set(p2_alleles) & set(ch_alleles)
    p1_violation = len(p1_shared) == 0
    p2_violation = len(p2_shared) == 0
    # px_violation - there's a violation, but we can't isolate it to one parent
    # e.g. P1 = A/A, P2 = A/A, CH = A/T
    px_violation = (len(p1_shared & p2_shared) == 1) and (p1_shared == p2_shared) and (len(set(ch_alleles)) == 2)
    return p1_violation or p2_violation or px_violation

def cands(al, phased, state):
    '''
    Get the set of canditate alleles:
    If this variant is phased, it's the allele corresponding to this phasing state.
    If unphased, then it's both alleles.
    '''
    if phased:
        return set([al[state]])
    else:
        return set(al)

def compat(al1, phased1, state1, al2, phased2, state2):
    '''
    Given the parent and child genotypes, and which alleles are supposedly being inherited,
    see if the actual haplotypes are consistent.
    '''
    r = len(cands(al1,phased1,state1) & cands(al2,phased2,state2)) > 0
    return r

def infer_trio_pattern(trio_data):
    ''' Match transmitted haplotypes and child phasing.'''

    n = len(trio_data)

    score = np.zeros((n, len(states)), dtype=np.int32)
    traceback = np.zeros((n, len(states)), dtype=np.int8)
    error_track = np.zeros((n, len(states)), dtype=np.int32)

    for i in range(n):
        # Phase sets
        curr_trio_row = TrioRow(trio_data[i])

        (ps_p1, ps_p2, ps_ch) = curr_trio_row.phase_sets
        (p1_alleles, p2_alleles, ch_alleles) = curr_trio_row.alleles
        (p1_phased, p2_phased, ch_phased) = curr_trio_row.phased

        for s_current in states:
            (p1_trans, p2_trans, ch_phasing) = state_idx[s_current]

            # State mask - describes the set of error modes in the data given this state
            state_mask = 0

            # State penalty - describes the matching score of the alleles in this state
            state_penalty = 0

            if is_mendelian_violation(p1_alleles, p2_alleles, ch_alleles):
                # keep track of this separately
                state_mask |= MENDELIAN_VIOLATION

            if p1_alleles[p1_trans] != ch_alleles[ch_phasing]:
                state_penalty += MATCH_ERROR

                # Annotate likely reason for error
                if compat(p1_alleles, p1_phased, p1_trans, ch_alleles, ch_phased, ch_phasing):
                    state_mask |= P1_UNPHASED if not p1_phased else 0
                    state_mask |= CH_UNPHASED if not ch_phased else 0
                else:
                    state_mask |= P1_ERROR

            if p2_alleles[p2_trans] != ch_alleles[1 - ch_phasing]:
                state_penalty += MATCH_ERROR

                # Annotate likely reason for error
                if compat(p2_alleles, p2_phased, p2_trans, ch_alleles, ch_phased, 1-ch_phasing):
                    state_mask |= P2_UNPHASED if not p2_phased else 0
                    state_mask |= CH_UNPHASED if not ch_phased else 0
                else:
                    state_mask |= P2_ERROR

            if i == 0:
                # initial state (no transitions to calculate)
                score[i, s_current] = state_penalty
                error_track[i, s_current] = state_mask

            else:
                best_score = 9999999
                best_prev_state = -1

                prev_trio_row = TrioRow(trio_data[i-1])

                for s_prev in states:
                    (prev_p1_trans, prev_p2_trans, prev_ch_phasing) = state_idx[s_prev]
                    (prev_ps_p1, prev_ps_p2, prev_ps_ch) = prev_trio_row.phase_sets

                    prev_score = score[i-1, s_prev]
                    trans_penalty = 0
                    trans_mask = state_mask

                    # Phasing switches -- only incurred if (a) we're in the same phase set
                    # and (b) the current variant is phased
                    if p1_phased and prev_ps_p1 == ps_p1 and p1_trans != prev_p1_trans:
                        trans_penalty += SWITCH_ERROR
                        trans_mask |= P1_LONG_SWITCH

                    if p2_phased and prev_ps_p2 == ps_p2 and p2_trans != prev_p2_trans:
                        trans_penalty += SWITCH_ERROR
                        trans_mask |= P2_LONG_SWITCH

                    if ch_phased and prev_ps_ch == ps_ch and ch_phasing != prev_ch_phasing:
                        trans_penalty += SWITCH_ERROR
                        trans_mask |= CH_LONG_SWITCH

                    total_penalty = prev_score + state_penalty + trans_penalty

                    if total_penalty < best_score:
                        best_trans_mask = trans_mask
                        best_score = total_penalty
                        best_prev_state = s_prev

                # Save the state
                traceback[i, s_current] = best_prev_state
                score[i, s_current] = best_score
                error_track[i, s_current] = best_trans_mask

    # Traceback to find the optimal matching
    best_final_score = 99999
    best_final_state = 0

    for s in states:
        if score[n-1, s] < best_final_score:
            best_final_state = s
            best_final_score = score[n-1, s]

    best_state = [0] * n
    best_state[n - 1] = best_final_state

    for i in range(n-2,-1,-1):
        prev_state = best_state[i+1]
        best_state[i] = traceback[i+1, prev_state]

    error_vec = [get_error_tags(error_track[i, best_state[i]]) for i in range(n)]
    best_score = [score[i, best_state[i]] for i in range(n)]
    state_vec = [state_idx[s] for s in best_state]

    return (state_vec, error_vec, best_score)

class Call:
    def __init__(self, current_phase_set, record):
        self.chrom = tk_io.get_record_chrom(record)
        self.pos = tk_io.get_record_pos(record)
        self.key = (self.chrom, self.pos)

        self.ref = tk_io.get_record_ref(record)
        self.filters = tk_io.get_record_passes_filters(record)

        alt_alleles = tk_io.get_record_alt_alleles(record)
        all_alleles = [self.ref] + alt_alleles

        (genotype, self.phased) = tk_io.get_record_genotype_phased(record)

        # always set homozygotes as phased
        if genotype[0] == genotype[1]:
            self.phased = True

        # note -- if there are two alts, this will just pick one.
        self.phase_set = current_phase_set
        self.hap = (all_alleles[genotype[0]], all_alleles[genotype[1]])

        self.record = record

# TODO - many of these methods are shared with analyze_snpindel_calls. they should be consolidated.

def variant_tuples(variant_iter):
    """
    Return iterator of Call objects. Keep track of the last phase set seen,
    so that we can assign unphased variants to them to ensure contiguity of phase blocks.
    """
    current_phase_set = 0
    for record in variant_iter:
        phase_set = tk_io.get_record_phase_set(record)
        if phase_set is not None:
            current_phase_set = phase_set
        yield Call(current_phase_set, record)

def mk_var_tuples(fn, targets, locus):
    reader = tk_io.VariantFileReader(fn)
    var_iter = tk_io.get_variant_iterator_pos(reader, targets, locus)
    var_tuples = variant_tuples(var_iter)
    return var_tuples

def is_unfiltered(filter_field):
    return (filter_field is None or (filter_field is not None and (len(filter_field) == 1 and filter_field[0] == "PASS")) or filter_field == [])

def arg_min(vals, k):
    min_key = None
    min_val = None
    min_idx = -1

    for (idx, v) in enumerate(vals):
        if v is None:
            continue
        key = k(v)
        if min_idx == -1:
            min_val = v
            min_key = key
            min_idx = idx
        elif key < min_key:
            min_val = v
            min_key = key
            min_idx = idx

    return (min_idx, min_val, min_key)

def mark_mins(vals, k):
    (min_idx, min_val, min_key) = arg_min(vals, k)
    at_min = [v is not None and k(v) == min_key for v in vals]
    return at_min

def parallel_iter(_iters, key):
    ''' Iterate through sorted i1 and i2 simulateneously. Emit (v1, v2) for items
        that compare equal. Emit (v1, None) for value v1 from i1 that doesn't have
        an equal item in i2. Emit (None, v2) for value v2 from i2 that doesn't have
        and equal items in i1.  Comparisons are made with key(v).  Assumes
        that i1 and i2 are sorted with respect to key(). Also assumes that
        a sigle iterator does not contain items equal to each other.'''

    iters = [iter(x) for x in _iters]
    n = len(iters)
    cache_vals = [None] * n
    it_done = [False] * n

    for i in xrange(n):
        try:
            cache_vals[i] = iters[i].next()
        except StopIteration:
            it_done[i] = True

    while any(cache_vals):
        # Find current values equal to min key
        mins = mark_mins(cache_vals, key)

        # yield set at min
        yield [v if at_min else None for (v, at_min) in zip(cache_vals, mins)]

        # get next value from emitted iterators
        for (i, at_min) in enumerate(mins):
            if at_min:
                cache_vals[i] = None
                if not it_done[i]:
                    try:
                        cache_vals[i] = iters[i].next()
                    except StopIteration:
                        it_done[i] = True
