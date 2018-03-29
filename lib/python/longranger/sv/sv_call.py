#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for computing fragment overlaps and calling large scale structural variants.
#
import numpy as np
from tenkit.stats import robust_divide
import longranger.sv.io as tk_sv_io
from longranger.sv.constants import MIN_CORRECT_HAP_FRAC
from collections import namedtuple
from enum import Enum

DOWNSTREAM = '+'
UPSTREAM = '-'

class Zygosity(Enum):
    hom = 'HOM'
    het = 'HET'
    unknown = '.'

class Haplotype(Enum):
    hap0 = '0'
    hap1 = '1'
    unknown = '.'

_SvType = namedtuple('SvType', ['name', 'orientation'])
UNK = _SvType('UNK', ('.', '.'))
DEL = _SvType('DEL', ('.', '.'))
INV = _SvType('INV', ('.', '.'))
DUP = _SvType('DUP', ('.', '.'))
INV_DUP_OUT = _SvType('INVDUPOUT', ('.', '.'))
INV_DUP_IN = _SvType('INVDUPIN', ('.', '.'))
TRANS_00 = _SvType('DISTAL', ('+', '+'))
TRANS_ee = _SvType('DISTAL', ('-', '-'))
TRANS_0e = _SvType('DISTAL', ('+', '-'))
TRANS_e0 = _SvType('DISTAL', ('-', '+'))

PROXIMAL_SV_TYPES = [DEL, INV, DUP]
DISTAL_SV_TYPES = [TRANS_ee, TRANS_00, TRANS_e0, TRANS_0e]


class SvCall(object):
    def __init__(self, sv_type, break1=(None, None), break2=(None, None),
                 chrom1='', chrom2='', qual=0, name=None,
                 hap1=None, hap2=None, zygosity=Zygosity.unknown, info=None):
        self.sv_type = sv_type
        if not isinstance(break1, tuple):
            assert isinstance(break1, int)
            self.break1 = (break1, break1 + 1)
        else:
            self.break1 = break1
        if not isinstance(break2, tuple):
            assert isinstance(break2, int)
            self.break2 = (break2, break2 + 1)
        else:
            self.break2 = break2
        self.qual = qual
        self.hap1 = hap1
        self.hap2 = hap2
        self.zygosity = zygosity
        self.chrom1 = chrom1
        self.chrom2 = chrom2
        self.name = name
        self.info = info
        self.add_info('ORIENT', sv_type.orientation[0] + sv_type.orientation[1])


    @classmethod
    def from_em_results(cls, c1, c2, ps1, ps2, max_lrs,
                        max_locus, sv_type, zygosity, max_hap, support,
                        posterior_probs):

        if zygosity == Zygosity.unknown or zygosity == Zygosity.hom:
            max_hap = (Haplotype.unknown, Haplotype.unknown)

        no_sv_max, sv_max, het_sv_max = max_lrs
        lr = sv_max - no_sv_max if max_hap is None else het_sv_max - no_sv_max

        hap_probs1, hap_probs2, _ = posterior_probs

        if sv_max >= het_sv_max:
            max_hap = (Haplotype.unknown, Haplotype.unknown)
            sv_support = int(np.sum(support > 0))
            hap_allelic_frac = robust_divide(np.sum(support > 0),
                                             np.sum(support > 0) + np.sum(support < 0))
            correct_hap_frac = None
        else:
            hap_probs_sv = np.minimum(hap_probs1[:, int(max_hap[0].value)],
                                      hap_probs2[:, int(max_hap[1].value)])
            
            # Support on the best haplotype (weighted by phasing probabilities)
            sv_support = int(np.ceil(np.sum(hap_probs_sv[support > 0])))

            sv_no_support = int(np.ceil(np.sum(hap_probs_sv[support < 0])))
            # Fraction of assigned haplotype that supports the event
            hap_allelic_frac = robust_divide(sv_support, sv_support + sv_no_support)

            # Fraction of support coming from assigned haplotypes (weighted by
            # phasing probabilities)
            correct_hap_frac = robust_divide(sv_support, int(np.sum(support > 0)))
            if correct_hap_frac < MIN_CORRECT_HAP_FRAC:
                max_hap = (Haplotype.unknown, Haplotype.unknown)
                
        break1 = (max_locus[0], max_locus[0] + 1)
        break2 = (max_locus[1], max_locus[1] + 1)

        info = {}
        info['PS1'] = ps1
        info['PS2'] = ps2
        info['HAP_ALLELIC_FRAC'] = hap_allelic_frac
        info['FRAC_HAP_SUPPORT'] = correct_hap_frac
        info['ALLELIC_FRAC'] = np.mean(support > 0)
        info['LR'] = lr

        return cls(sv_type, break1=break1, break2=break2,
                   chrom1=c1, chrom2=c2, qual=sv_support, zygosity=zygosity,
                   hap1=max_hap[0], hap2=max_hap[1], info=info)

    @staticmethod
    def svs_to_dataframe(sv_call_list):
        quals = []
        chroms1 = []
        starts1 = []
        stops1 = []
        chroms2 = []
        starts2 = []
        stops2 = []
        infos = []
        names = []

        for idx, sv_call in enumerate(sv_call_list):
            if sv_call.name is None:
                names.append(str(idx))
            else:
                names.append(sv_call.name)
            quals.append(sv_call.qual)
            chroms1.append(sv_call.chrom1)
            chroms2.append(sv_call.chrom2)
            starts1.append(sv_call.break1[0])
            stops1.append(sv_call.break1[1] + 1)
            starts2.append(sv_call.break2[0])
            stops2.append(sv_call.break2[1] + 1)

            info_keys = sv_call.get_info_keys()
            info_vals = [sv_call.get_info(key) for key in info_keys]
            info_keys.append('TYPE')
            info_vals.append(sv_call.sv_type.name)
            info_keys.append('HAPS')
            info_vals.append(sv_call.hap1.value + ',' + sv_call.hap2.value)
            info_keys.append('ZS')
            info_vals.append(sv_call.zygosity.value)
            infos.append(tk_sv_io.update_info('.', info_keys, info_vals))

        out_df = tk_sv_io.create_sv_df(chroms1, starts1, stops1,
                                       chroms2, starts2, stops2,
                                       names, quals, info_strs=infos)
        return out_df


    def hap_prob_str(self):
        return ','.join(['{:.3f}'.format(p) for p in self.hap_probs])


    def __str__(self):
        format_str = '{} ({}:{}, {}:{}), LR = {}'
        values = [self.sv_type, self.chrom1, self.break1[0],
                  self.chrom2, self.break2[0], self.qual]
        if not self.hap1 is None and not self.hap2 is None:
            format_str = format_str + ' ({}, {})'
            values.extend([self.hap1, self.hap2])
        return format_str.format(*values)


    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.qual < other.qual

    def __le__(self, other):
        return self.qual <= other.qual

    def __gt__(self, other):
        return self.qual > other.qual

    def __ge__(self, other):
        return self.qual >= other.qual

    def add_info(self, key, value):
        self.info[key] = value

    def get_info_keys(self):
        return self.info.keys()

    def get_info(self, key):
        return self.info[key]
