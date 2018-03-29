#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for reading/writing/manipulating BEDPE files from SV-calling.
#

import numpy as np
import pandas as pd
import martian
import gzip
import vcf.model
import vcf.parser
from collections import OrderedDict
import tenkit.reference as tk_reference
import tenkit.bio_io as tk_io
import cStringIO
import subprocess
from itertools import repeat

BEDPE_HEADER = '''# Fields are documented in the BEDPE section at http://support.10xgenomics.com
'''

def result_list_to_df(res):
    chroms1 = [r.chrom1 for r in res]
    chroms2 = [r.chrom2 for r in res]
    starts1 = [r.break1[0] for r in res]
    stops1 = [r.break1[1] for r in res]
    starts2 = [r.break2[0] for r in res]
    stops2 = [r.break2[1] for r in res]
    quals = [r.lr for r in res]
    infos = [r.info for r in res]
    names = np.arange(len(res))
    return create_sv_df(chroms1, starts1, stops1,
                        chroms2, starts2, stops2,
                        names, quals, infos)

def create_sv_df(chroms1, starts1, stops1, chroms2, starts2, stops2,
                 names, quals, info_strs='.', filter_strs='.'):
    out_df = pd.DataFrame({'chrom1':chroms1, 'start1':starts1, 'stop1':stops1,
                           'chrom2':chroms2, 'start2':starts2, 'stop2':stops2,
                           'name':names, 'qual':quals,
                           'strand1':'.', 'strand2':'.',
                           'filters':filter_strs, 'info':info_strs},
                          index=np.arange(len(names)))
    return out_df


def drop_info(calls, info_list):
    out_infos = []
    for _, row in calls.iterrows():
        out_infos.append(update_info(row.info, info_list, repeat(None, len(info_list))))
    calls['info'] = out_infos
    return calls


def update_info(filters, new_fields, new_values):
    '''None will remove the info field'''

    filter_map = {}
    filters = filters.strip(';').split(';')
    for f in filters:
        fields = f.split('=')
        if len(fields) == 2:
            filter_map[fields[0]] = fields[1]
    for f, v in zip(new_fields, new_values):
        if v == '':
            filter_map[f] = '.'
        else:
            filter_map[f] = v
    if len(filter_map) == 0:
        info = '.'
    else:
        info = ';'.join([str(f) + '=' + str(filter_map[f]) for f in sorted(filter_map.keys()) if not filter_map[f] is None])
    return info


def update_filters(filters, name, value, min_val=-np.inf, max_val=np.inf):
    if min_val is None:
        min_val = -np.inf
    if max_val is None:
        max_val = np.inf
    filter_set = set([f for f in set(filters.strip('.;').split(';')) if len(f.strip()) > 0])
    if value is None or value < min_val or value > max_val:
        filter_set.add(name)
    if len(filter_set) == 0:
        new_filters = '.'
    else:
        new_filters = ';'.join(sorted(filter_set)).strip(';')
    return new_filters


def extract_sv_info(filter_str, filter_names):
    filter_fields = filter_str.split(';')
    filter_vals = ['' for n in filter_names]
    for f in filter_fields:
        sub_filter_fields = f.split('=')
        if len(sub_filter_fields) < 2:
            continue
        for idx, name in enumerate(filter_names):
            if name == sub_filter_fields[0]:
                filter_vals[idx] = sub_filter_fields[1]
                if filter_vals[idx] == 'None':
                    filter_vals[idx] = None
    return filter_vals


def get_sv_zygosity(info_str):
    zs = extract_sv_info(info_str, ['ZS'])[0]
    if zs is None or len(zs) == 0:
        return '.'
    return zs

def extract_sv_hap(info_str):
    haps = extract_sv_info(info_str, ['HAPS'])[0]
    if haps is None or len(haps) == 0:
        return ('.', '.')

    return haps

def get_sv_type(info_str):
    val = extract_sv_info(info_str, ['TYPE'])[0]
    if len(val) == 0:
        return 'UNK'
    return val

def get_npairs(info):
    val = extract_sv_info(info, ['NPAIRS'])[0]
    if len(val) == 0:
        return 0
    return int(val)

def get_nsplit(info):
    val = extract_sv_info(info, ['NSPLIT'])[0]
    if len(val) == 0:
        return 0
    return int(val)

def get_phase_set(info, break_num=1):
    val = extract_sv_info(info, ['PS' + str(break_num)])[0]
    if val is None or len(val) == 0:
        return '.'
    return val

def get_break_orientation(info):
    # orientation is a string +-, --, -+, or ++
    # by convention, we use .. to denote missing or unknown
    val = extract_sv_info(info, ['ORIENT'])[0]
    if len(val) == 0:
        return '..'
    return val

def get_break_orientation_arr(info):
    # orientation is a string +-, --, -+, or ++
    # by convention, we use .. to denote missing or unknown
    val = extract_sv_info(info, ['ORIENT'])[0]
    arr = np.zeros((2, ), dtype=np.int)
    if len(val) == 0:
        return arr
    arr[0] = 1 if val[0] == '+' else -1
    arr[1] = 1 if val[1] == '+' else -1
    return arr

def get_tier(info):
    val = extract_sv_info(info, ['TIER'])[0]
    if len(val) == 0:
        return 0
    return int(val)

def get_nbcs1(info):
    val = extract_sv_info(info, ['NBCS1'])[0]
    if len(val) == 0:
        return None
    return int(val)

def get_nbcs2(info):
    val = extract_sv_info(info, ['NBCS2'])[0]
    if len(val) == 0:
        return None
    return int(val)

def get_noov(info):
    val = extract_sv_info(info, ['NOOV'])[0]
    if len(val) == 0:
        return None
    return int(val)

def get_sv_df_dists(sv_df):
    dist_arr = np.array(sv_df['start2'] - sv_df['stop1'], dtype=np.float)
    dist_arr[np.array(sv_df['chrom1'] != sv_df['chrom2'])] = np.inf
    return dist_arr


def check_sv_names(bedpe_df):
    if bedpe_df.shape[1] < 7:
        martian.throw('Ground truth BEDPE must have a name column.')
    names = set(bedpe_df['name'])
    if len(names) != len(bedpe_df):
        martian.throw('Names in ground truth BEDPE must be unique.')


def read_sv_bedpe_to_df(bedpe):
    col_names = ['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2',
                 'name', 'qual', 'strand1', 'strand2', 'filters', 'info']
    if bedpe is None:
        return pd.DataFrame(columns=col_names)

    try:
        df = pd.read_csv(bedpe, sep='\t', header=None,
                         index_col=None, comment='#')
    except ValueError:
        return pd.DataFrame(columns=col_names)
    if df.shape[1] < 6:
        martian.throw('BEDPE file must have at least 6 columns')
    ncols = min(len(col_names), df.shape[1])
    df = df.iloc[:, 0:ncols]
    df.columns = col_names[0:ncols]
    df[['chrom1', 'chrom2']] = df[['chrom1', 'chrom2']].astype(object)
    df[['start1', 'stop1', 'start2', 'stop2']] = df[['start1', 'stop1', 'start2', 'stop2']].astype(int)
    if not 'name' in df.columns:
        df['name'] = np.arange((len(df),))
    if not 'qual' in df.columns:
        df['qual'] = 1
    if not 'strand1' in df.columns:
        df['strand1'] = '.'
    if not 'strand2' in df.columns:
        df['strand2'] = '.'
    if not 'filters' in df.columns:
        df['filters'] = '.'
    if not 'info' in df.columns:
        df['info'] = '.'

    return df


def write_sv_df_to_bedpe(in_df, out_bedpe):
    col_names = ['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2',
                 'name', 'qual', 'strand1', 'strand2', 'filters', 'info']
    if in_df is None:
        out_df = pd.DataFrame(columns=col_names)
    else:
        if not 'name' in in_df.columns:
            in_df['name'] = [str(i) for i in range(len(in_df))]
        if not 'qual' in in_df.columns:
            in_df['qual'] = 0
        for f in ['strand1', 'strand2', 'filters', 'info']:
            if not f in in_df.columns:
                in_df[f] = '.'

        out_df = in_df[col_names]
        out_df[['start1', 'stop1', 'start2', 'stop2', 'qual']] = out_df[['start1', 'stop1', 'start2', 'stop2', 'qual']].astype(int)
        out_df = out_df.sort(['qual', 'chrom1', 'start1', 'stop1',
                              'chrom2', 'start2', 'stop2'],
                             ascending=[0, 1, 1, 1, 1, 1, 1])

    with open(out_bedpe, 'w') as out_fn:
        out_fn.write(BEDPE_HEADER + '#')
        out_df.to_csv(out_fn, sep='\t', header=True, index=False)


def switch_sv_breaks(sv_df, sel_cols):
    """Switches break1 and break2 in the sv dataframe.
    DOES NOT switch any info fields."""

    chroms1 = np.array([s for s in sv_df.chrom1])
    starts1 = np.array(sv_df.start1)
    stops1 = np.array(sv_df.stop1)
    chroms2 = np.array([s for s in sv_df.chrom2])
    starts2 = np.array(sv_df.start2)
    stops2 = np.array(sv_df.stop2)

    chroms1[sel_cols] = chroms2[sel_cols]
    chroms2[sel_cols] = np.array([s for s in sv_df.chrom1])[sel_cols]
    starts1[sel_cols] = starts2[sel_cols]
    starts2[sel_cols] = np.array(sv_df.start1)[sel_cols]
    stops1[sel_cols] = stops2[sel_cols]
    stops2[sel_cols] = np.array(sv_df.stop1)[sel_cols]
    new_df = sv_df.copy()

    new_df.chrom1 = chroms1
    new_df.chrom2 = chroms2
    new_df.start1 = starts1
    new_df.start2 = starts2
    new_df.stop1 = stops1
    new_df.stop2 = stops2
    return new_df


def extract_sv_filters(filter_str):
    """Converts the FILTERS field of a BEDPE row into a list of filters"""
    filters = [filt for filt in filter_str.strip('.;').split(';') if filt]
    return filters


def get_vcf_gt(bedpe_hap, bedpe_zs):
    """Converts the HAPS info field into a VCF commpatible genotype call"""
    if bedpe_zs == 'HOM':
        return '1|1'
    elif bedpe_zs == '.':
        return './.'
    elif bedpe_hap == '.':
        return '0/1'
    return '0|1' if bedpe_hap == '1' else '1|0'


class BedpeInfoField(object):
    """A BEDPE INFO field"""
    def __init__(self, name, bedpe_name, field_type, desc):
        self.name = name
        self.bedpe_name = bedpe_name
        self.field_type = field_type
        self.description = desc

    def vcf_repr(self):
        return vcf.parser._Info(self.name, 1, self.field_type, self.description, None, None)

    def extract_from_bedpe_row(self, row):
        val = extract_sv_info(row, [self.bedpe_name])[0]
        if val is None or val == '':
            return '.'
        return val


class BedpeFilterField(object):
    """A BEDPE filter field"""
    def __init__(self, name, bedpe_name, desc):
        self.name = name
        self.bedpe_name = bedpe_name
        self.description = desc

    def vcf_repr(self):
        return vcf.parser._Filter(self.name, self.description)


def write_sv_vcf_template(out_vcf, reference_path, sample_name):
    """Writes a basic header for a VCF containing SVs to out_vcf.

    out_vcf should be an open stream.
    """
    out_vcf.write('##fileformat=VCFv4.2\n')
    out_vcf.write('##reference=' + reference_path + '\n')
    out_vcf.write('##phasing=none\n')
    out_vcf.write('##ALT=<ID=UNK,Description="SV of unknown type">\n')
    out_vcf.write('##ALT=<ID=DEL,Description="Deletion">\n')
    out_vcf.write('##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">\n')
    out_vcf.write('##ALT=<ID=INV,Description="Inversion">\n')
    out_vcf.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                             'FILTER', 'INFO', 'FORMAT', sample_name]) + '\n')


def create_breakend_records(name, chrom1, pos, cipos, chrom2, end_pos, ciend,
                            qual, bedpe_info_str, out_filters,
                            format_field_defs, sample_idx, sample_name, fasta):
    """Get a list of VCF breakend records from a BEDPE row"""

    data_instantiator = vcf.model.make_calldata_tuple([field.id for field in format_field_defs])
    format_str = ':'.join(field.id for field in format_field_defs)

    sv_type = get_sv_type(bedpe_info_str)
    orient_str = get_break_orientation(bedpe_info_str)

    ps1 = get_phase_set(bedpe_info_str, 1)
    ps2 = get_phase_set(bedpe_info_str, 2)
    haps = extract_sv_hap(bedpe_info_str)
    zs = extract_sv_info(bedpe_info_str, ['ZS'])[0]

    uncertain_orient = False
    if sv_type == 'DEL':
        orient = (False, True)
    elif sv_type == 'INV':
        orient = (False, False)
    elif sv_type == 'DUP':
        orient = (True, False)
    else:
        sv_type = 'BND'
        if orient_str == '..':
            orient = (True, True)
            uncertain_orient = True
        else:
            orient = (True if orient_str[0] == '+' else False,
                      True if orient_str[1] == '+' else False)

    ref1 = get_fasta_base(fasta, chrom1, pos)
    ref2 = get_fasta_base(fasta, chrom2, end_pos)

    alt1 = [vcf.model._Breakend(chrom2, end_pos + 1, orient[0], orient[1], ref1, True)]
    alt2 = [vcf.model._Breakend(chrom1, pos + 1, orient[1], orient[0], ref2, True)]

    info1 = OrderedDict()
    info1['SVTYPE'] = 'BND' # Has to be BND for breakend records...
    info1['SVTYPE2'] = sv_type
    info1['MATEID'] = name + '_2'
    info1['CIPOS'] = cipos
    info1['PS'] = ps1
    info1['EVENT'] = name

    info2 = OrderedDict()
    info2['SVTYPE'] = 'BND'
    info2['SVTYPE2'] = sv_type
    info2['MATEID'] = name + '_1'
    info2['CIPOS'] = ciend
    info2['PS'] = ps2
    info2['EVENT'] = name

    if uncertain_orient:
        info1['IMPRECISE_DIR'] = True
        info2['IMPRECISE_DIR'] = True

    gt1 = get_vcf_gt(haps[0], zs)
    gt2 = get_vcf_gt(haps[1], zs)

    record1 = vcf.model._Record(chrom1, pos + 1, name + '_1',
                                ref1, alt1, qual, out_filters,
                                info1, format_str, sample_idx,
                                samples=[vcf.model._Call(None, sample_name, data_instantiator(gt1))])
    record2 = vcf.model._Record(chrom2, end_pos + 1, name + '_2',
                                ref2, alt2, qual, out_filters,
                                info2, format_str, sample_idx,
                                samples=[vcf.model._Call(None, sample_name, data_instantiator(gt2))])

    if sv_type == 'INV':
        # For inversions described as breakends we need to add two extra events.
        pos += 1
        end_pos += 1

        ref3 = get_fasta_base(fasta, chrom1, max(0, pos - 1))
        ref4 = get_fasta_base(fasta, chrom2, max(0, end_pos - 1))

        alt3 = [vcf.model._Breakend(chrom2, end_pos + 1, True, True, ref3, True)]
        alt4 = [vcf.model._Breakend(chrom1, pos + 1, True, True, ref4, True)]

        info3 = info1.copy()
        info3['MATEID'] = name + '_4'
        info3['CIPOS'] = cipos

        info4 = info2.copy()
        info4['MATEID'] = name + '_3'
        info4['CIPOS'] = ciend

        record3 = vcf.model._Record(chrom1, pos + 1, name + '_3',
                                    ref3, alt3, qual, out_filters,
                                    info3, format_str, sample_idx,
                                    samples=[vcf.model._Call(None, sample_name, data_instantiator(gt1))])
        record4 = vcf.model._Record(chrom2, end_pos + 1, name + '_4',
                                    ref4, alt4, qual, out_filters,
                                    info4, format_str, sample_idx,
                                    samples=[vcf.model._Call(None, sample_name, data_instantiator(gt2))])

        return [record1, record2, record3, record4]

    return [record1, record2]


def get_pos_confidence_intervals(start, stop):
    middle = (start + stop) / 2
    conf_interval = [start - middle, stop - middle]
    return middle, conf_interval


def get_fasta_base(fasta, chrom, pos):
    if str(chrom) in fasta and 0 <= pos < len(fasta[str(chrom)]):
        return fasta[str(chrom)][pos].upper()
    return '.'


def bedpe_to_vcf(bedpe_file, vcf_file, sample_name, source_str, reference_path):
    """Converts a BEDPE file to VCF
    source_str: the string to add in the VCF header as the new source
    """

    # INFO field definitions that are VCF-specific (not from BEDPE)
    info_field_defs = []
    info_field_defs.append(vcf.parser._Info('END', 1, 'Integer', 'Confidence interval around END for imprecise variants', None, None))
    info_field_defs.append(vcf.parser._Info('CIEND', 2, 'Integer', 'Confidence interval around END for imprecise variants', None, None))
    info_field_defs.append(vcf.parser._Info('CIPOS', 2, 'Integer', 'Confidence interval around POS for imprecise variants', None, None))
    info_field_defs.append(vcf.parser._Info('IMPRECISE_DIR', 0, 'Flag', 'Imprecise orientation of breakend', None, None))
    info_field_defs.append(vcf.parser._Info('POS', 1, 'Integer', 'Start position of the variant described in this record', None, None))
    info_field_defs.append(vcf.parser._Info('SVTYPE', 1, 'String', 'Type of structural variant', None, None))
    info_field_defs.append(vcf.parser._Info('SVTYPE2', 1, 'String', 'Type of structural variant implied by breakends', None, None))
    info_field_defs.append(vcf.parser._Info('SVLEN', 1, 'Integer', 'Difference in length between REF and ALT alleles', None, None))
    # We won't just copy the phase set from the BEDPE because each BEDPE row
    # has PS1 and PS2, the phase sets associated with both ends.
    info_field_defs.append(vcf.parser._Info('PS', 1, 'Integer', 'Phase set assignment', None, None))
    info_field_defs.append(vcf.parser._Info('MATEID', '.', 'String', 'ID of mate breakends', None, None))
    info_field_defs.append(vcf.parser._Info('EVENT', 1, 'String', 'ID of event associated to breakend', None, None))

    # More INFO field definitions that will be just copied from the BEDPE
    bedpe_fields = []
    bedpe_fields.append(BedpeInfoField('HAP_ALLELIC_FRAC', 'HAP_ALLELIC_FRAC',
                                       'Float', 'Haplotype allelic fraction'))
    bedpe_fields.append(BedpeInfoField('ALLELIC_FRAC', 'ALLELIC_FRAC',
                                       'Float', 'Allelic fraction'))
    bedpe_fields.append(BedpeInfoField('PAIRS', 'NPAIRS',
                                       'Integer', 'Number of supporting read pairs'))
    bedpe_fields.append(BedpeInfoField('SPLIT', 'NSPLIT',
                                       'Integer', 'Number of supporting split reads'))
    bedpe_fields.append(BedpeInfoField('WildCov', 'WildCov',
                                       'Float', 'Number of barcodes in the wild haplotype'))
    bedpe_fields.append(BedpeInfoField('MolTtl', 'MolTtl',
                                       'Integer', 'Total number of barcodes in the region of the deletion'))
    bedpe_fields.append(BedpeInfoField('MolTtlNoR', 'MolTtlNoR',
                                       'Integer', 'Number of barcodes without covering reads in the region of the deletion'))
    bedpe_fields.append(BedpeInfoField('MolDel', 'MolDel',
                                       'Integer', 'Total number of barcodes in deleted haplotype'))
    bedpe_fields.append(BedpeInfoField('MolDelNoR', 'MolDelNoR',
                                       'Integer', 'Number of barcodes without covering reads in deleted haplotype'))
    bedpe_fields.append(BedpeInfoField('MolWild', 'MolWild',
                                       'Integer', 'Total number of barcodes in wild haplotype'))
    bedpe_fields.append(BedpeInfoField('MolWildNoR', 'MolWildNoR',
                                       'Integer', 'Number of barcodes without covering reads in wild haplotype'))
    bedpe_fields.append(BedpeInfoField('PVAL', 'PVAL',
                                       'Float', 'p-value for the called event'))

    bedpe_fields.append(BedpeInfoField('BGSize', 'BGSize',
                                       'Integer', 'Size of the background region'))
    bedpe_fields.append(BedpeInfoField('BGImparityPval', 'BGImparityPval',
                                       'Float', 'The p-value for the hypothesis that the two haplotypes in the background region have different barcode coverage'))
    bedpe_fields.append(BedpeInfoField('BGTtlRCnt', 'BGTtlRCnt',
                                       'Integer', 'Total number of barcodes in the background region'))
    bedpe_fields.append(BedpeInfoField('BGHP1RCnt', 'BGHP1RCnt',
                                       'Integer', 'Number of barcodes on haplotype 1 in the background region'))
    bedpe_fields.append(BedpeInfoField('BGHP2RCnt', 'BGHP2RCnt',
                                       'Integer', 'Number of barcodes on haplotype 2 in the background region'))
    bedpe_fields.append(BedpeInfoField('BGBaseCov', 'BGBaseCov',
                                       'Float', 'Base level coverage in the background region'))
    bedpe_fields.append(BedpeInfoField('BGPhaseFrac', 'BGPhaseFrac',
                                       'Float', 'Fraction of reads phased'))
    bedpe_fields.append(BedpeInfoField('NRead', 'NRead',
                                       'Float', 'Number of barcodes covering the homo deletion region'))
    bedpe_fields.append(BedpeInfoField('SOURCE', 'SOURCE',
                                       'String', 'Method from which the event is called'))


    info_field_defs.extend(bf.vcf_repr() for bf in bedpe_fields)

    # FORMAT field definitions
    format_field_defs = []
    format_field_defs.append(vcf.parser._Format('GT', 1, 'String', 'Genotype'))

    # FILTER definitions copied from the BEDPE.
    bedpe_filters = []
    bedpe_filters.append(BedpeFilterField('BLACK_DIST', 'BLACK_DIST', 'Too close to the SV blacklist'))
    bedpe_filters.append(BedpeFilterField('BLACK_FRAC', 'BLACK_FRAC', 'Too much overlap with the SV blacklist'))
    bedpe_filters.append(BedpeFilterField('SEG_DUP', 'SEG_DUP', 'Breakpoints overlap copies of the same segmental duplication'))
    bedpe_filters.append(BedpeFilterField('LOWQ', 'LOWQ', 'Low confidence call'))
    bedpe_filters.append(BedpeFilterField('LARGE_PVALUE', 'LARGE_PVALUE', 'Insignificant pvalue'))
    bedpe_filters.append(BedpeFilterField('LOW_WILD_COV', 'LOW_WILD_COV',
                                          'Too few reads on the wild type haplotype and hence an unliable pvalue calculation'))
    bedpe_filters.append(BedpeFilterField('IMPARITY_IN_HAP_COV', 'IMPARITY_IN_HAP_COV',
                                          'In the background region, one haplotype has significanlty more reads than the other haplotype'))
    bedpe_filters.append(BedpeFilterField('LOW_COV_REGION', 'LOW_COV_REGION',
                                          'The background region has significantly low barcode coverage'))
    bedpe_filters.append(BedpeFilterField('LOW_PHASING_REGION', 'LOW_PHASING_REGION',
                                          'Low fraction of reads phased in the background region'))
    bedpe_filters.append(BedpeFilterField('CNV_SEG_DUP', 'SEGDUP',
                                          'Called events overlap with segdup tracks with sequence identity is 98% or higher'))
    bedpe_filters.append(BedpeFilterField('LowCov', 'LowCov',
                                          'Called events overlap with low exome target regions'))
    bedpe_filters.append(BedpeFilterField('HighContamination', 'HighContamination',
                                          'Too many potential contaminated barcodes in a supposedly homo deletion region'))

    filter_defs = [bf.vcf_repr() for bf in bedpe_filters]

    fasta_path = tk_reference.get_fasta(reference_path)

    sample_idx = {sample_name: 0}

    # Write a template VCF into a stream.
    # We need this template to feed into the vcf Writer.
    tmp_template_stream = cStringIO.StringIO()
    write_sv_vcf_template(tmp_template_stream, fasta_path, sample_name)
    # Not sure why, but trying to just read from the first stream won't work.
    template_stream = cStringIO.StringIO(tmp_template_stream.getvalue())
    tmp_template_stream.close()

    fasta = tk_reference.open_reference(reference_path)
    data_instantiator = vcf.model.make_calldata_tuple([field.id for field in format_field_defs])
    format_str = ':'.join(field.id for field in format_field_defs)

    bedpe_df = read_sv_bedpe_to_df(bedpe_file)

    compressed = vcf_file.endswith('.gz')
    out_fn = gzip.open(vcf_file, 'wb') if compressed else open(vcf_file, 'w')

    vcf_writer = tk_io.VariantFileWriter(out_fn, template_file=template_stream,
                                         new_info_fields=info_field_defs,
                                         new_format_fields=format_field_defs,
                                         new_filters=filter_defs,
                                         new_source=source_str)
    template_stream.close()
    all_records = []
    for _, row in bedpe_df.iterrows():
        sv_type = get_sv_type(row.info)
        if sv_type is None:
            sv_type = '.'

        # The returned positions will be 0-based
        pos, cipos = get_pos_confidence_intervals(row.start1, row.stop1 - 1)
        end_pos, ciend = get_pos_confidence_intervals(row.start2, row.stop2 -1)

        sv_len = end_pos - pos
        if sv_type == 'DEL':
            sv_len = -sv_len

        info_fields = OrderedDict()
        info_fields['SVTYPE'] = sv_type

        out_filters = extract_sv_filters(row.filters)

        ps1 = get_phase_set(row.info, 1)
        ps2 = get_phase_set(row.info, 2)
        if not ps1:
            ps1 = '.'
        if not ps2:
            ps2 = '.'
        haps = extract_sv_hap(row.info)
        zs = extract_sv_info(row.info, ['ZS'])[0]

        if sv_type != 'DISTAL' and row.chrom1 == row.chrom2 and ps1 == ps2:
            if sv_type == 'DUP':
                alt = [vcf.model._SV('DUP:TANDEM')]
            else:
                alt = [vcf.model._SV(sv_type)]
            info_fields['SVLEN'] = sv_len
            info_fields['CIPOS'] = cipos
            info_fields['CIEND'] = ciend
            info_fields['END'] = end_pos
            info_fields['PS'] = ps1

            gt = get_vcf_gt(haps[0], zs)

            samples = [vcf.model._Call(None, sample_name, data_instantiator(gt))]

            # According to the specification:
            # "Strings must include the base before the event (which must be reflected
            # in the POS field), unless the event occurs at position 1 on the contig in
            # which case it must include the base after the event".
            #
            # the event starts at position pos (0-based) so the position before the 
            # event (1-based) is pos. We'll write pos in the VCF (which is 1-based), and
            # fetch the (0-based) position pos - 1
            ref = get_fasta_base(fasta, row.chrom1, max(0, pos - 1))

            records = [vcf.model._Record(row.chrom1, pos, row['name'],
                                         ref, alt, row['qual'], out_filters,
                                         info_fields, format_str,
                                         sample_idx, samples=samples)]
        else:
            records = create_breakend_records(str(row['name']), row.chrom1, pos, cipos,
                                              row.chrom2, end_pos, ciend, row['qual'],
                                              row['info'], out_filters, format_field_defs,
                                              sample_idx, sample_name, fasta)

        # All records generated from the same row of the BEDPE will share
        # the same values of the BEDPE info fields.
        for record in records:
            for field in bedpe_fields:
                record.INFO[field.name] = field.extract_from_bedpe_row(row['info'])
            all_records.append(record)


    all_records.sort(key=lambda x: (x.CHROM, x.POS))
    for r in all_records:
        vcf_writer.write_record(r)

    out_fn.close()


def index_sv_vcf(vcf):
    ''' bgzip, and tabix a VCF.'''
    subprocess.check_call("bgzip -c {0} > {0}.gz".format(vcf), shell=True)
    subprocess.check_call("tabix -p vcf {0}.gz".format(vcf), shell=True)


