#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Filter variants based on whether the reads' barcodes supporting alleles are consistent with
#  the barcodes phasing information

import tenkit
import tenkit.constants
import tenkit.bio_io as tk_io
import pysam
import tenkit.bam as tk_bam
import tenkit.tabix as tk_tabix
import os.path
import math
import tenkit.chunk_utils as tk_chunks

__MRO__ = """
stage BARCODE_AWARE_FILTER_SNPINDELS(
    in  bam    bam_file,
    in  string sex,
    in  string restrict_locus,
    in  string vc_precalled,
    in  string vc_mode,
    in  tsv.gz fragment_phasing,
    in  vcf.gz variants,
    out vcf.gz,
    src py     "stages/snpindels/barcode_aware_filter_snpindels",
) split using (
    in  string locus,
)
"""

def split(args):
    input_bam = tk_bam.create_bam_infile(args.bam_file)
    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    if args.restrict_locus is None:
        loci = tk_chunks.chunk_by_locus(chroms, chrom_lengths, tenkit.constants.PARALLEL_LOCUS_SIZE)
    else:
        loci = [{'locus': args.restrict_locus}]

    return {'chunks': loci}

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    tk_io.combine_vcfs(outs.default.strip(".gz"), [c.default.strip(".gz") for c in chunk_outs if os.path.isfile(c.default.strip(".gz"))])

def main(args, outs):
    vc_mode, _, _, _ = tk_io.get_vc_mode(args.vc_precalled, args.vc_mode)

    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    chrom = str(chrom)

    if chrom in ['chrM', 'MT', 'M'] or (args.sex.lower() in ["f", "female"] and chrom in ["chrY", "Y"]):
        return

    fragment_barcode_info = pysam.Tabixfile(args.fragment_phasing)
    AH_0_BH_0 = ('AH_0_BH_0','1','Integer','Number of barcodes that have been called as supporting haplotype 0 which are on reads that have support for the allele which has been phased as haplotype 0')
    AH_1_BH_1 = ('AH_1_BH_1','1', 'Integer', 'Number of barcodes that have been called as supporting haplotype 1 which are on reads that have support for the allele which has been phased as haplotype 1')
    AH_0_BH_1 = ('AH_0_BH_1','1', 'Integer', 'Number of barcodes that have been called as supporting haplotype 0 which are on reads that have support for the allele which has been phased as haplotype 1')
    AH_1_BH_0 = ('AH_1_BH_0','1', 'Integer', 'Number of barcodes that have been called as supporting haplotype 1 which are on reads that have support for the allele which has been phased as haplotype 0')
    BX_HAP_OR = ('BX_HAP_OR', '1', 'Float', "Barcode aware haplotype filtering score (log odds ratio currently)")
    BARCODE_AWARE_FILTER = [("BARCODE_AWARE_FILTER", "Uses haplotype information from the fragments and the alleles to filter some variants that are not consistent with haplotype (ie variants should have most of their allele haplotype 0 alleles coming from barcodes whose fragments are haplotype 0 etc)")]
    extra_fields = [AH_0_BH_0, AH_1_BH_1, AH_0_BH_1, AH_1_BH_0, BX_HAP_OR]
    input_variants = tk_io.VariantFileReader(args.variants)
    with open(outs.default.strip(".gz"),'w') as output_file:
        output_variants = tk_io.VariantFileWriter(output_file, template_file=open(args.variants, 'r'), new_info_fields=extra_fields, new_filters = BARCODE_AWARE_FILTER)
        variant_iterator = tk_io.get_variant_iterator_pos(input_variants, None, args.locus)
        for record in variant_iterator:
            sample = record.samples[0]
            ref = tk_io.get_record_ref(record)
            alt_alleles = tk_io.get_record_alt_alleles(record)

            if not tk_io.get_record_passes_filters(record):
                output_variants.write_record(record)
                continue
            if len(sample.gt_alleles) > 1:
                genotype_1 = int(sample.gt_alleles[0])
                genotype_2 = int(sample.gt_alleles[1])
                if genotype_1 == genotype_2:
                    output_variants.write_record(record)
                    continue #homozygous, can't filter this way
            else:
                output_variants.write_record(record)
                continue #homozygous, can't filter this way


            chrom = tk_io.get_record_chrom(record)
            if not chrom == "chrM":
                variant_barcode_info = load_variant_barcode_phasing_info(record, fragment_barcode_info)
                if not barcode_aware_filter(record, variant_barcode_info):
                    if record.FILTER is None:
                        record.FILTER = []
                    if tk_io.get_var_type(ref, alt_alleles[0]) == "S" and ((vc_mode == 'call') or (vc_mode == "precalled_plus" and "TENX" in record.INFO)):
                        record.FILTER.append("BARCODE_AWARE_FILTER")
            output_variants.write_record(record)

def load_variant_barcode_phasing_info(record, fragment_barcode_info):
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    end = pos + tk_io.get_record_max_length(record)
    barcode_info = {}
    sample = record.samples[0]
    phase_set = int(get_data(sample.data, "PS", -1))
    for line in tk_tabix.tabix_safe_fetch(fragment_barcode_info, chrom, pos, end + 1):
        info = line.strip("\n").split("\t")
        barcode = info[6]
        frag_phase_set = int(info[3])
        if frag_phase_set != phase_set and phase_set != -1:
            continue
        assert(not barcode in barcode_info)
        barcode_info[barcode] = (float(info[7]), float(info[8]), float(info[9]))
    return barcode_info


def get_data(data, field, default):
    v = getattr(data, field, default)
    if v is None:
        return default
    return v

def get_bcs(bc_string):
    bcs = bc_string.split(";")
    if bcs == ['']:
        bcs = []
    return bcs

def barcode_aware_filter(record, phase_set_barcodes):
    sample = record.samples[0]
    genotype_0 = -1
    genotype_1 = -1
    if len(sample.gt_alleles) > 1:
        genotype_0 = int(sample.gt_alleles[0])
        genotype_1 = int(sample.gt_alleles[1])
        if genotype_0 == genotype_1:
            return True #homozygous, can't filter this way
    else:
        return True #homozygous, can't filter this way

    barcodes = tk_io.get_record_barcodes(record)
    def get_bc(pos):
        if barcodes is None:
            return []

        if len(barcodes) > pos:
            return barcodes[pos]
        else:
            return []

    bcs_for_allele_hap0 = get_bc(genotype_0)
    bcs_for_allele_hap1 = get_bc(genotype_1)

    haplotype_0_bcs_supporting_allele_hap0 = 0
    haplotype_0_bcs_supporting_allele_hap1 = 0
    haplotype_1_bcs_supporting_allele_hap0 = 0
    haplotype_1_bcs_supporting_allele_hap1 = 0

    for barcode_for_allele_hap0 in bcs_for_allele_hap0:
        barcode_for_allele_hap0_id = barcode_for_allele_hap0.split("_")[0]
        if barcode_for_allele_hap0_id in phase_set_barcodes:
            barcode_for_allele_hap0_frag_hap = phase_set_barcodes.get(barcode_for_allele_hap0_id)
            (bx_a0_hap_0, bx_a0_hap_1, bx_a0_hap_mixed) = barcode_for_allele_hap0_frag_hap
            if bx_a0_hap_0 > 0.999:
                haplotype_0_bcs_supporting_allele_hap0 += 1
            elif bx_a0_hap_1 > 0.999:
                haplotype_1_bcs_supporting_allele_hap0 += 1
    for barcode_for_allele_hap1 in bcs_for_allele_hap1:
        barcode_for_allele_hap1_id = barcode_for_allele_hap1.split("_")[0]
        if barcode_for_allele_hap1_id in phase_set_barcodes:
            barcode_for_allele_hap1_frag_hap = phase_set_barcodes.get(barcode_for_allele_hap1_id)
            (bx_a1_hap_0, bx_a1_hap_1, bx_a1_hap_mixed) = barcode_for_allele_hap1_frag_hap
            if bx_a1_hap_0 > 0.999:
                haplotype_0_bcs_supporting_allele_hap1 += 1
            elif bx_a1_hap_1 > 0.999:
                haplotype_1_bcs_supporting_allele_hap1 += 1
    score = score_barcode_haplotype_support(sample.phased, haplotype_0_bcs_supporting_allele_hap0, haplotype_1_bcs_supporting_allele_hap1, haplotype_0_bcs_supporting_allele_hap1, haplotype_1_bcs_supporting_allele_hap0)
    record.INFO['AH_0_BH_0'] = haplotype_0_bcs_supporting_allele_hap0
    record.INFO['AH_1_BH_1'] = haplotype_1_bcs_supporting_allele_hap1
    record.INFO['AH_1_BH_0'] = haplotype_0_bcs_supporting_allele_hap1
    record.INFO['AH_0_BH_1'] = haplotype_1_bcs_supporting_allele_hap0
    record.INFO['BX_HAP_OR'] = score
    return score >= 0

def score_barcode_haplotype_support(phased, h00, h11, h01, h10):
    p = 0.001
    h_good = h00+h11
    h_bad = h01+h10
    h_total = h_good+h_bad
    log_odds = h_good * math.log(1-p) + h_bad * math.log(p) - h_total * math.log(0.5)
    if phased == 0:
        log_odds = max(log_odds, h_bad * math.log(1-p) + h_good * math.log(p) - h_total * math.log(0.5))
    return log_odds




