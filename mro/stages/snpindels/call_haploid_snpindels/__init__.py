#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Call variants on each haplotype
#
import tenkit.bam as tk_bam
import tenkit.constants
import tenkit.bio_io as tk_io
import tenkit.chunk_utils as tk_chunks
import tenkit.reference as tk_reference
import subprocess
import tenkit.reference as tk_ref
import tenkit.log_subprocess
import longranger.variants
import longranger.variant_caller as vc

__MRO__ = """
stage CALL_HAPLOID_SNPINDELS(
    in vcf.gz putative_variants,
    in bam    input,
    out string locus,
    out vcf,
) split using (
    in string locus,
)
"""

def split(args):

    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)
    bam_chunk_size_gb = 3.0
    loci = tk_chunks.get_sized_bam_chunks(args.input, bam_chunk_size_gb, contig_whitelist=primary_contigs) 

    _, caller, _, _ = tk_io.get_vc_mode(args.precalled, args.mode)
    chunks = []
    gb = 4
    threads = 1
    if caller == "gatk":
        gb = 6
        threads = 2
    for chunk in loci:
        chunk.update({"__mem_gb": gb, "__threads": threads})
        chunks.append(chunk)
    return {'chunks': chunks}


def main(args, outs):
    outs.chunk_locus = args.locus
    mode, _, _, _ = tk_io.get_vc_mode(args.precalled, args.mode)
    if mode == "precalled":
        return
    else:
        reference_pyfasta = tk_ref.open_reference(args.reference_path)
        bam_in = tenkit.bam.create_bam_infile(args.input)
        vc_mode, variant_caller, precalled_file, gatk_path = tk_io.get_vc_mode(args.precalled, args.mode)
        (hap1_vcf, bam1) = call_haploid(1, bam_in, args.locus, args.reference_path, variant_caller, gatk_path, args.__mem_gb)
        (hap2_vcf, bam2) = call_haploid(2, bam_in, args.locus, args.reference_path, variant_caller, gatk_path, args.__mem_gb)
        merge_naive(hap1_vcf, hap2_vcf, args.locus, outs.default, bam_in, bam1, bam2, reference_pyfasta, args.putative_variants)

def merge_naive(vcf1, vcf2, locus, output_filename, bam, bam1, bam2, reference_pyfasta, putative_vcf):
    vcf_hap1 = tenkit.bio_io.VariantFileReader(vcf1)
    vcf_hap2 = tenkit.bio_io.VariantFileReader(vcf2)
    with open(output_filename,'w') as synthetic_diploid:
        HAPLOCALLED = ('HAPLOCALLED','1','Integer','called by the haplotype caller in which each haplotype is separated and variant calling is done with ploidy=1', None, None)
        new_info_fields = [HAPLOCALLED]
        vcf_out1 = tenkit.bio_io.VariantFileWriter(synthetic_diploid, template_file = open(putative_vcf), new_info_fields=new_info_fields)
        for (var1, var2) in pair_iter(vcf_hap1.record_getter(), vcf_hap2.record_getter()):
            if var1 == None:
                tenkit.bio_io.set_record_genotype_phased(var2, [0,1], True)
                var2.INFO['AF'] = 0.5
                var2.INFO['AC'] = 1
                var2.INFO['AN'] = 2
                filter_variant(var2, bam, reference_pyfasta)
                vcf_out1.write_record(var2)
            elif var2 == None:
                tenkit.bio_io.set_record_genotype_phased(var1, [1,0], True)
                var1.INFO['AF'] = 0.5
                var1.INFO['AC'] = 1
                var1.INFO['AN'] = 2
                filter_variant(var1, bam, reference_pyfasta)
                vcf_out1.write_record(var1)
            elif tenkit.bio_io.get_record_alt_alleles(var1)[0] != tenkit.bio_io.get_record_alt_alleles(var2)[0]:
                filter_variant(var1, bam, reference_pyfasta)
                filter_variant(var2, bam, reference_pyfasta)
                var = var1
                if tk_io.get_record_passes_filters(var1):
                    if tk_io.get_record_passes_filters(var2):
                        tenkit.bio_io.set_record_genotype_phased(var1, [1,2], True)
                        var1.ALT = [var1.ALT[0], var2.ALT[0]]
                        var1.INFO['AF'] = 0.5
                        var1.INFO['AC'] = 2
                        var1.INFO['AN'] = 2
                    else:
                        tenkit.bio_io.set_record_genotype_phased(var1, [0,1], True)
                        var1.INFO['AF'] = 0.5
                        var1.INFO['AC'] = 1
                        var1.INFO['AN'] = 2
                elif tk_io.get_record_passes_filters(var2):
                    tenkit.bio_io.set_record_genotype_phased(var2, [1,0], True)
                    var = var2
                    var.INFO['AF'] = 0.5
                    var.INFO['AC'] = 1
                    var.INFO['AN'] = 2
                else:
                    tenkit.bio_io.set_record_genotype_phased(var1, [1,2], True)
                    var1.ALT = [var1.ALT[0], var2.ALT[0]]
                    var1.INFO['AF'] = 0.5
                    var1.INFO['AC'] = 2
                    var1.INFO['AN'] = 2
                vcf_out1.write_record(var)
            else:
                var1.QUAL = (var1.QUAL if var1.QUAL else 0) + (var2.QUAL if var2.QUAL else 0)
                tenkit.bio_io.set_record_genotype_phased(var1, [1,1], True)
                var1.INFO['AC'] = 2
                var1.INFO['AN'] = 1
                vcf_out1.write_record(var1)

def get_phase_set(record, bam):
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    for read in bam.fetch(chrom, pos-1, pos+1):
        if dict(read.tags).get('PS') is not None:
            return dict(read.tags).get('PS')
    return None

def filter_variant(var, bam, reference_pyfasta):
    if tk_io.get_record_qual(var) < 50:
        tk_io.set_record_filters(var, ['10X_QUAL_FILTER'])
        return
    chrom = tk_io.get_record_chrom(var)
    pos = tk_io.get_record_pos(var)
    ref = tk_io.get_record_ref(var)
    alts = tk_io.get_record_alt_alleles(var)
    (counts, _, _, _, _, _) = tk_bam.get_allele_read_info(chrom, pos, ref, alts, 30, 30, 30, 45, bam, reference_pyfasta)
    if float(counts[1]) < 2 or float(counts[1])/float(counts[0] + counts[1]) < 0.15:
        tk_io.set_record_filters(var, ['10X_ALLELE_FRACTION_FILTER'])
def pair_iter(i1, i2):
    v1 = None
    v2 = None

    while True:

        if v1 is None:
            try:
                v1 = i1.next()
            except StopIteration:
                if v2 is not None:
                    yield (None, v2)
                for x2 in i2:
                    yield (None, x2)
                break
        if v2 is None:
            try:
                v2 = i2.next()
            except StopIteration:
                if v1 is not None:
                    yield (v1, None)
                for x1 in i1:
                    yield (x1, None)
                break
        k1 = tk_io.get_record_pos(v1)
        k2 = tk_io.get_record_pos(v2)
        if k1 == k2:
            yield (v1, v2)
            v1 = None
            v2 = None
        elif k1 < k2:
            yield (v1, None)
            v1 = None
        else:
            yield (None, v2)
            v2 = None

def call_haploid(haplotype, bam, locus, reference_path, variant_caller, gatk_path, mem_gb):
    bam_name = "hap"+str(haplotype)+".bam"
    haploid_bam, _ = tenkit.bam.create_bam_outfile(bam_name,None, None, template = bam)
    (chrom, start, stop) = tk_io.get_locus_info(locus)
    for read in bam.fetch(chrom, start, stop):
        readhap = dict(read.tags).get('HP')
        if readhap != None and int(readhap) == haplotype:
            haploid_bam.write(read)
    haploid_bam.close()
    tk_bam.index(bam_name)
    tmp_vcf_name = "tmp_hap"+str(haplotype)+".vcf"
    vcf_name = "hap"+str(haplotype)+".vcf"

    fasta_path = tk_ref.get_fasta(reference_path)
    vc.run_variant_caller(variant_caller, gatk_path, mem_gb, fasta_path, bam_name, tmp_vcf_name, haploid_mode=True)

    longranger.variants.canonicalize(tmp_vcf_name, vcf_name)
    tenkit.tabix.index_vcf(vcf_name)
    bam_in = tk_bam.create_bam_infile(bam_name)
    return (vcf_name+".gz", bam_in)

def join(args, outs, chunk_defs, chunk_outs):
    mode, _, precalled, _ = tk_io.get_vc_mode(args.precalled, args.mode)
    if mode == "precalled":
        subprocess.check_call(['cp', args.putative_variants, outs.precalled+".gz"])
        outs.precalled = outs.precalled + ".gz"
        subprocess.check_call(['cp', args.putative_variants+".tbi", outs.precalled+".tbi"])
    else:
        outs.precalled = None

    outs.default = [x.default for x in chunk_outs]
    outs.chunk_locus = [chunk.chunk_locus for chunk in chunk_outs]



