#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Populates INFO fields of VCF file with statistics useful for filtering
#  also strips fields from precalled VCF that are not in VCF_WHITE_LIST_INFO_FIELDS
#
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
from tenkit.constants import VARIANT_CALL_FILTER
import tenkit.log_subprocess
from tenkit.constants import VCF_WHITE_LIST_INFO_FIELDS
from longranger.variants import canonicalize
import tenkit.reference
import os
import os.path
import vcf
import numpy
import martian
import math

__MRO__ = """
stage POPULATE_INFO_FIELDS(
    in  string vc_precalled,
    in  string variant_mode,
    in  vcf    input,
    in  string reference_path,
    in  bam    bam,
    in  int    min_mapq_attach_bc,
    out vcf,
    src py     "stages/snpindels/populate_info",
) split using (
    in  vcf    chunk_input,
)
"""

def split(args):
    vcf_list = args.input
    precalled_list = args.vc_precalled
    locus_list = args.chunk_locus
    if type(vcf_list) is not list:
        vcf_list = [vcf_list]
    if locus_list is None:
        locus_list = [None for x in vcf_list]
    if precalled_list is None:
        precalled_list = [None for x in vcf_list]
    if type(precalled_list) is not list:
        precalled_list = [precalled_list]

    chunk_defs = []
    for vcf_file, precalled, locus in zip(vcf_list, precalled_list, locus_list):
        chunk_defs.append({'chunk_input':vcf_file, 'precalled_chunk':precalled, 'precalled_no_chunk': (type(args.vc_precalled) is not list), 'locus': locus, '__mem_gb': 4})
    return {'chunks': chunk_defs}

def main(args, outs):
    def real_file(f):
        return f is not None and os.path.isfile(f) and os.path.getsize(f) > 0
    precalled = False
    called = False
    if real_file(args.precalled_chunk):
        precalled = True
        if args.precalled_no_chunk:
            return
    elif real_file(args.chunk_input):
        called = True
    if called:
        input_vcf = "canonicalized.vcf"
        canonicalize(args.chunk_input, input_vcf)
    elif precalled:
        canonicalize(args.precalled_chunk, "canonicalized.vcf")
        input_vcf = outs.default+"tmp.vcf"
        vcf_reader = tk_io.VariantFileReader("canonicalized.vcf")
        with open(input_vcf, "w") as file_write:
            output = tk_io.VariantFileWriter(file_write, template_file=open("canonicalized.vcf"))
            for record in vcf_reader.record_getter():
                (bad, message) = tk_io.is_record_bad_variant(record)
                if bad:
                    try:
                        record_strified = str(record)
                    except:
                        record_strified = "<error displaying record>"
                    martian.exit("error on vcf record: "+record_strified+", "+message)
                if record.QUAL == None:
                    record.QUAL = 0
                info_fields = {key: record.INFO.get(key) for key in VCF_WHITE_LIST_INFO_FIELDS if key in record.INFO}
                record.INFO = info_fields
                sample_call = tk_io.get_record_sample_call(record)
                data = sample_call.data
                data_dict = data._asdict()
                if "GT" in data_dict:
                    new_sample_vals = [data_dict["GT"]]
                    new_format = "GT"
                    new_fields = ["GT"]
                else:
                    new_sample_vals = ["./."]
                    new_format = "GT"
                    new_fields = ["GT"]
                data_instantiator = vcf.model.make_calldata_tuple(new_fields)
                data = data_instantiator(*new_sample_vals)
                sample_call.data = data
                record.samples[0] = sample_call
                record.FORMAT = new_format
                output.write_record(record)
    else:
        outs.default = None
        return

    reference_pyfasta = tenkit.reference.open_reference(args.reference_path)

    chunk_record_size = 10000

    vcf_reader = tk_io.VariantFileReader(input_vcf)

    # Check whether to write GL field
    v = vcf.Reader(open(input_vcf))
    write_gl = v.formats.has_key('GL')

    bam = tk_bam.create_bam_infile(args.bam)
    outfile = outs.default
    if not precalled:
        outfile = outs.default+"tmp2.vcf"
    with open(outfile, "w") as file_write:
        new_source = "10X/pipelines/stages/snpindels/attach_bcs_snpindels %s" % martian.get_pipelines_version()
        new_formats = [("BX", ".", "String", "Barcodes and Associated Qual-Scores Supporting Alleles")]
        AO = ('AO','.','Integer','Alternate allele observed count', None, None)
        RO = ('RO','1','Integer','Reference allele observed count', None, None)
        post_homopolymer_counts_field = ('POSTHPC','.', 'Integer', 'Postvariant homopolymer count', None, None)
        post_homopolymer_base_field = ('POSTHPB','.', 'Character','Postvariant homopolymer base', None, None)
        post_dinucleotide_base_field = ('POSTDNB','.', 'String', 'Post variant dinucleotide repeat sequence', None, None)
        post_dinucleotide_counts_field = ('POSTDNC','.','Integer','Post variant dinucleotide repeat count', None, None)
        post_trinucleotide_base_field = ('POSTTNB','.','String','Post variant trinucleotide repeat sequence', None, None)
        post_trinucleotide_counts_field = ('POSTTNC','.','Integer','Post variant trinucleotide repeat count', None, None)
        mean_map_alt = ('MUMAP_ALT', '.', 'Float', 'Mean mapping scores of alt alleles', None, None)
        mean_map_ref = ('MUMAP_REF', '1', 'Float', 'Mean mapping score of ref allele', None, None)
        rescued = ('RESCUED','.', 'Integer','How many reads were rescued via cross barcode mapq correction', None, None)
        not_rescued = ('NOT_RESCUED','.','Integer', 'How many reads were not rescued via cross barcode mapq correction', None, None)
        mean_molecule_difference = ("MMD", '.', 'Float', 'Mean molecule divergence from reference per read', None, None)
        haplocalled = ("HAPLOCALLED", "1", "Integer", "1 for variants that were called after phasing via splitting the bam into its component haplotypes and calling variants in haploid mode", None, None)
        extra_fields = [post_homopolymer_counts_field, post_homopolymer_base_field, mean_map_ref, mean_map_alt, AO, RO, mean_molecule_difference, rescued, not_rescued,
                        post_dinucleotide_base_field, post_dinucleotide_counts_field, post_trinucleotide_base_field, post_trinucleotide_counts_field, haplocalled]
        output = tk_io.VariantFileWriter(file_write, template_file=open(input_vcf), new_source=new_source, new_info_fields=extra_fields,
                                  new_format_fields=new_formats)
        record_list = []
        for record in vcf_reader.record_getter():
            if len(record_list) == chunk_record_size:
                for record_out in record_list:
                    if tk_io.get_record_passes_filters(record_out):
                        populate_fields(record_out, bam, reference_pyfasta, args)
                for record_out in record_list:
                    output.write_record(record_out)
                record_list = []

            record_list.append(record)

        for record_out in record_list:
            if tk_io.get_record_passes_filters(record_out):
                populate_fields(record_out, bam, reference_pyfasta, args)
        for record_out in record_list:
            output.write_record(record_out)

    if not precalled and args.haploid_merge is not None:
        output_filename = outs.default+"tmp3.vcf"
        merge_haploid(outfile, args.haploid_merge, args.locus, output_filename, bam, reference_pyfasta, args, write_gl)
        outfile = output_filename

    # Filter variants
    if not precalled:
        last_file = outfile
        for filter in VARIANT_CALL_FILTER:
            with open(outs.default+"tmp.vcf",'w') as output_file:
                filt = VARIANT_CALL_FILTER[filter]
                tenkit.log_subprocess.check_call(['bcftools','filter','-O', 'v','--soft-filter',filter,'-e',filt,'-m','\'+\'', last_file],stdout=output_file)
            tenkit.log_subprocess.check_call(['mv', outs.default+"tmp.vcf", outs.default+"tmp2.vcf"])
            last_file = outs.default+"tmp2.vcf"

        tenkit.log_subprocess.check_call(['mv', last_file, outs.default])

def merge_haploid(novel_vcf, putative_vcf, locus, output_filename, bam, reference_pyfasta, args, add_gl = True):
    novel = tenkit.bio_io.VariantFileReader(novel_vcf)
    putative_variants = tenkit.bio_io.VariantFileReader(putative_vcf)
    (chrom, start, stop) = tk_io.get_locus_info(locus)
    with open(output_filename,'w') as output:
        output_vcf = tenkit.bio_io.VariantFileWriter(output, template_file = open(putative_vcf))
        for (novel, putative) in pair_iter(novel.record_getter(), putative_variants.record_getter(fetch_chrom=chrom, fetch_start = start, fetch_end=stop)):
            if putative is not None and (tk_io.get_record_passes_filters(putative) or tk_io.get_record_filters(putative) != ['10X_QUAL_FILTER']):
                putative.INFO['HAPLOCALLED'] = 0
                output_vcf.write_record(putative)
            elif novel is not None:
                tk_io.set_record_phase_set(novel, get_phase_set(novel, bam))
                tk_io.set_record_phase_qual(novel, 25)
                tk_io.set_record_junction_qual(novel, 25)
                populate_fields(novel, bam, reference_pyfasta, args)
                novel.INFO['HAPLOCALLED'] = 1
                if add_gl:
                    tk_io.set_record_genotype_likelihoods(novel, calculate_psuedo_genotype_likelihoods(novel))
                if novel.QUAL is not None: # freebayes ploidy 1 gives some variants with '.' as the qual. These are extremely low quality not worth even tracking
                    output_vcf.write_record(novel)
            else:
                putative.INFO['HAPLOCALLED'] = 0
                output_vcf.write_record(putative)
def calculate_psuedo_genotype_likelihoods(rec):
    ao = rec.INFO['AO'][0]
    ro = rec.INFO['RO']
    novar = ao*(-10)
    het = ao*(-3) + ro*(-3)
    hom = ro*(-10)
    best = max(max(novar,het),hom)
    return str(novar-best)+','+str(het-best)+','+str(hom-best)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    mode, _, _,_ = tk_io.get_vc_mode(args.vc_precalled, args.variant_mode)
    if mode=="precalled" and args.vc_precalled is not None and type(args.vc_precalled) is not list:
        outs.default = args.vc_precalled
        return
    # List of inputs
    files = [f.default for f in chunk_outs if os.path.isfile(f.default)]
    if len(files) == 0:
        outs.default = None
        return
    cat_filename = outs.default[:-3]
    tk_io.combine_vcfs(cat_filename, files)

def populate_fields(record, bam, reference_pyfasta, args):
    alleles = tk_io.get_record_alt_alleles(record)
    ref = tk_io.get_record_ref(record)
    post_homopolymer_counts = []
    post_homopolymer_bases = []
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    ref = tk_io.get_record_ref(record)
    post_homopolymer_counts = []
    post_homopolymer_bases = []
    post_dinucleotide_counts = []
    post_dinucleotide_bases = []
    post_trinucleotide_counts = []
    post_trinucleotide_bases = []
    for allele in alleles:
        variant_length = tk_io.get_allele_length(ref, allele)
        if variant_length != 0:
            post_hp_c, post_hp_b = populate_repeat_info(record, bam, variant_length, reference_pyfasta, 1)
            post_dn_c, post_dn_b = populate_repeat_info(record, bam, variant_length, reference_pyfasta, 2)
            post_tn_c, post_tn_b = populate_repeat_info(record, bam, variant_length, reference_pyfasta, 3)
            post_homopolymer_counts.append(post_hp_c)
            post_homopolymer_bases.append(post_hp_b)
            post_dinucleotide_counts.append(post_dn_c)
            post_dinucleotide_bases.append(post_dn_b)
            post_trinucleotide_counts.append(post_tn_c)
            post_trinucleotide_bases.append(post_tn_b)
    if len(post_homopolymer_counts) != 0:
        record.INFO['POSTHPC'] = post_homopolymer_counts
        record.INFO['POSTHPB'] = post_homopolymer_bases
        record.INFO['POSTDNC'] = post_dinucleotide_counts
        record.INFO['POSTDNB'] = post_dinucleotide_bases
        record.INFO['POSTTNC'] = post_trinucleotide_counts
        record.INFO['POSTTNB'] = post_trinucleotide_bases

    (counts, mean_mapqs, bc_qual_string, molecule_differences, AS, rescue) = tk_bam.get_allele_read_info(chrom, pos, ref, alleles, 30, -1, args.min_mapq_attach_bc, args.default_indel_qual, bam, reference_pyfasta)

    tk_io.set_record_barcodes(record, bc_qual_string)
    record.INFO['MMD'] = numpy.mean(molecule_differences[1])
    if math.isnan(record.INFO['MMD']):
        record.INFO['MMD'] = -1
    record.INFO['MUMAP_REF'] = mean_mapqs[0]
    record.INFO['MUMAP_ALT'] = mean_mapqs[1:]
    record.INFO['RO'] = counts[0]
    record.INFO['AO'] = counts[1:]
    record.INFO['RESCUED'] = numpy.sum(numpy.sum(x) for x in rescue)
    record.INFO['NOT_RESCUED'] = numpy.sum([y for y in [numpy.sum([1-z for z in x]) for x in rescue]])

def get_phase_set(record, bam):
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    for read in bam.fetch(chrom, pos-1, pos+1):
        if dict(read.tags).get('PS') is not None:
            return dict(read.tags).get('PS')
    return None
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

def populate_repeat_info(record, bam, variant_length, reference_pyfasta, length):
    post_poly_count = 0
    post_poly_base = None
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    lastBase = None
    gap = min(30, len(reference_pyfasta[chrom])-pos-1)
    #sequence = {x: tk_bam.get_base_counts_at_locus(chrom, pos + x, bam) for x in range(0 , gap + max(-variant_length,1))}
    sequence = reference_pyfasta[chrom][(pos+1):(pos+gap+1)].upper()
    #from the base after the indel to the end of the gap
    for base in range(0, gap, length):
        if lastBase is None:
            post_poly_count = 1
            post_poly_base = sequence[base:base+length]
            lastBase = post_poly_base
        elif lastBase is not None:
            if lastBase == sequence[base:base+length]:
                post_poly_count += 1
            else:
                break
        else:
            break
    return post_poly_count, post_poly_base
