#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Call variants from a BAM file
#
import os.path
import tenkit.log_subprocess
import tenkit.constants
from tenkit.exceptions import NotSupportedException
import tenkit.bio_io as tk_io
import tenkit.reference as tk_reference
import martian
from tenkit.regions import Regions
import tenkit.tabix as tk_tabix
import tenkit.chunk_utils as tk_chunks
import longranger.variant_caller as vc

__MRO__ = """
stage CALL_SNPINDELS(
    in  string vc_precalled,
    in  string variant_mode    "'freebayes' or 'gatk:/path/to/GenomeAnalysisTK.jar'",
    in  string restrict_locus,
    in  bed    targets_file,
    in  bam    input           "sorted and indexed bam file",
    in  string reference_path  "path to reference genome",
    in  bed    high_coverage_excluded_bed,
    out vcf                    "output vcf",
    out vcf    precalled,
    src py     "stages/snpindels/call_snpindels",
) split using (
    in  string locus,
)
"""

def split(args):
    vc_mode, variant_caller, precalled_filename, gatk_path = tk_io.get_vc_mode(args.vc_precalled, args.variant_mode)
    precalled_file = None
    if vc_mode == "precalled" or vc_mode == "precalled_plus":
        mem_gb = 8
        threads = 1
        precalled_file  = martian.make_path("precalled_vcf.vcf")
        tenkit.log_subprocess.check_call(['cp', precalled_filename, precalled_file])
        tk_tabix.index_vcf(precalled_file)
        precalled_file = precalled_file+".gz"
    if vc_mode != "precalled":
        if variant_caller == 'freebayes':
            mem_gb = 5
            threads = 1
        elif variant_caller == "gatk":
            mem_gb = 8
            threads = 2
            # make sure the gatk jar file exists
            if gatk_path is None:
                martian.throw("variant_caller 'gatk' selected, must supply path to gatk jar file -- e.g. \"gatk:/path/to/GenomeAnalysisTK.jar\"")

            gatk_loc = gatk_path
            if not (os.path.exists(gatk_loc)):
                martian.throw("variant_caller 'gatk' selected, gatk jar file does not exist: %s" % gatk_loc)
        else:
            raise NotSupportedException('Variant caller not supported: ' + vc_mode)

    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)
    bam_chunk_size_gb = 3.0

    if args.restrict_locus is None:
        loci = tk_chunks.get_sized_bam_chunks(args.input, bam_chunk_size_gb, contig_whitelist=primary_contigs, 
            extra_args = {'__mem_gb': mem_gb, '__threads': threads, 'split_input': precalled_file})
    else:
        loci = [{'locus': args.restrict_locus}]

    return {'chunks': loci}

def join(args, outs, chunk_defs, chunk_outs):
    outs.default = [chunk.default for chunk in chunk_outs]
    outs.precalled = [chunk.precalled for chunk in chunk_outs]

def get_coverage_regions(args):
    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    regions = Regions(tk_io.get_target_regions(open(args.high_coverage_excluded_bed)).get(chrom))
    if regions == None:
        regions = Regions()
    return regions

def main(args, outs):
    vc_mode, variant_caller, precalled_file, gatk_path = tk_io.get_vc_mode(args.vc_precalled, args.variant_mode)
    locus = args.locus
    (chrom, start, stop) = tk_io.get_locus_info(locus)
    fasta_path = tk_reference.get_fasta(args.reference_path)

    bedfile = outs.default+".bed"
    regions = Regions()
    if args.targets_file is not None:
        for (chrom, start, end) in tk_io.get_bed_iterator(args.targets_file, args.locus):
            regions.add_region((start, end))
    else:
        (chrom, start, stop) = tk_io.get_locus_info(args.locus)
        regions.add_region((start, stop))
    coverage_regions = None
    if (vc_mode != "precalled") and args.high_coverage_excluded_bed is not None:
        coverage_regions = get_coverage_regions(args)
        regions = regions.intersect(coverage_regions)

    bed_length = 0
    with open(bedfile, 'w') as bed_writer:
        for region in regions.get_region_list():
            (start, end) = region
            bed_writer.write(chrom+"\t"+str(start)+"\t"+str(end)+"\n")
            bed_length += 1
    if vc_mode == "precalled" or vc_mode == "precalled_plus":
        outs.default = None
        precalled_vars_path = args.split_input
        vcf = tk_io.VariantFileReader(precalled_vars_path)
        with open(outs.precalled, "w") as file_write:
            output = tk_io.VariantFileWriter(file_write, template_file=open(precalled_vars_path))
            variant_iter = tk_io.get_variant_iterator_pos(vcf, bedfile, args.locus)
            for record in variant_iter:
                output.write_record(record)
    if not (vc_mode == "precalled"):
        outs.precalled = None
        primary_contigs = tk_reference.load_primary_contigs(args.reference_path)
        if bed_length > 0 and chrom in primary_contigs:
            vc.run_variant_caller(variant_caller, gatk_path, args.__mem_gb, fasta_path, args.input, outs.default, bedfile)
