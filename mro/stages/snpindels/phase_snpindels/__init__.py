#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Attach read barcode information to variants
#
import tenkit.log_subprocess
import tenkit
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.tabix as tk_tabix
import tenkit.stats as tk_stats
import tenkit.constants
import tenkit.reference
import os
import pysam
from phaser import Phaser
import stitcher
import martian
import random
import numpy as np

__MRO__ = """
stage PHASE_SNPINDELS(
    in  string reference_path,
    in  string sex,
    in  string restrict_locus,
    in  bam    bam_file,
    in  h5     fragments,
    in  vcf.gz input,
    in  float  bc_mix_prob,
    in  float  min_var_hap_conf,
    in  float  min_junction_hap_conf,
    in  int    hap_block_size,
    in  int    hap_block_buffer_size,
    in  int    max_reassign_rounds,
    in  int    chunk_stitching_overlap,
    out vcf.gz,
    out tsv.gz fragment_phasing,
    src py     "stages/snpindels/phase_snpindels",
) split using (
    in  string locus,
)
"""

def split(args):
    input_bam = tk_bam.create_bam_infile(args.bam_file)

    if args.sex is None or args.sex.lower() in ['m', 'male']:
        remove_chroms = ['chrX', 'chrY', 'chrM', 'X', 'Y', 'MT', 'M']
    elif args.sex.lower() in ['f', 'female']:
        remove_chroms = ['chrY', 'chrM', 'Y', 'MT', 'M']
    else:
        martian.throw("Unrecognized sex: %s" % args.sex)

    primary_contigs = tenkit.reference.load_primary_contigs(args.reference_path)

    # estimate density of het snps and barcodes
    primary_contig_lengths = [(chrom, length) for (chrom, length) in zip(input_bam.references, input_bam.lengths) if chrom in primary_contigs]
    (frac_het_snps, bcs_per_het_snp, het_rate) = smooth_sample_bcs_and_het_snps(args.input, primary_contig_lengths)
    martian.log_info("Fraction of SNPs that are het: %f, BCs per SNP: %f, Hets per bp: %f" % (frac_het_snps, bcs_per_het_snp, het_rate))

    # Set up dynamic chunk sizes to make smaller chunks for highly het organisms
    het_rate = max(min(0.05, het_rate), 0.0001)
    parallel_locus_size = min(tenkit.constants.PARALLEL_LOCUS_SIZE, tk_stats.robust_divide(60000, het_rate))

    if args.restrict_locus is None:
        loci = tk_bam.generate_tiling_windows(input_bam, parallel_locus_size, overlap=args.chunk_stitching_overlap)
    else:
        loci = [args.restrict_locus]

    chunks = []
    for (idx, locus) in enumerate(loci):
        (chrom, start, end) = tk_io.get_locus_info(locus)

        if args.fragments is None or chrom not in primary_contigs or chrom in remove_chroms:
            mem_gb = 3
            martian.log_info("Chunk %d: No phasing, requesting %d GB" % (idx, mem_gb))
            chunk = {'locus': locus, '__mem_gb': mem_gb, 'do_phasing': False}
        else:
            est_het_snps = round(frac_het_snps * count_records(args.input, chrom, start, end))
            est_gb = np.ceil(0.5 + 400.0/1e9 * est_het_snps * bcs_per_het_snp) # empirical memory usage
            min_gb = 4
            if np.isnan(est_gb):
                mem_gb = min_gb
            else:
                mem_gb = max(min_gb, int(est_gb))
            martian.log_info("Chunk %d: Estimated %f het SNPs, requesting %f GB" % (idx, est_het_snps, mem_gb))
            chunk = {'locus': locus, '__mem_gb': mem_gb, 'do_phasing': True}

        chunks.append(chunk)

    return {'chunks': chunks, 'join': {'__mem_gb': 16, '__threads': 4}}

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()
    frag_files = [c.fragment_phasing.strip('.gz') for c in chunk_outs]
    vcf_files = [c.default.strip('.gz') for c in chunk_outs]
    loci = [tk_io.get_locus_info(c.locus) for c in chunk_defs]

    out_frag_base = outs.fragment_phasing.strip('.gz')
    out_vcf_base = outs.default.strip('.gz')

    # stitch phase blocks
    (stitched_chrom_vcfs, stitched_chrom_frags) = stitcher.multi_join_parallel(frag_files, vcf_files, loci, args.__threads)


    # combine the chromosome-level outputs
    combine_frags(out_frag_base, stitched_chrom_frags)
    tk_io.combine_vcfs(out_vcf_base, stitched_chrom_vcfs)
    if args.vc_precalled is not None:
        outs.vc_precalled = outs.default
    else:
        outs.vc_precalled = None

    # final indexing
    pysam.tabix_index(out_frag_base, seq_col=0, start_col=1, end_col=2)

def combine_frags(out_frags, in_frags_list):
    tmp_fn = out_frags + ".tmp"

    # cat files to tmp
    with open(tmp_fn, 'w') as fout:
        for fn in in_frags_list:
            with open(fn) as fin:
                for l in fin:
                    fout.write(l)

    # write header to final file
    with open(out_frags, 'w') as fout:
        header = ['chrom', 'frag_start', 'frag_end', 'phase_set', 'ps_start', 'ps_end', 'bc', 'h0', 'h1', 'hmix', 'reads', 'molecule_id']
        fout.write("#" + "\t".join(header) + "\n")

    # append sorted data to final file
    # limit sort to 4G to stay under the memory consumption limit
    tenkit.log_subprocess.check_call('sort -S 4G -k1,1 -k2,2n ' + tmp_fn  + ' >> ' + out_frags, shell=True)
    os.remove(tmp_fn)

def pass_variants(vfr, vfw, chrom, start, stop, strip_phasing_info=True):
    ''' Pass PyVCF records through, optionally removing the existing phasing info '''
    for record in vfr.record_getter(fetch_chrom=chrom, fetch_start=start, fetch_end=stop):
        if strip_phasing_info:
            try:
                (gt, phased) = tk_io.get_record_genotype_phased(record)
                tk_io.set_record_genotype_phased(record, gt, False)
            except:
                pass # This is for ./. or other malformed genotypes
            tk_io.set_record_phase_set(record, -1)
            tk_io.set_record_phase_qual(record, 0)
            tk_io.set_record_junction_qual(record, 0)
        vfw.write_record(record)

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    input_vfr = tk_io.VariantFileReader(args.input)

    bc_mix_prob = args.bc_mix_prob
    min_var_hap_conf = args.min_var_hap_conf
    min_junction_hap_conf = args.min_junction_hap_conf
    hap_block_size = args.hap_block_size
    hap_block_buffer_size = args.hap_block_buffer_size
    max_reassign_rounds = args.max_reassign_rounds
    chrom, start, stop = tk_io.get_locus_info(args.locus)

    output_file = open(outs.default.strip('.gz'), 'w')
    fragment_output_file = open(outs.fragment_phasing.strip('.gz'), 'w')
    vc_mode, _, _, _ = tk_io.get_vc_mode(args.vc_precalled, args.vc_mode)

    # Add the component name and the version of the phasing code
    new_source = "10X/pipelines/stages/snpindels/phase_snpindels %s" % martian.get_pipelines_version()
    new_filters = [("10X_PHASING_INCONSISTENT", "Uses haplotype information from the fragments and the alleles to filter some variants that are not consistent with phasing."),
                   ("10X_HOMOPOLYMER_UNPHASED_INSERTION", "Unphased insertions in homopolymer regions tend to be false positives")]
    new_formats = [("PS", 1, "Integer", "ID of Phase Set for Variant"),
                   ("PQ", 1, "Integer", "Phred QV indicating probability at this variant is incorrectly phased"),
                   ("JQ", 1, "Integer", "Phred QV indicating probability of a phasing switch error in gap prior to this variant"),
                   ]
    vfw = tk_io.VariantFileWriter(output_file, template_file=open(args.input), new_source=new_source,
                                new_format_fields=new_formats, new_filters = new_filters)
    if args.do_phasing:
        phaser = Phaser(input_vfr, args.fragments, chrom, start, stop, bc_mix_prob, min_junction_hap_conf, min_var_hap_conf, hap_block_buffer_size, hap_block_size, max_reassign_rounds, vc_mode)
        phaser.call_haps(vfw, fragment_output_file)
    else:
        pass_variants(input_vfr, vfw, chrom, start, stop, strip_phasing_info=True)
    output_file.close()
    fragment_output_file.close()

    tk_tabix.sort_unique_tabix_vcf(outs.default.strip('.gz'))

## Methods for dynamic memory estimation
# TODO move these to tenkit if they're generally useful

def count_records(vcf, chrom, start, end):
    tbx = pysam.Tabixfile(vcf)

    if chrom not in tbx.contigs:
        return 0

    count = 0
    for x in tbx.fetch(chrom, start, end):
        count += 1
    return count

def sample_by_locus(vcf, locus):
    (chrom, start, end) = locus
    recs = tk_io.VariantFileReader(vcf).record_getter(fetch_chrom = chrom, fetch_start = start, fetch_end = end)
    total_recs = 0
    het_snp_recs = 0
    het_snp_bcs = 0
    for rec in recs:
        total_recs += 1
        if rec.var_type == 'snp' and not tk_io.get_record_homozygous(rec):
            het_snp_recs += 1
            bcs_per_hap = tk_io.get_record_barcodes(rec)
            het_snp_bcs += sum([len(hap) for hap in bcs_per_hap])
    return (total_recs, het_snp_recs, het_snp_bcs)

def sample_bcs_and_het_snps(vcf, contigs):
    block_size = 1e6
    num_samples = 10
    total_recs = 0
    het_snp_recs = 0
    het_snp_bcs = 0
    total_size = 0
    for i in xrange(num_samples):
        (chrom, length) = random.choice(contigs)
        start = random.randint(0, max(0, length - block_size))
        end = min(start + block_size, length)
        (tot, snp, bcs) = sample_by_locus(vcf, (chrom, start, end))
        total_recs += tot
        het_snp_recs += snp
        het_snp_bcs += bcs
        total_size += (end-start)
    return (tk_stats.robust_divide(het_snp_recs, total_recs), tk_stats.robust_divide(het_snp_bcs, het_snp_recs), tk_stats.robust_divide(het_snp_recs, total_size))

def smooth_sample_bcs_and_het_snps(vcf, contigs):
    # take the mean of means, to get a smoother estimate.
    num_reps = 5
    mean_bcs = []
    mean_het_snps = []
    mean_het_rate = []
    for i in xrange(num_reps):
        (het_snps, bcs, het_rate) = sample_bcs_and_het_snps(vcf, contigs)
        mean_bcs.append(bcs)
        mean_het_snps.append(het_snps)
        mean_het_rate.append(het_rate)
    return (np.mean(mean_het_snps), np.mean(mean_bcs), np.mean(mean_het_rate))
