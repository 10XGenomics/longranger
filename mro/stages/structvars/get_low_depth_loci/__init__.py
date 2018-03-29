#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Finds loci with decreased depth.
#

import tenkit.pandas as pd
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import numpy as np
import cPickle
from itertools import groupby
from tenkit.chunk_utils import generate_chrom_loci, pack_loci
from tenkit.constants import PARALLEL_LOCUS_SIZE
import tenkit.hdf5 as tk_hdf5
import tenkit.reference as tk_reference

__MRO__ = """
stage GET_LOW_DEPTH_LOCI(
    in  bam    possorted_bam,
    in  string reference_path,
    in  bed    targets,
    in  h5     hap_coverage,
    in  int    bin_size,
    in  int    min_len,
    out pickle loci,
    out tsv    cov_summary,
    src py  "stages/structvars/get_low_depth_loci",
) split using (
    in  string[] loci,
)
"""


def split(args):
    assert(args.min_len >= 2 * args.bin_size)

    # We only need the bam to get the chromosome names and lengths.
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    if args.targets is None:
        target_regions = None
    else:
        with open(args.targets, 'r') as f:
            target_regions = tk_io.get_target_regions(f)

    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)

    all_loci = []
    for (chrom_name, chrom_size) in zip(in_bam.references, in_bam.lengths):
        if not chrom_name in primary_contigs:
            continue
        # The chunks will overlap by min_len. This will ensure that we don't
        # miss any regions of low depth that span chunk boundaries.
        new_loci = generate_chrom_loci(target_regions, chrom_name, chrom_size,
            PARALLEL_LOCUS_SIZE / 2, overlap = args.min_len)
        all_loci.extend(new_loci)
    in_bam.close()

    # Group loci
    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci, '__mem_gb':16} for loci in locus_sets]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    out_loci = []
    summary_df = None
    for chunk in chunk_outs:
        with open(chunk.loci, 'r') as f:
            out_loci.extend(cPickle.load(f))
        summary_df = pd.concat([summary_df,
            pd.read_csv(chunk.cov_summary, sep = '\t', header = 0, index_col = None)],
            ignore_index = True)

    # There might be some overlapping loci due to the overlap between chunks
    # but we hope that subsequent stages will deal with this.
    with open(outs.loci, 'w') as f:
        cPickle.dump(out_loci, f)
    summary_df.to_csv(outs.cov_summary, sep = '\t', header = True, index = False)

def main(args, outs):
    reader = tk_hdf5.DataFrameReader(args.hap_coverage)
    sel_cols = ['cov_q30_hap0', 'cov_q30_hap1', 'cov_q30_hap2']
    ext_cols = list(sel_cols)
    ext_cols.append('total_cov')

    out_loci = []
    summary_df = None
    for (chrom, start, stop) in (tk_io.get_locus_info(l) for l in args.loci):
        cov = reader.query((chrom, start, stop))
        cov['bin'] = np.array(cov['pos'] / args.bin_size, dtype = np.int)
        cov['total_cov'] = cov[sel_cols].sum(axis = 1)
        mean_cov = np.mean(cov['total_cov'])
        summary_df = pd.concat([summary_df,
            pd.DataFrame({'chrom':chrom, 'start':start, 'stop':stop, 'mean_cov':mean_cov}, index = [0])],
            ignore_index = True)
        # Remove very small phase sets. These tend to be single-SNP phase sets
        # and can result from erroneous SNPs.
        cov = cov.groupby('phase_set').filter(lambda x: len(x) > 1000)
        sum_df = cov.groupby(['bin', 'phase_set'])[ext_cols].mean().reset_index()
        sum_df['low'] = sum_df.total_cov < 0.8 * mean_cov
        sum_df['low_hap0'] = np.logical_and(sum_df.total_cov < mean_cov,
            sum_df.cov_q30_hap0 < 0.8 * sum_df.cov_q30_hap1)
        sum_df['low_hap1'] = np.logical_and(sum_df.total_cov < mean_cov,
            sum_df.cov_q30_hap1 < 0.8 * sum_df.cov_q30_hap0)

        if not sum_df.empty:
            any_low = np.logical_or(sum_df.low,
                                    np.logical_or(sum_df.low_hap1, sum_df.low_hap0))

            bins = np.array(sum_df['bin'])
            bins = np.concatenate([bins, [np.max(bins) + 1]])
            pos = 0
            # Get runs of 0s and 1s in any_low
            for bit, group in groupby(any_low):
                group_size = len(list(group))
                group_start = bins[pos] * args.bin_size
                group_stop = bins[pos + group_size] * args.bin_size
                region_len = group_stop - group_start
                if bit and region_len >= args.min_len:
                    out_loci.append((chrom, max(0, group_start - args.bin_size),
                                     group_start + args.bin_size,
                                     chrom, max(0, group_stop - args.bin_size),
                                     group_stop + args.bin_size))
                pos += group_size

    with open(outs.loci, 'w') as f:
        cPickle.dump(out_loci, f)

    summary_df.to_csv(outs.cov_summary, sep = '\t', header = True, index = False)
