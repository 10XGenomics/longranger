#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Coallate metrics for monitoring LandLord coverage spikes
#

import tenkit.pandas as p
import tenkit.hdf5
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.stats as tk_stats
import tenkit.reference
import tenkit.chunk_utils as tk_chunks
import numpy as np

__MRO__ = """
stage COMPUTE_COVERAGE_SPIKES(
    in  h5     coverage,
    in  csv    cov_hist,
    in  bam    bam_infile,
    in  string reference_path,
    out bed    spikes,
) split using (
    in  string locus,
)
"""


def split(args):
    input_bam = tk_bam.create_bam_infile(args.bam_infile)

    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    cov_hist = p.read_csv(args.cov_hist)
    weighted_count = cov_hist.counts[1:] * cov_hist.coverage[1:]
    mean_pos_cov = tk_stats.robust_divide(weighted_count.sum(), cov_hist.counts[1:].sum())

    primary_contigs = tenkit.reference.load_primary_contigs(args.reference_path) - {'chrM', 'chrY'}

    loci = tk_chunks.chunk_by_locus(chroms, chrom_lengths, tenkit.constants.PARALLEL_LOCUS_SIZE * 2,
        contig_whitelist = primary_contigs, extra_args = {'mean': mean_pos_cov})

    # Handle empty case
    if len(loci) == 0:
        loci = [{'locus': None, 'mean': None}]

    return {'chunks': loci, 'join': {'__mem_gb': 12.0}}

def main(args, outs):
    if args.locus is None:
        outs.spikes = None
        return

    if args.mean is None or (not args.mean):
        outs.spikes = None
        return

    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    chrom = str(chrom)
    cov = tenkit.hdf5.read_data_frame_indexed(args.coverage, [(chrom, start, stop)], query_cols=['pos', 'coverage'])
    #cov.coverage = tk_stats.robust_divide(cov.coverage, args.mean)
    cov.coverage /= max(1.0, args.mean) # args.mean should not be zero unless due to other issue, like missing args.cov_hist
    cov = cov[cov.coverage > 10]
    spikes_pos = cov["pos"].values
    breaks = list(np.where(np.diff(spikes_pos) !=1)[0]+1)

    if len(breaks) > 0:
        starts = [ spikes_pos[b] for b in [0]+breaks]
        ends = [ spikes_pos[b-1]+1 for b in breaks+[len(spikes_pos)]]
    else:
        starts = []
        ends = []

    with open(outs.spikes, "w") as fout:
        for s, e in zip(starts, ends):
            fout.write(chrom+"\t"+str(s)+"\t"+str(e)+"\n")

def join(args, outs, chunk_defs, chunk_outs):
    in_files = [out.spikes for out in chunk_outs if out.spikes is not None]

    with open(outs.spikes, "w") as fout:
        for chunk_in in in_files:
            with open(chunk_in) as fin:
                for line in fin:
                    fout.write(line)
