#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Finds pairs of genomic windows with significant barcode overlaps.
#

import numpy as np
import tenkit.pandas as pd
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
from tenkit.constants import PARALLEL_LOCUS_SIZE
from tenkit.chunk_utils import generate_chrom_loci, pack_loci
import tenkit.hdf5 as tk_hdf5
from tenkit.coverage import get_hap_coverage

__MRO__ = """
stage GET_HAP_COVERAGE(
    in  bam possorted_bam,
    in  bed targets,
    in  h5  phase_set_h5,
    out h5  hap_coverage,
    src py  "stages/structvars/get_hap_coverage",
) split using (
    in  string[] loci,
)
"""

COV_QUALS = [30, 60]


def split(args):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    # Load pull-down targets
    if args.targets is None:
        target_regions = None
    else:
        with open(args.targets, 'r') as f:
            target_regions = tk_io.get_target_regions(f)

    all_loci = []
    for (chrom_name, chrom_size) in zip(in_bam.references, in_bam.lengths):
        all_loci.extend(generate_chrom_loci(target_regions, chrom_name, chrom_size, PARALLEL_LOCUS_SIZE))
    in_bam.close()

    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci, '__mem_gb':16} for loci in locus_sets]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    in_files = [out.hap_coverage for out in chunk_outs]
    tk_hdf5.combine_data_frame_files(outs.hap_coverage, in_files)

    # tabix index the coverage file
    tk_hdf5.create_tabix_index(outs.hap_coverage, 'chrom', 'pos', 'pos')


def main(args, outs):

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    for (chrom, start, stop) in (tk_io.get_locus_info(l) for l in args.loci):
        cov_df = get_hap_coverage(in_bam, args.phase_set_h5, chrom, start, stop, cov_quals=COV_QUALS)
        tk_hdf5.append_data_frame(outs.hap_coverage, cov_df)

    in_bam.close()
