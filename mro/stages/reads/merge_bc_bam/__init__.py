#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make a barcoded-sorted bam file from a bunch of sorted buckets.
#

import subprocess
import os

__MRO__ = """
stage MERGE_BC_BAM(
    in  bam      bc_sorted_bams,
    out bam      bc_sorted_bam,
    src py       "stages/reads/merge_bc_bam",
) split using (
    in string dummy,
)
"""

def split(args):
   return { 'chunks': [], 'join': {'__mem_gb': 6, '__threads': 4}}

def main(args, outs):
    pass

def join(args, outs, chunk_defs, chunk_outs):
    filtered_inputs = [a for a in args.bc_sorted_bams if os.path.isfile(a)]
    #tk_bam.concatenate(outs.bc_sorted_bam, filtered_inputs)
    args = ['samtools', 'merge', '-@', str(args.__threads), '-t', 'BX', '-n', outs.bc_sorted_bam]
    args.extend(filtered_inputs)
    subprocess.check_call(args)
