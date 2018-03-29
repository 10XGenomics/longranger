#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Copy BAM file to output
#
import shutil

__MRO__ = """
stage COPY_BAM(
    in  bam     input,
    in  bam.bai input_index,
    out bam     possorted_bam,
    out bam.bai possorted_bam_index,
    src py      "stages/reads/copy_bam",
)
"""

def main(args, outs):

    shutil.copyfile(args.input, outs.possorted_bam)

    outs.possorted_bam_index = outs.possorted_bam + ".bai"
    shutil.copyfile(args.input_index, outs.possorted_bam_index)
