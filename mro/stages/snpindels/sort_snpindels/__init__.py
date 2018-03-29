#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort VCF file
#
import os
import os.path
import tenkit.tabix as tk_tabix
import tenkit.bio_io as tk_io
import martian

__MRO__ = """
stage SORT_SNPINDELS(
    in  vcf    input    "input vcf file",
    out vcf.gz          "output compressed vcf",
    src py     "stages/snpindels/sort_snpindels",
)
"""

def main(args, outs):
    main_sort_variants(args, outs)

def main_sort_variants(args, outs):
    if args.input is None or args.input == [None]:
        outs.default = None
    else:
        outs.coerce_strings()

        # List of inputs
        sort_filename = outs.default[0:(len(outs.default)-3)]
        if type(args.input) == type([]):
            files = [f for f in args.input if os.path.isfile(f)]
            if len(files) == 0:
                outs.defaut = None
                return
            cat_filename = outs.default[:-3]
            tk_io.combine_vcfs(cat_filename, args.input)

        # Single input
        else:
            if not os.path.exists(args.input):
                outs.default = None
                return
            tk_tabix.sort_vcf(args.input, sort_filename)
            tk_tabix.index_vcf(sort_filename)
