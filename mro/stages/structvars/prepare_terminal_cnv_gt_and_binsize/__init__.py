#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

__MRO__ = """
stage PREPARE_TERMINAL_CNV_GT_AND_BINSIZE(
    in  bedpe   gt_variant,
    in  int     min_cnv_len,
    in  string  reference_path,
    out bedpe   gt_variant_large_cnv,
)
"""

def main(args, outs):
    if args.gt_variant and args.min_cnv_len :
        with open(args.gt_variant) as gt_in, open(outs.gt_variant_large_cnv, "w") as gt_out:
            for line in gt_in.readlines():
                if line.startswith("#"):
                    gt_out.write(line)
                else:
                    infos  = line.strip().split("\t")
                    chrom1 = infos[0]
                    chrom2 = infos[3]
                    start1 = int(infos[1])
                    end1   = int(infos[4])
                    extra  = infos[11].split(";")
                    extra_fields = { x.split("=")[0]:x.split("=")[1] for x in extra }
                    type_ = extra_fields.get("TYPE", None)

                    #print chrom1, chrom2, start1, end1, extra, type_
                    if chrom1 == chrom2 and (end1 - start1) >= int(args.min_cnv_len) and type_ in ["DEL", "DUP"]:
                        gt_out.write(line)

    if args.gt_variant is None:
        outs.gt_variant_large_cnv = None

