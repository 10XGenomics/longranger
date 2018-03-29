#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# convert fragments.h5 file into csv with selected columns
#

import tenkit.hdf5

__MRO__ = """
stage CONVERT_HDF5_TO_CSV(
    in  h5       fragments_h5,
    out csv      fragments,
    src py       "stages/convert_hdf5_to_csv",
) split using () 
"""

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    if not args.fragments_h5:
        outs.fragments = None
        return

    data = tenkit.hdf5.read_data_frame(args.fragments_h5, ["molecule_id","start_pos", "end_pos", "chrom", "bc"])
    data.to_csv(outs.fragments)
