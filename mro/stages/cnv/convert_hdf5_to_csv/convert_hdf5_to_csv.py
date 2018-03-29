# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.hdf5
#import pandas as p
import sys

data = tenkit.hdf5.read_data_frame(sys.argv[1], ["molecule_id","start_pos", "end_pos", "chrom"])

data.to_csv(sys.argv[1].replace("h5","csv"))
