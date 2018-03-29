#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make sparse matrices of bcs X genomic windows
#

import tenkit.pandas as pd
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.hdf5 as tk_hdf5
import numpy as np

__MRO__ = """
stage GET_READS(
    in  bam  possorted_bam,
    in  int  mapq,
    in  h5   coverage,
    out h5   reads,
    src py   "stages/structvars/get_reads",
) split using (
    in string chrom,
    in int    start_pos,
    in int    end_pos,
)
"""

CHUNK_SIZE = 50000000

def split(args):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    chroms = in_bam.references
    lengths = in_bam.lengths

    chunks = []
    for chrom, length in zip(chroms, lengths):
        nchunks = int(np.ceil(length / float(CHUNK_SIZE)))
        for c in range(nchunks):
            chunk_start = c * CHUNK_SIZE
            chunk_stop = min((c + 1) * CHUNK_SIZE, length)
            chunks.append({'chrom':chrom, 'start_pos':chunk_start, 
                'end_pos':chunk_stop, '__mem_gb':16})

    return {'chunks': chunks}


def main(args, outs):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    chrom = args.chrom

    poses = []
    mol_qs = []
    bcs = []

    for read in in_bam.fetch(str(chrom), int(args.start_pos), int(args.end_pos)):
        if not read.is_secondary and not read.is_duplicate and read.is_read1 and \
            not read.is_unmapped and read.mapq >= args.mapq:
                poses.append(read.pos)
                mol_qs.append(tk_io.get_read_molecule_conf(read))
                bcs.append(tk_io.get_read_barcode(read))
    ret_df = pd.DataFrame({'chrom':chrom, 'pos':poses, 'bc':bcs, 'mol_qual':mol_qs})
                
    if len(ret_df) > 0:
        start_pos = poses[0]
        end_pos = poses[-1]
        cov_df = tk_hdf5.read_data_frame_indexed(args.coverage, [(chrom, start_pos, end_pos + 1)])

        # Boolean array with length equal to the range of positions in ret_df
        on_target = np.zeros((end_pos - start_pos + 1,), dtype = np.bool)
        on_target[cov_df.pos - start_pos] = True
        cum_on_target = np.cumsum(on_target)

        ret_df['on_target'] = on_target[ret_df.pos - start_pos]
        # Note that fraction is set to 1 if there are no bases
        frac_on_target = np.ones((len(ret_df),)) * on_target[0]
        for i, p in enumerate(poses):
            if i > 0:
                nbp = float(p - poses[i - 1] - 1)
                if nbp > 0:
                    frac_on_target[i] = (cum_on_target[p - start_pos] -
                                         cum_on_target[poses[i - 1] - start_pos] -
                                         int(on_target[p - start_pos])) / nbp
                else:
                    frac_on_target[i] = float(on_target[p - start_pos])
        ret_df['frac_on_target'] = frac_on_target
    else:
        ret_df['on_target'] = False
        ret_df['frac_on_target'] = 0.0
    tk_hdf5.write_data_frame(outs.reads, ret_df)


def join(args, outs, chunk_defs, chunk_outs):
    in_files = [out.reads for out in chunk_outs]
    tk_hdf5.combine_data_frame_files(outs.reads, in_files)
    tk_hdf5.create_tabix_index(outs.reads, 'chrom', 'pos', 'pos')
