#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Finds pairs of genomic windows with significant barcode overlaps.
#

import os.path
import numpy as np
import tenkit.pandas as pd
import tenkit.reference as tk_reference
import tenkit.bam as tk_bam
import longranger.sv.io as tk_sv_io
from tenkit.constants import PARALLEL_LOCUS_SIZE
from tenkit.chunk_utils import generate_chrom_loci, pack_loci
import subprocess

import martian

__MRO__ = """
stage GET_DEL_CANDIDATES(
    in  bam    possorted_bam,
    in  string sex,
    in  string reference_path,
    in  string barcode_whitelist,
    in  int    min_mapq,
    in  int    min_del_len,
    in  int    max_del_len,
    in  int    min_bad_region,
    in  float  transition_prob,
    in  float  het_read_prob,
    out bedpe  del_candidates,
    src py     "stages/structvars/get_del_candidates",
) split using (
    in  string[] loci,
)
"""

def split(args):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)

    if not args.reference_path is None and os.path.exists(tk_reference.get_primary_contigs(args.reference_path)):
        with open(tk_reference.get_primary_contigs(args.reference_path), 'r') as f:
            primary_contigs = set([line.strip() for line in f.readlines()])
    else:
        # Default is to include all contigs
        primary_contigs = set(in_bam.references)

    all_loci = []
    for (chrom_name, chrom_size) in zip(in_bam.references, in_bam.lengths):
        if chrom_name in primary_contigs and not (args.sex in ['f', 'female'] and chrom_name in ['Y', 'chrY']):
            all_loci.extend(generate_chrom_loci(None, chrom_name, chrom_size, PARALLEL_LOCUS_SIZE))
    in_bam.close()

    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci, '__mem_gb':8} for loci in locus_sets]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    out_df = None
    for chunk in chunk_outs:
        tmp_df = tk_sv_io.read_sv_bedpe_to_df(chunk.del_candidates)
        out_df = pd.concat([out_df, tmp_df], ignore_index=True)

    out_df['name'] = np.arange(len(out_df))
    tk_sv_io.write_sv_df_to_bedpe(out_df, outs.del_candidates)



def main(args, outs):

    if args.barcode_whitelist is None:
        # write empty dataframe
        tk_sv_io.write_sv_df_to_bedpe(None, outs.del_candidates)
        martian.log_info('Data seem un-barcoded. No deletion candidates will be computed.')
        return


    if True:
        '''
          pvc [options] call-one <fasta> <bam> <locus>
          pvc [options] call-bed -o <out> <fasta> <bam> <bed> [<which>]
          pvc [options] bam-svs <out> <bam>
          pvc [options] cands <bam> <locus> <out>
          pvc (-h | --help)
          pvc --version

        Options:
          --min-size=<m>       Mininum event size
          --min-kmer-obs=<k>   Minimum number of kmer observations
          -h --help            Show this screen.
          --version            Show version.
          --trace              Trace logging
          -d --debug           Debug logging
        '''

        for locus in args.loci:
            tmp_file = "tmp.bedpe"

            min_detect_size = 25
            pvc_args = ['pvc', '--het-read-prob=%f' % args.het_read_prob, '--min-size=%d' % min_detect_size, 'cands', args.possorted_bam, locus, tmp_file]
            subprocess.check_call(pvc_args)
            subprocess.check_call('cat %s >> %s' % (tmp_file, outs.del_candidates), shell=True)
