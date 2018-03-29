#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Summarize genome coverage from BAM file
#

#import json
import tenkit.pandas as p
import tenkit.hdf5
import tenkit.coverage
import tenkit.reference
import tenkit.stats as tk_stats
import tenkit.bio_io as tk_io
import h5py
#import tenkit.bio_io
#import tenkit.bam as tk_bam
#import tenkit.dict_utils as tk_dict_utils
#import tenkit.regions as tk_regions
#from tenkit.constants import COVERAGE_TRIM_TAIL

#import martian

__MRO__ = """
stage REPORT_COVERAGE(
    in  string   reference_path,
    in  map      baits_file_map,
    in  h5       coverage,
    out csv      bait_csv,
    src py       "stages/reporter/report_baits",
) split using (
    in  string   tag,
    in  string   bait_file,
)
"""


def split(args):

    if args.baits_file_map is None:
        baits_path_map = None
    else:
        baits_path_map = args.baits_file_map

    # Bail if it's not a targeted run
    if baits_path_map is None:
        return {'chunks':[{'tag': None, 'bait_file': None}]}

    # Load pull-down targets
    chunk_defs =  []
    # loop over bait files
    for (tag, bait_file) in baits_path_map.iteritems():
        chunk_defs.append({'tag':tag, 'bait_file':bait_file})

    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    # Combine the coverage hdf5 files
    frame =  p.DataFrame()
    list_ = []
    if args.baits_file_map and outs.bait_csv:
        in_files = [out.bait_csv  for (cdef, out) in zip(chunk_defs, chunk_outs)]
        for file_ in in_files:
            df = p.read_csv(file_, index_col=None, header=0)
            list_.append(df)
        frame = p.concat(list_)

        # write csv
        frame.to_csv(outs.bait_csv)
    else:
        outs.target_coverage = None


def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    if args.coverage is None or args.bait_file is None:
        outs.bait_csv = None
        return

    f = h5py.File(args.coverage)
    has_basic = ('coverage_deduped' in f) and ('mapq30_coverage_deduped' in f)
    has_subsampling = ('coverage_subsampled' in f) and  ('coverage_deduped_subsampled' in f)
    f.close()
    if not has_basic: return

    fasta = tenkit.reference.open_reference(args.reference_path)

    df = p.DataFrame()
    coverage_reader = tenkit.hdf5.DataFrameReader(args.coverage)

    #regs = tenkit.bio_io.get_target_regions_dict(args.bait_file)
    #for chrom in regs:
    #    for start, end in regs[chrom]:

    bedIte = tk_io.get_bed_iterator(args.bait_file)
    for chrom, start, end in bedIte:
        if has_subsampling:
            coverage = coverage_reader.query([chrom, start, end], query_cols=['coverage_deduped','coverage_deduped_subsampled', 'coverage_subsampled','mapq30_coverage_deduped'] ,coords = False)
            mean_cov = coverage.mean()
            gc= get_gc(chrom, (start,end), fasta)
            #df = df.append({'name': 'Zed', 'age': 9, 'height': 2}, ignore_index=True)
            df = df.append({
                'chrom':chrom,
                'start':start,
                'end':end,
                'tag':args.tag,
                'coverage_deduped':mean_cov['coverage_deduped'],
                'coverage_deduped_subsampled':mean_cov['coverage_deduped_subsampled'],
                'coverage_subsampled':mean_cov['coverage_subsampled'],
                'mapq30_coverage_deduped':mean_cov['mapq30_coverage_deduped'],
                'gc':gc
                }, ignore_index=True)
        else:
            coverage = coverage_reader.query([chrom, start, end], query_cols=['coverage_deduped', 'mapq30_coverage_deduped'] ,coords = False)
            mean_cov = coverage.mean()
            gc= get_gc(chrom, (start,end), fasta)
            #df = df.append({'name': 'Zed', 'age': 9, 'height': 2}, ignore_index=True)
            df = df.append({
                'chrom':chrom,
                'start':start,
                'end':end,
                'tag':args.tag,
                'coverage_deduped':mean_cov['coverage_deduped'],
                'mapq30_coverage_deduped':mean_cov['mapq30_coverage_deduped'],
                'gc':gc
                }, ignore_index=True)


    df.to_csv(outs.bait_csv)

def get_gc(chrom, interval, fasta):
    (start, end) = interval
    seq = fasta[chrom][start:end]
    isGC = [x in ['G','g','C','c'] for x in seq if x not in ['N', 'n']]
    return tk_stats.robust_divide(sum(isGC)*1.0, len(isGC))
