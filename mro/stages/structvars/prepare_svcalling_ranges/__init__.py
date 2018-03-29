#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import os.path
import numpy as np
import cPickle
import tenkit.pandas as pd
import tenkit.bam as tk_bam
import tenkit.hdf5 as tk_hdf5
import tenkit.regions as tk_regions
import tenkit.reference as tk_reference
import martian

__MRO__ = """
stage PREPARE_SVCALLING_RANGES(
    in  bam possorted_bam,
    in  string reference_path,
    in  h5  coverage,
    in  csv coverage_csv,
    in  int min_region_len,
    in  int region_pad,
    in  string sex,
    in  int  max_bc_coverage,
    out pickle ranges,
    src py    "stages/structvars/prepare_svcalling_ranges",
) split using (
    in  string chrom,
    in  int    size,
)
"""

MAX_BASEPAIRS = 100000000
DEFAULT_MAX_COV = 300

def split(args):
    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    chroms = in_bam.references
    chrom_lens = in_bam.lengths

    if not args.reference_path is None and os.path.exists(tk_reference.get_primary_contigs(args.reference_path)):
        with open(tk_reference.get_primary_contigs(args.reference_path), 'r') as f:
            primary_contigs = set([line.strip() for line in f.readlines()])
    else:
        # Default is to include all contigs
        primary_contigs = set(chroms)

    chunks = []
    for c, s in zip(chroms, chrom_lens):
        if not c in primary_contigs:
            continue
        if args.sex is not None and args.sex.lower() in ['f', 'female'] and c in ['Y', 'chrY']:
            continue
        chunks.append({'chrom':c, 'size':s, '__mem_gb':8})
    return {'chunks':chunks, 'join':{'__mem_gb':8}}


def main(args, outs):

    if not os.path.isfile(args.coverage_csv):
        max_cov = DEFAULT_MAX_COV

    cov_summary_df = pd.read_csv(args.coverage_csv, index_col = 0)
    if len(cov_summary_df) == 0:
        max_cov = DEFAULT_MAX_COV
    max_cov = cov_summary_df.iloc[-1].coverage

    if args.max_bc_coverage is None:
        max_bc_coverage = np.inf
    else:
        max_bc_coverage = args.max_bc_coverage
        
    ranges = []
    nchunks = int(np.ceil(args.size / float(MAX_BASEPAIRS)))
    for i in range(nchunks):
        start = i * MAX_BASEPAIRS
        # chunks will be slightly overlapping
        stop = min(args.size, (i + 1) * MAX_BASEPAIRS + args.min_region_len)

        # Do a test read to see if the file contains the bc_coverage column. This adds a few seconds...
        cov_df = tk_hdf5.read_data_frame_indexed(args.coverage, queries = [(args.chrom, 0, 1)])
        if not 'bc_coverage' in cov_df.columns:
            martian.log_info('Data seem un-barcoded. SV-calling ranges will be empty.')
        else:
            cov_df = tk_hdf5.read_data_frame_indexed(args.coverage,
                                                     queries = [(args.chrom, start, stop)],
                                                     query_cols = ['bc_coverage', 'coverage_deduped'])
            cov_df = np.array(cov_df[np.logical_and(cov_df.bc_coverage < max_bc_coverage,
                                                    cov_df.coverage_deduped < max_cov)].pos, dtype = np.int)
            new_ranges = get_good_ranges(cov_df, args.min_region_len, args.region_pad)
            ranges.extend(new_ranges)

    with open(outs.ranges, 'w') as f:
        cPickle.dump({'chrom':args.chrom, 'ranges':ranges}, f)


def join(args, outs, chunk_defs, chunk_outs):
    out_ranges = {}
    for chunk in chunk_outs:
        with open(chunk.ranges, 'r') as f:
            ranges = cPickle.load(f)
        if len(ranges['ranges']) > 0:
            new_ranges = [(a, b) for (a, b) in ranges['ranges'] if b - a > args.min_region_len]
            # This will merge overlapping regions
            out_ranges[ranges['chrom']] = tk_regions.Regions(new_ranges)
    with open(outs.ranges, 'w') as f:
        cPickle.dump(out_ranges, f)


def get_good_ranges(good_poses, min_gap = 100, pad = 0):
    """Converts a list of positions to a list of ranges.
    Args:
    - min_gap: minimum length of output region.
    Return value:
    A list of tuples (start, stop) in BED-like format (start inclusive,
    end exclusive).
    """

    if len(good_poses) == 0:
        return []

    start = good_poses[0]
    good_regions = []
    for i, p in enumerate(good_poses):
        # If the gap is too long, emit the last good region.
        if i > 0 and p - good_poses[i - 1] - 1 > 0:
            good_regions.append((max(0, start - pad),
                                 max(0, good_poses[i - 1] + 1 + pad)))
            start = p
    good_regions.append((max(0, start - pad),
                         max(0, good_poses[-1] + 1 + pad)))
    good_regions = [r for r in good_regions if r[1] - r[0] >= min_gap]
    return good_regions


def get_good_ranges_old(bad_poses, chrom_len, min_gap = 100):
    """Gets a list of "bad" positions on the chromosome and returns a list of "good" ranges.
    Args:
    - bad_poses: array of bad positions.
    - chrom_len: Length of chromosome.
    - min_gap: minimum length of a "bad region". Smaller bad regions will be ignored.
    This is also the minimum length of a good region. A good region with length < min_gap
    surrounded by bad regions will be merged with the surrounding bad regions.
    - win_around_bad: extra buffer to add around the bad regions

    Return value:
    A list of tuples (start, stop) with the good regions on the chromosome. These are BED-like
    (start pos inclusive, stop pos exclusive).
    """

    # Get the ranges of bad regions
    bad_regions = []
    if len(bad_poses) > 0:
        start = bad_poses[0]
        stop = start
        for i, p in enumerate(bad_poses):
            if p > chrom_len:
                break
            if p - stop <= min_gap:
                # Keep merging into the previous bad region
                stop = p
            else:
                if stop - start + 1 >= min_gap:
                    # bad region long enough to not ignore
                    bad_regions.append((start, stop + 1))
                start = p
                stop = p
        if stop - start + 1 >= min_gap:
            bad_regions.append((start, stop + 1))

    good_regions = []
    if len(bad_regions) > 0:
        if bad_regions[0][0] > 0:
            # Chromosome doesn't start with a bad region
            good_regions.append((0, bad_regions[0][0]))
        if len(bad_regions) > 1:
            for i, b in enumerate(bad_regions[1:]):
                # bad_regions[i] will be the region before b
                # Add to the good regions the region between the end of the previous
                # bad region and the beginning of b.
                # No need to test if the middle region is long enough
                good_regions.append((bad_regions[i][1], b[0]))

        if chrom_len - bad_regions[-1][1] + 1 >= min_gap:
            good_regions.append((bad_regions[-1][1], chrom_len))
    elif chrom_len >= min_gap:
        good_regions = [(0, chrom_len)]
    return good_regions
