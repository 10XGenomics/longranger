#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# extract statistics about single partition behaviour
import tenkit.pandas as p
import numpy as np
import scipy.stats
import itertools
import os
import os.path
import math

import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.stats as tk_stats
import tenkit.seq as tk_seq
from tenkit.constants import FRAG_LEN_HIST_BIN_SIZE
from tenkit.constants import FRAG_LEN_HIST_BIN_SIZE_FOR_LW_MODE
import tenkit.hdf5
import tenkit.safe_json

__MRO__ = """
stage REPORT_SINGLE_PARTITION(
    in  bam    input,
    in  string barcode_whitelist,
    in  bed    targets_file,
    in  int    read_link_distance,
    out json   single_partition,
    out json   fragment_size,
    out h5     fragments,
    out h5     barcodes,
    out json   barcode_histogram,
    src py     "stages/reporter/report_single_partition",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
"""

NUM_BCS_LOADING_ESTIMATE = 1000
MAPQ_THRESHOLD = 30

def chunk_split_func(r):
    return tk_io.get_read_barcode(r)

def split(args):
    if args.input is None or args.barcode_whitelist is None:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}

    # Some R&D bc sets have very small diversity -- don't run on them
    barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)
    if len(barcode_whitelist) < 100:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}

    min_chunks = 5
    if len(barcode_whitelist) > 1e6:
        min_chunks = 10

    bam_in = tk_bam.create_bam_infile(args.input)
    chunks = tk_bam.chunk_bam_records(bam_in, chunk_split_func, chunk_size_gb=8.0, min_chunks=min_chunks)
    for c in chunks:
        c['__mem_gb'] = 7.0

    return {'chunks': chunks, 'join': {'__mem_gb': 32.0}}

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    main_report_single_partition(args, outs)


def high_level_stats(prefix, fragment_df, bc_df):
    stats = {}

    if fragment_df is not None:
        stats["num_bcs"] = len(bc_df)
        stats["mean_reads_per_barcode"] = bc_df.bc_num_reads.mean()
        stats["stddev_reads_per_barcode"] = bc_df.bc_num_reads.std()
        stats["cv_reads_per_barcode"] = stats["stddev_reads_per_barcode"]/stats["mean_reads_per_barcode"]

        stats['mean_reads_per_fragment'] = bc_df.bc_mean_reads_per_fragment.mean()
        stats['mean_reads_per_fragment_pf'] = np.average(bc_df.bc_mean_reads_per_fragment.values, weights=bc_df.bc_num_fragments)
        stats['stddev_reads_per_fragment'] = bc_df.bc_mean_reads_per_fragment.std()
        stats['linked_fragment_fraction'] = bc_df.bc_linked_fragment_fraction.mean()
        stats['linked_read_fraction'] = bc_df.bc_linked_read_fraction.mean()
        stats['n50_reads_per_fragment'] = tk_stats.N50(fragment_df.num_reads.values)

        stats["mean_fragment_length"] = fragment_df.est_len.mean()
        stats["stddev_fragment_length"] = fragment_df.est_len.std()
        stats["median_fragment_length"] = fragment_df.est_len.median()

        stats["fraction_fragments_gt_100kb"] = (fragment_df.est_len > 100000).mean()
        stats["fraction_fragments_gt_20kb"] = (fragment_df.est_len > 20000).mean()
        stats["fraction_fragments_lt_5kb"] = (fragment_df.est_len < 5000).mean()

        est_len = fragment_df.est_len.values

        def weighted_avg_and_std(values, weights):
            """
            Return the weighted average and standard deviation.

            values, weights -- Numpy ndarrays with the same shape.
            """
            average = np.average(values, weights=weights)
            variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
            return (average, math.sqrt(variance))

        (lw_mean, lw_std) = weighted_avg_and_std(est_len, est_len)
        
        fragment_df['len_bin'] = np.floor_divide(fragment_df.est_len.values, FRAG_LEN_HIST_BIN_SIZE_FOR_LW_MODE).astype(int) * FRAG_LEN_HIST_BIN_SIZE_FOR_LW_MODE
        multi_read_frags = fragment_df[fragment_df.num_reads > 1]
        len_bins = multi_read_frags.groupby(['len_bin']).apply(len)
        lw_len_bins = len_bins * multi_read_frags.groupby(['len_bin']).est_len.apply(np.average)                                                                                                          
        lw_scale_factor = tk_stats.robust_divide(float(sum(multi_read_frags.est_len)), float(len(multi_read_frags.est_len)))                                                                                                       
        lw_len_mode = None
        if len(lw_len_bins)>0:
            lw_len_hist = {k:(tk_stats.robust_divide(float(v), lw_scale_factor)) for (k,v) in lw_len_bins.iteritems()}                                                                                                                 
            lw_len_mode = max(lw_len_hist, key=lambda key: lw_len_hist[key]) + (FRAG_LEN_HIST_BIN_SIZE_FOR_LW_MODE/2)
        del multi_read_frags
        stats["lw_mode_of_binned_fragment_length"] = lw_len_mode
        
        stats["lw_mean_fragment_length"] = lw_mean
        stats["lw_stddev_fragment_length"] = lw_std
        stats["lw_median_fragment_length"] = tk_stats.N50(est_len)

        stats["lw_fraction_fragments_gt_100kb"] = float(est_len[est_len > 100000].sum()) / est_len.sum()
        stats["lw_fraction_fragments_gt_20kb"] = float(est_len[est_len > 20000].sum()) / est_len.sum()
        stats["lw_fraction_fragments_lt_5kb"] =  float(est_len[est_len < 5000].sum()) / est_len.sum()

        stats["mean_fragments_per_barcode"] = bc_df.bc_num_fragments.mean()
        stats["stddev_fragments_per_barcode"] = bc_df.bc_num_fragments.std()
        stats["cv_fragments_per_barcode"] = stats["stddev_fragments_per_barcode"]/stats["mean_fragments_per_barcode"]

        stats["max_fragments_per_barcode"] = bc_df.bc_num_fragments.max()

    else:
        stats["num_bcs"] = None
        stats["mean_reads_per_barcode"] = None
        stats["stddev_reads_per_barcode"] = None
        stats["cv_reads_per_barcode"] = None

        stats['mean_reads_per_fragment'] = None
        stats['mean_reads_per_fragment_pf'] = None
        stats['stddev_reads_per_fragment'] = None
        stats['linked_fragment_fraction'] = None
        stats['linked_read_fraction'] = None
        stats['n50_reads_per_fragment'] = None

        stats["mean_fragment_length"] = None
        stats["stddev_fragment_length"] = None
        stats["median_fragment_length"] = None

        stats["fraction_fragments_gt_100kb"] = None
        stats["fraction_fragments_gt_20kb"] = None
        stats["fraction_fragments_lt_5kb"] = None

        stats["lw_mode_of_binned_fragment_length"] = None

        stats["lw_mean_fragment_length"] = None
        stats["lw_stddev_fragment_length"] = None
        stats["lw_median_fragment_length"] = None

        stats["lw_fraction_fragments_gt_100kb"] = None
        stats["lw_fraction_fragments_gt_20kb"] = None
        stats["lw_fraction_fragments_lt_5kb"] = None

        stats["mean_fragments_per_barcode"] = None
        stats["stddev_fragments_per_barcode"] = None
        stats["cv_fragments_per_barcode"] = None

        stats["max_fragments_per_barcode"] = None

    final_stats = { (prefix + "_" + key):value for (key, value) in stats.iteritems() }
    return final_stats

def partition_chroms(chroms, lengths, max_length=200e6):
    '''
    Partition a list of chromosomes into sets of roughly equal length.
    '''
    allofem = [] # list of lists
    current_chroms = []
    current_length = 0
    for chrom, length in zip(chroms, lengths):
        # each group should have at least one, even if it's over max_length
        current_chroms.append(chrom)
        current_length += length
        if current_length >= max_length:
            allofem.append(current_chroms)
            current_chroms = []
            current_length = 0
    if len(current_chroms) > 0:
        allofem.append(current_chroms)
    return allofem

def join(args, outs, chunk_defs, chunk_outs):
    bam_in = tk_bam.create_bam_infile(args.input)
    # Combine fragment h5 files
    in_files = [out.fragments for out in chunk_outs if out.fragments and os.path.exists(out.fragments)]
    nfrags = 0

    chrom_partition = partition_chroms(bam_in.references, bam_in.lengths)

    if len(in_files) > 0:
        readers = [tenkit.hdf5.DataFrameReader(f) for f in in_files]

        for chrom_set in chrom_partition:
            chunks = [r.query_chroms(chrom_set) for r in readers]
            chunk = p.concat(chunks)
            chunk.sort(['chrom', 'start_pos', 'bc'], inplace=True)

            # Always save the BC as categorical
            chunk['bc'] = chunk['bc'].astype('category')

            chunk['molecule_id'] = np.arange(len(chunk), dtype=np.int32) + nfrags
            nfrags += len(chunk)
            tenkit.hdf5.append_data_frame(outs.fragments, chunk)

        for r in readers:
            r.close()

        tenkit.hdf5.create_tabix_index(outs.fragments, 'chrom', 'start_pos', 'end_pos')

    else:
        outs.fragments = None

    # Combine BC h5 files
    in_files = [out.barcodes for out in chunk_outs if out.barcodes and os.path.exists(out.barcodes)]
    if len(in_files) > 0:
        tenkit.hdf5.combine_data_frame_files(outs.barcodes, in_files)
    else:
        outs.barcodes = None

    summary = {}

    # Compute high-level BC summary metrics
    # Load BC data
    if outs.barcodes:
        bc_df = tenkit.hdf5.read_data_frame(outs.barcodes)
        fragment_df = tenkit.hdf5.read_data_frame(outs.fragments, query_cols=['bc', 'num_reads', 'est_len', 'chrom', 'start_pos'])

        bc_df.sort('bc_num_reads', inplace=True)

        # bin the bc counts and write a json histogram file
        n_reads = bc_df.bc_num_reads.values
        max_val = np.percentile(n_reads, 99.99) * 1.3
        min_val = n_reads.min()
        num_bins = 400
        step = math.ceil((max_val - min_val)/num_bins)
        bins = np.arange(min_val, max_val, step)
        (hist, edges) = np.histogram(n_reads, bins=bins)
        bc_count_hist = {int(edges[i]):hist[i] for i in range(len(bins)-1)}

        # Summarize properties of n50 and n90 BC set
        bc_df['cum_reads'] = np.cumsum(bc_df.bc_num_reads)
        n50_read_thresh = sum(bc_df.bc_num_reads) * 0.5
        n50_bcs = bc_df[bc_df.cum_reads > n50_read_thresh]
        n50_fra = fragment_df[fragment_df.bc.isin(n50_bcs.bc)]
        n50_stats = high_level_stats("n50", n50_fra, n50_bcs)
        del n50_fra

        n90_read_thresh = sum(bc_df.bc_num_reads) * 0.1
        n90_bcs = bc_df[bc_df.cum_reads > n90_read_thresh]
        n90_fra = fragment_df[fragment_df.bc.isin(n90_bcs.bc)]
        n90_stats = high_level_stats("n90", n90_fra, n90_bcs)
        del n90_fra

        for (k,v) in n50_stats.iteritems():
            summary[k] = v

        for (k,v) in n90_stats.iteritems():
            summary[k] = v

        # Generate a fragment length histogram
        fragment_df['len_bin'] = np.floor_divide(fragment_df.est_len.values, FRAG_LEN_HIST_BIN_SIZE).astype(int) * FRAG_LEN_HIST_BIN_SIZE

        multi_read_frags = fragment_df[fragment_df.num_reads > 1]
        len_bins = multi_read_frags.groupby(['len_bin']).apply(len)
        del multi_read_frags

        len_hist = {k:v for (k,v) in len_bins.iteritems()}

        # Write fragment length hist to json
        with open(outs.fragment_size, 'w') as fragment_size_file:
            tenkit.safe_json.dump_numpy(len_hist, fragment_size_file)

        # Estimate total DNA per partition by looking at hottest 1000 GEMs or GEMs w/ bc_mean_reads_per_fragment > 2, whichever is fewer
        hot_bcs = bc_df[np.logical_and(bc_df.bc_mean_reads_per_fragment > 2.0, bc_df.bc_num_reads > 25)]
        hot_bcs.sort('bc_mean_reads_per_fragment', inplace=True)
        if len(hot_bcs) > 50:
            hot_bcs = hot_bcs[-NUM_BCS_LOADING_ESTIMATE:]
            summary['estimated_dna_per_partition'] = round(scipy.stats.tmean(hot_bcs.bc_est_len, scipy.percentile(hot_bcs.bc_est_len, (1,99))))
        else:
            summary['estimated_dna_per_partition'] = None

        # Read-based effective diversity
        reads = bc_df.bc_num_reads.values
        sum_sq = (reads**2.0).sum()
        effective_diversity = tk_stats.robust_divide((reads.sum()**2.0), float(sum_sq))
        summary['effective_diversity_reads'] = effective_diversity

        # Fragment-based effective diversity
        fragments = bc_df.bc_num_fragments.values
        sum_sq = (fragments**2.0).sum()
        effective_diversity = tk_stats.robust_divide((fragments.sum()**2.0), float(sum_sq))
        summary['effective_diversity_fragments'] = effective_diversity

    else:
        # No fragment_size file emitted
        outs.fragment_size = None

        n50_stats = high_level_stats("n50", None, None)
        n90_stats = high_level_stats("n90", None, None)

        for (k,v) in n50_stats.iteritems():
            summary[k] = v

        for (k,v) in n90_stats.iteritems():
            summary[k] = v

        bc_count_hist = {}

        summary['estimated_dna_per_partition'] = None
        summary['effective_diversity_reads'] = None
        summary['effective_diversity_fragments'] = None

    with open(outs.barcode_histogram, 'w') as barcode_hist_file:
        tenkit.safe_json.dump_numpy(bc_count_hist, barcode_hist_file)

    # Write summary to json
    with open(outs.single_partition, 'w') as summary_file:
        tenkit.safe_json.dump_numpy(summary, summary_file, pretty=True)

def read_has_barcode(r):
    bc = tk_io.get_read_barcode(r)
    if bc is None:
        return False
    else:
        return True

def limit_groupby(gen, key_func, max_items):
    ''' A replacement for itertools.groupby that has a limit on the maximum group size.
        Elements beyond max_items are dropped '''
    group = []
    group_key = None

    for item in gen:
        new_key = key_func(item)

        if group_key is None:
            group_key = new_key
            group = [item]

        elif new_key == group_key and len(group) < max_items:
            group.append(item)

        else:
            yield (group_key, group)
            group_key = new_key
            group = [item]

    if group is not None:
        yield (group_key, group)



def main_report_single_partition(args, outs):
    # Bail out if there no valid barcodes
    if args.barcode_whitelist is None or args.input is None:
        outs.fragments = None
        return

    bam_in = tk_bam.create_bam_infile(args.input)

    if args.targets_file is None:
        target_regions = None
    else:
        target_regions = tk_io.get_target_regions(open(args.targets_file))

    # Bail out if we're on a small genome
    if sum(bam_in.lengths) < 3e6:
        outs.fragments = None
        return

    bam_chunk = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    # Skip reads without a barcode
    bam_chunk_filt = itertools.ifilter(read_has_barcode, bam_chunk)
    bc_read_iter = itertools.groupby(bam_chunk_filt, lambda x: tk_io.get_read_barcode(x))

    bc_data = (summarize_barcode(bc, list(reads), args.read_link_distance, bam_in.references, target_regions) for (bc, reads) in bc_read_iter)
    bc_data_filt = (x for x in bc_data if x is not None)

    frag_tbl, bc_tbl = make_output_dataframes(bc_data_filt)
    if frag_tbl is not None:
        # Sort and index fragment table, so that we can combine the fragments files per-chromosome to reduce memory consumption
        frag_tbl.sort(['chrom', 'start_pos'], inplace=True)
        tenkit.hdf5.write_data_frame(outs.fragments, frag_tbl)
        tenkit.hdf5.create_tabix_index(outs.fragments, 'chrom', 'start_pos', 'end_pos')
    if bc_tbl is not None:
        tenkit.hdf5.write_data_frame(outs.barcodes, bc_tbl)


def create_fragments(reads, max_distance):
    '''
    Scan the position sorted list of reads from a barcode, and estimate what
    genome fragments were in the partition
    '''

    frags = []
    current_frag = []
    last_q30 = 0

    for r in reads:
        if current_frag == []:
            # Never start a fragment at a low quality read
            if r.mapq >= MAPQ_THRESHOLD:
                current_frag.append(r)
                last_q30 = 0
        else:
            last_read = current_frag[-1]
            if r.tid != last_read.tid or r.pos - last_read.pos > max_distance:
                # Save this fragment and prepare for next
                # Trim the trailing low quality reads from the fragment
                frags.append([current_frag[f] for f in range(last_q30 + 1)])
                last_q30 = 0
                if r.mapq >= MAPQ_THRESHOLD:
                    current_frag = [r]
                else:
                    current_frag = []
            else:
                current_frag.append(r)
                if r.mapq >= MAPQ_THRESHOLD:
                    last_q30 = len(current_frag) - 1

    if current_frag != []:
        frags.append([current_frag[f] for f in range(last_q30 + 1)])

    return frags


def summarize_fragment(reads, references, targets):
    '''
    Summarize stats about reads from a single fragment.
    Assume reads are position sorted and from a single chrom
    '''

    n_se = len([r for r in reads])
    n = len(set(r.qname for r in reads))
    start = reads[0].pos
    end = max(reads[-1].aend, reads[-1].pos)

    stats = {}
    stats['num_reads'] = n
    stats['num_reads_se'] = n_se
    stats['chrom'] = references[reads[0].tid]
    stats['start_pos'] = start
    stats['end_pos'] = end

    # Observed length and adjusted length
    l = end - start

    # Note -- this is the estimator for the interval length when both ends are unknown
    # it doesn't work when there is only 1 observation
    if n > 1:
        adj_len = l * (n+1) / (n-1)
    else:
        # For single-read fragment, just assign a fragment length of 1000 for observed and adjusted
        adj_len = 1000
        l = 1000

    # Check if there are any MAPQ=30 reads in the fragment
    confident_mapping = any(r.mapq >= MAPQ_THRESHOLD for r in reads)

    # If targets are present, check if we intersect any
    if targets is not None:
        if targets.has_key(stats['chrom']):
            on_target = targets[stats['chrom']].overlaps_region(start, end)
        else:
            on_target = False
    else:
        on_target = True

    # Only include on-target fragments with at least 1 confidently mapped read
    if not (confident_mapping and on_target):
        return None

    stats['obs_len'] = l
    stats['est_len'] = adj_len
    return stats

def is_read_and_mate_unmapped(read):
    if read.is_paired:
        return read.is_unmapped and read.mate_is_unmapped
    else:
        return read.is_unmapped

def summarize_barcode(bc, all_reads, max_distance, references, targets):
    reads = filter(lambda r: not (r.is_duplicate or r.is_unmapped), all_reads)
    reads_fully_unmapped = filter(is_read_and_mate_unmapped, all_reads)

    npairs = len(set(r.qname for r in reads))
    npairs_unmapped = len(set(r.qname for r in reads_fully_unmapped))

    # SORT_BY_BC uses read.qname as a secondary sort.  Resort reads by position
    reads.sort(key=lambda r: (r.tid, r.pos))

    fragments = create_fragments(reads, max_distance)
    frag_stats = [summarize_fragment(f, references, targets) for f in fragments]
    frag_stats = [f for f in frag_stats if f is not None]

    # Don't emit 0-fragment BCs
    if len(frag_stats) == 0:
        return None

    stats = {}
    stats['bc'] = bc
    stats['bc_num_fragments'] = len(frag_stats)
    stats['bc_num_reads'] = npairs
    stats['bc_num_unmapped_reads'] = npairs_unmapped
    stats['bc_mean_reads_per_fragment'] = np.mean([s['num_reads'] for s in frag_stats])

    # Fraction of reads in fragments with >1 read
    linked_read_pairs = set()
    for frag in fragments:
        rds = set(r.qname for r in frag)
        if len(rds) > 1:
            linked_read_pairs.update(rds)

    stats['bc_linked_read_fraction'] = float(len(linked_read_pairs)) / npairs

    # Fraction of fragments with >1 read
    stats['bc_linked_fragment_fraction'] = float(sum(1.0 for frag in fragments if len(frag) > 1)) / len(fragments)
    stats['bc_est_len'] = sum(f['est_len'] for f in frag_stats)

    return (stats, frag_stats)

def make_output_dataframes(bcs_frags_in):
    fragments = []
    bcs = []

    bc_dfs = []
    fragment_dfs = []

    for (bc_stats, frags) in bcs_frags_in:
        # Denormalize selected bc columns into the fragments dataframe
        for (k,v) in bc_stats.items():
            if k in ['bc', 'bc_num_reads', 'bc_mean_reads_per_fragment', 'bc_est_len', 'bc_num_unmapped_reads']:
                for frag in frags:
                    frag[k] = v

        fragments.extend(frags)
        bcs.append(bc_stats)

        if len(fragments) > 2e6:
            (frag_df, bc_df) = make_df_chunk(fragments, bcs)

            fragment_dfs.append(frag_df)
            fragments = []

            bc_dfs.append(bc_df)
            bcs = []

    (frag_df, bc_df) = make_df_chunk(fragments, bcs)
    fragment_dfs.append(frag_df)
    bc_dfs.append(bc_df)

    frag_dfs = [ x for x in fragment_dfs if x is not None ]
    bc_dfs = [ x for x in bc_dfs if x is not None ]

    if len(bc_dfs) > 0:
        frag_df = p.concat(frag_dfs)
        bc_df = p.concat(bc_dfs)
    else:
        frag_df = None
        bc_df = None

    return (frag_df, bc_df)


def make_df_chunk(fragments, bcs):
    # No BC results -- will write an empty file
    if len(bcs) == 0:
        return (None,None)
    else:
        fragment_df = p.DataFrame(fragments)
        bc_df = p.DataFrame(bcs)

        # Set good types for the fragment data frame to reduce size
        fragment_df.start_pos = fragment_df.start_pos.astype(np.int32)
        fragment_df.end_pos   = fragment_df.end_pos.astype(np.int32)
        fragment_df.obs_len   = fragment_df.obs_len.astype(np.int32)
        fragment_df.est_len   = fragment_df.est_len.astype(np.int32)
        fragment_df.num_reads = fragment_df.num_reads.astype(np.int32)
        fragment_df.num_reads_se = fragment_df.num_reads_se.astype(np.int32)
        fragment_df.bc_num_reads = fragment_df.bc_num_reads.astype(np.int32)
        fragment_df.bc_est_len   = fragment_df.bc_est_len.astype(np.int32)
        fragment_df.bc_mean_reads_per_fragment = fragment_df.bc_mean_reads_per_fragment.astype(np.float32)

        return (fragment_df, bc_df)
