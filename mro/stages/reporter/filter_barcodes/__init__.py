#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# filter barcodes that may have coalesced into one gem (on the basis of genomic overlap)
import tenkit.pandas as p
import numpy as np

import tenkit.hdf5
import tenkit.safe_json
import coalescence

import martian

__MRO__ = """
stage FILTER_BARCODES(
    in  bam    input,
    in  h5     fragments,
    in  h5     barcodes,
    out tsv    barcode_blacklist,
    out json   filter_barcodes_results,
    src py     "stages/reporter/filter_barcodes",
) split using (
)
"""

def split(args):
    return {'chunks': [{}], 'join': {'__mem_gb': 16}}


def main(args, outs):
    pass

def join(args, outs, chunk_defs, chunk_outs):

    summary = {}
    # Compute high-level BC summary metrics
    # Load BC data
    if args.barcodes:
        bc_df = tenkit.hdf5.read_data_frame(args.barcodes)
        fragment_df = tenkit.hdf5.read_data_frame(args.fragments, query_cols=['bc', 'chrom', 'start_pos'])

        bc_df.sort('bc_num_reads', inplace=True)
        bc_df['cum_reads'] = np.cumsum(bc_df.bc_num_reads)

        # Measure coalescence rate on all BCs that could conceivably be used
        # to call SVs - i.e. ignore BCs that contribute the cumulative bottom 1% of reads
        n99_read_thresh = sum(bc_df.bc_num_reads) * 0.01
        n99_bcs = bc_df[bc_df.cum_reads > n99_read_thresh]
        martian.log_info("number of bcs to screen for coalescence: %d" % len(n99_bcs))
        martian.log_info("subsetting fragments to use")

        if len(n99_bcs) > 1:
            selected_frags = fragment_df[fragment_df.bc.isin(n99_bcs.bc)]
            del fragment_df
            martian.log_info("Doing coalescence calculation")
            coa_calc = coalescence.BcSimilarity(selected_frags, set(n99_bcs.bc), args.input)
            coa_bc_tbl = coa_calc.coalescence_analysis()

            # Also add barcodes that are extreme outliers in the number of fragments observed
            med_frags_per_bc = n99_bcs.bc_num_fragments.median()
            high_quantile = n99_bcs.bc_num_fragments.quantile(0.98)
            bc_num_fragments_threshold = max(med_frags_per_bc*5.0, high_quantile)

            med_reads_per_bc = n99_bcs.bc_num_reads.median()
            high_quantile = n99_bcs.bc_num_reads.quantile(0.98)
            bc_num_reads_threshold = max(med_reads_per_bc*5.0, high_quantile)

            overloaded_bcs = n99_bcs[(n99_bcs.bc_num_fragments > bc_num_fragments_threshold) | (n99_bcs.bc_num_reads > bc_num_reads_threshold)]
            summary['fract_bcs_overloaded'] = float(len(overloaded_bcs)) / len(n99_bcs)

            # Remove bcs that are already in the blacklist
            nr_overloaded_bcs = overloaded_bcs[~overloaded_bcs.bc.isin(coa_bc_tbl.bc)]

            # Add overloaded bcs to blacklist
            overloaded_bc_tbl = p.DataFrame({'bc': nr_overloaded_bcs.bc, 'cluster_id': -1, 'cluster_size': -1})

            # Write barcode blacklist
            bad_bc_tbl = p.concat([coa_bc_tbl, overloaded_bc_tbl])
            bad_bc_tbl.to_csv(outs.barcode_blacklist, sep="\t", index=False)

            # Compute coalescence stats
            summary['fract_bcs_in_clusters_all'] = float(len(coa_bc_tbl)) / len(n99_bcs)
            summary['fract_bcs_in_clusters_eq_2'] = float((coa_bc_tbl.cluster_size == 2).sum()) / len(n99_bcs)
            summary['fract_bcs_in_clusters_gt_2'] = float((coa_bc_tbl.cluster_size > 2).sum()) / len(n99_bcs)
            summary['num_clusters_gt_8'] = (coa_bc_tbl.cluster_size > 8).sum()

            # Compute stats ignoring clusters of Hamming distance 2
            hd2_clusters = []
            for cluster in coa_bc_tbl.groupby('cluster_id'):
                if all_within_hamming_distance(cluster[1].bc.values, 2):
                    hd2_clusters.append(cluster[0])

            coa_tbl_no_hd2 = coa_bc_tbl[~coa_bc_tbl.cluster_id.isin(hd2_clusters)]
            summary['fract_bcs_in_clusters_all_no_hd2'] = float(len(coa_tbl_no_hd2)) / len(n99_bcs)
            summary['fract_bcs_in_clusters_eq_2_no_hd2'] = float((coa_tbl_no_hd2.cluster_size == 2).sum()) / len(n99_bcs)
            summary['fract_bcs_in_clusters_gt_2_no_hd2'] = float((coa_tbl_no_hd2.cluster_size > 2).sum()) / len(n99_bcs)

        else:
            empty_df = p.DataFrame({'bc':[], 'cluster_id':[], 'cluster_size':[]})
            empty_df.to_csv(outs.barcode_blacklist, sep="\t", index=False)

            # null coalescence stats
            summary['fract_bcs_overloaded'] = None
            summary['fract_bcs_in_clusters_all'] = None
            summary['fract_bcs_in_clusters_eq_2'] = None
            summary['fract_bcs_in_clusters_gt_2'] = None
            summary['num_clusters_gt_8'] = None
            summary['fract_bcs_in_clusters_all_no_hd2'] = None
            summary['fract_bcs_in_clusters_eq_2_no_hd2'] = None
            summary['fract_bcs_in_clusters_gt_2_no_hd2'] = None

    else:
        outs.barcode_blacklist = None
        summary['fract_bcs_overloaded'] = None
        summary['fract_bcs_in_clusters_all'] = None
        summary['fract_bcs_in_clusters_eq_2'] = None
        summary['fract_bcs_in_clusters_gt_2'] = None
        summary['num_clusters_gt_8'] = None
        summary['fract_bcs_in_clusters_all_no_hd2'] = None
        summary['fract_bcs_in_clusters_eq_2_no_hd2'] = None
        summary['fract_bcs_in_clusters_gt_2_no_hd2'] = None


    # Write summary to json
    with open(outs.filter_barcodes_results, 'w') as results_file:
        tenkit.safe_json.dump_numpy(summary, results_file, pretty=True)

def all_within_hamming_distance(bcs, dist):
    for i in range(0, len(bcs)):
        for j in range(i+1, len(bcs)):
            if not within_hamming_distance(bcs[i], bcs[j], dist):
                return False
    return True

def within_hamming_distance(bc1, bc2, dist):
    assert(len(bc1) == len(bc2))
    diffs = 0
    for c1, c2 in zip(bc1, bc2):
        if c1 != c2:
            diffs += 1
        if diffs > dist:
            return False
    return True
