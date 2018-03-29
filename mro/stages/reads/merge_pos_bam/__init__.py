#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make a position-sorted bam file from a bunch of sorted buckets.
#

import json
import heapq
import os.path
import tenkit.bam as tk_bam
from sets import Set
import tenkit.reference as tk_reference
from tenkit.constants import BARCODE_LOCATION
from tenkit.constants import FRAG_LEN_HIST_BIN_SIZE
from collections import defaultdict
import tenkit.pandas as p
import numpy as np
import scipy.stats
import os
import math
from tenkit import safe_json
import tenkit.stats as tk_stats
import tenkit.seq as tk_seq
import tenkit.hdf5
import tenkit.safe_json
import subprocess
__MRO__ = """
stage MERGE_POS_BAM(
    in  json     position_chunks,
    in  string   reference_path,
    in  string   barcode_whitelist,
    in  json     barcode_count,
    in  bed      targets_file,
    in  int      bc_max_error_allowed,
    in  float    bc_pseudo_count,
    in  bool     bc_use_mapping,
    in  int      bc_mapq,
    in  bool     bc_allow_indel,
    in  bool     frag_no_merging,
    in  int      frag_mapq,
    in  float    frag_pval,
    in  int      frag_freq,
    in  int      frag_min_num_reads_off_target,
    out bam      pos_sorted_bam,
    out bam.bai  pos_sorted_bam_index,
    out json     single_partition,
    out json     fragment_size,
    out h5       fragments,
    out h5       barcodes,
    out json     barcode_histogram,
    out json     summary,
    src py       "stages/reads/merge_pos_bam",
) split using (
    in  string[] position_chunk,
    in  string tid,
    in  bool     unmapped_chunk,
)
"""

NUM_BCS_LOADING_ESTIMATE = 1000
MAPQ_THRESHOLD = 30

def split(args):
    chunk_defs = []
    with open(args.position_chunks) as position_chunks_file:
        position_chunks = json.load(position_chunks_file)

    tid2bams = defaultdict(list)
    for (key,val) in position_chunks.iteritems():
        for b in val:
            tid = os.path.basename(b).split('-')[0]
            tid2bams[tid].append(b)

    for tid, position_chunk in sorted(tid2bams.iteritems()):
        #if tid not in ["000020"]: continue
        chunk_defs.append({'position_chunk':position_chunk, 'tid':tid, "__mem_gb": 8, '__threads':2})
    return {'chunks':chunk_defs, 'join': {'__mem_gb': 16}}

def main(args, outs):
    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)
    primary_contigs_file = outs.summary+"_primary_contigs.txt"
    with open(outs.summary+"_primary_contigs.txt", "w") as fout:
        for c in primary_contigs:
            fout.write(c+"\n")

    samtools_merge_args = ['samtools','merge', '-s', '0', '-c', '-p', '-@', '4', outs.pos_sorted_bam+"tmp0.bam"]
    samtools_merge_args.extend(args.position_chunk)
    subprocess.check_call(samtools_merge_args)
    #subprocess.check_call(['samtools','index',outs.pos_sorted_bam+"tmp0.bam"])
    #tk_bam.merge(outs.pos_sorted_bam+"tmp0.bam", args.position_chunk, threads=4)

    ## update barcode and doing fragmentation
    tool = "report_single_partition"
    cmd = tool
    cmd += (" "+os.path.join(BARCODE_LOCATION, "4M-with-alts-february-2016"+".txt") )
    cmd += (" "+ args.barcode_count )
    cmd += (" " + outs.pos_sorted_bam+"tmp0.bam")
    cmd += (" " + primary_contigs_file)
    cmd += (" " + outs.summary+"_left7")
    cmd += (" " + outs.summary+"_right9")
    cmd += (" " + outs.summary+"_lectable")
    cmd += (" " + outs.summary+"_rectable")
    cmd += (" " + outs.pos_sorted_bam)
    cmd += (" --out "+outs.fragments+"_"+args.tid+".csv")

    if args.targets_file:
        targets_new = outs.summary+"_targets.bed"
        simplify_centromere_to_bed(args.targets_file, targets_new, idx = [0, 1, 2])
        cmd += " --target "+ str(targets_new)
    centromere_file = tk_reference.get_centromere_regions(args.reference_path)
    if os.path.exists(centromere_file):
        centromere_file_new = outs.summary+"_centromere.bed"
        simplify_centromere_to_bed(centromere_file, centromere_file_new)
        cmd += " --centromere "+str(centromere_file_new)

    cmd += " --maxerr "+str(args.bc_max_error_allowed)
    cmd += " --pseudo "+str(args.bc_pseudo_count)
    if args.bc_use_mapping:  cmd += " --mappingbc"
    cmd += (" --mapqbc " + str(args.bc_mapq))

    if args.frag_no_merging: cmd += " --nofragmerg"
    cmd += ("  --mapq " + str(args.frag_mapq))
    cmd += ("  --pval " + str(args.frag_pval))
    cmd += ("  --freq " + str(args.frag_freq))

    if args.bc_allow_indel:
        cmd += "  --allowindel"


    print cmd
    subprocess.check_call(cmd, shell=True)

def load_centromere_file(fn):
    centromeres = {}
    if not os.path.exists(fn):
        return centromeres

    with open(fn) as regions_file:
        for line in regions_file:
            if line.startswith('#'):
                continue
            elif line.startswith("CEN"):
                tokens = line.split("\t")
                if len(tokens) < 4:
                    continue
                chrom = tokens[1]
                start = int(tokens[2])
                end = int(tokens[3])
                centromeres[chrom] = (start, end)
    return centromeres


def join(args, outs, chunk_defs, chunk_outs):
    args_dict ={}
    args_dict["bc_allow_indel"]=args.bc_allow_indel
    args_dict["bc_max_error_allowed"]=args.bc_max_error_allowed
    args_dict["bc_pseudo_count"]=args.bc_pseudo_count
    args_dict["bc_use_mapping"]=args.bc_use_mapping
    args_dict["bc_mapq"]=args.bc_mapq
    args_dict["frag_no_merging"]=args.frag_no_merging
    args_dict["frag_mapq"]=args.frag_mapq
    args_dict["frag_pval"]=args.frag_pval
    args_dict["frag_freq"]=args.frag_freq
    fsummary = open(outs.summary, "w")
    fsummary.write(safe_json.safe_jsonify(args_dict))
    fsummary.close()

    tk_bam.concatenate(out_file_name=outs.pos_sorted_bam, all_in_file_names=[chunk.pos_sorted_bam for chunk in chunk_outs])
    tk_bam.index(outs.pos_sorted_bam)
    outs.pos_sorted_bam_index = outs.pos_sorted_bam + '.bai'

    bam_in = tk_bam.create_bam_infile(outs.pos_sorted_bam)
    chroms = bam_in.references
    barcode_whitelist = list(tk_seq.load_barcode_whitelist(args.barcode_whitelist))
    barcode_whitelist.sort()

    # Combine fragment csv files into a single h5 file
    in_csv_files = [co.fragments+"_"+cd.tid+".csv" for (cd, co)
        in zip(chunk_defs, chunk_outs) if os.path.exists(co.fragments+"_"+cd.tid+".csv")]


    nfrags = 0
    if len(in_csv_files) > 0:
        bc_num_frags = defaultdict(int)
        bc_num_reads = defaultdict(int)
        bc_num_single_reads = defaultdict(int)
        bc_num_lens = defaultdict(int)

        temp_csv_barcodes = outs.barcodes+"_temp.csv"
        nfrags = 0

        for f in in_csv_files:
            # TODO - sequentially append to fragments.h5 file to keep memory under control
            # - handle multiple GEM groups properly.
            # ensure the chroms column has string /categorical type in hdf5
            # - same fixes for barcodes.h5 file
            # handle 0-length outputs -- does that result in None file outs?
            frag_in = p.read_csv(f, names=["tid", "start_pos", "end_pos", "bc_id", "num_reads"])
            frag_in["obs_len"] = frag_in.end_pos - frag_in.start_pos
            frag_in[frag_in.num_reads <= 1].obs_len = 1000

            frag_in["est_len"] = np.maximum(1, frag_in["obs_len"] * (frag_in.num_reads + 1) / np.maximum(1, frag_in.num_reads - 1)).astype("int")
            frag_in[frag_in.num_reads <= 1].est_len = 1000
            
            barcode_seqs = []
            molecule_ids = []
    
            for (i, row) in frag_in.iterrows():

                bc_num_frags[row.bc_id] += 1
                bc_num_reads[row.bc_id] += row.num_reads
                bc_num_lens[row.bc_id] += row.est_len
                    
                bc_wl_id = int(row.bc_id) % len(barcode_whitelist)
                gg = int(row.bc_id) / len(barcode_whitelist) + 1
                barcode_seq = "%s-%d" % (barcode_whitelist[bc_wl_id], gg)
                barcode_seqs.append(barcode_seq)
                molecule_ids.append(nfrags)

                nfrags += 1

            frag_in["bc"] = p.Categorical(barcode_seqs)
            frag_in["chrom"] = p.Categorical.from_codes(frag_in.tid, chroms)
            frag_in["molecule_id"] = molecule_ids
            del frag_in["tid"]
            del frag_in["bc_id"]

            if len(frag_in) > 0:
                tenkit.hdf5.append_data_frame(outs.fragments, frag_in)


        with open(temp_csv_barcodes, "w") as csv_out:
            csv_out.write("bc,bc_est_len,bc_linked_read_fraction,bc_linked_fragment_fraction,bc_mean_reads_per_fragment,bc_num_fragments,bc_num_reads\n")
            for bc_id in range(len(barcode_whitelist)):
                bc = barcode_whitelist[bc_id]+"-1"
                if bc_id in bc_num_frags:
                    bc_est_len = bc_num_lens[bc_id]
                    bc_linked_read_fraction = 1.0 - bc_num_single_reads[bc_id]*1.0/bc_num_reads[bc_id]
                    bc_linked_fragment_fraction = 1.0 - bc_num_single_reads[bc_id]*1.0/bc_num_frags[bc_id]
                    bc_mean_reads_per_fragment = bc_num_reads[bc_id]*1.0/bc_num_frags[bc_id]
                    csv_out.write("%s,%d,%f,%f,%f,%d,%d\n" % (bc, bc_est_len, bc_linked_read_fraction, bc_linked_fragment_fraction, bc_mean_reads_per_fragment, bc_num_frags[bc_id], bc_num_reads[bc_id]))


        if nfrags == 0:
            outs.fragments = None
            outs.barcodes = None

        else:
            tenkit.hdf5.create_tabix_index(outs.fragments, 'chrom', 'start_pos', 'end_pos')

            df_barcodes = p.read_csv(temp_csv_barcodes)
            tenkit.hdf5.append_data_frame(outs.barcodes, df_barcodes)

    else:
        outs.fragments = None
        outs.barcodes= None

    summary =  {}
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





def update_mapqs(bamfilename, outfile, reference_path, centromeres):
    bam = tk_bam.create_bam_infile(bamfilename)
    bam_out, _ = tk_bam.create_bam_outfile(outfile, None, None, template=bam)
    variant_heap = []
    variant_map = {}
    read_heap = []
    primary_contigs = tk_reference.load_primary_contigs(reference_path)
    for read in bam:
        tags = [(key, value) for (key, value) in dict(read.tags).iteritems()]
        tags.append(('OM',int(read.mapq)))
        read.tags = tags
        chrom = bam.references[read.tid]
        if chrom not in primary_contigs or (chrom in centromeres and read.pos >= centromeres[chrom][0] and read.pos <= centromeres[chrom][1]):
            read.tags = [(key, value) for (key, value) in read.tags if key != 'AC' and key != 'XC']
            bam_out.write(read)
            continue
        add_variant_counts(read, variant_heap, variant_map)
        heapq.heappush(read_heap, (read.pos, read))
        update_updatable(read_heap, read.pos, variant_heap, variant_map, bam_out)
    update_updatable(read_heap, 500000000, variant_heap, variant_map, bam_out, empty_me = True)
    bam_out.close()


def add_variant_counts(read, variant_heap, variant_map):
    for variant in get_variants_best(read):
        variant_id = get_variant_id(variant, read)
        if variant_id in variant_map:
            variant_map[variant_id] += 1
        else:
            variant_map[variant_id] = 1
            heapq.heappush(variant_heap, variant_id)

def update_updatable(read_heap, pos, variant_heap, variant_map, bam_out, empty_me = False):
    while len(read_heap) > 0 and ((pos - read_heap[0][0] > 300 or empty_me) or len(read_heap) > 3000):
        _, read = heapq.heappop(read_heap)
        v_best = Set()
        v_best_full = {}
        for variant_in_best in get_variants_best(read):
            bases = bases_used(variant_in_best)
            v_best.add(bases)
            v_best_full[bases] = variant_in_best
        v_second_best = Set()
        for variant_in_second_best in get_variants_second_best(read):
            v_second_best.add(bases_used(variant_in_second_best))
        diff1 = v_best - v_second_best
        best_candidates = [v_best_full[b] for b in diff1]
        for variant_in_best in best_candidates:
            count = variant_map[get_variant_id(variant_in_best, read)]
            if count > 1:
                penalty = get_penalty(variant_in_best)
                modified_penalty = float(penalty)/float(count)
                AS = float(dict(read.tags).get('AS'))
                calculated_xs = AS - (read.mapq/10.0) # use this instead of the actual XS as this has already taken into account multiple alternate alignments
                new_mapq = (AS - penalty + modified_penalty - calculated_xs) * 10.0
                new_mapq = int(min(new_mapq, 60.0))
                read.mapq = new_mapq
        read.tags = [(key, value) for (key, value) in read.tags if key != 'AC' and key != 'XC']
        bam_out.write(read)
    while len(variant_heap) > 0 and (pos - variant_heap[0][0] > 400 or empty_me):
        variant_map.pop(variant_heap[0], None)
        heapq.heappop(variant_heap)

def get_penalty(variant):
    varlength = variant[2]
    if varlength == 1:
        return -2
    else:
        return -3

def get_variants_best(read):
    return get_variants(read, "AC")

def get_variants_second_best(read):
    return get_variants(read, "XC")

def get_variants(read, tag):
    var_string = dict(read.tags).get(tag)
    if var_string is None:
        return []
    else:
        variants = var_string.split(";")
        ret = []
        for v in variants:
            if len(v) > 1:
                info = v.split(',')
                if len(info) < 3:
                    print v
                    print tag
                ret.append((int(info[0]), int(info[1]), int(info[2]))) # format is pos, start base in read, length
    return ret

def bases_used(variant):
    return (variant[1],variant[2]) # returning start base in read, length

def get_variant_id(variant,read):
    return (variant[0], variant[2], read.seq[variant[1]: variant[1] + max(0, variant[2])]) # returning pos in reference, length, bases (empty string for deletion)


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

def simplify_centromere_to_bed(old, new, idx = [1, 2, 3]):
    with open(new, "w") as fnew, open(old) as fold:
        #fnew.write("chrom\tstart\tend\n")
        for line in fold:
            if line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track') or line.startswith('#'):
                continue
            items = line.split("\t")
            fnew.write("%s\t%s\t%s\n" % (items[idx[0]], items[idx[1]], items[idx[2]]))

