#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Coallate basic information about reads from a BAM file
#
import math
import numpy as np
import tenkit.pandas as p
import cPickle as pickle
import pysam
import shutil

import tenkit.safe_json
import tenkit.fasta as tk_fasta
import tenkit.bam as tk_bam
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import tenkit.bio_io as tk_io
import tenkit.hdf5
import tenkit.summary_manager as tk_summary
from tenkit.constants import (HIGH_CONF_MAPQ, BC_QUAL_CUTOFFS, INSERT_MAPQ_CUTOFFS, TARGET_MAPQ_CUTOFFS,
                              MAX_INSERT_SIZE, MAX_TARGET_DIST, MODERATE_CONF_MAPQ, READ_MATE_FAR_DIST,
                              READ_MATE_CHIM_TOO_CLOSE_DIST)

import martian


__MRO__ = """
stage REPORT_BASIC(
    in  bam    input,
    in  bam    input_pos,
    in  bed    targets_file,
    in  string barcode_whitelist,
    in  json   lot_info,
    out json   summary,
    out h5     barcode_counts,
    out h5     barcode_counts_mapped,
    out json   insert_sizes,
    out json   target_dists,
    out json   mapq_counts,
    out txt    misc_sm,
    out txt    qual_sms,
    out json   lot_info,
    src py     "stages/reporter/report_basic",
) split using (
    in  string chunk_start,
    in  string chunk_end,
    in  int    n_chunks,
)
"""

def split(args):
    bam = pysam.Samfile(args.input, check_sq=False)

    min_chunks = 1
    if args.barcode_whitelist is not None:
        barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)
        if len(barcode_whitelist) > 1e6:
            min_chunks = 4

    # Split to ensure read pairs always go together
    chunks = tk_bam.chunk_bam_records(bam, lambda x: x.qname, min_chunks=min_chunks)
    for dict in chunks:
        dict['n_chunks'] = len(chunks)
        dict['__mem_gb'] = 3
    return {'chunks': chunks, 'join': {'__mem_gb': 8}}

def main(args, outs):
    main_report_basic(args, outs)

def combine_summary_managers(outs):
    sm = None
    for out in outs:
        tmp_sm = tk_summary.SummaryManager.load( out )
        if sm is None:
            sm = tmp_sm
        else:
            sm = sm.combine( tmp_sm )
    return sm

def combine_nested_summary_managers(outs):
    dict = {}
    for out in outs:
        tmp_dict = None
        with open(out, 'rb') as out_handle:
            tmp_dict = pickle.load( out_handle )
        for key in tmp_dict:
            tmp_sm = tmp_dict[key]
            if not(key in dict):
                dict[key] = tmp_sm
            else:
                dict[key] = dict[key].combine( tmp_sm )
    return dict


def join(args, outs, chunk_defs, chunk_outs):
    ''' join the various outputs created by report_basic '''
    chunk_outs = list(chunk_outs)

    martian.log_info("combining misc summary managers")
    misc_sm_outs = [x.misc_sm for x in chunk_outs]
    misc_sm = combine_summary_managers( misc_sm_outs )

    martian.log_info("combining nested summary managers")
    qual_sms_outs = [x.qual_sms for x in chunk_outs]
    qual_sms = combine_nested_summary_managers( qual_sms_outs )

    martian.log_info("computing summary metrics")
    compute_summary_metrics( misc_sm, qual_sms)

    if args.barcode_whitelist:
        barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)
    else:
        barcode_whitelist = None

    # barcode hdf5
    if outs.barcode_counts or outs.barcode_counts_mapped:
        bc_table, bc_table_mapped = summarize_barcode_data(misc_sm, qual_sms, barcode_whitelist)
    if outs.barcode_counts:
        tenkit.hdf5.write_data_frame(outs.barcode_counts, bc_table)
    if outs.barcode_counts_mapped:
        tenkit.hdf5.write_data_frame(outs.barcode_counts_mapped, bc_table_mapped)
    
    # insert sizes output
    insert_size_dists = {}
    for qual in INSERT_MAPQ_CUTOFFS:
        insert_size_dists[qual] = qual_sms['insert_size_dists'].get_summarizer(qual).dict
    insert_sizes_output_file = open(outs.insert_sizes, 'w')
    insert_sizes_output_file.write(tenkit.safe_json.safe_jsonify(insert_size_dists) + '\n')
    insert_sizes_output_file.close()

    # target distances
    nearest_targ_dists = {}
    for qual in TARGET_MAPQ_CUTOFFS:
        nearest_targ_dists[qual] = qual_sms['nearest_targ_dists'].get_summarizer(qual).dict
    target_dists_output_file = open(outs.target_dists, 'w')
    target_dists_output_file.write(tenkit.safe_json.safe_jsonify(nearest_targ_dists))
    target_dists_output_file.close()

    # overall summary metrics
    summary_output_file = open(outs.summary, 'w')
    summary_output_file.write(tenkit.safe_json.safe_jsonify(misc_sm.get_summarizer('metrics').dict, pretty=True))
    summary_output_file.close()

    # mapq counts
    mapq_output_file = open(outs.mapq_counts, 'w')
    mapq_output_file.write(tenkit.safe_json.safe_jsonify(misc_sm.get_summarizer('mapq_counts').dict))
    mapq_output_file.close()

    # lot info - just copy the input file and return the copy as an output.
    # this is because the file comes from TRIM_READS originally, and we want to eliminate
    # dependencies on that stage's outputs (which are large) as early as possible.
    # TODO this should no longer be necessary when we have parameter-specific VDR
    if args.lot_info is not None:
        shutil.copyfile(args.lot_info, outs.lot_info)
    else:
        outs.lot_info = None

def main_report_basic(args, outs):
    bam_in = pysam.Samfile(args.input, check_sq=False)
    targets_filename = args.targets_file
    references = bam_in.references

    if args.input_pos is not None:
        bam_in_pos = tk_bam.create_bam_infile(args.input_pos)
        n_mapped = bam_in_pos.mapped
        n_chunk = math.ceil( n_mapped / args.n_chunks )
        bam_in_pos.close()
    else:
        n_mapped = 0
        n_chunk = 0

    if targets_filename is None or targets_filename == '':
        target_regions = None
    else:
        targets_file = open(targets_filename, 'r')
        target_regions = tk_io.get_target_regions(targets_file)

    if args.barcode_whitelist:
        barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)
    else:
        barcode_whitelist = None

    bam_slice = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))

    # do basic counting
    misc_sm, qual_sms = \
            compute_basic_stats(bam_slice,
                    target_regions,
                    n_chunk,
                    references,
                    barcode_whitelist)


    misc_sm.save( outs.misc_sm )
    with open(outs.qual_sms, 'wb') as out_handle:
        pickle.dump( qual_sms, out_handle )


def compute_basic_stats(bam_in, target_regions, n_mapped, references,
        barcode_whitelist=None):

    # Things we will compute from the read info in 'metrics':
    # num_reads, num_unmapped, num_single_mapped, num_conf_mapped,
    # num_pos_chimeras, num_far_chimeras, num_same_dir_chimeras,
    # num_outward_dir_chimeras, bad_bc_count, total_bases, mapped_bases, ...

    # For computing duplicate rate:
    # duplicate_read_count, non_duplicate_read_count
    metrics = tk_summary.DictionaryDistribution()

    insert_sizes_hist = tk_summary.QuantileEstimator( 0.01, n_mapped )
    r1_len_hist = tk_summary.QuantileEstimator( 0.01, n_mapped )
    r2_len_hist = tk_summary.QuantileEstimator( 0.01, n_mapped )

    mapq_counts = tk_summary.DictionaryDistribution()

    # Split by min QV - just used for JSON output
    # since theses are stores as nested dictionaries, each has their own summarizer
    qual_sms = {}

    # Counting all bc reads
    raw_bc_counts = tk_summary.DictionaryDistribution()
    processed_bc_counts = tk_summary.DictionaryDistribution()
    raw_bc_counts_mapped = tk_summary.DictionaryDistribution()
    
    # Mean sample index and BC quality
    mean_si_qual = tk_summary.BasicSummary()
    mean_bc_qual = tk_summary.BasicSummary()

    misc_sm = tk_summary.SummaryManager()
    misc_sm.add_summarizer( metrics, 'metrics' )
    misc_sm.add_summarizer( insert_sizes_hist, 'insert_sizes_hist' )
    misc_sm.add_summarizer( r1_len_hist, 'r1_len_hist')
    misc_sm.add_summarizer( r2_len_hist, 'r2_len_hist')
    misc_sm.add_summarizer( mapq_counts, 'mapq_counts' )
    misc_sm.add_summarizer( raw_bc_counts, 'raw_bc_counts' )
    misc_sm.add_summarizer( raw_bc_counts_mapped, 'raw_bc_counts_mapped' )
    misc_sm.add_summarizer( processed_bc_counts, 'processed_bc_counts' )
    misc_sm.add_summarizer( mean_si_qual, 'mean_si_qual' )
    misc_sm.add_summarizer( mean_bc_qual, 'mean_bc_qual' )

    insert_size_dists = qual_sms['insert_size_dists'] = tk_summary.SummaryManager()
    for mapq_cutoff in INSERT_MAPQ_CUTOFFS:
        insert_size_dists.add_summarizer( tk_summary.DictionaryDistribution(), mapq_cutoff )

    nearest_targ_dists = qual_sms['nearest_targ_dists'] = tk_summary.SummaryManager()
    for mapq_cutoff in TARGET_MAPQ_CUTOFFS:
        nearest_targ_dists.add_summarizer( tk_summary.DictionaryDistribution(), mapq_cutoff )
    nearest_targ_dists.set_attrib('empty', True)

    for info in read_pair_info_iter(bam_in, target_regions, references):

        metrics.add('num_reads', 2) # NEED ACCOUNT FOR SINGLE END!!!
        metrics.add('mapped_bases', info['mapped_bases'])
        metrics.add('soft_clipped_bases', info['soft_clipped_bases'])

        r1_seq_len = info['r1_seq_len']
        r2_seq_len = info['r2_seq_len']

        r1_mapped = info['r1_mapped']
        r2_mapped = info['r2_mapped']

        r1_mapq = info['r1_mapq']
        r2_mapq = info['r2_mapq']

        if r1_mapq >= HIGH_CONF_MAPQ:
            metrics.add('num_conf_mapped')

        if r2_mapq >= HIGH_CONF_MAPQ:
            metrics.add('num_conf_mapped')

        mapq_counts.add(r1_mapq)
        mapq_counts.add(r2_mapq)
        if info['sample_index_mean_qual'] is not None:
            mean_si_qual.add(info['sample_index_mean_qual'])

        # Internal runs might not have 10X bc
        if info['10X_bc_mean_qual'] is not None:
            mean_bc_qual.add(info['10X_bc_mean_qual'])

        if r1_mapq >= MODERATE_CONF_MAPQ and r2_mapq >= MODERATE_CONF_MAPQ and (not(info['r1_chrom'] == info['r2_chrom']) or abs(info['r1_pos'] - info['r2_pos']) > READ_MATE_CHIM_TOO_CLOSE_DIST):
            metrics.add('num_pos_chimeras')
            if not(info['r1_chrom'] == info['r2_chrom']):
                metrics.add('num_far_chimeras')
                metrics.add('num_inter_chimeras')
            else:
                if abs(info['r1_pos'] - info['r2_pos']) > READ_MATE_FAR_DIST:
                    metrics.add('num_far_chimeras')
                elif info['r1_direct'] == info['r2_direct']:
                    metrics.add('num_same_dir_chimeras')
                elif info['r1_pos'] < info['r2_pos']:
                    if not(info['r1_direct']) and info['r2_direct']:
                        metrics.add('num_outward_dir_chimeras')
                elif info['r2_pos'] < info['r1_pos']:
                    if not(info['r2_direct']) and info['r1_direct']:
                        metrics.add('num_outward_dir_chimeras')

        if not(r1_seq_len is None):
            metrics.add('total_bases', r1_seq_len)

        if not(r2_seq_len is None):
            metrics.add('total_bases', r2_seq_len)

        if not(r1_mapped and r2_mapped):
            metrics.add('num_unmapped', 2)

            if r1_mapped or r2_mapped:
                metrics.add('num_single_mapped', 2)

        else:
            if info['r1_chrom'] == info['r2_chrom']:

                # Summarize read length
                r1_len_hist.add(info['r1_seq_len'])
                r2_len_hist.add(info['r2_seq_len'])

                # Ignore insert sizes that are huge (and likely chimeras) to
                # avoid skewing the insert size distribution
                insert_size = info['insert_length']
                if insert_size is not None:
                    if insert_size <= MAX_INSERT_SIZE:
                        insert_sizes_hist.add(insert_size)

                    # Clip insert size to reasonable number to prevent blow-up of dictionary entries
                    insert_size = min(insert_size, MAX_INSERT_SIZE)

                    min_mapq = min(r1_mapq, r2_mapq)
                    for mapq_cutoff in INSERT_MAPQ_CUTOFFS:

                        if min_mapq >= mapq_cutoff:
                            dists = insert_size_dists.get_summarizer(mapq_cutoff)
                            dists.add(insert_size)

                target_dist = info['read_pair_targ_dist']
                if not(target_dist is None):
                    target_dist = min(target_dist, MAX_TARGET_DIST)
                    for mapq_cutoff in TARGET_MAPQ_CUTOFFS:
                        if min_mapq >= mapq_cutoff:
                            dists = nearest_targ_dists.get_summarizer(mapq_cutoff)
                            dists.add(target_dist)
                            if nearest_targ_dists.get_attrib('empty'):
                                nearest_targ_dists.set_attrib('empty', False)

                # Count duplicate and non-duplicate read pairs
                if info['is_duplicate']:
                    metrics.add('duplicate_read_count')
                else:
                    metrics.add('non_duplicate_read_count')

        # Base quality scores
        r1_q20_bases = info['r1_q20_bases']
        r1_q30_bases = info['r1_q30_bases']

        if not(r1_q20_bases is None or r1_q30_bases is None or r1_seq_len is None or r1_seq_len == 0):
            metrics.add('r1_q20_bases', r1_q20_bases)
            metrics.add('r1_q30_bases', r1_q30_bases)
            metrics.add('r1_tot_bases', r1_seq_len + info['10X_raw_bc_len'] + len(info['r1_trimmed_qual']))
            # include lengths of barcode + trimmed bases in r1_tot_bases to match their inclusion in QV metrics

        r2_q20_bases = info['r2_q20_bases']
        r2_q30_bases = info['r2_q30_bases']

        if not(r2_q20_bases is None or r2_q30_bases is None or r2_seq_len is None or r2_seq_len == 0):
            metrics.add('r2_q20_bases', r2_q20_bases)
            metrics.add('r2_q30_bases', r2_q30_bases)
            metrics.add('r2_tot_bases', r2_seq_len)

        si_q20_bases = info['si_q20_bases']
        si_q30_bases = info['si_q30_bases']
        si_tot_bases = info['sample_index_len']

        if not(si_q20_bases is None or si_q30_bases is None or si_tot_bases is None or si_tot_bases == 0):
            metrics.add('si_q20_bases', si_q20_bases)
            metrics.add('si_q30_bases', si_q30_bases)
            metrics.add('si_tot_bases', si_tot_bases)

        bc_q20_bases = info['bc_q20_bases']
        bc_q30_bases = info['bc_q30_bases']
        bc_tot_bases = info['10X_raw_bc_len']

        if not(bc_q20_bases is None or bc_q30_bases is None or bc_tot_bases is None or bc_tot_bases == 0):
            metrics.add('bc_q20_bases', bc_q20_bases)
            metrics.add('bc_q30_bases', bc_q30_bases)
            metrics.add('bc_tot_bases', bc_tot_bases)

        # If the raw_bc is None, then barcodes weren't attached
        # Don't count barcode info when using degenerate barcodes
        if barcode_whitelist is not None:
            raw_bc = info['10X_raw_bc']
            if raw_bc is not None and raw_bc in barcode_whitelist:
                raw_bc_counts.add(raw_bc)
                if info['r1_mapped'] or info['r2_mapped']:
                    raw_bc_counts_mapped.add(raw_bc)  # include for barcode_counts_mapped.h5 output
            proc_bc = info['10X_called_bc']
            if proc_bc:
                processed_bc_counts.add(proc_bc)
            elif raw_bc:
                metrics.add('bad_bc_count')

        r1_contains_N = info['r1_contains_N']
        r2_contains_N = info['r2_contains_N']

        if r1_contains_N:
            metrics.add('r1_contains_N', 1)
        if r2_contains_N:
            metrics.add('r2_contains_N', 1)

    return (misc_sm, qual_sms)

def summarize_barcode_data(misc_sm, qual_sms, barcode_whitelist):
    processed_bc_counts = misc_sm.get_summarizer('processed_bc_counts').dict
    metrics = misc_sm.get_summarizer('metrics')
    raw_bc_counts = misc_sm.get_summarizer('raw_bc_counts').dict
    raw_bc_counts_mapped = misc_sm.get_summarizer('raw_bc_counts_mapped').dict

    # Compute high-level stats on the barcodes. Only emit bc metrics if we have a whitelist and have barcodes
    # attached
    if barcode_whitelist and len(raw_bc_counts) > 0:

        # What fraction of the whitelist did we see?
        # Don't trust that the whitelist used during attach_bcs is the same as the one used here
        # Must intersect with the raw_bc sequences.  The processed_bc keys have the GEM group prepended.
        observed_good = set(raw_bc_counts.keys()).intersection(barcode_whitelist)

        metrics['fraction_bcs_observed'] = \
                min(1.0, tk_stats.robust_divide(float(len(observed_good)), len(barcode_whitelist)))

        # What fraction of clusters had a correct barcode
        total_good_bc_observations = float(sum(processed_bc_counts.values()))
        metrics['correct_bc_rate'] = \
                tk_stats.robust_divide(total_good_bc_observations,
                        (total_good_bc_observations + metrics['bad_bc_count']))

        # 'Effective diversity' of barcodes
        sum_sq = sum(( v**2 for v in processed_bc_counts.values()))
        effective_diversity = tk_stats.robust_divide((total_good_bc_observations**2.0), float(sum_sq))
        metrics['effective_diversity'] = effective_diversity

        def fraction_within_f(counts, f):
            med_counts = np.median(counts)
            counts_in_range = np.logical_and(counts >= med_counts/f, counts <= med_counts*f)
            return np.mean(counts_in_range)

        # fraction of barcodes with abundance within 2x of median
        count_array = np.array(processed_bc_counts.values())
        metrics['fraction_bc_within_2x_median'] = fraction_within_f(count_array, 2.0)
        metrics['fraction_bc_within_1.15x_median'] = fraction_within_f(count_array, 1.15)

        metrics['barcode_count_cv'] = np.std(count_array) / np.mean(count_array)


        raw_count_vect = list(raw_bc_counts.items())
        raw_bc_list = [ x[0] for x in raw_count_vect ]
        raw_bc_count = np.array([ x[1] for x in raw_count_vect ], dtype=np.int32)
        is_valid_bc = [ bc in barcode_whitelist for bc in raw_bc_list ]

        bc_table_cols = { 'bc_sequence': raw_bc_list, 'count': raw_bc_count, 'valid_bc': is_valid_bc}
        bc_table = p.DataFrame(bc_table_cols)

        raw_count_vect_mapped = list(raw_bc_counts_mapped.items())
        raw_bc_list_mapped = [ x[0] for x in raw_count_vect_mapped ]
        raw_bc_count_mapped = np.array([ x[1] for x in raw_count_vect_mapped ], dtype=np.int32)
        is_valid_bc_mapped = [ bc in barcode_whitelist for bc in raw_bc_list_mapped ]
        
        bc_table_cols_mapped = { 'bc_sequence': raw_bc_list_mapped, 'count': raw_bc_count_mapped, 'valid_bc': is_valid_bc_mapped}
        bc_table_mapped = p.DataFrame(bc_table_cols_mapped)
        
        
    else:
        metrics['fraction_bcs_observed'] = None
        metrics['correct_bc_rate'] = None
        metrics['effective_diversity'] = None
        metrics['fraction_bc_within_2x_median'] = None

        dummy_bc_table = { 'bc_sequence': np.array([], dtype=np.object),
                           'count': np.array([], dtype=np.int32),
                           'valid_bc': np.array([], dtype=np.bool) }
        for minqv in BC_QUAL_CUTOFFS:
            name = "mapped_minqv_%d_count" % minqv
            dummy_bc_table[name] = np.array([], dtype=np.int32)

            name = "unmapped_minqv_%d_count" % minqv
            dummy_bc_table[name] = np.array([], dtype=np.int32)

        bc_table = p.DataFrame(dummy_bc_table)
        bc_table_mapped = p.DataFrame(dummy_bc_table)
    
    return bc_table, bc_table_mapped

def compute_summary_metrics(misc_sm, qual_sms):
    ''' called in join step - extract summary metrics from combined summarizer objects '''
    metrics = misc_sm.get_summarizer('metrics')
    insert_hist = misc_sm.get_summarizer('insert_sizes_hist')
    r1_len_hist = misc_sm.get_summarizer('r1_len_hist')
    r2_len_hist = misc_sm.get_summarizer('r2_len_hist')

    nearest_targ_dists = qual_sms['nearest_targ_dists']
    print 'nearest_targ_dists',nearest_targ_dists.get_attrib('empty')

    unmapped_fract = float(metrics['num_unmapped']) / float(metrics['num_reads'])
    single_mapped_fract = float(metrics['num_single_mapped']) / float(metrics['num_reads'])

    metrics['unmapped_fract'] = unmapped_fract
    metrics['single_mapped_fract'] = single_mapped_fract

    metrics['clipped_base_fract'] = tk_stats.robust_divide(float(metrics['soft_clipped_bases']), metrics['total_bases'])
    metrics['mapped_base_fract'] = tk_stats.robust_divide(float(metrics['mapped_bases']), metrics['total_bases'])

    metrics['r1_q20_bases_fract'] = tk_stats.robust_divide(float(metrics['r1_q20_bases']), metrics['r1_tot_bases'])
    metrics['r1_q30_bases_fract'] = tk_stats.robust_divide(float(metrics['r1_q30_bases']), metrics['r1_tot_bases'])

    metrics['r2_q20_bases_fract'] = tk_stats.robust_divide(float(metrics['r2_q20_bases']), metrics['r2_tot_bases'])
    metrics['r2_q30_bases_fract'] = tk_stats.robust_divide(float(metrics['r2_q30_bases']), metrics['r2_tot_bases'])

    metrics['si_q20_bases_fract'] = tk_stats.robust_divide(float(metrics['si_q20_bases']), metrics['si_tot_bases'])
    metrics['si_q30_bases_fract'] = tk_stats.robust_divide(float(metrics['si_q30_bases']), metrics['si_tot_bases'])

    metrics['bc_q20_bases_fract'] = tk_stats.robust_divide(float(metrics['bc_q20_bases']), metrics['bc_tot_bases'])
    metrics['bc_q30_bases_fract'] = tk_stats.robust_divide(float(metrics['bc_q30_bases']), metrics['bc_tot_bases'])

    median_insert_size = insert_hist.quantile(0.5)
    iqr_insert_size = insert_hist.quantile(0.75) - insert_hist.quantile(0.25)

    metrics['median_insert_size'] = median_insert_size
    metrics['iqr_insert_size'] = iqr_insert_size

    metrics['r1_len'] = r1_len_hist.quantile(0.5)
    metrics['r2_len'] = r2_len_hist.quantile(0.5)

    conf_map_fract = float(metrics['num_conf_mapped'])/float(metrics['num_reads'])
    metrics['conf_map_fract'] = conf_map_fract

    metrics['mean_si_qual'] = misc_sm.get_summarizer('mean_si_qual').mean()
    metrics['mean_bc_qual'] = misc_sm.get_summarizer('mean_bc_qual').mean()

    metrics['reads_containing_N_fract'] = tk_stats.robust_divide(metrics['r1_contains_N'] + metrics['r2_contains_N'], metrics['num_reads'])

    # Target distance summary
    if not nearest_targ_dists.get_attrib('empty'):
        targ_dists = nearest_targ_dists.get_summarizer(TARGET_MAPQ_CUTOFFS[-1]).dict
        total_reads = sum(targ_dists.values())
        within_200bp = sum(v for (k,v) in targ_dists.iteritems() if k <= 200)
        fraction_fragments_on_target = float(targ_dists.get(0,0)) / total_reads
        fraction_fragments_within_200bp_target = float(within_200bp) / total_reads
    else:
        fraction_fragments_on_target = None
        fraction_fragments_within_200bp_target = None

    metrics['fraction_fragments_on_target'] = fraction_fragments_on_target
    metrics['fraction_fragments_within_200bp_target'] = fraction_fragments_within_200bp_target

    # Target distance summary -- using all aligned reads (the above only uses Q60 reads)
    if not nearest_targ_dists.get_attrib('empty'):
        targ_dists = nearest_targ_dists.get_summarizer(TARGET_MAPQ_CUTOFFS[0]).dict
        total_reads = sum(targ_dists.values())
        within_200bp = sum(v for (k,v) in targ_dists.iteritems() if k <= 200)
        if 0 in targ_dists:
            fraction_fragments_on_target = float(targ_dists[0]) / total_reads
        else:
            fraction_fragments_on_target = None
        fraction_fragments_within_200bp_target = float(within_200bp) / total_reads
    else:
        fraction_fragments_on_target = None
        fraction_fragments_within_200bp_target = None

    metrics['fraction_fragments_on_target_mapq0'] = fraction_fragments_on_target
    metrics['fraction_fragments_within_200bp_target_mapq0'] = fraction_fragments_within_200bp_target

    metrics['far_chimera_rate'] = tk_stats.robust_divide(
            metrics['num_far_chimeras'], metrics['num_pos_chimeras'])
    metrics['inter_chimera_rate'] = tk_stats.robust_divide(
            metrics['num_inter_chimeras'], metrics['num_pos_chimeras'])
    metrics['same_dir_chimera_rate'] = tk_stats.robust_divide(
            metrics['num_same_dir_chimeras'], metrics['num_pos_chimeras'])
    metrics['outward_dir_chimera_rate'] = tk_stats.robust_divide(
            metrics['num_outward_dir_chimeras'], metrics['num_pos_chimeras'])

    return metrics

# DOESNT SUPPORT SINGLE END SEQUENCING YET!
def read_pair_info_iter(bam_in, target_regions, references):

    while(True):
        try:
            keep_going = True
            while(keep_going):
                r1 = bam_in.next()
                if not(r1.is_secondary):
                    keep_going = False

            keep_going = True
            while(keep_going):
                r2 = bam_in.next()
                if not(r2.is_secondary):
                    keep_going = False
        except StopIteration:
            break

        name1 = r1.qname.split()[0]
        name2 = r2.qname.split()[0]
        assert(name1 == name2)

        if r1.is_read1:
            read1 = r1
            read2 = r2
        else:
            read1 = r2
            read2 = r1

        if r1.tid == r2.tid and r2.pos < r1.pos:
            pos2 = r1
            pos1 = r2
        else:
            pos1 = r1
            pos2 = r2

        read_pair_info = {}
        read_pair_info['name'] = name1
        read_pair_info['is_duplicate'] = read1.is_duplicate or read2.is_duplicate

        read_pair_info['r1_mapped'] = not(read1.is_unmapped)
        read_pair_info['r2_mapped'] = not(read2.is_unmapped)

        read_pair_info['r1_mapq'] = read1.mapq
        read_pair_info['r2_mapq'] = read2.mapq

        if read1.tid == -1:
            ref1 = ''
        else:
            ref1 = references[read1.tid]

        read_pair_info['r1_chrom'] = ref1
        read_pair_info['r1_pos'] = read1.pos
        read_pair_info['r1_direct'] = not(read1.is_reverse)
        read_pair_info['r1_min_qual'] = tk_fasta.get_min_qual(read1.qual)
        read_pair_info['r1_seq_len'] = len(read1.seq) if read1.seq is not None else 0
        if read1.seq is not None and 'N' in read1.seq:
            read_pair_info['r1_contains_N'] = True
        else:
            read_pair_info['r1_contains_N'] = False

        if read2.tid == -1:
            ref2 = ''
        else:
            ref2 = references[read2.tid]

        read_pair_info['r2_chrom'] = ref2
        read_pair_info['r2_pos'] = read2.pos
        read_pair_info['r2_direct'] = not(read2.is_reverse)
        read_pair_info['r2_min_qual'] = tk_fasta.get_min_qual(read2.qual)
        read_pair_info['r2_seq_len'] = len(read2.seq) if read2.seq is not None else 0
        if read2.seq is not None and 'N' in read2.seq:
            read_pair_info['r2_contains_N'] = True
        else:
            read_pair_info['r2_contains_N'] = False

        read_raw_bc = tk_io.get_read_raw_barcode(read1)
        if read_raw_bc is None:
            read_raw_bc = ''
        read_bc = tk_io.get_read_barcode(read1)
        if read_bc is None:
            read_bc = ''
        read_bc_qual = tk_io.get_read_barcode_qual(read1)
        read_sample_index = tk_io.get_read_sample_index(read1)
        if read_sample_index is None:
            read_sample_index = ''
        read_sample_index_qual = tk_io.get_read_sample_index_qual(read1)

        read_pair_info['10X_raw_bc_len'] = len(read_raw_bc)
        read_pair_info['sample_index_len'] = len(read_sample_index)

        read_pair_info['10X_raw_bc'] = read_raw_bc
        read_pair_info['10X_called_bc'] = read_bc
        read_pair_info['sample_index'] = read_sample_index

        read_pair_info['10X_bc_min_qual'] = tk_fasta.get_min_qual(read_bc_qual)
        read_pair_info['10X_bc_mean_qual'] = tk_fasta.get_mean_qual(read_bc_qual)

        if read_sample_index:
            read_pair_info['sample_index_min_qual'] = tk_fasta.get_min_qual(read_sample_index_qual)
            read_pair_info['sample_index_mean_qual'] = tk_fasta.get_mean_qual(read_sample_index_qual)
        else:
            read_pair_info['sample_index_min_qual'] = None
            read_pair_info['sample_index_mean_qual'] = None

        if not(ref1 == ref2) or target_regions is None:
            read_pair_targ_dist = None
        else:
            left_pos = pos1.pos
            right_pos = pos2.pos + (len(pos2.seq) if pos2.seq is not None else 0)
            regions = target_regions.get(ref1, None)

            if regions is None:
                # No targets on this chrom -- distance is large positive number
                read_pair_targ_dist = 100000000
            else:
                read_pair_targ_dist = get_read_regions_dist(left_pos, right_pos, regions)

        read_pair_info['read_pair_targ_dist'] = read_pair_targ_dist

        r1_mapped, r1_soft_clipped = cigar_count_bases(read1)
        r2_mapped, r2_soft_clipped = cigar_count_bases(read2)

        read_pair_info['mapped_bases'] = r1_mapped + r2_mapped
        read_pair_info['soft_clipped_bases'] = r1_soft_clipped + r2_soft_clipped

        if (not read1.is_unmapped and not read2.is_unmapped) and read1.tid == read2.tid:
            read_pair_info['insert_length'] = pos2.aend - pos1.pos
        else:
            read_pair_info['insert_length'] = None

        # track QVs for trimmed bases
        r1_trimmed_qual = ''
        if read1.has_tag('QT'):
          r1_trimmed_qual = read1.get_tag('QT')
        read_pair_info['r1_trimmed_qual'] = r1_trimmed_qual

        qual1 = read1.qual
        if qual1 is None:
            qual1 = np.array([])

        qual2 = read2.qual
        if qual2 is None:
            qual2 = np.array([])

        read_pair_info['r1_q20_bases'] = tk_fasta.get_bases_qual(qual1, 20)
        ## also include barcode and trimmed bases in QV counts for read1
        if not read_bc_qual == None:
          read_pair_info['r1_q20_bases'] += tk_fasta.get_bases_qual(read_bc_qual, 20)
        read_pair_info['r1_q20_bases'] += tk_fasta.get_bases_qual(r1_trimmed_qual, 20)

        read_pair_info['r1_q30_bases'] = tk_fasta.get_bases_qual(qual1, 30)
        ## also include barcode and trimmed bases in QV counts for read1
        if not read_bc_qual == None:
          read_pair_info['r1_q30_bases'] += tk_fasta.get_bases_qual(read_bc_qual, 30)
        read_pair_info['r1_q30_bases'] += tk_fasta.get_bases_qual(r1_trimmed_qual, 30)

        read_pair_info['r2_q20_bases'] = tk_fasta.get_bases_qual(qual2, 20)
        read_pair_info['r2_q30_bases'] = tk_fasta.get_bases_qual(qual2, 30)

        read_pair_info['si_q20_bases'] = tk_fasta.get_bases_qual(read_sample_index_qual, 20)
        read_pair_info['si_q30_bases'] = tk_fasta.get_bases_qual(read_sample_index_qual, 30)

        read_pair_info['bc_q20_bases'] = tk_fasta.get_bases_qual(read_bc_qual, 20)
        read_pair_info['bc_q30_bases'] = tk_fasta.get_bases_qual(read_bc_qual, 30)

        yield read_pair_info


def cigar_count_bases(read):
    ''' Total clipped bases of mapped reads. Unmapped reads not counted '''
    mapped = 0
    soft_clipped = 0

    if read.is_unmapped:
        return (0, 0)

    for (code, bp) in read.cigar:
        if code == 4:
            soft_clipped += bp
        if code == 0 or code == 1:
            mapped += bp


    return (mapped, soft_clipped)


# ASSUMES READ IS PERFECTLY ALIGNED!
def get_read_regions_dist(left_pos, right_pos, regions):
    left_start, left_end, left_internal = regions.get_closest_region(left_pos)
    right_start, right_end, right_internal = regions.get_closest_region(right_pos)

    left_sign = math.copysign(1, left_pos - left_start)
    right_sign = math.copysign(1, right_pos - right_start)
    if left_internal or right_internal or not(left_sign == right_sign):
        return 0

    return min(min(abs(left_pos - left_start), abs(left_pos - left_end)), min(abs(right_pos - right_start), abs(right_pos - right_end)))
