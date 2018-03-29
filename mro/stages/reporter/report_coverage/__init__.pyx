#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Summarize genome coverage from BAM file
#

import json
import random
import os
import tenkit.pandas as p
import numpy as np
cimport numpy as np
import math
from libc cimport math
import itertools
import subprocess
import tenkit.safe_json
import tenkit.hdf5
import tenkit.coverage
import tenkit.stats as tk_stats
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.dict_utils as tk_dict_utils
import tenkit.regions as tk_regions
from tenkit.constants import COVERAGE_TRIM_TAIL
from tenkit.chunk_utils import generate_chrom_loci, pack_loci, adjust_start
import tenkit.reference

import martian

__MRO__ = """
stage REPORT_COVERAGE(
    in  bam    input,
    in  string reference_path,
    in  bed    targets_file,
    in  bed    confident_regions,
    in  h5     fragments,
    in  string restrict_locus,
    in  map    regions_of_interest,
    out json   summary,
    out csv    coverage_csv,
    out h5     coverage,
    out h5     coverage_off_target,
    out h5     target_coverage,
    out csv    target_csv,
    src py     "stages/reporter/report_coverage",
) split using (
    in  string[] loci,
    in  float  subsample_rate,
)
"""

def get_subsample_coverage(genome, targets):
    try:
        # If kitten is available it will tell us what the subsampled coverage should be for this genome.
        import kitten.constants
        if kitten.constants.SUBSAMPLE_COVERAGE.has_key(genome):
            return kitten.constants.SUBSAMPLE_COVERAGE[genome][targets is not None]
        else:
            return None

    except ImportError:
        # If kitten isn't available (e.g. for customer builds) then we turn off subsampling by setting a very high coverage
        return None

def split(args):
    bam_in = tk_bam.create_bam_infile(args.input)
    ref_len = sum(bam_in.lengths)
    mem_gb = max(4, int(math.ceil(6.0 * float(ref_len) / 3e9)))

    if args.targets_file is None:
        targets_path = None
    else:
        targets_path = args.targets_file

    mean_depth = tenkit.coverage.estimate_mean_coverage(targets_path, bam_in)
    martian.log_info("Estimated mean coverage: %f" % mean_depth)

    # Kitten has constants for the subsample target coverage depending on genome and whether we are targetting
    genome = tenkit.reference.get_genome(args.reference_path)
    ss_cov = get_subsample_coverage(genome, args.targets_file)

    if ss_cov is None or mean_depth < ss_cov:
        subsample_rate = None
    else:
        subsample_rate = ss_cov/mean_depth

    # Load pull-down targets
    if args.targets_file is None:
        target_regions = None
    else:
        with open(args.targets_file, 'r') as t_file:
            target_regions = tk_io.get_target_regions(t_file)

    if args.restrict_locus is None:
        all_loci = []
        for (chrom_name, chrom_size) in zip(bam_in.references, bam_in.lengths):
            all_loci.extend(generate_chrom_loci(target_regions, chrom_name, chrom_size, tenkit.constants.PARALLEL_LOCUS_SIZE))

        locus_sets = pack_loci(all_loci)
    else:
        locus_sets = [[args.restrict_locus]]

    chunk_defs = [{'loci': loci, 'subsample_rate': subsample_rate, '__mem_gb': mem_gb} for loci in locus_sets]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    chunk_outs = list(chunk_outs)
    cat_args = ['cat']
    cat_args.extend([chunk_out.high_coverage_excluded_bed for chunk_out in chunk_outs])
    with open(outs.high_coverage_excluded_bed + "tmp.bed",'w') as bed_out:
        subprocess.check_call(cat_args,stdout=bed_out)
    with open(outs.high_coverage_excluded_bed, 'w') as bed_out:
        subprocess.check_call(['bedtools','merge','-i',outs.high_coverage_excluded_bed+"tmp.bed"],stdout=bed_out)

    top_level_keys = ['summary_depth_info_deduped_subsampled',
                      'summary_depth_info',
                      'summary_depth_info_subsampled',
                      'summary_depth_info_deduped',
                      'summary_depth_info_mapq30_deduped',

                      'confident_summary_depth_info_deduped_subsampled',
                      'confident_summary_depth_info',
                      'confident_summary_depth_info_subsampled',
                      'confident_summary_depth_info_deduped',
                      'confident_summary_depth_info_mapq30_deduped',

                      'summary_bc_depth_info',
                      'confident_summary_bc_depth_info',
                      'target_info',
                      'target_info_subsampled',
                      'mapq30_info']
    top_level_depths = [1] * len(top_level_keys)
    output_dict = {}
    for key in top_level_keys:
        output_dict[key] = {}

    total_male = []
    male_coverage = []
    autosomal_total = []
    autosomal_coverage = []
    for (chunk_def, chunk_out) in zip(chunk_defs, chunk_outs):
        in_file = open(chunk_out.summary)
        top_input_dict = json.loads(in_file.readline())
        if 'male_coverage' in top_input_dict:
            male_coverage.append(top_input_dict['male_coverage'])
            total_male.append(top_input_dict['male_loci'])
        if 'autosomal_coverage' in top_input_dict:
            autosomal_coverage.append(top_input_dict['autosomal_coverage'])
            autosomal_total.append(top_input_dict['autosomal_loci'])

        for (key, depth) in zip(top_level_keys, top_level_depths):
            input_dict = top_input_dict.get(key, {})
            if input_dict is None:
                input_dict = {}

            output_dict[key] = tk_dict_utils.add_dicts(output_dict[key], input_dict, depth)

    output_dict['detected_sex'] = None

    # Detect sex if possible
    if len(male_coverage) > 0 and len(autosomal_coverage) > 0:
        male_coverage = np.ndarray(shape=(len(male_coverage),),buffer = np.array(male_coverage),dtype=float)
        total_male = np.ndarray(shape=(len(total_male),),buffer = np.array(total_male),dtype=float)
        male_coverage = np.average(male_coverage, weights = total_male)
        output_dict['male_coverage'] = male_coverage
        autosomal_coverage = np.ndarray(shape=(len(autosomal_coverage),),buffer = np.array(autosomal_coverage),dtype=float)
        autosomal_total = np.ndarray(shape=(len(autosomal_total),),buffer = np.array(autosomal_total),dtype=float)
        autosomal_coverage = np.average(autosomal_coverage, weights = autosomal_total)
        output_dict['autosomal_coverage'] = autosomal_coverage
        # assuming autosomes are diploid, and sex chromosomes are haploid,
        # male_copy = male_coverage / (autosomal_coverage / 2)
        male_copy = tk_stats.robust_divide(male_coverage, tk_stats.robust_divide(autosomal_coverage, 2.0))
        output_dict['male_chr_copy'] = round(male_copy) if not math.isnan(male_copy) else float('nan')
        if args.sex is None:
            if output_dict['male_chr_copy'] >= 1:
                outs.sex = "male"
            else:
                outs.sex = "female"
        if output_dict.get('male_chr_copy') is not None and output_dict.get('male_chr_copy') >= 1:
            output_dict['detected_sex'] = 'male'
        else:
            output_dict['detected_sex'] = 'female'

    if args.sex is not None:
        outs.sex = args.sex
    else:
        outs.sex = output_dict['detected_sex']

    # Check that we generated the right number of points
    # start with the non-ambiguous regions, ignoring chroms that are not in the BAM header
    unambiguous = tenkit.reference.get_unambiguous_regions(args.reference_path)
    bam_in = tk_bam.create_bam_infile(args.input)
    chroms = set(bam_in.references).intersection(unambiguous)
    # TODO all this multi-chromosome region arithmetic should probably be handled by tenkit.regions
    regions = {chrom: unambiguous[chrom] for chrom in chroms}

    # factor in restrict_locus
    if args.restrict_locus is not None:
        (chrom, start, end) = tk_io.get_locus_info(args.restrict_locus)
        restrict = {chrom: tk_regions.Regions(regions=[(start, end)])}
        chroms = chroms.intersection(restrict)
        regions = {chrom: regions[chrom].intersect(restrict[chrom]) for chrom in chroms}

    # factor in targets
    if args.targets_file is not None:
        with open(args.targets_file, 'r') as t_file:
            targets = tk_io.get_target_regions(t_file)
        chroms = chroms.intersection(targets)
        regions = {chrom: regions[chrom].intersect(targets[chrom]) for chrom in chroms}

    # Generate a csv file of depths
    # figure out the 99.9th percentile of coverage, and set that as the top bin
    cov_pairs = sorted(list((int(cov), num_sites) for (cov, num_sites) in output_dict['summary_depth_info'].iteritems()))
    cov_level = np.array([int(cov) for (cov, num_sites) in cov_pairs])
    num_sites = np.array([num_sites for (cov, num_sites) in cov_pairs])
    cum_sites = np.cumsum(num_sites)

    top_indexes = np.where(cum_sites >= num_sites.sum() * 0.995)[0]
    if len(top_indexes) > 0:
        top_index = top_indexes[0]
    else:
        top_index = len(cov_level) - 1

    top_coverage = cov_level[top_index]

    # Generate a data frame with cov levels from 0 to top_coverage
    coverage_bins = np.arange(top_coverage+1)
    cov_df = p.DataFrame({"coverage":coverage_bins})
    subsample_cov = get_subsample_coverage(tenkit.reference.get_genome(args.reference_path), args.targets_file)

    if subsample_cov is None:
        column_names = {
                'confident_summary_depth_info_deduped' : 'counts' }
        confident_keys = ('confident_summary_depth_info_deduped',)
    else:
        column_names = {
                'confident_summary_depth_info_deduped' : 'counts',
                'confident_summary_depth_info_deduped_subsampled': 'counts_subsampled_%dx' % round(subsample_cov)}
        confident_keys = ('confident_summary_depth_info_deduped_subsampled', 'confident_summary_depth_info_deduped')

    for top_level_summary in top_level_keys:
        if top_level_summary in confident_keys:

            cov_counts = np.zeros(top_coverage+1)
            for (cov_level, count) in output_dict[top_level_summary].iteritems():
                cov_level = int(cov_level)
                if cov_level <= top_coverage:
                    cov_counts[cov_level] = count
                else:
                    cov_counts[top_coverage] += count

            column_name = column_names[top_level_summary]
            cov_df[column_name] = cov_counts

            cum_column_name = 'cumulative_' + column_name
            cov_df[cum_column_name] = np.cumsum(cov_counts)

    cov_df.to_csv(outs.coverage_csv)

    # Combine the coverage hdf5 files
    cov_h5 = outs.coverage

    if cov_h5:
        in_files = [out.coverage for out in chunk_outs]
        tenkit.hdf5.combine_data_frame_files(cov_h5, in_files)

        # tabix index the coverage file
        tenkit.hdf5.create_tabix_index(cov_h5, 'chrom', 'pos', 'pos')

    if args.targets_file and outs.target_coverage:
        in_files = [out.target_coverage for out in chunk_outs]
        tenkit.hdf5.combine_data_frame_files(outs.target_coverage, in_files)

        # write csv
        df = tenkit.hdf5.read_data_frame(outs.target_coverage)
        df.to_csv(outs.target_csv)
    else:
        outs.target_coverage = None

    # combine the off-target coverage for exome
    if args.targets_file:
        in_files = [out.coverage_off_target for out in chunk_outs]
        tenkit.hdf5.combine_data_frame_files(outs.coverage_off_target, in_files)
        tenkit.hdf5.create_tabix_index(outs.coverage_off_target, 'chrom', 'pos', 'pos')
    else:
        outs.coverage_off_target = None

    # Write final summary
    output_file = open(outs.summary, 'w')
    output_file.write(tenkit.safe_json.safe_jsonify(output_dict))

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    bam_in = tk_bam.create_bam_infile(args.input)
    subsample_rate = args.subsample_rate

    output_file = open(outs.summary, 'w')

    # Set random seed to get deterministic subsampling
    random.seed(0)

    # Load pull-down targets
    if args.targets_file is None:
        target_regions_genome = None
    else:
        with open(args.targets_file, 'r') as t_file:
            target_regions_genome = tk_io.get_target_regions(t_file)

    def get_target_regions(chrom):
        if target_regions_genome is None:
            return None
        # An empty array means we are doing targetting but there are no targets on this chromosome
        return target_regions_genome.get(chrom, tk_regions.Regions([]))

    # Load genome confident regions file
    if args.confident_regions is None:
        conf_regions_genome = None
    else:
        with open(args.confident_regions, 'r') as c_file:
            conf_regions_genome = tk_io.get_target_regions(c_file)

    def get_conf_regions(chrom):
        if conf_regions_genome == None:
            return tk_regions.Regions([(0,10**10)])

        return tk_regions.Regions(conf_regions_genome.get(chrom, []))

    # Get non-N bases of genome
    unambiguous_regions_genome = tenkit.reference.get_unambiguous_regions(args.reference_path)
    def get_unambiguous_regions(chrom):
        return unambiguous_regions_genome.get(chrom, tk_regions.Regions([]))

    # Load additional regions of interest, if specified
    misc_regions_all_chroms = {}
    if args.regions_of_interest is not None:
        for (name, region) in args.regions_of_interest.iteritems():
            if os.path.exists(region):
                with open(region) as r:
                    misc_regions_all_chroms[name] = tk_io.get_target_regions(r)

    def get_misc_regions(chrom):
        misc_regions = {}
        for (name, regions) in misc_regions_all_chroms.iteritems():
            if chrom in regions:
                misc_regions[name] = regions[chrom]
            else:
                misc_regions[name] = tk_regions.Regions()
        return misc_regions

    # Adjust start_pos to come after target regions if they intersect on
    chunk_size = int(1e6)
    chunks = []

    for (chrom, start, end) in (tk_io.get_locus_info(l) for l in args.loci):
        starts = [adjust_start(chrom, s, target_regions_genome) for s in range(start, end, chunk_size)]
        starts[0] = start
        ends = starts[1:] + [end]
        _chunks = zip(itertools.repeat(chrom), starts, ends)
        chunks.extend(_chunks)

    target_info_acc = {}
    target_info_ss_acc = {}

    hist_acc = {}
    hist_dd_acc = {}
    hist_mq30_dd_acc = {}
    hist_ss_acc = {}
    hist_ss_dd_acc = {}

    conf_hist_acc = {}
    conf_hist_dd_acc = {}
    conf_hist_mq30_dd_acc = {}
    conf_hist_ss_acc = {}
    conf_hist_ss_dd_acc = {}

    bc_depth_acc = {}
    bc_conf_depth_acc = {}

    if args.fragments is not None:
        frag_reader = tenkit.hdf5.DataFrameReader(args.fragments)
    else:
        frag_reader = None
    male_coverage = []
    autosomal_coverage = []
    male_total = []
    male_chromosomes = tenkit.reference.load_male_chromosomes(args.reference_path)
    if male_chromosomes is None:
        male_chromosomes = set()
    autosomal_chromosomes = tenkit.reference.load_autosomal_chromosomes(args.reference_path)
    if autosomal_chromosomes is None:
        autosomal_chromosomes = set()
    autosomal_total = []

    for (chrom, start, end) in chunks:
        target_regions = get_target_regions(chrom)
        conf_regions = get_conf_regions(chrom)
        unambiguous = get_unambiguous_regions(chrom)
        misc_regions = get_misc_regions(chrom)

        if frag_reader is not None:
            # Get fragment coverage data
            fragments = frag_reader.query((chrom, start, end), query_cols=['chrom', 'start_pos', 'end_pos'])
            frag_tuples = zip(fragments.start_pos.values, fragments.end_pos.values)

            r = get_bc_depth_info(frag_tuples, chrom, start, end, target_regions, conf_regions, unambiguous)
            bc_depth_df, bc_summary_depth, bc_conf_summary_depth = r
            bc_depth_acc = tk_dict_utils.add_dicts(bc_depth_acc, bc_summary_depth, 1)
            bc_conf_depth_acc = tk_dict_utils.add_dicts(bc_conf_depth_acc, bc_conf_summary_depth, 1)

        def filter_bam(chrom, start, end, read_filter=lambda x: True):
            return tk_bam.filter_bam(bam_in.fetch(chrom, start, end), remove_unmapped=True, read_filter=read_filter)

        # Filter secondary reads
        # Filter unmapped reads - they can get through because their mate is mapped
        all_reads = filter_bam(chrom, start, end)

        # All non-secondary reads
        r = get_depth_info(all_reads, chrom, start, end, target_regions, conf_regions, unambiguous, outs.high_coverage_excluded_bed, args.high_coverage_threshold, misc_regions)
        depth_df, summary_depth_info, conf_summary_depth_info, target_info, target_cov_df = r

        # Add in BC coverage
        if frag_reader is not None:
            depth_df['bc_coverage'] = bc_depth_df['bc_coverage']
        target_info_acc = tk_dict_utils.add_dicts(target_info_acc, target_info, 1)
        hist_acc = tk_dict_utils.add_dicts(hist_acc, summary_depth_info, 1)
        conf_hist_acc = tk_dict_utils.add_dicts(conf_hist_acc, conf_summary_depth_info, 1)

        # Deduped reads
        deduped_reads = filter_bam(chrom, start, end, read_filter=lambda x: not x.is_duplicate)
        r_deduped = get_depth_info(deduped_reads, chrom, start, end, target_regions, conf_regions, unambiguous, None, args.high_coverage_threshold, None)
        dd_df, dd_summary, conf_dd_summary, dd_target_info, dd_target_cov = r_deduped

        hist_dd_acc = tk_dict_utils.add_dicts(hist_dd_acc, dd_summary, 1)
        conf_hist_dd_acc = tk_dict_utils.add_dicts(conf_hist_dd_acc, conf_dd_summary, 1)
        depth_df['coverage_deduped'] = dd_df['coverage']

        # Deduped, mapq30 reads
        mq30_deduped_reads = filter_bam(chrom, start, end, read_filter=lambda x: (not x.is_duplicate) and x.mapq >= 30)
        r_mq30_deduped = get_depth_info(mq30_deduped_reads, chrom, start, end, target_regions, conf_regions, unambiguous, None, args.high_coverage_threshold, None)
        mq_dd_df, mq_dd_summary, conf_mq_dd_summary, _, _ = r_mq30_deduped

        hist_mq30_dd_acc = tk_dict_utils.add_dicts(hist_mq30_dd_acc, mq_dd_summary, 1)
        conf_hist_mq30_dd_acc = tk_dict_utils.add_dicts(conf_hist_mq30_dd_acc, conf_mq_dd_summary, 1)
        depth_df['mapq30_coverage_deduped'] = mq_dd_df['coverage']

        if subsample_rate and subsample_rate <= 1.0:
            # Subsampled reads
            subsampled_reads = filter_bam(chrom, start, end, read_filter=lambda x: random.random() < subsample_rate)
            r_ss = get_depth_info(subsampled_reads, chrom, start, end, target_regions, conf_regions, unambiguous, None, args.high_coverage_threshold, None)
            ss_df, ss_summary, conf_ss_summary, ss_target_info, ss_target_cov = r_ss

            target_info_ss_acc = tk_dict_utils.add_dicts(target_info_ss_acc, ss_target_info, 1)
            hist_ss_acc = tk_dict_utils.add_dicts(hist_ss_acc, ss_summary, 1)
            conf_hist_ss_acc = tk_dict_utils.add_dicts(conf_hist_ss_acc, conf_ss_summary, 1)

            # Subsampled deduped reads
            ss_dd_reads = filter_bam(chrom, start, end, read_filter=lambda x: random.random() < subsample_rate and not x.is_duplicate)
            r_deduped_ss = get_depth_info(ss_dd_reads, chrom, start, end, target_regions, conf_regions, unambiguous, None, args.high_coverage_threshold, None)
            ss_dd_df, ss_dd_summary, conf_ss_dd_summary, ss_dd_target_info, ss_dd_target_cov = r_deduped_ss

            hist_ss_dd_acc = tk_dict_utils.add_dicts(hist_ss_dd_acc, ss_dd_summary, 1)
            conf_hist_ss_dd_acc = tk_dict_utils.add_dicts(conf_hist_ss_dd_acc, conf_ss_dd_summary, 1)

            depth_df['coverage_subsampled'] = ss_df['coverage']
            depth_df['coverage_deduped_subsampled'] = ss_dd_df['coverage']

        # Write this chunk of the coverage track
        # For exome, write two versions (with and without off-target regions),
        # Using the smaller one as the main track

        if outs.coverage:
            if args.targets_file:
                depth_df_targeted = depth_df[(depth_df.target_regions == 1)]
                sub_df = depth_df_targeted
                tenkit.hdf5.append_data_frame(outs.coverage, depth_df_targeted)
                tenkit.hdf5.append_data_frame(outs.coverage_off_target, depth_df)
            else:
                sub_df = depth_df
                tenkit.hdf5.append_data_frame(outs.coverage, depth_df)
                outs.coverage_off_target = None
        if chrom in male_chromosomes:
            cov_male = np.mean(sub_df['coverage'])
            if not math.isnan(cov_male):
                male_coverage.append(cov_male)
                male_total.append(len(sub_df))

        if chrom in autosomal_chromosomes:
            cov_autosomal = np.mean(sub_df['coverage'])
            if not math.isnan(cov_autosomal):
                autosomal_coverage.append(cov_autosomal)
                autosomal_total.append(len(sub_df))

        # For exome, write the target-level coverage
        if outs.target_coverage and args.targets_file:
            tenkit.hdf5.append_data_frame(outs.target_coverage, target_cov_df)

    if len(male_coverage) > 0:
        male_coverage = np.average(male_coverage, weights = male_total)
    if len(autosomal_coverage) > 0:
        autosomal_coverage = np.average(autosomal_coverage, weights = autosomal_total)
    results = {}

    # On-target bases info
    results['target_info'] = target_info_acc
    results['target_info_subsampled'] = target_info_ss_acc

    results['pos_coverage_info'] = {}
    results['pos_coverage_info_deduped'] = {}

    if frag_reader is not None:
        results['summary_bc_depth_info'] = bc_depth_acc
        results['confident_summary_bc_depth_info'] = bc_conf_depth_acc
    else:
        results['summary_bc_depth_info'] = {}
        results['confident_summary_bc_depth_info'] = {}
    if len(male_total) > 0:
        results['male_coverage'] =male_coverage
        results['male_loci'] = np.sum(male_total)
    if len(autosomal_total) > 0:
        results['autosomal_coverage'] = autosomal_coverage
        results['autosomal_loci'] = np.sum(autosomal_total)

    results['summary_depth_info'] = hist_acc
    results['summary_depth_info_deduped'] = hist_dd_acc
    results['summary_depth_info_mapq30_deduped'] = hist_mq30_dd_acc
    results['gc_vs_depth_info'] = {}
    results['gc_vs_depth_info_deduped'] = {}

    results['pos_coverage_info_subsampled'] = {}
    results['pos_coverage_info_deduped_subsampled'] = {}
    results['summary_depth_info_subsampled'] = hist_ss_acc
    results['summary_depth_info_deduped_subsampled'] = hist_ss_dd_acc
    results['gc_vs_depth_info_subsampled'] = {}
    results['gc_vs_depth_info_deduped_subsampled'] = {}

    results['confident_summary_depth_info'] = conf_hist_acc
    results['confident_summary_depth_info_deduped'] = conf_hist_dd_acc
    results['confident_summary_depth_info_mapq30_deduped'] = conf_hist_mq30_dd_acc

    results['confident_summary_depth_info_subsampled'] = conf_hist_ss_acc
    results['confident_summary_depth_info_deduped_subsampled'] = conf_hist_ss_dd_acc

    # MQ30 stats for each region of interest
    mapq30_rois = {}

    mapq30_10x = depth_df['mapq30_coverage_deduped'] > 10
    mapq30_20x = depth_df['mapq30_coverage_deduped'] > 20
    mapq30_rois['mapq30_bases_10x_coverage_all'] = mapq30_10x.sum()
    mapq30_rois['mapq30_bases_20x_coverage_all'] = mapq30_20x.sum()

    for roi in misc_regions_all_chroms:
        mask = (depth_df[roi] == 1)
        mapq30_rois['mapq30_bases_10x_coverage_' + roi] = (mapq30_10x & mask).sum()
        mapq30_rois['mapq30_bases_20x_coverage_' + roi] = (mapq30_20x & mask).sum()

    results['mapq30_info'] = mapq30_rois

    output_file.write(tenkit.safe_json.safe_jsonify(results))
    output_file.close()

def zcat(arrays, dtype=None):
    ''' Wrapper for np.concatenate that accepts 0-length input '''
    if len(arrays) == 0:
        return np.array([], dtype=dtype)
    else:
        return np.concatenate(arrays)

def region_mask(positions, targets):
    ''' Make a mask array parallel with the positions array indicating whether each position falls within
        the regions intervals'''

    mask = np.zeros(positions.shape, dtype=np.bool)
    if len(positions) == 0:
        return mask

    regs = targets.overlapping_regions(positions[0], positions[-1])

    starts = [s for (s,e) in regs]
    ends =   [e for (s,e) in regs]
    lbounds = np.searchsorted(positions, starts)
    rbounds = np.searchsorted(positions, ends)

    for (l,r) in zip(lbounds, rbounds):
        mask[l:r] = True

    return mask

def get_depth_info(read_iter, chrom, long cstart, long cend, targets, conf_regions, unambiguous, high_coverage_bed, int high_cov_threshold, misc_regions):
    cdef np.ndarray depths = np.zeros(cend-cstart, np.int32)
    cdef np.ndarray mapqv = np.zeros(cend-cstart, np.uint32)
    cdef long start, end, rstart, rend, index, start_region, end_region

    for read in read_iter:
        pos = read.pos
        rstart = max(pos, cstart)

        # Increment to the end of the window or the end of the
        # alignment, whichever comes first
        rend = min(read.aend, cend)
        depths[(rstart-cstart):(rend-cstart)] += 1
        mapqv[(rstart-cstart):(rend-cstart)] += read.mapq

    cdef bint in_region
    cdef int depth
    if high_coverage_bed is not None:
        with open(high_coverage_bed, 'a') as bed:
            in_region = False
            start_region = 0
            for (index, depth) in enumerate(depths):
                if not(in_region) and depth <= high_cov_threshold:
                    in_region = True
                    start_region = cstart + index
                if (in_region and depth > high_cov_threshold):
                    in_region = False
                    end_region = cstart + index
                    bed.write(chrom+"\t"+str(start_region)+"\t"+str(end_region)+"\n")
                if index == len(depths) - 1:
                    end_region = cstart + index + 2
                    bed.write(chrom+"\t"+str(start_region)+"\t"+str(end_region)+"\n")

    positions = np.arange(cstart, cend, dtype=np.int32)

    # Divide mapq accumulator by coverage, round, and store as byte
    with np.errstate(divide='ignore'):
        mapqv = mapqv / depths
    # store positions with 0 coverage as -1
    mapqv[np.isnan(mapqv)] = -1
    mapqv = np.round(mapqv).astype(np.int8)

    target_info = {}
    target_cov_df = None
    target_mask = None

    cdef list valid_targets
    if targets is not None:
        total_bases = depths.sum()

        # Measure the mean coverage of each target
        valid_targets = [(start, end) for (start, end) in targets if start >= cstart and end <= cend]
        starts = np.array([start for (start,end) in valid_targets], dtype=np.int32)
        ends = np.array([end for (start,end) in valid_targets], dtype=np.int32)
        means = np.array([depths[(start-cstart):(end-cstart)].mean() for (start,end) in valid_targets], dtype=np.float32)
        dpcvs = np.array([compute_dpcv(depths[(start-cstart):(end-cstart)]) for (start,end) in valid_targets], dtype=np.float32)

        target_cov_df = p.DataFrame({'chrom':chrom, 'start':starts, 'end':ends, 'mean': means, 'dpcv': dpcvs})
        target_cov_df['chrom'] = target_cov_df['chrom'].astype('category')

        target_mask = region_mask(positions, targets)

        on_target_bases = depths[target_mask].sum()
        target_info = { 'on_target_bases' : on_target_bases, 'total_bases' : total_bases }

    # for all depth calculations, ignore ambiguous bases and off-target bases
    depth_mask = region_mask(positions, unambiguous)
    if target_mask is not None:
        depth_mask = depth_mask & target_mask

    # overall summary depth info
    cdef dict summary_depth_info = {}
    cdef int d

    if depth_mask.sum() and depths.size > 0:
        depth_hist = np.bincount(depths[depth_mask])
        for d in xrange(depth_hist.size):
            if depth_hist[d] > 0:
                summary_depth_info[int(d)] = int(depth_hist[d])

    # conf regions summary depth info
    conf_mask = region_mask(positions, conf_regions) & depth_mask
    cdef dict conf_summary_depth_info = {}

    if conf_mask.sum() > 0 and depths.size > 0:
        depth_hist = np.bincount(depths[conf_mask])
        for d in xrange(depth_hist.size):
            if depth_hist[d] > 0:
                conf_summary_depth_info[int(d)] = int(depth_hist[d])

    depth_df = p.DataFrame({ "chrom": chrom, "pos": positions, "mapq": mapqv, "coverage": depths})
    depth_df['chrom'] = depth_df['chrom'].astype('category')

    if target_mask is not None:
        depth_df['target_regions'] = target_mask.astype('int8')

    if misc_regions is not None:
        for (name, regions) in misc_regions.iteritems():
            # note -- some of the names might be unicode (?), so convert to string before adding column to DF
            depth_df[str(name)] = region_mask(positions, regions).astype('int8') # int8 has better HDF5 compression?

    return depth_df, summary_depth_info, conf_summary_depth_info, target_info, target_cov_df

def get_bc_depth_info(fragment_iter, chrom, long cstart, long cend, targets, conf_regions, unambiguous):

    depths = np.zeros(cend-cstart, np.int32)
    positions = np.arange(cstart, cend, dtype=np.int32)
    target_mask = None

    cdef long pos, frag_end, rstart, rend
    for fragment in fragment_iter:
        (pos, frag_end) = fragment
        rstart = max(pos, cstart)

        # Increment to the end of the window or the end of the
        # alignment, whichever comes first
        rend = min(frag_end, cend)
        depths[(rstart-cstart):(rend-cstart)] += 1

    if targets is not None:
        # Measure the mean coverage of each target
        target_mask = region_mask(positions, targets)

    # for all depth calculations, ignore ambiguous bases and off-target bases
    depth_mask = region_mask(positions, unambiguous)
    if target_mask is not None:
        depth_mask = depth_mask & target_mask

    # overall summary depth info
    cdef dict summary_depth_info = {}

    cdef int d
    if depth_mask.sum() > 0 and depths.size > 0:
        depth_hist = np.bincount(depths[depth_mask])
        for d in xrange(depth_hist.size):
            if depth_hist[d] > 0:
                summary_depth_info[int(d)] = int(depth_hist[d])

    # conf regions summary depth info
    conf_mask = region_mask(positions, conf_regions) & depth_mask
    cdef dict conf_summary_depth_info = {}

    if conf_mask.sum() > 0 and depths.size > 0:
        depth_hist = np.bincount(depths[conf_mask])
        for d in xrange(depth_hist.size):
            if depth_hist[d] > 0:
                conf_summary_depth_info[int(d)] = int(depth_hist[d])

    depth_df = p.DataFrame({"bc_coverage": depths})

    if target_mask is not None:
        depth_df['target_regions'] = target_mask

    return depth_df, summary_depth_info, conf_summary_depth_info

def get_depth_positional_cv(depth_hist, trim_tail):
    total_count = sum(depth_hist.values())
    cutoff_count = total_count*trim_tail
    seen_count = 0
    for depth in sorted(depth_hist.iterkeys(), reverse=True):
        seen_count += depth_hist[depth]
        if seen_count >= cutoff_count:
            cutoff = depth
            break
    trimmed_hist = {x: y for (x, y) in depth_hist.iteritems() if x <= cutoff}
    mean_val, var_val = tk_stats.mean_var_from_counts(trimmed_hist)
    if mean_val > var_val:
        return float('NaN')
    return tk_stats.robust_divide(np.sqrt(var_val - mean_val), mean_val)

def compute_dpcv(cov_pts):
    if len(cov_pts) < 1000:
        return float('NaN')

    depth_hist = np.bincount(cov_pts)
    summary_depth_info = {}
    for d in xrange(depth_hist.size):
        if depth_hist[d] > 0:
            summary_depth_info[int(d)] = int(depth_hist[d])

    return get_depth_positional_cv(summary_depth_info, COVERAGE_TRIM_TAIL)
