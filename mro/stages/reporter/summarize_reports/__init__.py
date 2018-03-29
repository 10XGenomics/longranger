#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Summarize results of analytics pipeline
#
import copy
import json
import numpy
import math
import csv
import itertools
import os

import tenkit.bam
import tenkit.alarms
import tenkit.safe_json
import tenkit.hdf5
from tenkit.constants import COVERAGE_TRIM_TAIL, CUSTOMER_LEFT_TAIL_COVERAGE

import tenkit
import tenkit.stats as tk_stats
import tenkit.reference

import martian

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  string sample_id,
    in  string reference_path,
    in  bed    targets,
    in  int    trim_length,
    in  json   duplicate_summary,
    in  json   basic_results,
    in  h5     barcode_counts,
    in  json   filter_barcodes_results,
    in  json   coverage_results,
    in  h5     coverage_details,
    in  json   variant_results,
    in  json   sv_results,
    in  int    sv_min_call_qv_wgs,
    in  int    sv_min_call_qv_target,
    in  json   single_partition_results,
    in  json   length_mass_results,
    in  bam    bam_file,
    in  json   lot_info,
    in  json   downsample_info,
    out json   analysis_params,
    out json   summary,
    out csv    summary_cs,
    out json   alarms,
    out txt    alarms_summary,
    src py     "stages/reporter/summarize_reports",
)
"""

# Map of internal metrics keys to customer facing keys
CS_METRICS_KEY_MAP = [

    # Long Ranger version
    [ "longranger_version",             "longranger_version" ],

    # Instrument IDs
    [ "instrument_ids",                 "instrument_ids" ],

    # GEM Performance
    [ "inferred_number_gems",           "gems_detected" ],
    [ "mean_bp_per_bc",                 "mean_dna_per_gem" ],
    [ "correct_bc_rate",                "bc_on_whitelist" ],
    [ "mean_bc_qual",                   "bc_mean_qscore" ],
    [ "n90_n50_reads_per_fragment",     "n50_linked_reads_per_molecule" ],


    # Note: in LR 2.0 we reported loaded_mass_ng. This metric will double count denatured DNA, so we will end up telling the customer that they've
    # loaded too much DNA. In LR 2.1 we've switching to the 'corrected' metric which divides by 1.6 to account for the observed denaturation rate.
    #[ "loaded_mass_ng",                 "loaded_mass_ng" ],
    [ "corrected_loaded_mass_ng",       "corrected_loaded_mass_ng" ],

    # Phasing
    [ "fract_snps_phased",              "snps_phased" ],
    [ "fract_genes_lt_100kb_phased",    "genes_phased_lt_100kb" ],
    [ "longest_phase_block",            "longest_phase_block" ],
    [ "N50_phase_block",                "n50_phase_block" ],

    # Input DNA
    [ "inferred_lw_mean_length",        "molecule_length_mean" ],
    [ "n90_lw_stddev_fragment_length",  "molecule_length_stddev" ],

    # Sequencing
    [ "num_reads",                      "number_reads" ],
    [ "median_insert_size",             "median_insert_size" ],
    [ "mean_depth",                     "mean_depth" ],
    [ "zero_cov_fract",                 "zero_coverage" ],
    [ "unmapped_fract",                 "mapped_reads" ],
    [ "dup_fraction",                   "pcr_duplication" ],
    [ "fraction_on_target",             "on_target_bases" ],

    # Quality scores broken down by read
    [ "r1_q20_bases_fract",             "r1_q20_bases_fract" ],
    [ "r1_q30_bases_fract",             "r1_q30_bases_fract" ],
    [ "r2_q20_bases_fract",             "r2_q20_bases_fract" ],
    [ "r2_q30_bases_fract",             "r2_q30_bases_fract" ],
    [ "si_q20_bases_fract",             "si_q20_bases_fract" ],
    [ "si_q30_bases_fract",             "si_q30_bases_fract" ],
    [ "bc_q20_bases_fract",             "bc_q20_bases_fract" ],
    [ "bc_q30_bases_fract",             "bc_q30_bases_fract" ],
]

def main(args, outs):
    # Analysis parameters
    lead_trim = args.trim_length
    analysis_params = {}
    analysis_params['lead_trim'] = lead_trim
    analysis_params['analysis_version'] = martian.get_pipelines_version()
    analysis_params_output_file = open(outs.analysis_params, 'w')
    analysis_params_output_file.write(json.dumps(analysis_params))
    analysis_params_output_file.close()

    # Summary metrics
    summary_metrics = {}

    # Add longranger version to summary so we get it everywhere the summary.json goes
    # We'll also pipe it to the customer csv
    summary_metrics['longranger_version'] = martian.get_pipelines_version()

    basic_metrics = json.load(open(args.basic_results, 'r'))
    for (k,v) in basic_metrics.items():
        summary_metrics[k] = v

    # Get the set of instrument IDs observed in the BAM file
    if args.bam_file is not None:
        instrument_ids = get_instrument_ids(args.bam_file)
        summary_metrics['instrument_ids'] = ";".join(instrument_ids)
    else:
        summary_metrics['instrument_ids'] = ''

    # Copy over single_partition results
    sp_metrics = json.load(open(args.single_partition_results, 'r'))
    for (k,v) in sp_metrics.items():
        summary_metrics[k] = v

    # Load the duplicate summary results
    # only include the overall dup rate in the customer metrics
    dup_metrics = json.load(open(args.duplicate_summary, 'r'))

    key = 'full_use_bcs'
    dup_counts = dup_metrics[key]

    if dup_counts is None:
        key = 'full_ignore_bcs'
        dup_counts = dup_metrics[key]

    mean_tag = "mean_dup_rate"
    sd_tag = "sd_dup_rate"
    optical_tag = "optical_dup_rate"
    dup_frac_tag = "dup_fraction"

    if dup_counts:
        dd = { int(k):v for (k,v) in dup_counts.items() }
        n_dups = sum([ v*(k-1) for (k,v) in dd.items() if k > 1 ])
        n_non_dups = sum(dd.values())

        mean_dup_rate = tk_stats.robust_divide(float(n_dups + n_non_dups), n_non_dups)
        summary_metrics[mean_tag] = mean_dup_rate

        # Customer facing dup rate on 0 - 1 scale
        summary_metrics[dup_frac_tag] = (mean_dup_rate - 1.0) / mean_dup_rate

        optical_dup_count = dup_metrics['optical_' + key]['count']
        summary_metrics[optical_tag] = tk_stats.robust_divide(float(optical_dup_count), n_non_dups)

        sd_terms = [ (k-mean_dup_rate)**2.0 * v for (k,v) in dd.items() ]
        sd_dup_rate = math.sqrt(tk_stats.robust_divide(sum(sd_terms), sum(dd.values())))
        summary_metrics[sd_tag] = sd_dup_rate
    else:
        summary_metrics[dup_frac_tag] = 0.0
        summary_metrics[mean_tag] = 1.0
        summary_metrics[sd_tag] = 0.0
        summary_metrics[optical_tag] = 0.0

    # Load the bias results
    bias_results = json.load(open(args.coverage_results, 'r'))
    summary_depth_info = bias_results['summary_depth_info']
    mean_depth, median_depth, zero_cov_fract = get_depth_info(summary_depth_info)
    on_target_bases = get_on_target_bases(summary_depth_info)
    depth_positional_cv = get_depth_positional_cv(summary_depth_info, COVERAGE_TRIM_TAIL)

    summary_depth_info = bias_results['summary_depth_info_deduped']
    mean_depth_deduped, median_depth_deduped, garb = get_depth_info(summary_depth_info)
    depth_positional_cv_deduped = get_depth_positional_cv(summary_depth_info, COVERAGE_TRIM_TAIL)

    # low coverage tail for customers, based on deduped coverage profile
    summary_metrics['low_cov_' + str(CUSTOMER_LEFT_TAIL_COVERAGE)] = get_depth_tail_fract(summary_depth_info, CUSTOMER_LEFT_TAIL_COVERAGE, left_tail=True)

    if bias_results['target_info'] != {}:
        target_info = bias_results['target_info']
        summary_metrics['fraction_on_target'] = tk_stats.robust_divide(float(target_info['on_target_bases']), target_info['total_bases'])
    else:
        summary_metrics['fraction_on_target'] = None

    summary_metrics['detected_sex'] = bias_results.get('detected_sex')
    summary_metrics['mean_depth'] = mean_depth
    summary_metrics['male_chromosome_copies'] = bias_results.get('male_chromosome_copies')
    summary_metrics['median_depth'] = median_depth
    summary_metrics['mean_depth_deduped'] = mean_depth_deduped
    summary_metrics['median_depth_deduped'] = median_depth_deduped
    summary_metrics['on_target_bases'] = on_target_bases
    summary_metrics['depth_positional_cv'] = depth_positional_cv
    summary_metrics['depth_positional_cv_deduped'] = depth_positional_cv_deduped
    summary_metrics['zero_cov_fract'] = zero_cov_fract

    # Compute fraction of reads in high-coverage spikes
    cov_data = bias_results['summary_depth_info_deduped']
    _, conf_median, _ = get_depth_info(cov_data)
    conf_median = max(conf_median, 1)
    cov_variance = conf_median + (conf_median * depth_positional_cv_deduped)**2
    cov_sigma = math.sqrt(cov_variance)
    high_cutoff = conf_median + 5.0 * cov_sigma
    cov_data = {int(k):v for (k,v) in cov_data.iteritems()}
    total = sum(float(k*v) for (k,v) in cov_data.iteritems())
    outlier = sum( float(k*v) for (k,v) in cov_data.iteritems() if k > high_cutoff)
    summary_metrics['high_coverage_pileup_fraction'] = tk_stats.robust_divide(outlier, total)

    # Add metrics from variant_results
    if not(args.variant_results is None):
        with open(args.variant_results) as variant_results_file:
            variant_results = json.load(variant_results_file)
        summary_metrics.update(variant_results)

    # Copy of coalescence results
    coa_metrics = json.load(open(args.filter_barcodes_results))
    for (k,v) in coa_metrics.items():
        summary_metrics[k] = v

    if not(args.sv_results is None):
        with open(args.sv_results) as sv_results_file:
            sv_results = json.load(sv_results_file)
        summary_metrics.update(sv_results)

    if not(args.short_del_results is None):
        with open(args.short_del_results) as short_del_results_file:
            short_del_results = json.load(short_del_results_file)
        new_res = {}
        for k, v in short_del_results.iteritems():
            new_res['short_del_' + k] = v
        summary_metrics.update(new_res)

    # Length mass results
    # Only copy scalar results
    if args.length_mass_results is not None:
        with open(args.length_mass_results) as length_mass_file:
            lm_results = json.load(length_mass_file)
            for (k,v) in lm_results.iteritems():
                if type(v) == str or type(v) == int or type(v) == float or v is None:
                    summary_metrics[k] = v

    # Reference genome information
    summary_metrics['reference_name'] = reference_name = tenkit.reference.get_genome(args.reference_path)
    ref_fasta = tenkit.reference.open_reference(args.reference_path)
    summary_metrics['reference_contigs'] = reference_contigs = len(ref_fasta)
    summary_metrics['reference_bases'] = reference_bases = sum(len(ref_fasta[contig]) for contig in ref_fasta)
    martian.log_info("Reference: %s, %d contigs, %d bases" % (reference_name, reference_contigs, reference_bases))

    # Check for SV blacklist (only check if SV calling is enabled)
    summary_metrics['sv_blacklist_present'] = True
    if not(args.sv_results is None):
        if not tenkit.reference.is_tenx(args.reference_path):
            if not os.path.exists(tenkit.reference.get_sv_blacklist(args.reference_path)):
                summary_metrics['sv_blacklist_present'] = False
                martian.alarm("WARNING: Pipeline run without a region blacklist for SV calling. SV calls may contain many false positives due to problematic regions in the reference.")

    # Gelbead lot information
    if not(args.lot_info is None):
        with open(args.lot_info) as lot_info_file:
            lot_info_results = json.load(lot_info_file)
        summary_metrics.update(lot_info_results)

    # Downsampling information
    if not(args.downsample_info is None):
        with open(args.downsample_info) as downsample_info_file:
            downsample_info_results = json.load(downsample_info_file)
        summary_metrics.update(downsample_info_results)

    # Summary metrics are now finalized -- evaluate alarms
    # Select alarm file -- right now we always use the same one
    alarm_rules = tenkit.alarms.load_rules(args.targets)
    alarms = tenkit.alarms.evaluate_alarms(alarm_rules, summary_metrics)

    # Write alarm file
    with open(outs.alarms, 'w') as alarms_output_file:
        alarms_output_file.write(tenkit.safe_json.safe_jsonify(alarms))

    # Log alarms to martian
    with open(outs.alarms_summary, 'w') as al_summary_file:
        def wl(s):
            al_summary_file.write(s + "\n")

        wl("10X Genomics - Pipeline Run Details")
        wl("-" * 40)
        wl("Sample ID: %s" % args.sample_id)
        wl("Genome: %s" % tenkit.reference.get_genome(args.reference_path))
        wl("Reference Path: %s" % args.reference_path)
        wl("Targets file: %s" % args.targets)
        if alarms is not None and len(alarms) > 0:
            wl("")
            wl("Sequencing Metric Alarms:")
            wl("-" * 40)
            for alarm in alarms:
                wl("%s [%s] -- %s" % (alarm['level'], alarm['title'], alarm['message']))
        else:
            wl("")
            wl("No alarms raised.")

    summary_output_file = open(outs.summary, 'w')
    summary_output_file.write(tenkit.safe_json.safe_jsonify(summary_metrics, pretty=True))
    summary_output_file.close()

    # Generate CS summary metrics CSV
    sv_calls_metric = "num_calls"
    metrics_key_map = copy.deepcopy(CS_METRICS_KEY_MAP)
    metrics_key_map.append([sv_calls_metric, "large_sv_calls"])


    if args.targets is None:
        metrics_key_map.append(["short_del_calledDEL_num_calls", "short_deletion_calls"])
    else:
        metrics_key_map.append(["short_del_total_del_numPosCalls", "short_deletion_calls"])

    generate_summary_cs_csv(metrics_key_map, summary_metrics, outs.summary_cs)

def generate_summary_cs_csv(metrics_key_map, summary_metrics, csv_fname):
    # Extract and rename desired fields from internal summary
    summary_metrics_cs = {}
    for pair in metrics_key_map:
        if summary_metrics.has_key(pair[0]):
            summary_metrics_cs[pair[1]] = summary_metrics[pair[0]]

    if summary_metrics_cs.has_key("mapped_reads"):
        # Convert unmapped_fract to mapped_reads
        summary_metrics_cs["mapped_reads"] = - (summary_metrics_cs["mapped_reads"] - 1)

    with open(csv_fname, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=[ pair[1] for pair in metrics_key_map ])
        writer.writeheader()
        writer.writerow(summary_metrics_cs)

def get_on_target_bases(depths):
    on_target_bases = 0
    for depth, count in depths.iteritems():
        on_target_bases += int(depth)*int(count)
    return on_target_bases

def get_depth_info(info):
    fixed_info = {int(x): y for (x, y) in info.iteritems()}

    total_depth_counts = sum(fixed_info.values())
    median_depth = None
    sorted_depths = sorted(fixed_info.keys())
    seen_depth_count = 0
    mean_depth = 0.0
    for depth in sorted_depths:
        seen_depth_count += fixed_info[depth]
        mean_depth += float(depth*fixed_info[depth])/float(total_depth_counts)
        if seen_depth_count > total_depth_counts/2 and median_depth is None:
            median_depth = depth
    zero_cov_fract = float(fixed_info.get(0, 0.0))/float(total_depth_counts)

    return (mean_depth, median_depth, zero_cov_fract)

def get_depth_tail_fract(info, cutoff, left_tail=True):
    if info == {} or info is None:
        return float('NaN')

    fixed_info = {int(x): y for (x, y) in info.iteritems()}
    total_count = float(sum(fixed_info.values()))
    restricted_count = 0.0
    for depth, count in fixed_info.iteritems():
        if (left_tail and depth <= cutoff) or (not(left_tail) and depth >= cutoff):
            restricted_count += count
    return restricted_count/total_count

def get_depth_positional_cv(info, trim_tail):
    fixed_info = {int(x): y for (x, y) in info.iteritems()}
    total_count = sum(fixed_info.values())
    cutoff_count = total_count*trim_tail
    seen_count = 0
    for depth in sorted(fixed_info.iterkeys(), reverse=True):
        seen_count += fixed_info[depth]
        if seen_count >= cutoff_count:
            cutoff = depth
            break
    trimmed_info = {x: y for (x, y) in fixed_info.iteritems() if x <= cutoff}
    mean_val, var_val = tk_stats.mean_var_from_counts(trimmed_info)
    if mean_val > var_val:
        return float('NaN')
    return tk_stats.robust_divide(numpy.sqrt(var_val - mean_val), mean_val)

def get_instrument_ids(bam_fn):
    ''' Scan some reads and make list of observed instrument ids '''

    bam = tenkit.bam.create_bam_infile(bam_fn)

    instrument_ids = set()

    # Iterate over some mapped reads and record the observed instrument
    # IDs, which should be the first field in the qname
    for read in itertools.islice(bam.fetch(), 5000):
        qname_parts = read.qname.split(":")

        if len(qname_parts) > 0:
            instrument_ids.add(qname_parts[0])

    # If something went wrong and we're not seeing a small set
    # of ids, just return a sampling
    return list(instrument_ids)[:10]
