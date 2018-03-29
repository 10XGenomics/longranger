#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import os.path
import cPickle
import re
import json
import numpy as np
from collections import defaultdict
from itertools import product
from longranger.sv.constants import SV_DEFAULT_MAX_BC_COV
import tenkit.safe_json
from sklearn.metrics import average_precision_score
import tenkit.pandas as pd
import longranger.sv.utils as tk_sv_utils
import longranger.sv.io as tk_sv_io
import tenkit.stats as tk_stats
import tenkit.regions as tk_regions
import tenkit.reference as tk_reference
from tenkit.coverage import get_depth_info_json
import longranger.genomic_tracks as lr_gt

import martian

__MRO__ = """
stage ANALYZE_SV_CALLS(
    in  string     sample_id,
    in  bedpe      variants,
    in  bedpe      gt_variants,
    in  json       call_summary,
    in  json       coverage,
    in  bool       keep_filters,
    in  int        min_call_qv_wgs,
    in  int        min_call_qv_target,
    in  int        min_read_support,
    in  string     reference_path,
    in  bed        sv_blacklist_regions,
    in  bedpe      seg_dups,
    in  int        min_dist_from_black,
    in  float      max_frac_black,
    in  int        seg_dup_min_dist,
    in  int[]      detect_dists,
    in  float      min_rel_overlap,
    in  bed        targets,
    in  int[]      target_dists,
    in  int        min_sv_len,
    in  float      min_allelic_frac,
    in  bool       is_germline,
    in  float      max_bc_cov_factor,
    in  string     blacklist_mode,
    in  string     segdup_mode,
    out json       summary,
    out tsv        summary_tsv,
    out bedpe      sv_calls,
    out bedpe      sv_candidates,
    out vcf.gz     svs,
    out vcf.gz.tbi svs_index,
    out bedpe      feasible_gt,
    out tsv        call_tsv,
    src py         "stages/structvars/analyze_sv_calls",
) split using (
    in  int    start_idx,
    in  int    stop_idx,
)
"""

MAX_SENS_TIER = 1
MAX_PPV_TIER = 2

def split(args):
    sv_df = tk_sv_io.read_sv_bedpe_to_df(args.variants)
    gt_df = tk_sv_io.read_sv_bedpe_to_df(args.gt_variants)
    tk_sv_io.check_sv_names(gt_df)

    sv_df["name"] = ["call_%d" % idx for idx in range(len(sv_df))]

    variants_bedpe = os.path.join(os.getcwd(), "variants.bedpe")
    tk_sv_io.write_sv_df_to_bedpe(sv_df, variants_bedpe)

    nsvs = sv_df.shape[0]
    nbreaks_per_chunk = max(100, int(np.ceil(nsvs / 32.0))) # avoid overchunking
    nchunks = int(np.ceil(nsvs / float(nbreaks_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nbreaks_per_chunk
        chunk_end = min(nsvs, (i + 1) * nbreaks_per_chunk)
        chunk_defs.append({'renamed_variants': variants_bedpe, 'start_idx':chunk_start, 'stop_idx':chunk_end, '__mem_gb':12})

    if len(chunk_defs) == 0:
        chunk_defs = [{'renamed_variants': variants_bedpe, 'start_idx':0, 'stop_idx':0, '__mem_gb':12}]


    return {'chunks': chunk_defs, 'join': {'__mem_gb': 16}}


def main(args, outs):
    pred_df = tk_sv_io.read_sv_bedpe_to_df(args.renamed_variants)
    pred_df = tk_sv_utils.get_dataframe_loc(pred_df, list(range(args.start_idx, args.stop_idx)))

    if not args.gt_variants is None:
        true_df = prepare_gt(args)
        true_df.to_csv(outs.feasible_gt, index=False, header=True, sep='\t', na_rep='NaN')
    else:
        true_df = None

    #### Get matches between this chunk of the calls and the ground truth
    max_detect_dist = np.max(np.array(args.detect_dists))
    res = get_matches(pred_df, true_df, max_detect_dist, args.min_rel_overlap)
    pred_to_match, true_to_match, _ = res

    #### Apply filters
    if len(pred_df) > 0:
        # Loading all these files can take awhile. Don't do it if there are no SVs to analyze.

        # blacklist and segdups files can come from 3 places, in this order of precedence:
        # 1. mro argument sv_blacklist_regions
        # 2. <reference_path>/regions/sv_blacklist.bed (or segdups.bed)
        # 3. <tenkit install>/sv_data/<genome>/default_sv_blacklist.bed (accessed by tenkit.constants.find_sv_blacklist)

        if os.path.exists(tk_reference.get_sv_blacklist(args.reference_path)):
            blacklist_file = tk_reference.get_sv_blacklist(args.reference_path)
        else:
            blacklist_file = lr_gt.get_genomic_track(args.sv_blacklist_regions, args.blacklist_mode, args.reference_path, "default_blacklist.bed")

        # This will merge overlapping blacklist regions
        black_regions = tk_sv_utils.bed_to_region_map(blacklist_file, merge=True)
        # Match each region in black_regions to a set of entries from the bed
        # file that overlap it. This is done so we can output the names of
        # entries that were used to blacklist each sv.
        black_region_names = get_region_names(blacklist_file, black_regions)
        # compute the distance between the breakpoints and the blacklist
        # elements. Get the distance together with the names of the closest
        # blacklist elements.
        res = get_df_region_dist(pred_df, black_regions, black_region_names)
        black_dists1, black_dists2, _, _, black_names1, black_names2 = res

        if os.path.exists(tk_reference.get_segdups(args.reference_path)):
            seg_dups_file = tk_reference.get_segdups(args.reference_path)
        else:
            seg_dups_file =  lr_gt.get_genomic_track(args.seg_dups, args.segdup_mode, args.reference_path, "default_segdups.bedpe")

        # from call to matching seg dups
        seg_dup_calls, _, _ = tk_sv_utils.compare_breaks(pred_df, seg_dups_file, max_dist=args.seg_dup_min_dist)
        seg_dup_regions = tk_sv_utils.bedpe_to_region_map(seg_dups_file, merge=True)
        all_bad_regions = tk_sv_utils.merge_region_maps(black_regions, seg_dup_regions)
    else:
        black_dists1 = None
        black_dists2 = None
        black_names1 = None
        black_names2 = None
        seg_dup_calls = {}
        all_bad_regions = None

    pred_df, min_qv = add_filters(pred_df, pred_to_match,
                                  black_dists1, black_dists2, black_names1, black_names2,
                                  seg_dup_calls, all_bad_regions, args)

    with open(re.sub('.json', '.pickle', outs.summary), 'wb') as f:
        cPickle.dump(pred_to_match, f)
        cPickle.dump(true_to_match, f)
        cPickle.dump((pred_df, min_qv), f)


def merge_predictions(chunk_outs):
    join_pred_to_match = {}
    join_true_to_match = defaultdict(set)
    join_pred_df = None
    feasible_gt = None
    for ci, chunk in enumerate(chunk_outs):
        with open(re.sub('.json', '.pickle', chunk.summary)) as f:
            pred_to_match = cPickle.load(f)
            true_to_match = cPickle.load(f)
            pred_df, min_qv = cPickle.load(f)
        join_pred_df = pd.concat([join_pred_df, pred_df], ignore_index=True)
        for pred, matches in pred_to_match.iteritems():
            join_pred_to_match[pred] = matches
        for sv, matches in true_to_match.iteritems():
            join_true_to_match[sv] = join_true_to_match[sv].union(matches)
        if ci == 0 and os.path.exists(chunk.feasible_gt):
            feasible_gt = pd.read_csv(chunk.feasible_gt, header=0, index_col=None, sep='\t')
    return (join_pred_to_match, join_true_to_match, join_pred_df, feasible_gt, min_qv)


def join(args, outs, chunk_defs, chunk_outs):
    pred_to_match, _, pred_df, true_df, min_qv = merge_predictions(chunk_outs)

    # Change TRANS type to DISTAL. This change will only
    # affect the type reported not the names of the metrics.
    new_info = []
    for _, row in pred_df.iterrows():
        sv_type = tk_sv_io.get_sv_type(row.info)
        if sv_type == 'TRANS':
            sv_type = 'DISTAL'
        new_info.append(tk_sv_io.update_info(row.info, ['TYPE'], [sv_type]))
    pred_df['info'] = new_info

    if not true_df is None:
        true_df.to_csv(outs.feasible_gt, index=False, header=True, sep='\t', na_rep='NaN')

    ##### Write BEDPE/VCF outputs
    tk_sv_io.write_sv_df_to_bedpe(pred_df, outs.sv_candidates)
    source_str = '10X/pipelines/stages/analyze_sv_calls {}'.format(martian.get_pipelines_version())
    sample_id = 'sample' if args.sample_id is None else args.sample_id
    tk_sv_io.bedpe_to_vcf(outs.sv_candidates, outs.svs.strip('.gz'),
                          sample_id, source_str, args.reference_path)
    # this will sort and gzip
    tk_sv_io.index_sv_vcf(outs.svs.strip(".gz"))
    outs.svs_index = outs.svs + '.tbi'
    # delete the non-gzipped file
    os.remove(outs.svs.strip('.gz'))

    if not pred_df.empty:
        call_df = pred_df[np.logical_or(pred_df['filters'] == '.',
                                        pred_df['filters'] == "PASS")]
    else:
        call_df = None
    tk_sv_io.write_sv_df_to_bedpe(call_df, outs.sv_calls)

    # Annotate each call with the matching ground truth svs. The resulting
    # dataframe might have multiple rows for the same call if there were multiple
    # matching ground truth svs.
    martian.log_info("merging calls and gt")
    if not pred_df.empty:
        pred_df = merge_calls_and_gt(pred_df, true_df, pred_to_match)

    martian.log_info("writing call_tsv")
    pred_df.to_csv(outs.call_tsv, index=False, header=True, sep='\t', na_rep='NaN')

    pred_df = pred_df[np.logical_not(pd.isnull(pred_df['name']))]

    max_dists = sorted(np.array(args.detect_dists))

    gt_sv_types = get_all_sv_types(true_df)
    call_sv_types = get_all_sv_types(pred_df)

    if not true_df is None:
        # Use the default MAX_PPV_TIER unless this is greater than the maximum tier
        # present in the data.
        max_ppv_tier = min(MAX_PPV_TIER, np.max(true_df.tier))
        # Use the default unless this is smaller than the minimum tier present in
        # the data.
        max_sens_tier = max(MAX_SENS_TIER, np.min(true_df.tier))
    else:
        max_ppv_tier = 1
        max_sens_tier = 1

    tiers = [max_ppv_tier, max_sens_tier]

    # All combinations of filters in ground truth and call set
    if not args.targets is None and not args.target_dists is None:
        target_dists = list(sorted(np.array(args.target_dists, dtype = np.float)))
        target_dists.append(float('NaN'))
    else:
        target_dists = [float('NaN')]

    combs = product([0, 1, 2, None], target_dists, gt_sv_types, tiers,
                    [True, False], call_sv_types, max_dists)

    metrics = defaultdict(list)

    gt_filters = ['genic_breaks', 'target_dist', 'gt_sv_type', 'tier']
    call_filters = ['call_filtered', 'call_sv_type', 'match_dist']

    for (genic_breaks, tdist, gt_sv_type, tier, is_filtered, call_sv_type, dist) in combs:
        if gt_sv_type != 'NA' and call_sv_type != 'NA' and gt_sv_type != call_sv_type:
            continue

        metrics['genic_breaks'].append(genic_breaks)
        metrics['target_dist'].append(tdist)
        metrics['gt_sv_type'].append(gt_sv_type)
        metrics['tier'].append(tier)
        metrics['call_filtered'].append(is_filtered)
        metrics['call_sv_type'].append(call_sv_type)
        metrics['match_dist'].append(dist)

        if true_df is None:
            sel_true_df = None
        else:
            sel_true_df = true_df
            if gt_sv_type != 'NA':
                sel_true_df = sel_true_df[sel_true_df.sv_type == gt_sv_type]
            if not np.isnan(tdist):
                sel_true_df = sel_true_df[sel_true_df.targ_dist <= tdist]
            sel_true_df = sel_true_df[sel_true_df.tier <= tier]
            # Restrict to genic or non-genic or take everything if this is None.
            if not genic_breaks is None:
                sel_true_df = sel_true_df[sel_true_df.genic_breaks == genic_breaks]

            if len(sel_true_df) == 0:
                sel_true_df = None

        sel_pred_df = pred_df

        if is_filtered and not pred_df.empty:
            sel_pred_df = sel_pred_df[(sel_pred_df.filters == '.') | (sel_pred_df.filters == 'PASS')]
        if call_sv_type != 'NA' and not pred_df.empty:
            sel_pred_df = sel_pred_df[sel_pred_df.sv_type == call_sv_type]
        if not pred_df.empty and (args.min_rel_overlap is None or args.min_rel_overlap == 0):
            # Do not apply thi filter if the matching is done based on overlap.
            sel_pred_df = sel_pred_df[np.logical_or(np.isnan(sel_pred_df.match_dist),
                                                    sel_pred_df.match_dist <= dist)]

        add_metrics(sel_pred_df, sel_true_df, metrics)

    column_names = gt_filters
    column_names.extend(call_filters)
    other_names = set(metrics.keys()).difference(set(column_names))
    column_names.extend(other_names)

    metric_df = pd.DataFrame(metrics)
    metric_df = metric_df[column_names]

    martian.log_info("writing summary tsv")
    metric_df.to_csv(outs.summary_tsv, index=False, header=True, sep='\t', na_rep='NaN')

    short_metrics = get_short_metrics(metric_df, other_names, max_ppv_tier, max_sens_tier, args)

    if not args.call_summary is None:
        with open(args.call_summary, 'r') as in_summary_fn:
            in_summary = json.load(in_summary_fn)
            for key, val in in_summary.iteritems():
                short_metrics[key] = val

    short_metrics['min_qv'] = min_qv

    with open(outs.summary, 'w') as out_file:
        out_file.write(tenkit.safe_json.safe_jsonify(short_metrics, pretty=True))


def add_metrics(pred_df, true_df, metrics):
    if true_df is None:
        true_names = set([])
        true_feasible_names = set([])
    else:
        true_names = set(true_df['name'])
        true_feasible_names = set(true_df[true_df.feasible]['name'])

    if pred_df.empty:
        num_calls = 0
        metrics['num_calls'].append(0)
        metrics['frac_read_support'].append(float('NaN'))
        metrics['mean_break_res'].append(float('NaN'))
        metrics['median_break_res'].append(float('NaN'))
        uniq_pred_df = pred_df
    else:

        pred_df['break_res'] = np.maximum(pred_df.stop1 - pred_df.start1,
                                          pred_df.stop2 - pred_df.start2)

        # If there are multiple rows for the same call, pick the one with the minimum
        # distance to match. If there are multiple gt events with the same distance, pick
        # the one that is in the given subset of the ground truth. If there are still
        # multiple, pick the one that has the same type. If there are still more
        # then pick the one with the same orientation.
        def sel_fun(x):
            if len(x) == 1:
                # single match or no match (NaN match)
                return x.index[0]
            else:
                sel = np.array(x.match_dist == np.min(x.match_dist), dtype=np.bool)
                idx = np.array(x[sel].index)
                if len(idx) == 0:
                    return x.index[0]
                elif len(idx) == 1:
                    # single call with minimum matching distance
                    return idx[0]
                else:
                    sel2 = np.logical_and(sel, np.array([m in true_names for m in x.match], dtype=np.bool))
                    idx2 = np.array(x[sel2].index)
                    if len(idx2) == 0:
                        return idx[0]
                    elif len(idx2) == 1:
                        # single call with minimum distance and matching type
                        return idx2[0]
                    else:
                        sel3 = np.logical_and(sel2, x.sv_type == x.sv_type_gt)
                        idx3 = np.array(x[sel3].index)
                        if len(idx3) == 0:
                            return idx2[0]
                        elif len(idx3) == 1:
                            return idx3[0]
                        else:
                            sel4 = np.logical_and(sel3, x.orient == x.orient_gt)
                            idx4 = np.array(x[sel4].index)
                            if len(idx4) == 0:
                                return idx3[0]
                            return idx4[0]


        # sel_fun returns the indices instead of a dataframe. this allows us to
        # handle an empty list of indices bettes
        uniq_idx = np.array(pred_df.groupby('name').apply(sel_fun))
        uniq_pred_df = pred_df.loc[uniq_idx]

        num_calls = len(uniq_pred_df)
        metrics['num_calls'].append(num_calls)
        metrics['frac_read_support'].append(np.mean(uniq_pred_df.read_support > 0))
        metrics['mean_break_res'].append(np.mean(uniq_pred_df.break_res))
        metrics['median_break_res'].append(np.median(uniq_pred_df.break_res))

    if true_df is None:
        metrics['num_gt'].append(0)
        metrics['sensitivity'].append(float('NaN'))
        metrics['mean_dist_to_match'].append(float('NaN'))
        metrics['median_dist_to_match'].append(float('NaN'))
        metrics['25prc_dist_to_match'].append(float('NaN'))
        metrics['75prc_dist_to_match'].append(float('NaN'))
        metrics['PPV'].append(0)
        metrics['AUC'].append(0)
        metrics['num_type_gt'].append(0)
        metrics['type_accuracy'].append(float('NaN'))
        metrics['num_orient_gt'].append(0)
        metrics['orient_accuracy'].append(float('NaN'))
        metrics['num_orient_gt_missing'].append(float('NaN'))
        metrics['orient_accuracy_missing'].append(float('NaN'))
    elif pred_df.empty or pred_df is None:
        metrics['num_gt'].append(len(true_feasible_names))
        metrics['sensitivity'].append(0)
        metrics['mean_dist_to_match'].append(float('NaN'))
        metrics['median_dist_to_match'].append(float('NaN'))
        metrics['25prc_dist_to_match'].append(float('NaN'))
        metrics['75prc_dist_to_match'].append(float('NaN'))
        metrics['PPV'].append(0)
        metrics['AUC'].append(0)
        metrics['num_type_gt'].append(0)
        metrics['type_accuracy'].append(0)
        metrics['num_orient_gt'].append(0)
        metrics['orient_accuracy'].append(0)
        metrics['num_orient_gt_missing'].append(0)
        metrics['orient_accuracy_missing'].append(0)
    else:
        # Need to do this because pred_df might have matches that do not correspond
        # to this subset of the true set.
        is_good_match = np.array([not pd.isnull(m) and m in true_feasible_names for m in pred_df.match], dtype = np.bool)
        num_called = len(set(pred_df[is_good_match].match))
        metrics['num_gt'].append(len(true_feasible_names))
        metrics['sensitivity'].append(tk_stats.robust_divide(num_called, len(true_feasible_names)))

        new_dists = np.array(uniq_pred_df.match_dist, dtype = np.float)
        new_dists[np.array([not pd.isnull(m) and not m in true_names for m in uniq_pred_df.match], dtype=np.bool)] = float('NaN')
        uniq_pred_df['match_dist'] = new_dists
        metrics['mean_dist_to_match'].append(np.nanmean(uniq_pred_df.match_dist))
        metrics['median_dist_to_match'].append(np.nanmedian(uniq_pred_df.match_dist))
        metrics['25prc_dist_to_match'].append(np.nanpercentile(uniq_pred_df.match_dist, 25))
        metrics['75prc_dist_to_match'].append(np.nanpercentile(uniq_pred_df.match_dist, 75))
        pred_arr = np.logical_not(np.isnan(uniq_pred_df.match_dist))
        num_correct = np.sum(pred_arr)
        metrics['PPV'].append(tk_stats.robust_divide(num_correct, num_calls))

        if len(pred_arr) == 0:
            auc = float('NaN')
        elif np.all(pred_arr == 0):
            auc = 0
        elif np.all(pred_arr == 1):
            auc = 1
        else:
            auc = average_precision_score(pred_arr, uniq_pred_df.qual)
        metrics['AUC'].append(auc)

        # orientation
        orient_df = uniq_pred_df[np.logical_and(np.logical_not(np.isnan(uniq_pred_df.match_dist)),
                                                np.logical_and((uniq_pred_df.sv_type == 'TRANS') | (uniq_pred_df.sv_type == 'UNK'),
                                                               np.logical_and(uniq_pred_df.orient_gt != '..',
                                                                              uniq_pred_df.orient != '..')))]
        metrics['num_orient_gt'].append(len(orient_df))
        correct_orient = np.sum(orient_df.orient == orient_df.orient_gt)
        metrics['orient_accuracy'].append(tk_stats.robust_divide(correct_orient, len(orient_df)))

        # same as before but count cases where we did not make a prediction as false predictions
        orient_df = uniq_pred_df[np.logical_and(np.logical_not(np.isnan(uniq_pred_df.match_dist)),
                                                np.logical_and(uniq_pred_df.sv_type == 'TRANS',
                                                               uniq_pred_df.orient_gt != '..'))]
        metrics['num_orient_gt_missing'].append(len(orient_df))
        correct_orient = np.sum(orient_df.orient == orient_df.orient_gt)
        metrics['orient_accuracy_missing'].append(tk_stats.robust_divide(correct_orient, len(orient_df)))

        type_df = uniq_pred_df[np.logical_and(np.logical_not(np.isnan(uniq_pred_df.match_dist)),
                                              uniq_pred_df.sv_type_gt != 'UNK')]
        metrics['num_type_gt'].append(len(type_df))
        correct_type = np.sum(type_df.sv_type_gt == type_df.sv_type)
        metrics['type_accuracy'].append(tk_stats.robust_divide(correct_type, len(type_df)))


def get_short_metrics(metric_df, metric_names, max_ppv_tier, max_sens_tier, args):
    metric_df = metric_df[metric_df.call_filtered]
    max_detect_dist = np.max(np.array(args.detect_dists))
    metric_df = metric_df[metric_df.match_dist == max_detect_dist]
    metric_df = metric_df[np.logical_or(metric_df.call_sv_type == 'NA',
                                        metric_df.gt_sv_type == 'NA')]

    if args.targets is None:
        metric_df = metric_df[np.isnan(metric_df.target_dist)]

    # metrics independent of gt
    call_metrics = ['num_calls', 'frac_read_support', 'mean_break_res', 'median_break_res']

    # metrics independent of callset
    gt_metrics = ['num_gt']

    # metrics to report by sv-type of call
    type_dependent_metrics = ['PPV', 'type_accuracy', 'orient_accuracy']
    sens_metrics = ['sensitivity']

    metrics = {}
    for _, row in metric_df.iterrows():
        if row.tier != max_ppv_tier and row.tier != max_sens_tier:
            continue
        if row.call_sv_type == 'NA':
            # metrics for the entire call set (irrespective of type)
            for name in metric_names:
                if name in call_metrics:
                    key = '{}'.format(name)
                else:
                    # Only report sensitivity metrics and ground truth metrics for the
                    # sensitivity tier. Note that ground truth metrics are the only
                    # metrics that get reported for both tiers.
                    if (row.tier != max_sens_tier and name in sens_metrics) or \
                       (row.tier == max_sens_tier and max_sens_tier != max_ppv_tier and not name in sens_metrics and not name in gt_metrics):
                        continue
                    if pd.isnull(row.genic_breaks):
                        if args.targets is None or np.isnan(row.target_dist):
                            key = '{}_{}'.format(row.gt_sv_type, name)
                        else:
                            key = '{}_{}_feasible{:d}'.format(row.gt_sv_type, name, int(row.target_dist))
                    else:
                        prefix = {0:'non', 1:'semi', 2:''}[row.genic_breaks]
                        key = '{}_{}_{}'.format(row.gt_sv_type, name, prefix + 'genic')

                    # Only for ground truth related metrics, we add the "tier" info. For the
                    # rest of the metrics it's implied by the metric.
                    if name in gt_metrics:
                        # If the sensitivity and PPV tiers are the same, we need to explicitly add
                        # both keys (with equal values).
                        if max_sens_tier == max_ppv_tier:
                            metrics[key] = row[name]
                        key = '{}{}'.format('conf_' if row.tier == max_sens_tier else '', key)

                metrics[key] = row[name]
        else:
            for name in metric_names:
                if name in gt_metrics:
                    continue
                elif name in call_metrics:
                    key = 'called{}_{}'.format(row.call_sv_type, name)
                    metrics[key] = row[name]
                else:
                    if (row.tier != max_sens_tier and name in sens_metrics) or \
                       (row.tier == max_sens_tier and max_sens_tier != max_ppv_tier and not name in sens_metrics):
                        continue
                    if name in type_dependent_metrics and (args.targets is None or np.isnan(row.target_dist)):
                        key = 'called{}_{}'.format(row.call_sv_type, name)
                        metrics[key] = row[name]

    return metrics


def get_df_region_dist(sv_df, regions, region_names = None, use_orient = False):
    """Computes the distance between the breakpoints of the BEDPE
    (read as dataframe) and a set of regions.
    region_names is a dict (start, stop) -> name."""
    dists1 = np.inf * np.ones((len(sv_df), ))
    dists2 = np.inf * np.ones((len(sv_df), ))
    regions1 = [(None, None) for i in range(len(sv_df))]
    regions2 = [(None, None) for i in range(len(sv_df))]
    matched_names1 = ['.' for i in range(len(sv_df))]
    matched_names2 = ['.' for i in range(len(sv_df))]

    for i, (_, row) in enumerate(sv_df.iterrows()):
        if use_orient:
            sv_type = tk_sv_io.get_sv_type(row.info)
            if sv_type == 'DEL':
                orient = '-+'
            elif sv_type == 'DUP':
                orient = '+-'
            elif sv_type == 'INV':
                orient = '..'
            else:
                orient = tk_sv_io.get_break_orientation(row.info)
            if (orient == '..' and sv_type != 'INV') or sv_type == 'UNK' or sv_type == 'INS':
                continue
            orient1 = tk_regions.Dirs.from_str(orient[0])
            orient2 = tk_regions.Dirs.from_str(orient[1])
        else:
            orient1, orient2 = (None, None)
        chrom1, chrom2 = row.chrom1, row.chrom2
        if chrom1 in regions:
            s1, e1, d1 = regions[chrom1].get_closest_region_to_region(row.start1, row.stop1,
                                                                      direction=orient1)
            if not s1 is None:
                d1 = int(d1)
                regions1[i] = (s1, e1)
                if not region_names is None and (s1, e1) in region_names:
                    matched_names1[i] = ','.join(list(region_names[(s1, e1)]))
        else:
            d1 = np.inf
        if chrom2 in regions:
            s2, e2, d2 = regions[chrom2].get_closest_region_to_region(row.start2, row.stop2,
                                                                      direction=orient2)
            if not s2 is None:
                d2 = int(d2)
                regions2[i] = (s2, e2)
                if not region_names is None and (s2, e2) in region_names:
                    matched_names2[i] = ','.join(list(region_names[(s2, e2)]))
        else:
            d2 = np.inf
        dists1[i] = d1
        dists2[i] = d2
    return (dists1, dists2, regions1, regions2, matched_names1, matched_names2)


def get_region_names(infile, regions):
    """Match each region to all entries from infile that overlap it."""
    name_map = defaultdict(set)
    if infile is not None:
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track') or line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) > 3:
                    name = fields[3]
                    ov_regions = regions[fields[0]].overlapping_regions(int(fields[1]), int(fields[2]))
                    for r in ov_regions:
                        name_map[r].add(name)
    return name_map


def dist_to_breaks(pos, s, e):
    if pos >= s and pos <= e:
        return 0
    return pos - e if pos > e else s - pos


def get_matches(pred_df, gt_variants, max_dist, min_rel_overlap):
    cp_pred_df = pred_df.copy()

    # Replace breakpoints with their middle points
    cp_pred_df.start1 = np.array((pred_df.start1 + pred_df.stop1) / 2, dtype=np.int)
    cp_pred_df.stop1 = pred_df.start1 + 1
    cp_pred_df.start2 = np.array((pred_df.start2 + pred_df.stop2) / 2, dtype=np.int)
    cp_pred_df.stop2 = pred_df.start2 + 1

    if not min_rel_overlap is None and min_rel_overlap > 0:
        res = tk_sv_utils.overlap_breaks(cp_pred_df, gt_variants, min_rel_overlap)
    else:
        res = tk_sv_utils.compare_breaks(cp_pred_df, gt_variants, max_dist)

    return res


def add_filters(pred_df, pred_to_match,
                black_dists1, black_dists2, black_names1, black_names2,
                seg_dup_calls, all_bad_regions, args):

    if not args.targets is None:
        min_call_qv = args.min_call_qv_target
    else:
        min_call_qv = args.min_call_qv_wgs

    if args.coverage is None:
        # used for WGS
        max_bc_cov = SV_DEFAULT_MAX_BC_COV
        bc_mean_depth = 200
    else:
        # used for exome
        with open(args.coverage, 'r') as f:
            cov_res = json.load(f)
        bc_summary_depth_info = cov_res['summary_bc_depth_info']
        bc_mean_depth, _,  _ = get_depth_info_json(bc_summary_depth_info)
        max_bc_cov = args.max_bc_cov_factor * bc_mean_depth

    if args.keep_filters:
        filter_strs = [s for s in pred_df.filters]
    else:
        filter_strs = ['.' for i in range(len(pred_df))]

    info_strs = [s for s in pred_df['info']]
    rps = np.zeros((len(pred_df), ), dtype = np.int)

    def get_cov_frac(black_regions, chrom, start, stop):
        regions = tk_sv_utils.strictly_overlapping_regions(black_regions, chrom, start, stop)
        tot_black = np.sum([r[1] - r[0] for r in regions])
        tot_len = float(stop - start)
        black_frac = tk_stats.robust_divide(tot_black, tot_len)
        return black_frac

    for i, (_, row) in enumerate(pred_df.iterrows()):
        npairs = tk_sv_io.get_npairs(row['info'])
        nsplit = tk_sv_io.get_nsplit(row['info'])
        rps[i] = npairs + nsplit
        sv_type = tk_sv_io.get_sv_type(row['info'])
        name = row['name']
        qual = row.qual

        ####### Filtering for read-pair calls #######
        frac_on_hap = tk_sv_io.extract_sv_info(row.info, ['FRAC_HAP_SUPPORT'])[0]
        allelic_frac = tk_sv_io.extract_sv_info(row.info, ['HAP_ALLELIC_FRAC'])[0]
        if allelic_frac != '':
            allelic_frac = float(allelic_frac)

        if args.is_germline is None:
            if qual < min_call_qv:
                filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)
            if not args.min_allelic_frac is None and not frac_on_hap is None and \
               frac_on_hap != '' and float(frac_on_hap) < args.min_allelic_frac:
                filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)
            if not args.min_allelic_frac is None and allelic_frac != '' and \
               float(allelic_frac) < args.min_allelic_frac:
                filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)
        elif args.targets is None:
            if args.is_germline:
                martian.log_info('Mean barcode depth {}'.format(bc_mean_depth))

                min_call_qv = max(min_call_qv, bc_mean_depth / 10.0)

                martian.log_info('Support cutoff: {} barcodes'.format(min_call_qv))

                enough_bcs = qual >= min_call_qv
                is_good = allelic_frac > 0.8 or (sv_type == 'INV' and allelic_frac > 0.6)
                is_good = is_good and enough_bcs
                if not is_good:
                    filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)
            else:
                min_call_qv = max(min_call_qv, 4)
                is_good = allelic_frac > 0.6 and qual >= min_call_qv
                if not is_good:
                    filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)
        else:
            if args.is_germline:
                # Harder to get confident support in Exome
                min_call_qv = max(min_call_qv, bc_mean_depth / 10.0)
                martian.log_info('Support cutoff: {} barcodes'.format(min_call_qv))
                # Apply a very lenient filter on allelic fraction because lots of barcodes can be unphased
                is_good = qual >= min_call_qv and allelic_frac > 0.05
                af = tk_sv_io.extract_sv_info(row.info, ['ALLELIC_FRAC'])[0]
                if af != '':
                    af = float(af)
                is_good = is_good and af > 0.04
            else:
                min_call_qv = max(min_call_qv, 4)
                is_good = qual >= min_call_qv

            if not is_good:
                filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'LOWQ', None)

        if not black_dists1 is None:
            chrom1, chrom2 = row.chrom1, row.chrom2
            black_dist1, black_dist2 = black_dists1[i], black_dists2[i]

            if chrom1 == chrom2:
                if chrom1 in all_bad_regions:
                    black_frac = get_cov_frac(all_bad_regions, chrom1, row.stop1, row.start2)
                else:
                    black_frac = 0.0
            else:
                black_frac = float('NaN')
        else:
            black_dist1 = np.inf
            black_dist2 = np.inf
            black_frac = float('NaN')
        filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'BLACK_DIST', min(black_dist1, black_dist2), args.min_dist_from_black)
        filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'BLACK_FRAC', black_frac, 0, args.max_frac_black)

        bname1 = '.'
        bname2 = '.'
        if black_dist1 < args.min_dist_from_black or re.search('BLACK_FRAC', filter_strs[i]):
            bname1 = black_names1[i]
        if black_dist2 < args.min_dist_from_black or re.search('BLACK_FRAC', filter_strs[i]):
            bname2 = black_names2[i]

        if name in seg_dup_calls:
            filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'SEG_DUP', None)
            seg_dup_match = ','.join(list(seg_dup_calls[name]))
        else:
            seg_dup_match = '.'

        nbcs1 = tk_sv_io.get_nbcs1(row.info)
        nbcs2 = tk_sv_io.get_nbcs2(row.info)
        if not nbcs1 is None and not nbcs2 is None and (nbcs1 > max_bc_cov or nbcs2 > max_bc_cov):
            filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'HIGH_BC_COV', None)

        filter_strs[i] = tk_sv_io.update_filters(filter_strs[i], 'READ_SUPPORT',
            npairs + nsplit, min_val = args.min_read_support)

        match_str = ','.join([str(s) for s in pred_to_match.get(name, '.')])

        if not args.targets is None:
            # Disable orientation reporting in exome
            info_strs[i] = tk_sv_io.update_info(info_strs[i], ['ORIENT'], [None])

        info_strs[i] = tk_sv_io.update_info(info_strs[i],
            ['BLACK_DIST1', 'BLACK_DIST2', 'BLACK_FRAC', 'BLACK1', 'BLACK2', 'MATCHES', 'SEG_DUP'],
            [black_dist1, black_dist2, black_frac, bname1, bname2, match_str, seg_dup_match])

    pred_df['filters'] = filter_strs
    pred_df['info'] = info_strs
    pred_df['read_support'] = rps

    return pred_df, min_call_qv


def prepare_gt(args):
    if args.gt_variants is None:
        return None

    true_df = tk_sv_io.read_sv_bedpe_to_df(args.gt_variants)

    # Length of ground truth sv
    true_df['dist'] = tk_sv_io.get_sv_df_dists(true_df)

    sv_types = []
    orients = []
    tiers = []

    # Mark genic SVs
    is_genic1 = np.zeros((len(true_df),), dtype=np.int)
    is_genic2 = np.zeros((len(true_df),), dtype=np.int)
    gene_regions = tk_reference.load_gene_boundaries(args.reference_path,
                                                     protein_coding=False)

    for row_idx, (_, row) in enumerate(true_df.iterrows()):
        if not 'info' in true_df.columns:
            sv_types.append('UNK')
            orients.append('..')
            tiers.append(0)
        else:
            sv_type = tk_sv_io.get_sv_type(row.info)
            if sv_type == 'DISTAL':
                sv_type = 'TRANS'
            sv_types.append(sv_type)
            orients.append(tk_sv_io.get_break_orientation(row.info))
            tiers.append(tk_sv_io.get_tier(row.info))

        is_genic1[row_idx] = int(row.chrom1 in gene_regions and
                                 bool(gene_regions[row.chrom1].overlapping_regions(row.start1, row.stop1)))
        is_genic2[row_idx] = int(row.chrom2 in gene_regions and
                                 bool(gene_regions[row.chrom2].overlapping_regions(row.start2, row.stop2)))

    true_df['break1_genic'] = is_genic1
    true_df['break2_genic'] = is_genic2
    # number of breakpoints overlapping genes
    true_df['genic_breaks'] = is_genic1 + is_genic2

    # put all the un-tiered entries into the last tier
    tiers = np.array(tiers)
    if len(tiers) == 0:
        total_tiers = 0
    else:
        total_tiers = np.max(tiers)

    tiers[tiers == 0] = total_tiers + 1

    true_df['tier'] = tiers
    true_df['sv_type'] = sv_types
    true_df['orient'] = orients

    if not args.min_sv_len is None:
        # Select only intra-chromosomal or svs that have a minimum distance between breakpoints
        is_feasible = np.array(true_df['dist'] >= args.min_sv_len, dtype=np.bool)

    if not args.targets is None and not args.target_dists is None:
        target_regions = tk_sv_utils.bed_to_region_map(args.targets, merge=True)
        res = get_df_region_dist(true_df, target_regions, use_orient=True)
        targ_dists1, targ_dists2, targs1, targs2, _, _ = res

        new_starts1 = np.array(true_df.start1)
        new_stops1 = np.array(true_df.stop1)
        new_starts2 = np.array(true_df.start2)
        new_stops2 = np.array(true_df.stop2)

        for i, (t1, t2) in enumerate(zip(targs1, targs2)):
            if not t1[0] is None and not t2[0] is None:
                new_starts1[i], new_stops1[i] = t1
                new_starts2[i], new_stops2[i] = t2

        true_df['start1'] = new_starts1
        true_df['stop1'] = new_stops1
        true_df['start2'] = new_starts2
        true_df['stop2'] = new_stops2

        true_df['targ_dist'] = np.maximum(np.array(targ_dists1), np.array(targ_dists2))
    else:
        true_df['targ_dist'] = np.zeros((len(true_df),), dtype = np.int)

    true_df['feasible'] = is_feasible

    return true_df


def merge_calls_and_gt(call_df, gt_df, call_to_gt):

    if not gt_df is None:
        gt_df.index = gt_df['name']
    else:
        call_to_gt = {}

    out_call_df = None
    for _, row in call_df.iterrows():
        sv_type = tk_sv_io.get_sv_type(row.info)
        orient = tk_sv_io.get_break_orientation(row.info)
        row['orient'] = orient

        # revert sv type name from DISTAL to TRANS to match ground truth
        # conventions
        if sv_type == 'DISTAL':
            sv_type = 'TRANS'
        row['sv_type'] = sv_type

        matches = list(call_to_gt.get(row['name'], [None]))
        # One output row per match
        for m in matches:
            row['match'] = m
            if not m is None and not gt_df is None:
                x = gt_df.loc[m]
                row['match_dist'] = max(dist_to_breaks(int((row.start1 + row.stop1)/2), x.start1, x.stop1),
                                        dist_to_breaks(int((row.start2 + row.stop2)/2), x.start2, x.stop2))
            else:
                row['match_dist'] = float('NaN')

            out_call_df = pd.concat([out_call_df, pd.DataFrame([row])], ignore_index=True)

    if not gt_df is None:
        out_call_df = pd.merge(out_call_df, gt_df, left_on='match',
                               right_on='name', how='outer', suffixes=['', '_gt'])
        out_call_df.drop(['filters_gt', 'dist'], axis=1, inplace=True)
    out_call_df.sort('name', inplace = True)

    return out_call_df

def get_all_sv_types(df):
    if df is None or df.empty:
        t = []
    else:
        t = list(set(df['sv_type']))
    # This will be equivalent to "take all types"
    t.append('NA')
    return t
