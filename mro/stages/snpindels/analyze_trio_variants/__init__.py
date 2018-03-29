#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Analyze consistency of trio variant calls
#

import trio
import tenkit.bio_io as tk_io
import numpy as np
import itertools
import os
import tenkit.pandas as p
import tenkit.hdf5
import tenkit.constants
import os.path
import tenkit.reference
import tenkit.log_subprocess
import tenkit.safe_json
import tenkit.chunk_utils as tk_chunks

__MRO__ = """
stage ANALYZE_TRIO_VARIANTS(
    in  vcf.gz     p1_vcf,
    in  h5         p1_coverage,
    in  vcf.gz     p2_vcf,
    in  h5         p2_coverage,
    in  vcf.gz     ch_vcf,
    in  h5         ch_coverage,
    in  bed        targets_file,
    in  string     reference_path,
    in  string     restrict_locus,
    in  map        regions,
    out h5         variants,
    out tsv        state_intervals,
    out json       summary,
    src py         "stages/snpindels/analyze_trio_variants",
) split using (
    in  string     locus,
)
"""

def main(args, outs):

    trio_iters = [
        trio.mk_var_tuples(args.p1_vcf, args.targets_file, args.locus),
        trio.mk_var_tuples(args.p2_vcf, args.targets_file, args.locus),
        trio.mk_var_tuples(args.ch_vcf, args.targets_file, args.locus),
    ]

    matched_vars = trio.parallel_iter(trio_iters, lambda x: x.key)
    current_phase_set = [0,0,0] # keep track of P1,P2,CH phase sets

    trio_tuples = []
    pos_ref_info = []

    for v in matched_vars:
        first_call = next(x for x in v if x is not None) # get the first of P1,P2,CH that isn't None
        (chrom, pos, ref) = (first_call.chrom, first_call.pos, first_call.ref)

        # If all the calls are filtered, don't report this as a real variant
        unfiltered_calls = filter(lambda call: (call is not None) and trio.is_unfiltered(call.record.FILTER), v)
        if len(unfiltered_calls) == 0:
            continue

        trio_alleles = []
        for (idx, call) in enumerate(v):
            if call is None:
                # Missing variants for an individual get annotated as homozygous reference
                # tuple is (genotype, is_phased, phase_set, call)
                call_phasing = ((ref, ref), True, current_phase_set[idx], call)
            else:
                call_phasing = (call.hap, call.phased, call.phase_set, call)
                current_phase_set[idx] = call.phase_set

            trio_alleles.append(call_phasing)

        trio_tuples.append(trio_alleles)
        pos_ref_info.append((chrom, pos, ref))

    if len(trio_tuples) < 2:
        outs.variants = None
        outs.state_intervals = None
        return

    # Determine the best-fitting inheritance pattern
    state_vec, tag_vec, score_vec = trio.infer_trio_pattern(trio_tuples)

    # TODO refactor this with named tuples
    df_info = []
    for (tups, (chrom, pos, ref), state, tags, score) in zip(trio_tuples, pos_ref_info, state_vec, tag_vec, score_vec):
            trio_info = []
            for ind in range(3):
                # individual info: ((al0, al1), phased, ps, call), state)
                ind_info = (tups[ind], state[ind])
                trio_info.append(ind_info)

            # want: (trio_info, chrom, pos, ref, tags, score)
            info = (trio_info, chrom, pos, ref, tags, score)
            df_info.append(info)

    # Load coverage tracks so we can attach coverage data
    cov_locus = tk_io.get_locus_info(args.locus)

    def load_cov(fn):
        if fn is None:
            return None

        df = tenkit.hdf5.read_data_frame_indexed(fn,
                [cov_locus], query_cols=['pos', 'coverage_deduped', 'mapq'])
        return df

    cov_info = [
        load_cov(args.p1_coverage),
        load_cov(args.p2_coverage),
        load_cov(args.ch_coverage),
    ]

    # Make a bed file by segmenting the state vector into constant regions
    state_poses = zip(pos_ref_info, state_vec)
    state_intervals = []
    # note: this relies on the fact that groupby only groups contiguous regions
    for (state, poses) in itertools.groupby(state_poses, key=lambda x: x[1]):
        grp = list(poses)
        ((chrom, start, ref), state) = grp[0]
        ((_    , end  , _  ), _    ) = grp[-1]

        state_intervals.append((chrom, start, end, state[0], state[1], state[2]))

    state_intervals_df = p.DataFrame(state_intervals, columns=['chrom', 'start', 'end', 'p1', 'p2', 'ch'])
    state_intervals_df.to_csv(outs.state_intervals, sep="\t")

    variants_df = variant_match_to_df(df_info, cov_info, args.regions)
    tenkit.hdf5.write_data_frame(outs.variants, variants_df)

def split(args):
    ref_fasta = tenkit.reference.open_reference(args.reference_path)
    chroms = ref_fasta.keys()
    chrom_lengths = [len(ref_fasta[k]) for k in chroms]

    primary_contigs = tenkit.reference.load_primary_contigs(args.reference_path)

    if args.restrict_locus is None:
        loci = tk_chunks.chunk_by_locus(chroms, chrom_lengths, tenkit.constants.PARALLEL_LOCUS_SIZE,
            contig_whitelist = primary_contigs, extra_args = {'__mem_gb': 12.0})
    else:
        loci = [{'locus': args.restrict_locus}]

    return {'chunks': loci, 'join': {'__mem_gb': 12.0}}

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    in_files = [out.variants for out in chunk_outs if out.variants and os.path.exists(out.variants)]
    if len(in_files) > 0:
        tenkit.hdf5.combine_data_frame_files(outs.variants, in_files)
    else:
        outs.variants = None

    variant_df = tenkit.hdf5.read_data_frame(outs.variants)
    stats = compute_trio_stats(variant_df)

    with open(outs.summary, 'w') as summary_file:
        summary_file.write(tenkit.safe_json.safe_jsonify(stats, pretty=True))

    interval_tbls = [p.read_table(o.state_intervals) for o in chunk_outs if o.state_intervals is not None]
    all_intervals = p.concat(interval_tbls)
    all_intervals.to_csv(outs.state_intervals, sep="\t")

def compute_trio_stats(df):
    # TODO compute some stats
    # ideas: % errors, % long switches, % concordance with phase blocks, % mendelian violation
    # broken down by conf vs non-conf regions
    pass

def get_cov(pos, cov_df):
    if cov_df is not None:
        # NOTE -- dtype of query to searchsorted must match the input array
        # or you can get a ~10^5 fold slowdown.
        pos_vals = cov_df.pos.values
        idx = np.searchsorted(pos_vals, np.int32(pos))

        if idx < 0 or idx >= len(pos_vals):
            return (-1, -1, -1)

        if (pos_vals[idx] != pos):
            # Should only hit this for indels
            if abs(pos_vals[idx+1] - pos) < abs(pos_vals[idx] - pos):
                idx = idx + 1

        if hasattr(cov_df, 'mapq30_coverage'):
            cov_mapq30 = cov_df.mapq30_coverage[idx]
        else:
            cov_mapq30 = -1

        return (cov_df.coverage_deduped[idx], cov_mapq30, cov_df.mapq[idx])
    else:
        return (-1, -1, -1)

def variant_match_to_df(variant_match_iter, cov_info, regions_files):

    region_prefix = "bed_"
    error_prefix = "flag_"

    # load regions
    region_sets = {}
    for (reg_name, bed_fn) in regions_files.items():
        if bed_fn is not None:
            with open(bed_fn) as bed_file:
                regions = tk_io.get_target_regions(bed_file)
            region_sets[str(reg_name)] = regions

    # define basic fields
    record_type = [
            ('chrom', object),
            ('pos', np.int32),
            ('ref', object),

            ('any_filters', np.int8),
            ('any_flags', np.int8),

            ('p1_al0', object),
            ('p1_al1', object),
            ('p1_state', np.int8),
            ('p1_phased', np.int8),
            ('p1_phase_set', np.int64),
            ('p1_cov', np.int16),
            ('p1_cov_mapq30', np.int16),
            ('p1_qual', np.int16),
            ('p1_mapq', np.int16),
            ('p1_filters', object),
            ('p1_has_variant', np.int8),

            ('p2_al0', object),
            ('p2_al1', object),
            ('p2_state', np.int8),
            ('p2_phased', np.int8),
            ('p2_phase_set', np.int64),
            ('p2_cov', np.int16),
            ('p2_cov_mapq30', np.int16),
            ('p2_qual', np.int16),
            ('p2_mapq', np.int16),
            ('p2_filters', object),
            ('p2_has_variant', np.int8),

            ('ch_al0', object),
            ('ch_al1', object),
            ('ch_state', np.int8),
            ('ch_phased', np.int8),
            ('ch_phase_set', np.int64),
            ('ch_cov', np.int16),
            ('ch_cov_mapq30', np.int16),
            ('ch_qual', np.int16),
            ('ch_mapq', np.int16),
            ('ch_filters', object),
            ('ch_has_variant', np.int8),
        ]

    # add fields for region membership
    # prepend "bed_" so we can tell them apart from other columns
    for reg_name in region_sets:
        record_type.append((region_prefix + reg_name, np.int8))

    # add fields for concordance errors
    for (field, _) in trio.error_bits:
        record_type.append((error_prefix + field, np.int8))

    # write output in blocks of recarrays, to be more memory efficient
    bsize = 100 # number of records per block
    record_block = np.recarray((bsize,), dtype=record_type)
    counter = 0

    # set fields
    results = []
    for (trio_info, chrom, pos, ref, error_tags, score) in variant_match_iter:
        # base fields
        row = record_block[counter]
        row.chrom = chrom
        row.pos = pos - 1 # make pos zero-based to match analyze_snpindels
        row.ref = ref

        # set VCF fields for each sample in trio
        prefixes = ['p1', 'p2', 'ch']
        for (prefix, al_info, cov_df) in zip(prefixes, trio_info, cov_info):
            (((al0, al1), phased, ps, call), state) = al_info
            row[prefix+"_al0"] = al0
            row[prefix+"_al1"] = al1
            row[prefix+"_state"] = state
            row[prefix+"_phased"] = phased
            row[prefix+"_phase_set"] = ps

            if call is not None:
                record = call.record
                row[prefix+"_qual"] = record.QUAL if record.QUAL is not None else -1
                row[prefix+"_filters"] = str(record.FILTER)
                if len(record.FILTER) > 0:
                    row['any_filters'] = True

            (cov, cov_q20, mapq) = get_cov(pos, cov_df)
            row[prefix+"_cov"] = cov
            row[prefix+"_cov_mapq30"] = cov_q20
            row[prefix+"_mapq"] = mapq

            row[prefix+"_has_variant"] = al0 != ref or al1 != ref

        # set error tags
        if len(error_tags) > 0:
            row['any_flags'] = True
        for error_tag in error_tags:
            row[error_prefix + error_tag] = True

        # annotate regions membership
        for (reg_name, regions) in region_sets.items():
            if chrom in regions:
                if regions[chrom].contains_point(pos):
                    row[region_prefix + reg_name] = True

        # deal with blocks
        counter += 1
        if counter == bsize:
            counter = 0
            results.append(record_block)
            record_block = np.recarray((bsize,), dtype=record_type)

    results.append(record_block[:counter])
    df = p.concat([p.DataFrame(blk) for blk in results], ignore_index=True)
    return df
