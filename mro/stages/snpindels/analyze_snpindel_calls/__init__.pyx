#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Analyze variant TP/FPs
#
import os
import tenkit.pandas as tk_pd
import numpy as np
import bisect
import itertools
import tenkit.safe_json
import tenkit.stats as tk_stats
import tenkit.hdf5
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.constants
import os.path
import tenkit.reference
import math
import tenkit.log_subprocess
import gzip
import tenkit.chunk_utils as tk_chunks

__MRO__ = """
stage ANALYZE_SNPINDEL_CALLS(
    in  string     vc_precalled,
    in  string     variant_mode                  "'freebayes' or 'gatk:/path/to/GenomeAnalysisTK.jar'",
    in  bam        bam_file,
    in  vcf.gz     ground_truth                  "ground truth variants",
    in  vcf.gz     input                         "called variants",
    in  h5         coverage,
    in  bed        targets_file                  "file specifying targets to focus analysis on",
    in  string     reference_path,
    in  string     restrict_locus,
    in  map        regions_of_interest,
    in  int        long_switch_penalty_multiple,
    in  tsv.gz     fragment_phasing,
    in  bam        validation_bam,
    out vcf.gz     varcalls,
    out vcf.gz.tbi varcalls_index,
    out h5         variants,
    out h5         phase_blocks,
    out csv        gene_stats,
    out csv        variant_stats,
    out json       summary                       "results of analysis",
    src py         "stages/snpindels/analyze_snpindel_calls",
) split using (
    in  string     locus,
)
"""

# if this is in the regions_of_interest, use it to compute gene stats on just the core regions
CORE_TARGET_REGIONS = "core_exonic_regions"

def split(args):
    input_bam = tk_bam.create_bam_infile(args.bam_file)

    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    primary_contigs = tenkit.reference.load_primary_contigs(args.reference_path)

    if args.restrict_locus is None:
        # use a smaller chunk size to reduce memory footprint
        chunk_size = tk_chunks.get_parallel_locus_size(args.targets_file) / 2
        loci = tk_chunks.chunk_by_locus(chroms, chrom_lengths, chunk_size, contig_whitelist = primary_contigs)
    else:
        loci = [{'locus': args.restrict_locus, "__mem_gb": 5}]

    return {'chunks': loci, 'join': {'__mem_gb': 16}}


def load_gene_finder(args):
    genes_file = tenkit.reference.get_gene_boundaries(args.reference_path)
    if os.path.exists(genes_file):
        # Load chromosomes from fasta to support "chr" prepending
        chroms = tenkit.reference.open_reference(args.reference_path).keys()
        return GeneFinder(genes_file, chroms)
    else:
        return GeneFinder(None, None)

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    in_files = [out.variants for out in chunk_outs if out.variants and os.path.exists(out.variants)]
    tmp_variants = outs.variants + ".tmp"

    if len(in_files) > 0:
        tenkit.hdf5.combine_data_frame_files(tmp_variants, in_files)
    else:
        outs.variants = None

    if outs.variants is not None:

        gene_finder = load_gene_finder(args)
        if outs.variants is not None:
            variant_df = tenkit.hdf5.read_data_frame(tmp_variants)
            variant_df = compute_switches(variant_df, args.long_switch_penalty_multiple)

            if len(variant_df) != 0:
                tenkit.hdf5.write_data_frame(outs.variants, variant_df)
                phase_blocks = compute_phase_blocks(variant_df)

                if len(phase_blocks) != 0:
                    tenkit.hdf5.write_data_frame(outs.phase_blocks, phase_blocks)
                else:
                    outs.phase_blocks = None
                    phase_blocks = None

            else:
                outs.variants = None
                variant_df = None
        else:
            variant_df = None

        # before using the fragment phasing, convince ourselves that it's correct
        # TODO - reinstate validation of phase_blocks w/ more memory efficient method
        # current code can cause OOM when fragment_phasing file gets large
        #if args.fragment_phasing is not None:
        #    validation.validate_phase_blocks(args.input, args.fragment_phasing)

        if args.p1_genes_list is not None:
            p1_genes = {line.rstrip() for line in open(args.p1_genes_list)}
        else:
            p1_genes = {}

        is_targeted = args.targets_file is not None
        (stats, stats_table, gene_stats) = compute_stats(variant_df, phase_blocks, gene_finder, args.ground_truth, args.fragment_phasing, args.regions_of_interest, p1_genes, is_targeted)

        if stats_table is not None:
            stats_table.to_csv(outs.variant_stats)
        else:
            outs.variant_stats = None

        if gene_stats is not None:
            gene_stats.to_csv(outs.gene_stats)
        else:
            outs.gene_stats = None

        # log if we used pre-called variants or not
        vc_mode, _, _, _ = tk_io.get_vc_mode(args.vc_precalled, args.variant_mode)
        used_precalled = (vc_mode == "precalled") or (vc_mode == "precalled_plus")
        stats['used_precalled_variants'] = used_precalled

        with open(outs.summary, 'w') as summary_file:
            summary_file.write(tenkit.safe_json.safe_jsonify(stats, pretty=True))
    else:
        outs.phase_blocks = None
        outs.gene_stats = None
        outs.summary = None

    # need to limit the vcf line length to 65536 characters for htslib parsing,
    # will remove this when our htslib fix is in the main branch for a while
    #cap_barcode_list_to_2500(args.input, outs.varcalls)
    vcfs = [chunk_out.varcalls[:-3] for chunk_out in chunk_outs if os.path.isfile(chunk_out.varcalls[:-3])]
    tk_io.combine_vcfs(outs.varcalls[:-3], vcfs)
    outs.varcalls_index = outs.varcalls + ".tbi"

def cap_barcode_list_to_2500(infile, out, chrom, start, stop):
    vcf_in = tk_io.VariantFileReader(infile)

    with open(out, 'w') as vcf_out_file:
        vcf_out = tk_io.VariantFileWriter(vcf_out_file, template_file=open(infile))
        for record in vcf_in.record_getter(restrict_type=None, fetch_chrom=chrom, fetch_start=start, fetch_end=stop):
            barcode_list = tk_io.get_record_barcodes(record)
            if barcode_list is not None:
                num_barcodes = 0
                for barcodes in barcode_list:
                    num_barcodes += len(barcodes)
                if num_barcodes > 2500:
                    new_barcodes = []
                    for barcodes in barcode_list:
                        num_barcodes_to_remove = int(math.ceil(float(num_barcodes - 2500)*float(len(barcodes)/float(num_barcodes))))
                        new_barcodes.append(barcodes[:len(barcodes)-num_barcodes_to_remove])
                    tk_io.set_record_barcodes(record, new_barcodes)
            vcf_out.write_record(record)

def is_primitive(ref, alt):
    return len(ref) == 1 or len(alt) == 1

class GeneFinder:
    def __init__(self, genes_file, chroms):
        if genes_file is None:
            self.chrom_gene_starts = {}
            self.chrom_genes = {}
            return

        self.chrom_genes = self.load_genes(genes_file, chroms)

        self.chrom_gene_starts = {}
        for (chrom, genes) in self.chrom_genes.iteritems():
            self.chrom_gene_starts[chrom] = [x[0] for x in genes]


    def simplify_region_set(self, regs):
        ''' Take a set of regions, and simplify it to a set of disjoint regions.
            Overlapping regions are combined into a single region '''

        out_regs = []
        regs = sorted(regs)
        i = 0

        while i < len(regs):
            s,e,g = regs[i]
            next_i = i + 1
            for j in range(i+1, len(regs)):
                if e > regs[j][0]:
                    e = max(e, regs[j][1])
                    next_i = j + 1

            i = next_i
            out_regs.append((s,e,g))

        return out_regs


    def load_genes(self, genes_file_name, chroms):
        gene_start_ends = {}
        genes_file = gzip.open(genes_file_name)
        for line in genes_file:
            if line.startswith("#"):
                continue
            info = line.strip().split('\t')
            if info[2] != "gene":
                continue
            chrom = info[0]
            if chrom not in chroms and chrom[0:3] != 'chr' and "chr"+chrom in chroms:
                chrom = "chr"+chrom
            gene_start = int(info[3])
            gene_end = int(info[4])
            if len(info) > 8:
                gene_ids = {attr.split()[0]: attr.split()[1] for attr in info[8].split(";") if len(attr.split()) > 1}
                gene_name = gene_ids.get('gene_name')
                if gene_name != None:
                    gene_name = gene_name.translate(None,'"')

            else:
                # Must have gene name!
                continue
                # Give a default name if bed file doesn't have names
                #gene_name = chrom + "_" + str(gene_start)

            chrom_genes = gene_start_ends.setdefault(chrom, [])
            chrom_genes.append((gene_start, gene_end, gene_name))

        gene_results = {}
        for (chrom, gene_list) in gene_start_ends.iteritems():
            chrom_regs = []
            sorted_genes = sorted(set(gene_list), key=lambda x: x[2])

            for (grp, regs) in itertools.groupby(sorted_genes, key=lambda x: x[2]):
                simplified_regs = self.simplify_region_set(list(regs))
                chrom_regs.extend(simplified_regs)

            gene_results[chrom] = sorted(chrom_regs)

        return gene_results

    def get_genes_for_chrom(self, chrom):
        return self.chrom_genes.get(chrom, [])

    def get_gene(self, chrom, pos):
        if not self.chrom_gene_starts.has_key(chrom):
            return None
        starts = self.chrom_gene_starts[chrom]
        genes = self.chrom_genes[chrom]

        idx = bisect.bisect_right(starts, pos) - 1

        if idx < 0:
            return None

        gene = genes[idx]
        (start,end,name) = gene
        if pos >= start and pos < end:
            return name
        else:
            return None


def compute_switches(df, long_switch_multiple):
    ''' Mark phasing errors on variant data '''

    def mark_switches(blk):
        ''' Use a simple HMM to assign phasing errors as 'short' or 'long'
            with a tunable tradeoff between them. Following the CPT-Seq paper
            (http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3119.html)
            we weight long to short 5:1. '''

        short_score = -1
        long_score = -1 * long_switch_multiple

        parity = blk.parity.values
        n = len(parity)

        if n == 0:
            blk['short_switch'] = []
            blk['long_switch'] = []
            return blk

        score_mat = np.ones((n,2), dtype=np.int) * -999999999
        trace_mat = np.ones((n,2), dtype=np.int) * -999999999

        # Initial condition
        for p in [0,1]:
            score_mat[0,p] = 0 if p == parity[0] else short_score

        # Recursion
        for i in range(1, n):
            for p in [0,1]:
                for prev_p in [0,1]:
                    _prev = score_mat[i-1, prev_p]
                    _long = 0 if prev_p == p else long_score
                    _short = 0 if p == parity[i] else short_score
                    score = _prev + _long + _short

                    if score > score_mat[i, p]:
                        score_mat[i,p] = score
                        trace_mat[i,p] = prev_p

        # Trace back from best final state
        base_parity = np.zeros(n, dtype=np.int)
        base_parity[n-1] = 0 if score_mat[n-1, 0] > score_mat[n-1,1] else 1

        for i in range(n-2, -1, -1):
            base_parity[i] = trace_mat[i+1, base_parity[i+1]]

        prev_base_parity = np.concatenate([base_parity[0:1], base_parity[:-1]])
        short_switch = parity != base_parity
        long_switch = base_parity != prev_base_parity

        blk['short_switch'] = short_switch
        blk['long_switch'] = long_switch
        return blk

    # can_test_phase only includes phased hets that have phased ground truth available
    # and are marked
    df_test = df[df.can_test_phase>0]
    if len(df_test) > 0:
        phase_blocks = df_test.groupby(['gt_phase_set', 'phase_set', 'chrom'])
        df_switches = phase_blocks.apply(mark_switches)

        df['long_switch'] = df_switches.long_switch
        df['short_switch'] = df_switches.short_switch
        df['long_switch'] = df.long_switch.fillna(False)
        df['short_switch'] = df.short_switch.fillna(False)
    else:
        df['long_switch'] = False
        df['short_switch'] = False

    return df

def compute_phase_blocks(df):
    ''' Summarize phase blocks '''

    def summarize_phase_block(_block, row):

        n = len(_block)
        start = _block.pos.min()
        end = _block.pos.max()
        chrom = _block.chrom.iat[0]
        phase_set = _block.phase_set.iat[0]
        length = end-start

        # All the phasing stats are computed only on the hets
        block = _block[(_block.FILTER == 1) & (_block.phase_set>0) & (_block.obs_genotype1 != _block.obs_genotype2)]

        # Compute the probability that a pair of SNPs in the block are correctly phased
        n_hets = len(block)
        block.reset_index(inplace=True, drop=True)

        testable_rows = block[block.can_test_phase > 0]
        n_test = len(testable_rows)

        # Determine the number of correctly phased pairs within on phase block
        def correct_pairs(bl_rows):
            nt = len(bl_rows)
            parity0 = (bl_rows.parity == 0).sum()
            parity1 = (bl_rows.parity == 1).sum()
            correct_pairs = (parity0 * (parity0-1))/2 + (parity1 * (parity1-1))/2
            return tk_pd.DataFrame({'frac_correct': tk_stats.robust_divide(float(correct_pairs),  (nt * (nt-1) / 2)), 'n_test':nt}, index=[0])

        fract_phased = block.obs_phased.mean()

        if n_test < 2:
            pair_correct_rate = float('NaN')
        else:
            testable_groups = testable_rows.groupby(['gt_phase_set'])
            testable_block_summary = testable_groups.apply(correct_pairs)
            testable_block_summary = testable_block_summary[testable_block_summary.n_test > 1]
            pair_correct_rate = tk_stats.robust_divide((testable_block_summary.n_test.values * testable_block_summary.frac_correct.values).sum(),  testable_block_summary.n_test.sum())

        row.chrom = chrom
        row.phase_set = phase_set
        row.start = start
        row.end = end
        row.length = length
        row.num_variants = n
        row.num_hets = n_hets
        row.num_testable_variants = n_test
        row.pair_correct_rate = pair_correct_rate
        row.fract_phased = fract_phased

    # Only include hets in phase block stats
    df_test = df[np.logical_and(df.phase_set>0, df.obs_genotype1 != df.obs_genotype2)]

    record_type = [
            ('chrom', object),
            ('phase_set', np.int64),
            ('start', np.int32),
            ('end', np.int32),
            ('length', np.int32),
            ('num_variants', np.int32),
            ('num_hets', np.int32),
            ('num_testable_variants', np.int32),
            ('pair_correct_rate', np.float),
            ('fract_phased', np.float) ]

    phase_blocks = df_test.groupby(['chrom', 'phase_set'])

    # Pre-allocate a recarray for the phase block data to be memory efficient
    # Using a standard 'apply' will cause the creation of a large number of DataFrame objects
    # if there are many phase blocks -- this can cause the memory consumption to blow up.
    record_block = np.recarray((len(phase_blocks),), dtype=record_type)

    counter = 0
    for (phase_block_id, variants) in phase_blocks:
        row = record_block[counter]
        summarize_phase_block(variants, row)
        counter += 1

    phase_block_summary = tk_pd.DataFrame(record_block)
    return phase_block_summary


def compute_snpindel_stats(filter_condition, var_type, region_name, variant_df, ground_truth):
    unfiltered = filter_condition == "unfiltered"
    any_type = var_type == "any"
    type = ""
    if var_type == "snp":
        type = "S"
    elif var_type == "insertion":
        type = "I"
    elif var_type == "deletion":
        type = "D"
    elif var_type == "complex":
        type = "C"
    elif var_type == "insertion_lt5bp":
        type = "I"
    elif var_type == "deletion_lt5bp":
        type = "D"
    elif var_type == "insertion_5_to_50bp":
        type = "I"
    elif var_type == "deletion_5_to_50bp":
        type = "D"
    else:
        type = "any"
    suffix = "_" + filter_condition
    if not var_type == "any":
        suffix = suffix + "_" + var_type
    no_length_filter = not ("lt5bp" in var_type or "50bp" in var_type)
    min_length = 0
    max_length = 5000
    if "lt5bp" in var_type:
        min_length = 0
        max_length = 5
    elif "50bp" in var_type:
        min_length = 5
        max_length = 51
    region_suffix = ""
    ana_variants = variant_df
    if region_name is not None:
        region_suffix = "_"+str(region_name)
        ana_variants = variant_df[variant_df[str(region_name)] == 1]
    suffix += region_suffix

    # Setup stats
    reg_name = region_name if region_name is not None else "all"
    stats = tk_pd.DataFrame({ "suffix": suffix,  "filter": filter_condition, "variant_type": var_type, "region": reg_name }, index=[0])

    ana_variants = ana_variants[np.logical_or(ana_variants.variant_type == type, any_type)]
    ana_variants = ana_variants[np.logical_or(np.logical_and(ana_variants.variant_length.abs() < max_length, ana_variants.variant_length.abs() >= min_length), no_length_filter)]
    ana_variants['is_call'] = np.logical_and(ana_variants.in_obs == 1, np.logical_or(ana_variants.FILTER == 1, unfiltered))

    variant_metrics(stats, ana_variants, ground_truth)
    return stats


def variant_metrics(stats, filtered, ground_truth):
    ''' Stats about some set of variants '''

    tps = len(filtered[np.logical_and(filtered.is_call, filtered.in_gt == 1)])
    fns = len(filtered[np.logical_and(filtered.is_call == False, filtered.in_gt == 1)])

    # novel / FP calls:
    novel = filtered[np.logical_and(filtered.is_call, filtered.in_gt == 0)]
    fps = len(novel)

    filtered = filtered[filtered.is_call]
    stats['titv'] = tk_stats.robust_divide(float(filtered.ti.sum()), filtered.tv.sum())
    stats['phased'] = tk_stats.robust_divide(float(filtered.obs_phased.sum()), float(tps+fps))

    if ground_truth != None:
        stats['tps'] = tps
        stats['fps'] = fps
        stats['fns'] = fns
        stats['sensitivity'] = tk_stats.robust_divide(float(tps), float(tps + fns))
        stats['ppv'] = tk_stats.robust_divide(float(tps), float(tps + fps))

        # Validation rate of 'novel' calls
        novel_valid = (novel.validation_ao >= 2) & (novel.validation_ao / novel.validation_cov >= 0.1)
        novel_validatable = (novel.validation_cov >= 12) | novel_valid

        stats['novel_calls'] = fps
        stats['novel_validate_support'] = novel_valid.sum()
        stats['novel_validate_cov'] = novel_validatable.sum()
        stats['novel_validate_ppv'] = tk_stats.robust_divide(float(stats['novel_validate_support']), stats['novel_validate_cov'])
        stats['novel_titv'] = tk_stats.robust_divide(float(novel.ti.sum()), novel.tv.sum())

    else:
        nan = float('NaN')
        stats['tps'] = nan
        stats['fps'] = nan
        stats['fns'] = nan
        stats['sensitivity'] = nan
        stats['ppv'] = nan
        stats['novel_calls'] = nan
        stats['novel_validate_support'] = nan
        stats['novel_validate_cov'] = nan
        stats['novel_validate_ppv'] = nan
        stats['novel_titv'] = nan



def compute_stats(variant_df, phase_blocks, gene_finder, ground_truth, fragment_phasing, regions_of_interest, p1_genes, is_targeted):
    ''' Compute summary stats on phasing quality '''

    # Compute stats about called variants vs ground truth
    filter_states = ["unfiltered", "filtered"]
    var_types = ["any", "snp", "insertion", "deletion", "complex", "insertion_lt5bp", "deletion_lt5bp", "insertion_5_to_50bp", "deletion_5_to_50bp"]

    if regions_of_interest is not None:
        #region_filter = regions_of_interest
        region_filter = {}
        for (name,fn) in regions_of_interest.iteritems():
            if os.path.exists(fn):
                region_filter[name] = fn
    else:
        region_filter = {}

    # for exome, use the intersection of V6 with conf regions for evaluation
    if is_targeted and 'conf_regions' in region_filter and 'V6' in region_filter:
        variant_df['conf_regions'] = variant_df.V6 & variant_df.conf_regions

    region_filter[None] = None

    # Compute all snpindel stats cases:
    stat_rows = []
    for filter_condition in filter_states:
        for var_type in var_types:
            for (region_name,_) in region_filter.iteritems():
                stat = compute_snpindel_stats(filter_condition, var_type, region_name, variant_df, ground_truth)
                stat_rows.append(stat)

    stats_table = tk_pd.concat(stat_rows, ignore_index=True)

    # Copy over the metrics we want in the json summary.
    stats = {}
    stats_of_interest = stats_table[(stats_table['filter'] == "filtered") & (stats_table['region'] == 'conf_regions')]

    for (i, row) in stats_of_interest.iterrows():
        for metric in ['sensitivity', 'ppv', 'tps', 'fps', 'phased']:
            key = metric + row.suffix
            value = row[str(metric)]
            stats[key] = value

    # novel calls
    stats_of_interest = stats_table[(stats_table['filter'] == "filtered") & (stats_table['region'] == "degen_regions")]

    for (i, row) in stats_of_interest.iterrows():
        for metric in ['novel_calls', 'novel_validate_ppv']:
            key = metric + row.suffix
            value = row[str(metric)]
            stats[key] = value


    filtered_df = variant_df[np.logical_and(variant_df.in_obs == 1, variant_df.FILTER == 1)]
    stats['snps'] = len(filtered_df[filtered_df.variant_type == "S"])
    stats['insertions'] = len(filtered_df[filtered_df.variant_type == "I"])
    stats['deletions'] = len(filtered_df[filtered_df.variant_type == "D"])

    df_hets = filtered_df[filtered_df.obs_genotype1 != filtered_df.obs_genotype2]
    df_hets = df_hets.drop(['alt1', 'alt2', 'allele'], 1)
    df_hets.drop_duplicates(inplace=True)

    df_snps_testable = df_hets[np.logical_and(df_hets.variant_type == 'S', df_hets.can_test_phase>0)]
    if len(df_snps_testable) > 0:
        # Switch error rate on phased SNPs
        stats["snp_long_switch_error"] = df_snps_testable.long_switch.mean()
        stats["snp_short_switch_error"] = df_snps_testable.short_switch.mean()
    del df_snps_testable

    df_indels_testable = df_hets[np.logical_and(df_hets.variant_type.isin(['I', 'D']), df_hets.can_test_phase>0)]
    if len(df_indels_testable) > 0:
        # Switch error rate on phased INDELSs
        stats["indel_short_switch_error"] = df_indels_testable.short_switch.mean()
    del df_indels_testable


    df_testable = df_hets[df_hets.can_test_phase>0]
    if len(df_testable) > 0:
        # Switch error rate on phased SNPs
        stats["long_switch_error"] = df_testable.long_switch.mean()
        stats["short_switch_error"] = df_testable.short_switch.mean()
    del df_testable


    if fragment_phasing is not None:
        frag_mix_count = 0
        frag_conf_phased_count = 0
        total = 0
        fragment_phasing_file = gzip.open(fragment_phasing, 'rb')
        for line in fragment_phasing_file:
            if line.startswith("#"):
                continue
            info = line.strip("\n").split("\t")
            h0 = float(info[7])
            h1 = float(info[8])
            hm = float(info[9])
            total += 1
            if hm >= 0.95:
                frag_mix_count += 1
            elif h0 >= 0.995 or h1 >= 0.995:
               frag_conf_phased_count += 1

        stats['frac_fragment_with_het_mixed_phase'] = tk_stats.robust_divide(float(frag_mix_count), float(total))
        stats['frac_fragment_with_het_phased'] = tk_stats.robust_divide(float(frag_conf_phased_count), float(total))

    if phase_blocks is None or len(phase_blocks) == 0:
        stats["N50_phase_block"] = float('NaN')
        stats["N90_phase_block"] = float('NaN')
        stats["longest_phase_block"] = float('NaN')
        stats["mean_phase_block"] = float('NaN')
        stats["median_phase_block"] = float('NaN')
        stats["prob_variant_correct_in_block"] = float('NaN')

        stats["prob_phased_10k_gap"] = float('NaN')
        stats["prob_phased_20k_gap"] = float('NaN')
        stats["prob_phased_40k_gap"] = float('NaN')

    else:
        lengths = phase_blocks.length
        stats["N50_phase_block"] = int(tk_stats.N50(lengths.values))
        stats["N90_phase_block"] = int(tk_stats.NX(lengths.values, 0.9))
        stats["longest_phase_block"] = int(lengths.max())
        stats["mean_phase_block"] = float(lengths.mean())
        stats["median_phase_block"] = int(lengths.median())
        stats["prob_variant_correct_in_block"] = float(tk_stats.robust_divide((phase_blocks.pair_correct_rate * phase_blocks.num_testable_variants).sum(), phase_blocks.num_testable_variants.sum()))

        # Stats for probability of phasing over gaps of various lengths
        poses = df_hets.pos.values
        phase_set = df_hets.phase_set.values
        obs_phased = df_hets.obs_phased.values

        dist = poses[1:] - poses[:-1]
        phased = np.logical_and(phase_set[:-1] == phase_set[1:], np.logical_and(obs_phased[:-1], obs_phased[1:]))

        pair_df = tk_pd.DataFrame({'dist':dist, 'phased':phased})
        stats["prob_phased_10k_gap"] = pair_df[(pair_df.dist - 10000).abs() < 1000].phased.mean()
        stats["prob_phased_20k_gap"] = pair_df[(pair_df.dist - 20000).abs() < 2000].phased.mean()
        stats["prob_phased_40k_gap"] = pair_df[(pair_df.dist - 40000).abs() < 4000].phased.mean()


    df_snps = df_hets[df_hets.variant_type == 'S']
    df_indels = df_hets[df_hets.variant_type.isin(['I', 'D'])]

    # Fraction of SNPs with phase called
    stats["fract_snps_phased"] = df_snps.obs_phased.mean()
    stats["fract_indels_phased"] = df_indels.obs_phased.mean()

    # Fraction of SNPs w/ one or more barcode
    stats["fract_snps_barcode"] = np.logical_or(df_snps.bcs1 > 0, df_snps.bcs2 > 0).mean()
    stats["fract_indels_barcode"] = np.logical_or(df_indels.bcs1 > 0, df_indels.bcs2 > 0).mean()

    # Fraction w/ barcodes on both alleles
    stats["fract_snps_barcode_both_alleles"] = np.logical_and(df_snps.bcs1 > 0, df_snps.bcs2 > 0).mean()
    stats["fract_indels_barcode_both_alleles"] = np.logical_and(df_indels.bcs1 > 0, df_indels.bcs2 > 0).mean()

    def gene_phase_summary(rows, gene, chrom, start, end):
        '''summarize the phasing of a single gene'''

        length = end - start
        is_p1 = gene in p1_genes
        summary = {'gene': gene, 'chrom': chrom, 'start':start, 'end':end, 'length':length, 'is_p1':is_p1}

        full_regions_summary = gene_phase_summary_helper(rows)
        core_regions_summary = gene_phase_summary_helper(rows[rows[CORE_TARGET_REGIONS] == 1])

        summary.update(full_regions_summary)
        summary.update({"core_"+k: v for (k,v) in core_regions_summary.iteritems()})

        return tk_pd.DataFrame(summary, index=[0])

    def gene_phase_summary_helper(rows):
        # HACK: ignore dots for now due to bug in phase set labeling
        rows = rows[(rows.phase_set != -1)]

        n_snps = len(rows)
        rows.reset_index(inplace=True, drop=True)

        testable_rows = rows[rows.can_test_phase > 0]
        n_test = len(testable_rows)

        n_phase_blocks = len(rows.phase_set.unique())

        gene_phased = (n_phase_blocks == 1)
        completely_phased = (rows.obs_phased == 1).all() and (n_phase_blocks == 1)

        # Gene phased metrics don't apply to genes with <2 variants
        if n_snps == 1:
            completely_phased = float('NaN')
            gene_phased = float('NaN')

        if n_snps < 2:
            pair_phased_rate = float('NaN')
            num_phase_blocks = 1
        else:
            # Group the SNPs in the gene by their phase block
            phase_set_groups = rows.groupby(['phase_set'])
            in_gene_pairs_phased = phase_set_groups.apply(phased_pairs).sum()
            pair_phased_rate = float(in_gene_pairs_phased) / (n_snps * (n_snps-1) / 2)
            num_phase_blocks = len(phase_set_groups)

        if n_test < 2:
            pair_correct_rate = float('NaN')
        else:
            testable_groups = testable_rows.groupby(['phase_set'])
            in_gene_pairs_correct = testable_groups.apply(correct_pairs).sum()
            testable_block_sizes = testable_groups.apply(len)
            in_gene_pairs_testable = (testable_block_sizes * (testable_block_sizes - 1) / 2).sum()
            pair_correct_rate = tk_stats.robust_divide(in_gene_pairs_correct, in_gene_pairs_testable)

        return {
            'num_variants':n_snps, 'num_testable_variants': n_test, 'phased': gene_phased, 'completely_phased': completely_phased,
            'pair_phased_rate': pair_phased_rate, 'pair_correct_rate': pair_correct_rate, 'num_phase_blocks': num_phase_blocks
        }

    # Determine the number of correctly phased pairs within on phase block
    def correct_pairs(bl_rows):
        parity0 = (bl_rows.parity == 0).sum()
        parity1 = (bl_rows.parity == 1).sum()
        correct_pairs = (parity0 * (parity0-1))/2 + (parity1 * (parity1-1))/2
        return correct_pairs

    # Number of pairs that are phased
    def phased_pairs(bl_rows):
        n_phased = bl_rows.obs_phased.sum()
        return (n_phased * (n_phased - 1)) / 2

    # Pull out SNPs in genes
    if len(df_snps) > 0:
        gene_stats_rows = []

        # if CORE_TARGET_REGIONS was not included, add it as all zeros
        if CORE_TARGET_REGIONS not in df_snps.columns:
            df_snps[CORE_TARGET_REGIONS] = 0

        for (chrom, genes) in gene_finder.chrom_genes.iteritems():
            dfc = df_snps[(df_snps.chrom == chrom)]
            dfc.set_index('pos', inplace=True)
            for (start,end,gene) in genes:
                gene_rows = dfc.ix[start:end]
                gene_stats_rows.append(gene_phase_summary(gene_rows, gene, chrom, start, end))

        if len(gene_stats_rows) > 0:
            gene_stats = tk_pd.concat(gene_stats_rows)

            gt1_snp_genes = gene_stats[gene_stats.num_variants > 1]

            stats['fract_genes_phased'] = gt1_snp_genes.phased.mean()
            stats['num_genes_phaseable'] = len(gt1_snp_genes)

            short_genes = gt1_snp_genes[gt1_snp_genes.length < 250e3]
            stats['fract_genes_lt_250kb_phased'] = short_genes.phased.mean()
            stats['num_genes_lt_250kb_phaseable'] = len(short_genes)

            short_genes = gt1_snp_genes[gt1_snp_genes.length < 100e3]
            stats['fract_genes_lt_100kb_phased'] = short_genes.phased.mean()
            stats['num_genes_lt_100kb_phaseable'] = len(short_genes)

            stats['fract_genes_completely_phased'] = gt1_snp_genes.completely_phased.mean()

            stats['prob_snp_phased_in_gene'] = \
                tk_stats.robust_divide((gt1_snp_genes.pair_phased_rate.values * gt1_snp_genes.num_variants.values).sum(), gt1_snp_genes.num_variants.sum())

            gt1_snp_genes = gene_stats[np.logical_and(gene_stats.num_testable_variants > 1, np.logical_not(np.isnan(gene_stats.pair_correct_rate.values)))]
            stats['prob_snp_correct_in_gene'] = \
                tk_stats.robust_divide((gt1_snp_genes.pair_correct_rate.values * gt1_snp_genes.num_testable_variants.values).sum(), gt1_snp_genes.num_testable_variants.sum())
        else:
            stats['fract_genes_phased'] = float('NaN')
            stats['fract_genes_completely_phased'] = float('NaN')
            stats['prob_snp_phased_in_gene'] = float('NaN')
            stats['prob_snp_correct_in_gene'] = float('NaN')
            gene_stats = None
    else:
        stats['fract_genes_phased'] = float('NaN')
        stats['fract_genes_completely_phased'] = float('NaN')
        stats['prob_snp_phased_in_gene'] = float('NaN')
        stats['prob_snp_correct_in_gene'] = float('NaN')
        gene_stats = None

    return (stats, stats_table, gene_stats)

def main(args, outs):

    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    chrom = str(chrom)

    if args.vc_precalled is not None:
        args.input = args.vc_precalled

    cap_barcode_list_to_2500(args.input, outs.varcalls[:-3], chrom, start, stop)

    reference_pyfasta = tenkit.reference.open_reference(args.reference_path)

    if args.targets_file is None:
        targets = None
    else:
        targets = args.targets_file

    if args.ground_truth is not None:
        tups_gt = get_normalized_variant_tuples(args.ground_truth, outs, False, targets, args.locus, reference_pyfasta)
    else:
        tups_gt = []
    tups_obs = get_normalized_variant_tuples(args.input, outs, True, targets, args.locus, reference_pyfasta)
    if args.comparison_vcf is not None:
        tups_comp = get_normalized_variant_tuples(args.comparison_vcf, outs, False, targets, args.locus, reference_pyfasta)
    else:
        tups_comp = []
    var_pairs = trio_iter(iter(tups_gt), iter(tups_obs), iter(tups_comp), lambda x: x[0])
    if args.coverage is not None:
        # Deduped coverage data
        cov_df = tenkit.hdf5.read_data_frame_indexed(args.coverage, [(chrom,start-10000,stop+10000)], query_cols=['chrom', 'pos', 'coverage_deduped', 'bc_coverage'])
    else:
        cov_df = None

    gene_finder = load_gene_finder(args)

    bc_cov = None
    if cov_df is not None and 'bc_coverage' in cov_df.columns:
        bc_cov = cov_df.bc_coverage

    if args.validation_bam is not None:
        validation_bam = tk_bam.create_bam_infile(args.validation_bam)
    else:
        validation_bam = None
    sample_bam = tk_bam.create_bam_infile(args.bam_file)

    # Construct the variant data-frame
    if cov_df is not None:
        df = variant_pairs_to_df(cov_df.pos.values, cov_df.coverage_deduped.values, bc_cov, var_pairs, start, gene_finder, args.regions_of_interest, validation_bam, reference_pyfasta, sample_bam)
    else:
        df = variant_pairs_to_df(None, None, None, var_pairs, start, gene_finder, args.regions_of_interest, validation_bam, reference_pyfasta, sample_bam)

    # Write results
    if len(df) != 0:
        tenkit.hdf5.write_data_frame(outs.variants, df)

def get_normalized_variant_tuples(vcf, outs, keep_filtered, targets, locus, reference_pyfasta):
    vcf_r = tk_io.VariantFileReader(vcf)
    vcfname = vcf.split("/")[-1]
    with open(outs.variants+vcfname+"norm.vcf", 'w') as norm:
        normalizable_vcf = tk_io.VariantFileWriter(norm, template_file = open(vcf))
        with open(outs.variants+vcfname+"non_norm.vcf", 'w') as non_norm:
            non_normalizable_vcf = tk_io.VariantFileWriter(non_norm, template_file = open(vcf))
            variant_iter = tk_io.get_variant_iterator_pos(vcf_r, targets, locus)
            write_variants_by_normalizable(variant_iter, reference_pyfasta, normalizable_vcf, non_normalizable_vcf)
    with open(outs.variants+vcfname+"non_norm_primitives.vcf", 'w') as prim:
        tenkit.log_subprocess.check_call(['vcfallelicprimitives', '-t','VCFALLELICPRIMITIVE', outs.variants+vcfname+"non_norm.vcf"],stdout=prim)
    with open(outs.variants+vcfname+"non_norm_primitives_cleaned.vcf",'w') as clean:
        tenkit.log_subprocess.check_call(['bcftools','filter',outs.variants+vcfname+"non_norm_primitives.vcf"],stdout=clean)
    tk_io.combine_vcfs(str(outs.variants+vcfname+"final.vcf"),[outs.variants+vcfname+"non_norm_primitives_cleaned.vcf", outs.variants+vcfname+"norm.vcf"])
    vcf_r = tk_io.VariantFileReader(outs.variants+vcfname+"final.vcf.gz")
    tups = variant_tuples(tk_io.get_variant_iterator_exclusive(vcf_r, targets, None), None, reference_pyfasta)
    if keep_filtered:
        tups_list = {var_info: record for ((var_info), record) in tups}
    else:
        tups_list = {var_info: record for ((var_info), record) in tups if tk_io.get_record_passes_filters(record)}
    tups_list = [((var_info[0], var_info[1], var_info[2], var_info[3], var_info[4]), tups_list[var_info]) for var_info in tups_list]
    tups_list = sorted(tups_list, key = lambda x: x[0])
    #remove dups
    tups = []
    for index in range(len(tups_list)):
        if index < len(tups_list) - 1:
            if tups_list[index][0][0:4] != tups_list[index + 1][0][0:4]:
                tups.append(((tups_list[index][0][0], tups_list[index][0][1], tups_list[index][0][2], tups_list[index][0][3]), tups_list[index][1]))
        else:
            tups.append(((tups_list[index][0][0], tups_list[index][0][1], tups_list[index][0][2], tups_list[index][0][3]), tups_list[index][1]))
    return tups


def write_variants_by_normalizable(variant_iter, reference_pyfasta, normalizable_vcf, non_normalizable_vcf):
    for record in variant_iter:
        chrom = tk_io.get_record_chrom(record)
        pos = tk_io.get_record_pos(record)
        ref = tk_io.get_record_ref(record)
        alt_alleles = tk_io.get_record_alt_alleles(record)
        if len(alt_alleles) == 1:
            (new_pos, new_ref, new_alt) = tk_io.normalize_variant(pos, ref, alt_alleles[0], reference_pyfasta[chrom])
            record.POS = new_pos
            record.ALT = [new_alt]
            record.REF = new_ref
            if tk_io.get_var_type(new_ref, new_alt) != "C":
                normalizable_vcf.write_record(record)
            else:
                non_normalizable_vcf.write_record(record)
        else:
            v = [((chrom, pos, ref, allele), record) for allele in alt_alleles]
            v = sorted(v, key=lambda x: x[0])
            normalizable = True
            for x in v:
                ((x_chrom, x_pos, x_ref, x_alt), record) = x
                (new_pos, new_ref, new_alt) = tk_io.normalize_variant(x_pos, x_ref, x_alt, reference_pyfasta[chrom])
                if tk_io.get_var_type(new_ref, new_alt) == "C":
                    normalizable = False
            if normalizable:
                normalizable_vcf.write_record(record)
            else:
                non_normalizable_vcf.write_record(record)

def variant_tuples(variant_iter, targets_file_name, reference_pyfasta):
    """
    If targets file is specified then grabs variants only from targets region, otherwise
    gets all variants in file.
    Returns sets for: SNPs, INS, DELS, COMPLEX for both total and filtered variants
    """
    for record in variant_iter:
        chrom = tk_io.get_record_chrom(record)
        pos = tk_io.get_record_pos(record)
        ref = tk_io.get_record_ref(record)
        alt_alleles = tk_io.get_record_alt_alleles(record)
        filters = tk_io.get_record_passes_filters(record)
        sample = record.samples[0]
        all_alleles = [ref] + alt_alleles
        try:
            used_allele_idxs = [int(index) for index in sample.gt_alleles]
            used_alleles = [all_alleles[index] for index in used_allele_idxs] #dont want to count alleles that dont appear in genotype
        except:
            continue # this is for ./. genotypes
        if len(alt_alleles) == 1:
            if alt_alleles[0] in used_alleles:
                (new_pos, new_ref, new_alt) = tk_io.normalize_variant(pos, ref, alt_alleles[0], reference_pyfasta[chrom])
                yield ((chrom, new_pos - 1, new_ref, new_alt, filters), record)
        else:
            # yield one tuple per alt allele present in the genotype
            v = [((chrom, pos, ref, allele), record) for allele in alt_alleles]
            v = sorted(v, key=lambda x: x[0])
            for x in v:
                ((x_chrom, x_pos, x_ref, x_alt), record) = x
                if x_alt in used_alleles: # this is to exclude cases where the vcf has more alt alleles than it has genotypes
                    (new_pos, new_ref, new_alt) = tk_io.normalize_variant(x_pos, x_ref, x_alt, reference_pyfasta[chrom])
                    yield ((x_chrom, new_pos - 1, new_ref, new_alt, filters), record)

def trio_iter(i1, i2, i3, key):
    v1 = None
    v2 = None
    v3 = None

    while True:
        if v1 is None:
            try:
                v1 = i1.next()
            except StopIteration:
                pass
        if v2 is None:
            try:
                v2 = i2.next()
            except StopIteration:
                pass
        if v3 is None:
            try:
                v3 = i3.next()
            except StopIteration:
                pass
        if v1 is None and v2 is None and v3 is None:
            break
        if v1 is None:
            k1 = ("z",500000000) # need something to compare > than any chromosome name, is there a better way?
        else:
            k1 = key(v1)
        if v2 is None:
            k2 = ("z",500000000)
        else:
            k2 = key(v2)
        if v3 is None:
            k3 = ("z",500000000)
        else:
            k3 = key(v3)

        if k1 == k2 and k1 == k3:
            yield (v1, v2, v3)
            v1 = None
            v2 = None
            v3 = None
        elif k1 == k2:
            if k1 < k3:
                yield (v1, v2, None)
                v1 = None
                v2 = None
            else:
                yield (None, None, v3)
                v3 = None
        elif k1 == k3:
            if k1 < k2:
                yield (v1, None, v3)
                v1 = None
                v3 = None
            else:
                yield (None, v2, None)
                v2 = None
        elif k2 == k3:
            if k2 < k1:
                yield (None, v2, v3)
                v2 = None
                v3 = None
            else:
                yield (v1, None, None)
                v1 = None
        else:
            if k1 < k2 and k1 < k3:
                yield (v1, None, None)
                v1 = None
            elif k2 < k1 and k2 < k3:
                yield (None, v2, None)
                v2 = None
            elif k3 < k1 and k3 < k2:
                yield (None, None, v3)
                v3 = None
def pair_iter(i1, i2, key):
    ''' Iterate through sorted i1 and i2 simulateneously. Emit (v1, v2) for items
        that compare equal. Emit (v1, None) for value v1 from i1 that doesn't have
        an equal item in i2. Emit (None, v2) for value v2 from i2 that doesn't have
        and equal items in i1.  Comparisons are made with key(v).  Assumes
        that i1 and i2 are sorted with respect to key(). Also assumes that
        a sigle iterator does not contain items equal to each other.'''

    v1 = None
    v2 = None

    while True:

        if v1 is None:
            try:
                v1 = i1.next()
            except StopIteration:
                # return the rest of i2
                if v2 is not None:
                    yield (None, v2)
                for x2 in i2:
                    yield (None, x2)
                break

        if v2 is None:
            try:
                v2 = i2.next()
            except StopIteration:
                # return the rest of i2
                if v1 is not None:
                    yield (v1, None)
                for x1 in i1:
                    yield (x1, None)
                break

        k1 = key(v1)
        k2 = key(v2)
        if k1 == k2:
            yield (v1, v2)
            v1 = None
            v2 = None
        elif k1 < k2:
            yield (v1, None)
            v1 = None
        else:
            yield (None, v2)
            v2 = None

def is_unfiltered(filter_field):
    return (filter_field is None or (filter_field is not None and (len(filter_field) == 1 and filter_field[0] == "PASS")) or filter_field == [])

def populate_repeat_info(record, reference_pyfasta, length):
    post_poly_count = 0
    post_poly_base = None
    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    lastBase = None
    gap = min(30, len(reference_pyfasta[chrom])-pos-1)
    sequence = reference_pyfasta[chrom][(pos+1):(pos+gap+1)].upper()
    for base in range(0, gap, length):
        if lastBase is None:
            post_poly_count = 1
            post_poly_base = sequence[base:base+length]
            lastBase = post_poly_base
        elif lastBase is not None:
            if lastBase == sequence[base:base+length]:
                post_poly_count += 1
            else:
                break
        else:
            break
    return post_poly_count, post_poly_base

def variant_pairs_to_df(chrom_pos, chrom_cov, bc_cov, variant_pair_iter, offset, gene_finder, regions_of_interest, validation_bam, reference_pyfasta, sample_bam):
    '''Generate a row of summary information about each match between a called allele and ground-truth allele'''
    record_type = [
            ('variant_type', 'S1'),
            ('variant_length', np.int16),
            ('chrom', object),
            ('pos', np.int32),
            ('gene', object),
            ('ref', object),
            ('alt1', object),
            ('alt2', object),
            ('allele', object),
            ('ti', np.uint8),
            ('tv', np.uint8),
            ('cov', np.int16),
            ('bc_cov', np.int16),
            ('in_obs', np.uint8),
            ('qual', np.int16),
            ('dp', np.int16),
            ('obs_phased', np.uint8),
            ('obs_genotype1', np.uint8),
            ('obs_genotype2', np.uint8),
            ('bcs1', np.uint16),
            ('bcs2', np.uint16),
            ('phase_set', np.int64),
            ('pq', np.int16),
            ('jq', np.int16),
            ('dps', np.int16),
            ('in_gt', np.uint8),
            ('in_comparison', np.uint8),
            ('gt_phased', np.uint8),
            ('gt_genotype1', np.uint8),
            ('gt_genotype2', np.uint8),
            ('gt_phase_set', np.int64),
            ('parity', np.uint8),
            ('can_test_phase', np.uint8),
            ('ro',np.int16),
            ('ao1',np.int16),
            ('ao2',np.int16),
            ('FILTER',np.uint8),
            ('posthpc',np.uint8),
            ('posthpb',np.object),
            ('postdnc',np.uint8),
            ('postdnb',np.object),
            ('posttnc',np.uint8),
            ('posttnb',np.object),
            ('mumap_alt1', np.float32),
            ('mumap_alt2', np.float32),
            ('mumap_ref', np.float32),
            ('validation_cov',np.int32),
            ('validation_ao',np.int32),
            ('molecule_differences',np.float32),
            ('ASN',np.int32),
            ('ASP',np.int32),
            ('haplocalled', np.uint8),
            ('gc', np.float32),
            ('periodicity', np.float32),
            ('rescued', np.int32),
            ('not_rescued', np.int32)]

    def get_data(data, field, default):
        v = getattr(data, field, default)
        if v is None:
            return default

        return v

    regions = {}
    if regions_of_interest is not None:
        for (name,region) in regions_of_interest.iteritems():
            if os.path.exists(region):
                with open(region) as r:
                    regions[name] = tk_io.get_target_regions(r)
                record_type.append((str(name), np.uint8))

    results = []
    # Construct recarray
    bsize = 100
    record_block = np.recarray((bsize,), dtype=record_type)
    counter = 0

    for (gt_variant, obs_variant, compare_variant) in variant_pair_iter:
        row = record_block[counter]
        row.ti = 0
        row.tv = 0
        row.in_comparison = 0
        if compare_variant is not None:
            row.in_comparison = 1

        if obs_variant is not None:
            ((chrom, pos, ref, allele), obs_record) = obs_variant
            rec = obs_record
            row.ro = tk_io.get_record_ref_allele_count(rec) if tk_io.get_record_ref_allele_count(rec) is not None else -1
            # TODO understand why MMD is sometimes returned as a list.
            mmd = rec.INFO['MMD'] if 'MMD' in obs_record.INFO else -1.0
            if type(mmd) is list:
                if len(mmd) == 0:
                    mmd = 0
                else:
                    mmd = mmd[0]
            rescued = -1
            if "RESCUED" in obs_record.INFO:
                rescued = 0
                for x in obs_record.INFO["RESCUED"]:
                    try:
                        rescued += int(x)
                    except:
                        pass
            row.rescued = rescued
            row.not_rescued = int(float(obs_record.INFO["NOT_RESCUED"][0])) if "NOT_RESCUED" in obs_record.INFO else -1
            row.molecule_differences = mmd
            alt_allele_count = tk_io.get_record_alt_allele_counts(rec)
            row.ao1 = alt_allele_count[0] if alt_allele_count is not None else -1
            row.ao2 = alt_allele_count[1] if alt_allele_count is not None and len(alt_allele_count) > 1 else -1
            row.FILTER = 1 if is_unfiltered(tk_io.get_record_filters(obs_record)) else 0
            row.ASN = int(obs_record.INFO["AON"][0]) if "AON" in obs_record.INFO else -2
            row.ASP = int(obs_record.INFO["AOP"][0]) if "AOP" in obs_record.INFO else -2
            row.haplocalled = rec.INFO.get('HAPLOCALLED') if 'HAPLOCALLED' in rec.INFO else 0
        elif gt_variant is not None:
            ((chrom, pos, ref, allele), gt_record) = gt_variant
            alleles = [allele]
            (counts, mean_mapqs, bc_qual_string, molecule_differences, diffs, _) = tk_bam.get_allele_read_info(chrom, pos + 1, ref, alleles, 30, -1, 0, 0, sample_bam, reference_pyfasta)
            row.ro = counts[0]
            row.ao1 = counts[1]
            row.ao2 = -1
            row.FILTER = 0
            row.haplocalled = 0
            row.posthpc = 0
            row.posthpb = 'N'
            row.mumap_alt1 = mean_mapqs[1]
            row.mumap_ref = mean_mapqs[0]
            neg = 0
            pos = 0
            for x in diffs[1]:
                if x < 0:
                    neg += 1
                elif x > 0:
                    pos += 1
            row.ASN = neg
            row.ASP = pos
        else:
            ((chrom, pos, ref, allele), comp_record) = compare_variant
            rec = comp_record

        if (allele, ref) in [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]:
            row.ti = 1

        if (ref in ['A', 'G'] and allele in ['C', 'T']) or (ref in ['C', 'T'] and allele in ['A', 'G']):
            row.tv = 1


        if gt_variant is not None:
            ((chrom, pos, ref, allele), gt_record) = gt_variant
            rec = gt_record
        post_hp_c, post_hp_b = populate_repeat_info(rec, reference_pyfasta, 1)
        post_dn_c, post_dn_b = populate_repeat_info(rec, reference_pyfasta, 2)
        post_tn_c, post_tn_b = populate_repeat_info(rec, reference_pyfasta, 3)

        row.posthpc = post_hp_c
        row.posthpb = post_hp_b
        row.postdnc = post_dn_c
        row.postdnb = post_dn_b
        row.posttnc = post_tn_c
        row.posttnb = post_tn_b

        gc = 0.0
        tots = 0.0
        for i in range(max(0, pos - 15), min(len(reference_pyfasta[chrom])-1, pos + 15)):
            tots += 1
            if reference_pyfasta[chrom][i].upper() == 'G' or reference_pyfasta[chrom][i].upper() == 'C':
                gc += 1
        if tots >= 0:
            row.gc = gc / tots
        else:
            row.gc = -1

        periodicity = -1.0
        for period in range(1,4):
            seq1 = reference_pyfasta[chrom][max(0, pos - 15): min(len(reference_pyfasta[chrom])-1, pos + 15)].upper()
            seq2 = reference_pyfasta[chrom][max(0, pos - 15 + period): min(len(reference_pyfasta[chrom])-1, pos + 15 + period)].upper()
            match = 0.0
            tots = 0.0
            for (a,b) in zip(seq1, seq2):
                if a == b:
                    match += 1
                tots += 1
            if tots > 0:
                repitivity = match/tots
                if repitivity > periodicity:
                    periodicity = repitivity
        row.periodicity = periodicity

        row.variant_type = tk_io.get_var_type(ref, allele)
        row.variant_length = tk_io.get_allele_length(ref, allele)

        row.chrom = chrom
        alts = [allele]

        if validation_bam is not None: # and chrom in validation_bam.references:
            if validation_bam.references[0][0:3] != "chr":
                chrom = chrom[3:]
            counts, mapq_means, bc_qual_strings, _ ,_,_= tk_bam.get_allele_read_info(chrom,pos+1, ref, alts, 0, 0, 0, 0, validation_bam, reference_pyfasta, match = 1, mismatch = -3, gap_open = -1, gap_extend = -4)
            row.validation_cov = sum(counts)
            row.validation_ao = counts[1]
        else:
            row.validation_cov = -1
            row.validation_ao = -1

        for (name,region) in regions.iteritems():
            if chrom in region:
                if region[chrom].contains_point(pos):
                    row[str(name)] = 1

        row.pos = pos
        row.ref = ref
        alt_alleles = tk_io.get_record_alt_alleles(rec)
        row.alt1 = str(alt_alleles[0])
        row.alt2 = str(alt_alleles[1]) if len(rec.ALT) > 1 else 'N'

        gene = gene_finder.get_gene(chrom, rec.POS)
        if gene:
            row.gene = gene
        else:
            row.gene = ''

        row.allele = allele

        if chrom_cov is not None:
            # NOTE -- dtype of query to searchsorted must match the input array
            # or you can get a ~10^5 fold slowdown.
            idx = np.searchsorted(chrom_pos, np.int32(pos))

            if idx >= len(chrom_pos):
                idx = len(chrom_pos) - 2
            if (chrom_pos[idx] != pos):
                if abs(chrom_pos[idx+1] - pos) < abs(chrom_pos[idx] - pos):
                    idx = idx + 1
            row.cov = chrom_cov[idx]
        else:
            row.cov = -1


        if bc_cov is not None:
            # NOTE -- dtype of query to searchsorted must match the input array
            # or you can get a ~10^5 fold slowdown.
            idx = np.searchsorted(chrom_pos, np.int32(pos))

            if idx >= len(chrom_pos):
                idx = len(chrom_pos) - 2
            if (chrom_pos[idx] != pos):
                if abs(chrom_pos[idx+1] - pos) < abs(chrom_pos[idx] - pos):
                    idx = idx + 1
            row.bc_cov = bc_cov[idx]
        else:
            row.bc_cov = -1


        if obs_variant is not None:
            row.in_obs = 1
            try:
                # some 1k genomes vcfs have QUAL=nan
                row.qual = obs_record.QUAL
            except (ValueError, TypeError):
                row.qual = -1

            if obs_record.INFO.has_key('DP'):
                row.dp = obs_record.INFO['DP']
            else:
                row.dp = -1

            sample = obs_record.samples[0]
            row.obs_phased = int(sample.phased)

            # sometimes when a multisample vcf is made into a single sample vcf, some records get genotype ./. if they werent in the
            # sample remaining and pyvcf chokes on these sometimes. :/
            try:
                if len(sample.gt_alleles) > 1:
                    row.obs_genotype1 = int(sample.gt_alleles[0])
                    row.obs_genotype2 = int(sample.gt_alleles[1])
                else:
                    row.obs_genotype1 = int(sample.gt_alleles[0])
                    row.obs_genotype2 = int(sample.gt_alleles[0])
            except:
                pass

            def get_bcs(bc_string):
                bcs = bc_string.split(";")
                if bcs == ['']:
                    bcs = []
                return bcs

            alleles_in_order = [row.obs_genotype1, row.obs_genotype2]
            alleles_in_order.sort()

            bcs1 = []
            bcs2 = []
            bcs_data = get_data(sample.data, "BX", None)

            if bcs_data is not None:
                if len(bcs_data) > alleles_in_order[0]:
                    bcs1 = get_bcs(bcs_data[alleles_in_order[0]])
                if len(bcs_data) > alleles_in_order[1]:
                    bcs2 = get_bcs(bcs_data[alleles_in_order[1]])

            row.bcs1 = len(bcs1)
            row.bcs2 = len(bcs2)

            mumap_alt = tk_io.get_record_mean_mapq_alt(obs_record)
            row.mumap_alt1 = mumap_alt[0]
            row.mumap_ref = tk_io.get_record_mean_mapq_ref(obs_record)
            if len(mumap_alt) > 1:
                row.mumap_alt2 = mumap_alt[1]
            else:
                row.mumap_alt2 = -1

            row.phase_set = get_data(sample.data, "PS", -1)
            row.pq = get_data(sample.data, "PQ", -1)
            row.jq = get_data(sample.data, "JQ", -1)
            row.dps = get_data(sample.data, "DP", -1)
        if gt_variant is not None:
            row.in_gt = 1 if is_unfiltered(tk_io.get_record_filters(gt_record)) else 0
            sample = gt_record.samples[0]
            row.gt_phased = int(sample.phased)
            if len(sample.gt_alleles) > 1:
                row.gt_genotype1 = int(sample.gt_alleles[0])
                row.gt_genotype2 = int(sample.gt_alleles[1])
            else:
                row.gt_genotype1 = int(sample.gt_alleles[0])
                row.gt_genotype2 = int(sample.gt_alleles[0])

            row.gt_phase_set = get_data(sample.data, "PS", -1)

        # What constitutes a 'testable' variant:
        # heterozygous in observed and gt
        # phased in obs and gt
        # in a valid phase set
        # not filtered in obs or gt
        if obs_variant is not None and \
            gt_variant is not None and \
            row.gt_genotype1 != row.gt_genotype2 and \
            row.gt_phased and row.obs_phased and row.phase_set > 0 and \
            row.FILTER > 0 and row.in_gt > 0:
            obs_alleles = [tk_io.get_record_ref(obs_record)] + tk_io.get_record_alt_alleles(obs_record)
            gt_alleles = [tk_io.get_record_ref(gt_record)] + tk_io.get_record_alt_alleles(gt_record)
            if row.obs_genotype1 != 0:
                (obs_pos_g1, obs_new_ref_g1, obs_allele_genotype1) = tk_io.normalize_variant(tk_io.get_record_pos(obs_record), obs_alleles[0], obs_alleles[row.obs_genotype1], reference_pyfasta[chrom])
            else:
                obs_allele_genotype1 = 0
                obs_pos_g1 = pos
                obs_new_ref_g1 = ref
            if row.obs_genotype2 != 0:
                (obs_pos_g2, obs_new_ref_g2, obs_allele_genotype2) = tk_io.normalize_variant(tk_io.get_record_pos(obs_record), obs_alleles[0], obs_alleles[row.obs_genotype2], reference_pyfasta[chrom])
            else:
                obs_allele_genotype2 = 0
                obs_pos_g2 = pos
                obs_new_ref_g2 = ref
            if row.gt_genotype1 != 0:
                (gt_pos_g1, gt_new_ref_g1, gt_allele_genotype1) = tk_io.normalize_variant(tk_io.get_record_pos(gt_record), gt_alleles[0], gt_alleles[row.gt_genotype1], reference_pyfasta[chrom])
            else:
                gt_pos_g1 = pos
                gt_allele_genotype1 = 0
                gt_new_ref_g1 = ref
            if row.gt_genotype2 != 0:
                (gt_pos_g2, gt_new_ref_g2, gt_allele_genotype2) = tk_io.normalize_variant(tk_io.get_record_pos(gt_record), gt_alleles[0], gt_alleles[row.gt_genotype2], reference_pyfasta[chrom])
            else:
                gt_pos_g2 = pos
                gt_allele_genotype2 = 0
                gt_new_ref_g2 = ref
            if obs_allele_genotype1 == allele:
                row.parity = ((obs_allele_genotype1 == gt_allele_genotype1) & (obs_pos_g1 == gt_pos_g1) & (obs_new_ref_g1 == gt_new_ref_g1))
            elif obs_allele_genotype2 == allele:
                row.parity = ((obs_allele_genotype2 == gt_allele_genotype2) & (obs_pos_g2 == gt_pos_g2) & (obs_new_ref_g2 == gt_new_ref_g2))
            else:
                assert(False)
            row.can_test_phase = True
        else:
            row.parity = False
            row.can_test_phase = False

        counter += 1
        if counter == bsize:
            counter = 0
            results.append(record_block)
            record_block = np.recarray((bsize,), dtype=record_type)

    results.append(record_block[:counter])
    df = tk_pd.concat([tk_pd.DataFrame(blk) for blk in results], ignore_index=True)
    return df
