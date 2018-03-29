# imports
import tenkit.log_subprocess
import tenkit.bio_io as tk_io
import tenkit.tabix as tk_tabix
import itertools
import pandas
from multiprocessing import Pool
import traceback
import sys
from collections import namedtuple

"""
Some details on data structures:
* VARIANTS are represented as PyVCF Records, usually wrapped in a VariantFileReader instance
* FRAGMENTS are represented as Pandas DataFrames with the column names below
* LOCI are represented as named tuples of (chrom, start, end)
* PHASE BLOCKS are represented as named tuples of (id, locus)
"""
frag_col_names = ['chrom', 'frag_start', 'frag_end', 'phase_set', 'ps_start', 'ps_end', 'bc', 'h0', 'h1', 'hmix', 'reads', 'molecule_id']
Locus = namedtuple('Locus', ['chrom', 'start', 'end'])
PhaseBlock = namedtuple('PhaseBlock', ['id', 'locus'])

def un_bgzip(vcfgz):
    out = vcfgz.rstrip('.gz')
    tenkit.log_subprocess.check_call("bgzip -cd {0} > {1}".format(vcfgz, out), shell=True)

def get_variant_iterator(vfr, locus):
    """
    Wrapper to get around the fact that tk_io.get_variant_iterator() takes a locus string instead of a tuple.
    """
    locus_string = tk_io.create_locus_info(locus.chrom, locus.start, locus.end)
    return tk_io.get_variant_iterator_pos(vfr, None, locus_string)

def get_phase_block_for_pos(frags, chrom, pos):
    pbs = frags[['chrom', 'phase_set', 'ps_start', 'ps_end']].drop_duplicates()
    sliced = pbs[(pbs.chrom == chrom) & (pbs.ps_start <= pos) & (pbs.ps_end > pos)]
    if len(sliced) > 1:
        # This should never happen
        raise Exception("More than one phase set intersecting position {}".format(pos))
    elif len(sliced) == 0:
        # No phase block overlaps this position
        return None
    else:
        block = sliced.iloc[0]
        # NOTE: the block end is inclusive, but we want to keep the convention of all locus ends being exclusive, so add 1 here
        return PhaseBlock(block.phase_set, Locus(block.chrom, block.ps_start, block.ps_end + 1))

def lockstep_variant_iterator(vfr_left, vfr_right, shared_locus):
    """
    Traverse two copies of the same variants in lockstep, making sure
    we never get out of sync.
    """
    iter_left = get_variant_iterator(vfr_left, shared_locus)
    iter_right = get_variant_iterator(vfr_right, shared_locus)

    for (var_left, var_right) in zip(iter_left, iter_right):
        # keep these 1-indexed since they're only used for the error message
        pos_left = tk_io.get_record_pos(var_left)
        pos_right = tk_io.get_record_pos(var_right)
        if pos_left != pos_right:
            raise Exception("Variant positions are out of sync: {0}:{1}, {0}:{2}".format(shared_locus.chrom, pos_left, pos_right))
        yield (var_left, var_right)

def calc_match_scores(vfr_left, vfr_right, shared_locus):
    num_total_variants = 0
    num_considered_variants = 0
    num_matching_forward = 0
    num_matching_reverse = 0

    for (var_left, var_right) in lockstep_variant_iterator(vfr_left, vfr_right, shared_locus):
        num_total_variants += 1
        (gt_left, is_phased_left) = tk_io.get_record_genotype_phased(var_left)
        (gt_right, is_phased_right) = tk_io.get_record_genotype_phased(var_right)

        is_het_left = len(gt_left) > 1 and gt_left[0] != gt_left[1]
        is_het_right = len(gt_right) > 1 and gt_right[0] != gt_right[1]

        if is_phased_left and is_phased_right and is_het_left and is_het_right:
            gt_right_flipped = [gt_right[1], gt_right[0]]
            num_considered_variants += 1
            if gt_left == gt_right:
                num_matching_forward += 1
            elif gt_left == gt_right_flipped:
                num_matching_reverse += 1

    if num_considered_variants == 0:
        return (0, 0)
    else:
        return (1.0 * num_matching_forward / num_considered_variants, 1.0 * num_matching_reverse / num_considered_variants)

def combine_phasing(record_left, record_right, new_phase_set, flip_right=False):
    """
    1. if the variant is phased in both lists and they agree, keep it as-is
    2. if the variant is phased in both lists but they disagree, unphase it
    3. if the variant is phased in one list only, keep the phased version
    4. if the variant is phased in neither list, keep it as-is
    """
    if flip_right:
        flip_phasing(record_right)

    (gt_left, is_phased_left) = tk_io.get_record_genotype_phased(record_left)
    (gt_right, is_phased_right) = tk_io.get_record_genotype_phased(record_right)

    # start with left-hand copy by convention. the two copies should differ only in phasing.
    record_out = record_left
    if is_phased_left and is_phased_right and (gt_left != gt_right):
        unphase(record_out)
    elif (not is_phased_left) and is_phased_right:
        record_out = record_right

    tk_io.set_record_phase_set(record_out, new_phase_set)
    return record_out

def adjust_phasing(record, new_phase_set, flip=False):
    tk_io.set_record_phase_set(record, new_phase_set)
    if flip:
        flip_phasing(record)

def flip_phasing(record):
    (gt, phased) = tk_io.get_record_genotype_phased(record)
    if phased and len(gt) > 1:
        gt_flipped = [gt[1], gt[0]]
        tk_io.set_record_genotype_phased(record, gt_flipped, phased)

def unphase(record):
    (gt, phased) = tk_io.get_record_genotype_phased(record)
    tk_io.set_record_genotype_phased(record, gt, False)
    tk_io.set_record_phase_qual(record, 0)
    tk_io.set_record_junction_qual(record, 0)

def stitch_variant_iterator(vfr_left, vfr_right, chunk_locus_left, chunk_locus_right,
                    block_locus_left, block_locus_right, phase_set_left, flip_right=False):
    # first define some intervals.
    # we're assuming that chunk_locus_left / chunk_locus_right are the entire spans of the left / right chunks respectively,
    # whereas block_locus_left / block_locus_right are the spans of the overlapping blocks that are going to be merged together.

    # locus0 is from the start of the left-hand chunk to the start of the overlap. the records in here will not be changed.
    locus0 = Locus(chunk_locus_left.chrom, chunk_locus_left.start, chunk_locus_right.start)
    # locus1 spans the overlap. the records in here are duplicated and will need to be merged.
    locus1 = Locus(chunk_locus_left.chrom, chunk_locus_right.start, chunk_locus_left.end)
    # locus2 is from the end of the overlap to the end of the right-hand phase block.
    # in this region, we will need to overwrite the phase set and maybe flip the haplotypes
    locus2 = Locus(chunk_locus_left.chrom, chunk_locus_left.end, block_locus_right.end)
    # locus3 is from the end of the right-hand phase block to the end of the right-hand chunk. the records in here will not be changed.
    locus3 = Locus(chunk_locus_left.chrom, block_locus_right.end, chunk_locus_right.end)

    # 0. pass through the left-hand records
    for record_out in get_variant_iterator(vfr_left, locus0):
        yield record_out

    # 1. merge inside the overlap
    for (record_left, record_right) in lockstep_variant_iterator(vfr_left, vfr_right, locus1):
        record_out = combine_phasing(record_left, record_right, phase_set_left, flip_right=flip_right)
        yield record_out

    # 2. tweak the right hand records
    for record_out in get_variant_iterator(vfr_right, locus2):
        adjust_phasing(record_out, phase_set_left, flip=flip_right)
        yield record_out

    # 3. pass through the right-hand records
    for record_out in get_variant_iterator(vfr_right, locus3):
        yield record_out

def split_variant_iterator(vfr_left, vfr_right, new_locus_left, new_locus_right):
    # assert no overlap
    assert(new_locus_left.end <= new_locus_right.start)
    for record_out in get_variant_iterator(vfr_left, new_locus_left):
        yield record_out
    first_phase_set_right = None
    for record_out in get_variant_iterator(vfr_right, new_locus_right):
        if first_phase_set_right is None:
            first_phase_set_right = tk_io.get_record_pos(record_out) - 1
        # if we see a real phase set that's less than the new one,
        # then the block was truncated and should be updated
        current_ps = tk_io.get_record_phase_set(record_out)
        if current_ps > 0 and current_ps < first_phase_set_right:
            adjust_phasing(record_out, first_phase_set_right, flip=False)
        yield record_out

def join_loci(locus0, locus1):
    return Locus(locus0.chrom, locus0.start, locus1.end)

def flip_frag_phasing(frags):
    frags[['h1', 'h0']] = frags[['h0', 'h1']]

def unify_frag_rows(left_frags, right_frags, new_phase_set, new_start, new_end):
    """
    Given two dataframes of barcoded fragments from two overlapping phase blocks,
    take the 'union' of the rows to produce a unified set of fragments,
    with updated phase set info.
    """
    def calc_overlap(start0, end0, start1, end1):
        start = max(start0, start1)
        end = min(end0, end1)
        return max(end - start, 0)

    joined = pandas.concat([left_frags, right_frags]).reset_index()

    # calculate how much each fragment overlaps its phase block
    joined['overlaps'] = map(calc_overlap, joined.frag_start, joined.frag_end, joined.ps_start, joined.ps_end)

    # complicated part: sort by overlap size, then group by BC/frag start/frag end.
    # in each group, we want to keep only one row (the one with the biggest overlap),
    # which is always the first row since we sorted. then reset_index to flatten.
    fixed = joined.sort('overlaps', ascending=False).groupby(['bc', 'frag_start', 'frag_end']).first().reset_index()

    # finally, rewrite the phase set info and fix column order
    fixed[['phase_set', 'ps_start', 'ps_end']] = (new_phase_set, new_start, new_end)
    fixed = fixed[frag_col_names]
    return fixed

def stitch_frags(left_frags, right_frags, left_phase_set, right_phase_set,
                 left_block_start, right_block_end, flip_right=False):
    # 1. add everything leading up to the left-hand block
    output_head = left_frags[left_frags.phase_set < left_phase_set].reset_index()

    # 2. aggregate the fragments in the overlapping phase blocks. any other phase blocks in the overlap region
    # (i.e. trailing stuff on the left and leading stuff on the right) is just ignored.
    left_frags_in_block = left_frags[left_frags.phase_set == left_phase_set].reset_index()
    right_frags_in_block = right_frags[right_frags.phase_set == right_phase_set].reset_index()

    if flip_right:
        flip_frag_phasing(right_frags_in_block)

    output_body = unify_frag_rows(left_frags_in_block, right_frags_in_block, left_phase_set, left_block_start, right_block_end)

    # 3. add everything that comes after the right-hand block
    output_tail = right_frags[right_frags.phase_set > right_phase_set].reset_index()

    # 4. concatenate
    return pandas.concat([output_head, output_body, output_tail]).reset_index()

def split_frags(left_frags, right_frags, left_end_new, right_start_new):
    """
    Given two overlapping fragment sets, merge them by redefining the start / end as the specified points
    """
    # assert no overlap
    assert(left_end_new <= right_start_new)

    left_head = left_frags[(left_frags.ps_end < left_end_new) & (left_frags.frag_start < left_end_new)].reset_index()
    left_body = left_frags[(left_frags.ps_start <= left_end_new) & (left_frags.ps_end >= left_end_new) & (left_frags.frag_start <= left_end_new)].reset_index()
    right_body = right_frags[(right_frags.ps_start <= right_start_new) & (right_frags.ps_end >= right_start_new) & (right_frags.frag_end >= right_start_new)].reset_index()
    right_tail = right_frags[(right_frags.ps_start > right_start_new) & (right_frags.frag_end > right_start_new)].reset_index()

    left_body.ps_end = left_end_new
    right_body.ps_start = right_start_new
    right_body.phase_set = right_start_new
    # TODO should we update the haplotype probabilties?

    if len([x for x in [left_head, left_body, right_body, right_tail] if not x.empty]) == 0:
        return pandas.DataFrame({'chrom':[], 'frag_start':[], 'frag_end':[], 'phase_set':[], 'ps_start':[], 'ps_end':[], 'reads':[]})
    return pandas.concat([x for x in [left_head, left_body, right_body, right_tail] if not x.empty]).reset_index()

def stitch(vfr_left, vfr_right, vfw, frags_left, frags_right, frags_file_out, loc_left, loc_right, min_match_score = 0.95, min_overlap = 2e4):
    # calculate overlap between chunks
    chunk_overlap_locus = Locus(loc_left.chrom, loc_right.start, loc_left.end)

    # get blocks on borders
    bordering_block_left = get_phase_block_for_pos(frags_left, chunk_overlap_locus.chrom, chunk_overlap_locus.start)
    bordering_block_right = get_phase_block_for_pos(frags_right, chunk_overlap_locus.chrom, chunk_overlap_locus.end)

    record_iterator_out = []
    fragment_df_out = None
    stitched = False

    if bordering_block_left is not None and bordering_block_right is not None:
        # calculate overlap between the blocks of interest
        # TODO less tuple indexing?
        block_overlap_locus = Locus(bordering_block_left.locus.chrom, bordering_block_right.locus.start, bordering_block_left.locus.end)
        block_overlap_length = block_overlap_locus.end - block_overlap_locus.start
        if block_overlap_length >= min_overlap:
            (forward_score, reverse_score) = calc_match_scores(vfr_left, vfr_right, block_overlap_locus)
            if max(forward_score, reverse_score) >= min_match_score:
                # actually stitch
                flip_right = reverse_score > forward_score
                record_iterator_out = stitch_variant_iterator(vfr_left, vfr_right, loc_left, loc_right, bordering_block_left.locus, bordering_block_right.locus,
                    bordering_block_left.id, flip_right=flip_right)
                fragment_df_out = stitch_frags(frags_left, frags_right, bordering_block_left.id, bordering_block_right.id,
                    bordering_block_left.locus.start, bordering_block_right.locus.end, flip_right=flip_right)
                stitched = True

    if not stitched:
        # cut down the middle
        cut_point = chunk_overlap_locus.start + (chunk_overlap_locus.end - chunk_overlap_locus.start) / 2.0
        loc_left_adjusted = Locus(loc_left.chrom, loc_left.start, cut_point)
        loc_right_adjusted = Locus(loc_right.chrom, cut_point, loc_right.end)

        # see if we need to adjust the phase blocks
        block_to_truncate_left = get_phase_block_for_pos(frags_left, chunk_overlap_locus.chrom, cut_point)
        block_to_truncate_right = get_phase_block_for_pos(frags_right, chunk_overlap_locus.chrom, cut_point)

        if block_to_truncate_left is not None:
            search_locus_left = Locus(block_to_truncate_left.locus.chrom, block_to_truncate_left.locus.start, cut_point)
            left_end_new = get_closest_variant_pos([record for record in get_variant_iterator(vfr_left, search_locus_left)], cut_point, -1)
            if left_end_new is None:
                left_end_new = search_locus_left.start # no variants - stop at beginning of search
            else:
                left_end_new += 1 # since we want to include the last variant, but the locus is half-closed, we must add 1 here
            loc_left_adjusted = Locus(loc_left_adjusted.chrom, loc_left_adjusted.start, left_end_new)

        if block_to_truncate_right is not None:
            search_locus_right = Locus(block_to_truncate_right.locus.chrom, cut_point, block_to_truncate_right.locus.end)
            right_start_new = get_closest_variant_pos([record for record in get_variant_iterator(vfr_right, search_locus_right)], cut_point, 1)
            if right_start_new is None:
                right_start_new = search_locus_right.end # no variants - stop at end of search
            loc_right_adjusted = Locus(loc_right_adjusted.chrom, right_start_new, loc_right_adjusted.end)

        record_iterator_out = split_variant_iterator(vfr_left, vfr_right, loc_left_adjusted, loc_right_adjusted)
        fragment_df_out = split_frags(frags_left, frags_right, loc_left_adjusted.end, loc_right_adjusted.start)

    # make sure types are correct
    int_cols = ['frag_start', 'frag_end', 'phase_set', 'ps_start', 'ps_end', 'reads']
    fragment_df_out[int_cols] = fragment_df_out[int_cols].astype('int64')

    # make sure column + row order is correct
    if not fragment_df_out.empty:
        fragment_df_out = fragment_df_out[frag_col_names].sort(['chrom', 'frag_start', 'frag_end'])

    # write it out
    for record in record_iterator_out:
        vfw.write_record(record)

    fragment_df_out.to_csv(frags_file_out, sep='\t', header=False, index=False, mode='w')

def get_closest_variant_pos(variants, target_pos, direction, het_only=False):
    """
    Get closest variant to target, looking in specified direction (-1 = before target, 1 = after target)
    """
    # too lazy to implement binary search.
    for variant in variants[::direction]:
        pos = tk_io.get_record_pos(variant) - 1
        right_direction = (pos <= target_pos) if direction < 0 else (pos >= target_pos)
        if (gt_is_het(variant) or not het_only) and right_direction:
            return pos
    # edge case - no variants in that direction.
    return None

def gt_is_het(record):
    (gt, phased) = tk_io.get_record_genotype_phased(record)
    return len(gt) > 1 and gt[0] != gt[1]

### HIGH-LEVEL JOIN METHODS

def join_simple(vcf_left, vcf_right, vcf_out, ff_left, ff_right, ff_out, loc_left, loc_right):
    frags_left = pandas.read_table(ff_left, header=None, names=frag_col_names)
    frags_right = pandas.read_table(ff_right, header=None, names=frag_col_names)

    # make sure chrom is formatted as string categorical
    frags_left.chrom = frags_left.chrom.astype('str').astype('category')
    frags_right.chrom = frags_right.chrom.astype('str').astype('category')

    vfr_left = tk_io.VariantFileReader(vcf_left + ".gz")
    vfr_right = tk_io.VariantFileReader(vcf_right + ".gz")
    with open(vcf_left + ".gz", 'r') as template_in, open(vcf_out, 'w') as results_out:
        vfw = tk_io.VariantFileWriter(results_out, template_file=template_in)
        stitch(vfr_left, vfr_right, vfw, frags_left, frags_right, ff_out, loc_left, loc_right)
    tk_tabix.sort_unique_tabix_vcf(vcf_out)

def join_single_chrom(frags, vcfs, loci):
    merged_vcfs = map(lambda v: v.rstrip(".vcf") + ".out.vcf", vcfs)
    merged_frags = map(lambda f: f.rstrip(".tsv") + ".out.tsv", frags)
    merged_loci = map(lambda l: join_loci(loci[0], l), loci)

    # the first .out files are just copies, to make the iteration cleaner
    tenkit.log_subprocess.check_call("cp {} {}".format(vcfs[0] + ".gz", merged_vcfs[0] + ".gz"), shell=True)
    tenkit.log_subprocess.check_call("cp {} {}".format(vcfs[0] + ".gz.tbi", merged_vcfs[0] + ".gz.tbi"), shell=True)
    tenkit.log_subprocess.check_call("cp {} {}".format(frags[0], merged_frags[0]), shell=True)

    # merge right
    n = len(frags)
    if n == 1:
        # need an un-gzipped verzion for combine_vcfs to work later on.
        un_bgzip(merged_vcfs[0] + ".gz")
    else:
        for i in range(1, n):
            join_simple(merged_vcfs[i-1], vcfs[i], merged_vcfs[i], merged_frags[i-1], frags[i], merged_frags[i], merged_loci[i-1], loci[i])

    final_vcf = merged_vcfs[n-1]
    final_frags = merged_frags[n-1]

    return (final_vcf, final_frags)

def join_single_chrom_helper(arg_tuple):
    try:
        (loci, frags, vcfs) = arg_tuple
        (out_vcf, out_frags) = join_single_chrom(frags, vcfs, loci)
        return (out_vcf, out_frags)
    except:
        # raise an exception that encapsulates the traceback in pickle-able form,
        # so that multiprocessing can bubble it up
        raise Exception(get_stacktrace_string())

def get_stacktrace_string():
    # copied from Martian. we need this here so that we can peer inside the multiprocess
    # TODO - move this functionality into tenkit?
    etype, evalue, tb = sys.exc_info()
    stacktrace = ["Traceback (most recent call last):"]
    while tb:
        frame = tb.tb_frame
        filename, lineno, name, text = traceback.extract_tb(tb, limit=1)[0]
        stacktrace.append("  File '%s', line %d, in %s" % (filename ,lineno, name))
        if text:
            stacktrace.append("    %s" % text.strip())
        for key, value in frame.f_locals.items():
            try:
                stacktrace.append("        %s = %s" % (key, str(value)))
            except:
                pass
        tb = tb.tb_next
    stacktrace += [line.strip() for line in traceback.format_exception_only(etype, evalue)]
    return "\n".join(stacktrace)

def multi_join_parallel(all_frags, all_vcfs, all_loci, n_procs=20):
    # TODO - figure out how handle logging inside the process pool
    pool = Pool(n_procs)
    # convert loci to named representation
    all_loci = map(lambda x: Locus(x[0], x[1], x[2]), all_loci)
    # group by chrom
    grouped = itertools.groupby(zip(all_loci, all_frags, all_vcfs), lambda x: x[0].chrom)
    arg_tuples = [map(list, zip(*group)) for chrom, group in grouped]
    try:
        results = pool.map(join_single_chrom_helper, arg_tuples)
    finally:
        pool.close()
        pool.join()
    (stitched_chrom_vcfs, stitched_chrom_frags) = map(list, zip(*results))
    return (stitched_chrom_vcfs, stitched_chrom_frags)
