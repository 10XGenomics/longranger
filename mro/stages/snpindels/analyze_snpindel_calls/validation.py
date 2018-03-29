import tenkit.hdf5 as tk_h5
import tenkit.bio_io as tk_io
import pandas

def get_phase_blocks_from_vcf(vcf, filter_trivial=True, locus=None):
    if locus is not None:
        (loc_chrom, loc_start, loc_end) = locus
    else:
        (loc_chrom, loc_start, loc_end) = (None, None, None)

    records = tk_io.VariantFileReader(vcf).record_getter(fetch_chrom=loc_chrom, fetch_start=loc_start, fetch_end=loc_end)

    phase_set, chrom, start, end, length = [], [], [], [], []
    curr_set, curr_chrom, curr_start, curr_end = None, None, None, None

    def get_record_data(record):
        record_set = tk_io.get_record_phase_set(record)
        record_chrom = tk_io.get_record_chrom(record)
        record_pos = tk_io.get_record_pos(record) - 1
        return (record_set, record_chrom, record_pos)

    def write_curr():
        phase_set.append(curr_set)
        chrom.append(curr_chrom)
        start.append(curr_start)
        end.append(curr_end)
        length.append(curr_end - curr_start)

    for record in records:
        (next_set, next_chrom, next_pos) = get_record_data(record)
        if curr_set is None and next_set is not None:
            # initialize
            curr_set, curr_chrom, curr_start, curr_end = next_set, next_chrom, next_pos, next_pos
        elif next_set is not None and (curr_set != next_set or curr_chrom != next_chrom):
            # new set. write curr, then update all
            write_curr()
            curr_set, curr_chrom, curr_start, curr_end = next_set, next_chrom, next_pos, next_pos
        elif next_set is not None:
            # continuation of set. just update end
            curr_end = next_pos

    # write the last set
    if curr_set is not None:
        write_curr()

    pb = pandas.DataFrame({'chrom': chrom, 'phase_set': phase_set, 'start': start, 'end': end, 'length': length})
    if filter_trivial:
        pb = pb[(pb.length > 0) & (pb.phase_set != -1)]
    return pb[['chrom', 'phase_set', 'start', 'end', 'length']].sort(['chrom', 'phase_set']).reset_index(drop=True)

def get_phase_blocks_from_h5(phase_blocks_h5, filter_trivial=True, locus=None):
    pb = tk_h5.read_data_frame(phase_blocks_h5)
    if filter_trivial:
        pb = pb[(pb.num_hets > 0) & (pb.fract_phased > 0) & (pb.length > 0) & (pb.phase_set != -1)]
    if locus is not None:
        pb = pb[(pb.chrom == locus[0]) & (pb.start >= locus[1]) & (pb.end <= locus[2])]
    return pb[['chrom', 'phase_set', 'start', 'end', 'length']].sort(['chrom', 'phase_set']).reset_index(drop=True)

def get_phase_blocks_from_frags(frag_phasing_tsv, filter_trivial=True, locus=None):
    compress = 'gzip' if frag_phasing_tsv.endswith('.gz') else None
    fp = pandas.read_table(frag_phasing_tsv, compression=compress)

    # make sure chrom is formatted as string categorical
    fp['#chrom'] = fp['#chrom'].astype('str').astype('category')

    pb = fp[['#chrom', 'phase_set', 'ps_start', 'ps_end']].drop_duplicates()
    pb['ps_length'] = pb['ps_end'] - pb['ps_start']
    # rename columns to match phase_blocks.h5
    pb.columns = ['chrom', 'phase_set', 'start', 'end', 'length']
    if filter_trivial:
        pb = pb[(pb.length > 0) & (pb.phase_set != -1)]
    if locus is not None:
        pb = pb[(pb.chrom == locus[0]) & (pb.start >= locus[1]) & (pb.end <= locus[2])]
    return pb.sort(['chrom', 'phase_set']).reset_index(drop=True)

def count_bc_redundancies(fp):
    '''
    Count the number of cases where a barcode appears more than once per phase block.
    '''
    # TODO actually use this? it's surprisingly slow.
    fp_unique = fp[['#chrom', 'phase_set', 'bc']].drop_duplicates()
    return len(fp) - len(fp_unique)

def check_overlaps(pb, max_allowable_overlap=100):
    '''
    Find cases where adjacent phase blocks overlap.
    Ignore overlaps that are less than max_allowable_overlap
    (useful when working with canonicalized data).
    '''
    pbs = pb.sort(['chrom', 'start'])
    # create two dataframes, shifted by 1 from each other, to keep things vectorized
    left = pbs.iloc[:-1].reset_index()
    right = pbs.iloc[1:].reset_index()
    chroms_match = (left.chrom == right.chrom)
    overlaps = left.end - right.start
    discord = chroms_match & (overlaps > max_allowable_overlap)
    return (left[discord], right[discord])

def check_concordance(pb_left, pb_right, vcf, max_allowable_overhang_length=1000, max_allowable_overhang_variants=2, filter_trivial=True):
    '''
    Compare two sets of phase blocks - what's shared and what's different. 
    The max_allowable_overhang parameters allow one to ignore minor differences.
    '''
    # hack: change chrom names to numbers to get around pandas 0.15.2 bug when merging categoricals
    # (see https://github.com/pydata/pandas/issues/9426)
    all_chroms = pandas.Categorical(pb_left.chrom.append(pb_right.chrom).drop_duplicates())
    chrom_left_codes = pandas.Categorical(pb_left.chrom).set_categories(all_chroms).codes
    chrom_right_codes = pandas.Categorical(pb_right.chrom).set_categories(all_chroms).codes
    
    pb_left_tmp = pandas.DataFrame({'chrom': chrom_left_codes, 'phase_set': pb_left.phase_set, 'start_left': pb_left.start, 'end_left': pb_left.end})
    pb_right_tmp = pandas.DataFrame({'chrom': chrom_right_codes, 'phase_set': pb_right.phase_set, 'start_right': pb_right.start, 'end_right': pb_right.end})
    
    pb_merged = pandas.merge(pb_left_tmp, pb_right_tmp, on=['chrom', 'phase_set'], how='outer')

    if filter_trivial:
        pb_merged = pb_merged[(pb_merged.end_right - pb_merged.start_right > 1) | (pb_merged.end_left - pb_merged.start_left > 1)].reset_index()

    # map codes back to chrom names
    pb_chroms = pandas.Categorical.from_codes(pb_merged.chrom, all_chroms.categories)
    pb_merged.chrom = pb_chroms

    # split into shared/non-shared phase sets
    pb_merged_shared = pb_merged[pb_merged.start_left.notnull() & pb_merged.start_right.notnull()]
    pb_merged_left_only = pb_merged[pb_merged.start_right.isnull()]
    pb_merged_right_only = pb_merged[pb_merged.start_left.isnull()]

    # look at overlaps in the shared sets
    start_max = pb_merged_shared.start_left.combine(pb_merged_shared.start_right, max)
    end_min = pb_merged_shared.end_left.combine(pb_merged_shared.end_right, min)
    overlap = map(lambda x: max(x, 0), end_min - start_max)
    non_overlap = (pb_merged_shared.end_left - pb_merged_shared.start_left) + (pb_merged_shared.end_right - pb_merged_shared.start_right) - map(lambda x: 2*x, overlap)
    discord = (non_overlap > max_allowable_overhang_length)

    pb_merged_shared_concordant = pb_merged_shared[~discord]
    pb_merged_shared_discordant = pb_merged_shared[discord]

    # ignore cases with only a few overhang variants (probably due to off-by-one errors)
    if len(pb_merged_shared_discordant) > 0:
        vfr = tk_io.VariantFileReader(vcf)
        pb_merged_shared_discordant = pb_merged_shared_discordant[pb_merged_shared_discordant.apply(lambda x: too_many_overhang_variants(x, vfr, max_allowable_overhang_variants), axis=1)]

    return (pb_merged_shared_concordant, pb_merged_shared_discordant, pb_merged_left_only, pb_merged_right_only)

def too_many_overhang_variants(pb, vfr, max_allowable_overhang_variants):
    # find overhangs
    chrom = pb.chrom
    overhang0_start = min(pb.start_left, pb.start_right)
    overhang0_end = max(pb.start_left, pb.start_right)
    overhang1_start = min(pb.end_left, pb.end_right)
    overhang1_end = max(pb.end_left, pb.end_right)

    if overhang1_start <= overhang0_end:
        overhang_loci = [(chrom, overhang0_start, overhang1_end)]
    else:
        overhang_loci = [(chrom, overhang0_start, overhang0_end), (chrom, overhang1_start, overhang1_end)]

    badness = False
    for loc in overhang_loci:
        locus_string = tk_io.create_locus_info(*loc)
        variants = [record for record in tk_io.get_variant_iterator_pos(vfr, None, locus_string)]
        if len(variants) > max_allowable_overhang_variants:
            badness = True
    return badness

def validate_phase_blocks(vcf, frag_phasing, filter_trivial=True, locus=None):
    # validate phase blocks from two different sources (VCF, fragment phasing)
    # NB: specifying a locus that's not a whole chromosome may lead to some edge effects,
    # due to differences in the way vcf_blocks / frag_blocks are computed

    # TODO parallelize get_phase_blocks_from_vcf if no locus is specified? the whole VCF might take a while.
    vcf_blocks = get_phase_blocks_from_vcf(vcf, filter_trivial=filter_trivial, locus=locus)
    frag_blocks = get_phase_blocks_from_frags(frag_phasing, filter_trivial=filter_trivial, locus=locus)

    (vcf_overlapping_left, vcf_overlapping_right) = check_overlaps(vcf_blocks)
    (frag_overlapping_left, frag_overlapping_right) = check_overlaps(frag_blocks)
    (concordant, discordant, vcf_only, frags_only) = check_concordance(vcf_blocks, frag_blocks, vcf)

    # all these numbers should be zero if working with non-canonicalized data.
    try:
        assert(len(vcf_overlapping_left) == 0)
        assert(len(frag_overlapping_left) == 0)
        assert(len(discordant) == 0)
        assert(len(vcf_only) == 0)
        assert(len(frags_only) == 0)
    except AssertionError as e:
        e.args += ("The fragment_phasing disagrees with the phase blocks implied by the VCF", vcf, frag_phasing)
        raise