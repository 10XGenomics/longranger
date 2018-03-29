#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import sys
import numpy as np
import longranger.sv.io as tk_sv_io
import tenkit.stats as tk_stats
import longranger.sv.sv_call as sv_call
import longranger.sv.utils as tk_sv_utils
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.pandas as pd

from itertools import product, groupby, combinations
from collections import defaultdict

from tenkit.constants import READ_MATE_CHIM_TOO_CLOSE_DIST
from longranger.sv.constants import (MAX_FRAG_SIZE, MAX_STORED_READPAIRS,
                                     MAX_CLUSTER_READPAIRS, MIN_LOG_PROB,
                                     SPLIT_LEN)

# Supported types of SVs
TDUP_STR = 'DUP'
DEL_STR = 'DEL'
INV_STR = 'INV'
TRANS_STR = 'TRANS'
COMPLEX_STR = 'COMPLEX'
UNKNOWN_STR = 'UNK'
NORMAL_STR = 'NORMAL'
TRANS_STR = 'TRANS'
TRANS_FF_STR = 'TRANS_FF'
TRANS_FR_STR = 'TRANS_FR'
TRANS_RF_STR = 'TRANS_RF'
TRANS_RR_STR = 'TRANS_RR'
ALL_TRANS_STRS = ['TRANS', 'TRANS_FF', 'TRANS_FR', 'TRANS_RF', 'TRANS_RR']


def pos_overlaps(pos, region):
    return pos >= region[0] and pos < region[1]


def get_distant_readpairs(in_bam, chroms, starts, stops, max_insert = 500, min_mapq = 60):
    r1 = (starts[0], stops[0])
    r2 = (starts[1], stops[1])
    readpairs = get_readpairs(in_bam, chroms, starts, stops, max_insert = max_insert,
                              min_mapq = min_mapq)
    return [rp for rp in readpairs if ((pos_overlaps(rp.read1.pos, r1) and pos_overlaps(rp.read2.pos, r2)) or
                                       (pos_overlaps(rp.read1.pos, r2) and pos_overlaps(rp.read2.pos, r1)))]



def get_readpairs2(in_bam, chroms, starts, stops, max_insert = 500, min_mapq=60):
    """Gets all readpairs or split-reads in the specified set of regions"""
    qname_info = {}
    pair_info = []

    nreads = 0
    for chrom, start, stop in zip(chroms, starts, stops):
        for read in in_bam.fetch(str(chrom), start, stop):
            if read.is_duplicate or read.is_unmapped or read.mate_is_unmapped or read.mapq < min_mapq or read.rnext == -1:
                continue
            if tk_io.get_read_barcode(read) is None:
                continue
            
            l = qname_info.setdefault(read.qname, [])
            l.append(read)
            nreads += 1
            if nreads >= MAX_STORED_READPAIRS:
                continue

        if nreads >= MAX_STORED_READPAIRS:
            continue

    for (k, reads) in qname_info.iteritems():

        for (r1, r2) in combinations(reads, 2):
            if is_pair(r1, r2) or is_split(r1, r2):
                pair_info.append(ReadPair(r1, r2, max_insert))

        if len(pair_info) > MAX_STORED_READPAIRS:
            return pair_info

    return pair_info



def get_readpairs(in_bam, chroms, starts, stops, max_insert=500,
                  min_mapq=60, normal_only=False):
    """Gets all readpairs or split-reads in the specified set of regions"""
    pair_info = {}
    tmp_pair_info = {}

    for chrom, start, stop in zip(chroms, starts, stops):
        for read in in_bam.fetch(str(chrom), start, stop):
            if read.is_duplicate or read.is_unmapped or read.mate_is_unmapped or read.mapq < min_mapq or read.rnext == -1:
                continue
            if tk_io.get_read_barcode(read) is None:
                continue

            chrom1, pos1, chrom2, pos2 = in_bam.getrname(read.tid), read.pos, in_bam.getrname(read.rnext), read.mpos
            if (read.qname, chrom1, pos1) in tmp_pair_info:
                read1 = tmp_pair_info[(read.qname, chrom1, pos1)]
                if (is_pair(read1, read) or is_split(read1, read)):
                    pair_info[read.qname] = ReadPair(read1, read, max_insert)
                    if len(pair_info) > MAX_STORED_READPAIRS:
                        return pair_info.values()
                del tmp_pair_info[(read.qname, chrom1, pos1)]
            else:
                if not (read.qname, chrom2, pos2) in tmp_pair_info:
                    if not normal_only or (chrom2 == chrom1 and not read.is_reverse and read.mate_is_reverse and abs(pos2 - pos1) + len(read.seq) < max_insert):
                        tmp_pair_info[(read.qname, chrom2, pos2)] = read

    def get_additional_reads(c, s, e):
        for read in in_bam.fetch(str(c), s, e):
            if read.is_duplicate or read.is_unmapped or read.mate_is_unmapped or read.mapq < min_mapq or read.rnext == -1:
                continue
            if tk_io.get_read_barcode(read) is None:
                continue
            chrom1, pos1 = in_bam.getrname(read.tid), read.pos
            if (read.qname, chrom1, pos1) in tmp_pair_info:
                read1 = tmp_pair_info[(read.qname, chrom1, pos1)]
                if (is_pair(read1, read) or is_split(read1, read)):
                    pair_info[read.qname] = ReadPair(read1, read, max_insert)
                del tmp_pair_info[(read.qname, chrom1, pos1)]

    for chrom, start, stop in zip(chroms, starts, stops):
        get_additional_reads(chrom, max(0, start - max_insert), start)
        get_additional_reads(chrom, stop, stop + max_insert)

    return pair_info.values()


def get_discordant_loci(bam_file, chrom=None, starts=None, stops=None,
                        min_mapq=0, min_insert=0, max_insert=1000, max_merge_range=1000,
                        min_sv_len = 100, max_sv_len = np.inf,
                        ins_logsf_fun = None, min_lr_to_call = 0, min_reads_to_call = 0,
                        chimera_rate=1e-9, reads_as_qual = False):
    '''Gets discordant reads in a given st of loci.
    - chrom, start, stop: a region from which reads will be fetched.
    loci might be outside this region).
    - min/max_insert: insert sizes smaller or larger than this will be considered "discordant"
    - min_sv_len: if a read-pair / split-read suggests an SV shorter than this, it will be ignored.
    Return value:
    (df, read_counts, pair_info)
    - df: DataFrame with discordant loci
    - read_counts: dict with counts of discordant pairs and split reads
    '''
    in_bam = tk_bam.create_bam_infile(bam_file)
    tmp_pair_info = {}
    pair_info = {}
    df = None
    read_counts = {}
    read_counts['split'] = defaultdict(int)
    read_counts['pair'] = defaultdict(int)

    min_bad_split_dist = max(READ_MATE_CHIM_TOO_CLOSE_DIST, max(min_sv_len, max_insert))
    
    for start, stop in zip(starts, stops):
        for read in in_bam.fetch(str(chrom), start, stop):
            if read.mapq < min_mapq or read.is_unmapped or read.mate_is_unmapped or read.is_duplicate or read.rnext == -1:
                continue
            chrom1, pos1, chrom2, pos2 = in_bam.getrname(read.tid), read.pos, in_bam.getrname(read.rnext), read.mpos

            if max_sv_len < np.inf and chrom1 != chrom2:
                continue
            
            alignments = tk_io.get_read_chimeric_alignments(read)

            # This might not be accurate if read is not read2.
            ins = abs(pos1 - pos2) + len(read.seq)
            if (read.qname, chrom1, pos1) in tmp_pair_info or (read.qname, chrom1, pos1 - 1) in tmp_pair_info:
                # this is the discordant pair we were looking for
                if (read.qname, chrom1, pos1) in tmp_pair_info:
                    read1, _, _ = tmp_pair_info[(read.qname, chrom1, pos1)]
                    del tmp_pair_info[(read.qname, chrom1, pos1)]
                else:
                    read1, _, _ = tmp_pair_info[(read.qname, chrom1, pos1 - 1)]
                    del tmp_pair_info[(read.qname, chrom1, pos1 - 1)]
                if is_pair(read1, read) or is_split(read1, read):
                    rp = ReadPair(read1, read, max_insert)
                    pair_info[read.qname] = rp
                    if rp.is_split():
                        read_counts['split'][rp.sv_type] += 1
                    else:
                        read_counts['pair'][rp.sv_type] += 1

            elif not alignments is None and len(alignments) == 1:
                # do not consider reads that were split into >2 parts
                a = alignments[0]
                chrom2, pos2, split_mapq = a[0], a[1], a[4]
                split_dist = np.abs(pos1 - pos2)
                good_dist = (chrom1 != chrom2 and max_sv_len == np.inf) or \
                            (chrom1 == chrom2 and split_dist > min_bad_split_dist and split_dist < max_sv_len)
                if split_mapq >= min_mapq and good_dist:
                    if not (read.qname, chrom2, pos2) in tmp_pair_info and not read.qname in pair_info:
                        tmp_pair_info[(read.qname, chrom2, pos2)] = (read, SPLIT_LEN, True)
            elif ins < max_sv_len and (chrom1 != chrom2 or ins > max(min_sv_len, max_insert)):
                if not (read.qname, chrom2, pos2) in tmp_pair_info and not read.qname in pair_info:
                    tmp_pair_info[(read.qname, chrom2, pos2)] = (read, int(max_insert / 2.0), False)

            if len(pair_info) > MAX_STORED_READPAIRS:
                tmp_df = call_breaks_from_readpairs(in_bam, pair_info, ins_logsf_fun,
                    max_insert, max_merge_range, min_sv_len, min_reads_to_call, min_lr_to_call, chimera_rate)
                df = pd.concat([df, tmp_df], ignore_index = True)
                pair_info.clear()

        tmp_df = call_breaks_from_readpairs(in_bam, pair_info, ins_logsf_fun,
                                            max_insert, max_merge_range, min_sv_len, min_reads_to_call, min_lr_to_call, chimera_rate,
                                            reads_as_qual = reads_as_qual)
        df = pd.concat([df, tmp_df], ignore_index = True)
        pair_info.clear()

    in_bam.close()
    return df, read_counts, tmp_pair_info


def cluster_pair_info(in_bam, pair_info, max_merge_range):
    read_pair_loci = []
    for name, rp in pair_info.iteritems():
        chrom1, pos1, chrom2, pos2 = in_bam.getrname(rp.read1.tid), rp.read1.pos, in_bam.getrname(rp.read2.tid), rp.read2.pos
        if rp.is_split():
            ext = SPLIT_LEN
        else:
            ext = int(rp.max_insert / 2.0)
        if chrom1 == chrom2:
            ext = min(ext, int(0.5 * (pos2 - pos1)))
        read_pair_loci.append((chrom1, max(0, pos1 - ext), pos1 + ext, (name, 1)))
        read_pair_loci.append((chrom2, max(0, pos2 - ext), pos2 + ext, (name, 2)))
    _, mem_to_cluster, _ = tk_sv_utils.cluster_loci(read_pair_loci, 0, max_range = max_merge_range)

    # Add a link between clusters for every read-pair / split read.
    cluster_pairs = defaultdict(list)
    for name, rp in pair_info.iteritems():
        cluster_idx1 = mem_to_cluster[(name, 1)]
        cluster_idx2 = mem_to_cluster[(name, 2)]
        if cluster_idx1 != cluster_idx2:
            cluster_pairs[(cluster_idx1, cluster_idx2)].append(rp)
    return cluster_pairs


def call_breaks_from_clusters(in_bam, cluster_pairs, ins_logsf_fun, max_insert, max_merge_range,
                              min_sv_len, min_reads_to_call, min_lr_to_call, chimera_rate,
                              reads_as_qual = False):

    chroms1 = []
    starts1 = []
    stops1 = []
    chroms2 = []
    starts2 = []
    stops2 = []
    info_strs = []
    quals = []
    names = []

    for i, cluster in enumerate(cluster_pairs.values()):
        chrom1 = in_bam.getrname(cluster[0].read1.tid)
        chrom2 = in_bam.getrname(cluster[0].read2.tid)
        if len(cluster) > MAX_CLUSTER_READPAIRS or len(cluster) < min_reads_to_call:
            continue
        res_arr = get_sv_lr(cluster, ins_logsf_fun, max_insert, chimera_rate)
        for j, res in enumerate(res_arr):
            (lr, num_split, num_pairs, sv_len, support_range, sv_type, support_readpairs) = res
            if chrom1 != chrom2:
                sv_len = np.inf
            range1, range2 = support_range
            if num_split + num_pairs >= min_reads_to_call and lr >= min_lr_to_call and sv_len >= min_sv_len and not range1 is None and not range2 is None:
                chroms1.append(chrom1)
                chroms2.append(chrom2)
                starts1.append(int(range1[0]))
                stops1.append(int(range1[1]))
                starts2.append(int(range2[0]))
                stops2.append(int(range2[1]))
                sv_name = str(i) + '_' + str(j)
                names.append(sv_name)
                if reads_as_qual:
                    quals.append(num_split + num_pairs)
                else:
                    quals.append(lr)
                haps = ','.join([str(h) for h in tk_sv_utils.get_read_haps(support_readpairs)])
                bcs = set([tk_io.get_read_barcode(rp.read1) for rp in support_readpairs if not tk_io.get_read_barcode(rp.read1) is None])
                info_strs.append(tk_sv_io.update_info('.', ['NPAIRS', 'NSPLIT', 'NBCS', 'TYPE', 'HAP_READS', 'OVERLAP_SVS'],
                    [num_pairs, num_split, len(bcs), sv_type, haps, len(res_arr)]))

    df = tk_sv_io.create_sv_df(chroms1, starts1, stops1, chroms2, starts2, stops2,
        names, quals, info_strs = info_strs)
    return df


def call_breaks_from_readpairs(in_bam, pair_info, ins_logsf_fun, max_insert, max_merge_range,
                               min_sv_len, min_reads_to_call, min_lr_to_call, chimera_rate,
                               reads_as_qual = False):

    cluster_pairs = cluster_pair_info(in_bam, pair_info, max_merge_range)
    df = call_breaks_from_clusters(in_bam, cluster_pairs, ins_logsf_fun, max_insert, max_merge_range,
                                   min_sv_len, min_reads_to_call, min_lr_to_call, chimera_rate,
                                   reads_as_qual = reads_as_qual)
    return df


def get_clipped_loci(in_bam, chrom, start, stop, min_mapq = 60,
                     forward_only = False, reverse_only = False):
    supported_pos = defaultdict(int)
    for read in in_bam.fetch(str(chrom), start, stop):
        if read.is_duplicate or read.is_unmapped or read.mate_is_unmapped or read.mapq < min_mapq or read.rnext == -1:
            continue
        if tk_io.get_read_barcode(read) is None:
            continue
        if forward_only and read.is_reverse:
            continue
        if reverse_only and not read.is_reverse:
            continue
        if is_good_clipped(read):
            pos = get_clip_pos(read)
            supported_pos[pos] += 1
    return supported_pos


def is_split_consecutive(read1, read2):
    """Returns True if read1 and read2 are split reads mapped to (almost) consecutive positions."""
    # If the read was hard-clipped, then the length will be the length
    # of what was left. That's why we want to take the max here.
    qlen = max(read1.rlen, read2.rlen)

    # Number of matches or insertions. This is the number of bases of the
    # read that are consumed in the alignment. We want this to be close to
    # the length of the read.
    m1 = sum([p for (s, p) in read1.cigar if s == 0 or s == 1])
    m2 = sum([p for (s, p) in read2.cigar if s == 0 or s == 1])
    is_full_read_covered = abs(m1 + m2 - qlen) <= 2 * SPLIT_LEN

    start1 = map_qstart(read1)
    start2 = map_qstart(read2)
    if start1 < start2:
        are_consecutive = abs(start2 - (start1 + m1)) <= 2 * SPLIT_LEN
    else:
        are_consecutive = abs(start1 - (start2 + m2)) <= 2 * SPLIT_LEN
    return is_split(read1, read2) and is_full_read_covered and are_consecutive


def is_split(read1, read2):
    """Returns true iff read1/read2 are a pair of split reads.
    Only considers a pair of a primary and a secondary read as a split pair.
    """
    return read1.qname == read2.qname and read1.is_read1 == read2.is_read1 and read1.is_secondary != read2.is_secondary


def is_pair(read1, read2):
    """Returns true if read1 and read2 are a pair (of primary reads)."""
    return read1.qname == read2.qname and not read1.is_secondary and not read2.is_secondary and read1.is_read1 != read2.is_read1


def get_pair_break_dist(read, poses):
    """Distance from breakpoint to beginning of read."""
    if read.is_reverse:
        return read.pos + read.qlen - poses
    return poses - read.pos + map_qstart(read)


def get_clip_pos(read):
    if (read.is_reverse and map_qstart(read) == 0) or \
       (not read.is_reverse and map_qstart(read) > 0):
        return read.pos
    return read.aend


def is_good_clipped(read):
    ctuples = read.cigar
    if len(ctuples) != 2:
        return False
    if read.is_reverse:
        ctuples = ctuples[::-1]
    if ctuples[0][0] in [4, 5] and ctuples[1][0] == 0:
        return True
    if ctuples[1][0] in [4, 5] and ctuples[0][0] == 0:
        return True
    return False


def map_qstart(read):
    """Position within the read where mapping starts. Eg. if the first 4 bases of the read
    are soft-clipped, the position returned will be 5 (4th position, 0-based)."""
    ctuples = read.cigar
    if read.is_reverse:
        ctuples = ctuples[::-1]
    map_pos = 0
    for s, p in ctuples:
        if s == 0:
            break
        if s in [4, 5, 1]:
            map_pos += p
    return map_pos


def get_ranges(poses1, poses2, max_dist = 1):
    ranges = []
    if len(poses1) == 0:
        return []
    start1 = poses1[0]
    start2 = poses2[0]
    stop1 = poses1[0]
    stop2 = poses2[0]
    for (p1, p2) in zip(poses1, poses2):
        if (not overlaps((start1, stop1), (p1, p1), max_dist) or not overlaps((start2, stop2), (p2, p2), max_dist)):
            ranges.append(((start1, stop1 + 1), (start2, stop2 + 1)))
            start1 = p1
            start2 = p2
            stop1 = p1
            stop2 = p2
        start1 = min(start1, p1)
        start2 = min(start2, p2)
        stop1 = max(stop1, p1)
        stop2 = max(stop2, p2)
    ranges.append(((start1, stop1 + 1), (start2, stop2 + 1)))
    return ranges


def extend_range(r1, r2):
    return (min(r1[0], r2[0]), max(r1[1], r2[1]))

def default_range():
    return (1000000000000, 0)

def get_sv_lr(readpairs, ins_logsf_fun, max_ins, chimera_rate=1.0e-9):
    res = []

    sorted_rps = sorted(readpairs, key = lambda x:x.sv_type)
    for svt, readpair_iter in groupby(sorted_rps, lambda x:x.sv_type):
        if svt == NORMAL_STR:
            continue
        num_split = 0
        num_pairs = 0
        support_readpairs = []
        type_readpairs = list(readpair_iter)

        # Get the union of all breakpoint positions supported by any of the read pairs.
        r1 = default_range()
        r2 = default_range()
        for rp in type_readpairs:
            range1, range2 = rp.get_break_ranges(max_ins)
            r1 = extend_range(r1, range1)
            r2 = extend_range(r2, range2)

        space1 = max(1, int( float(r1[1] - r1[0]) / 2000.0))
        space2 = max(1, int( float(r2[1] - r2[0]) / 2000.0 ))
        amat, bmat = np.meshgrid(np.arange(r1[0], r1[1], space1, dtype=np.int32), np.arange(r2[0], r2[1], space2, dtype=np.int32))

        if amat.size > 0:
            a = np.array(amat.flat)
            b = np.array(bmat.flat)

            # Compute joint probabilities in the above ranges
            probs = np.zeros((len(a),))
            if isinstance(chimera_rate, dict):
                cr = chimera_rate.get(svt, chimera_rate.get(TRANS_STR, 1e-4))
            else:
                cr = chimera_rate
            for rp in type_readpairs:
                tmp_probs = rp.get_break_lr(a, b, max_ins, ins_logsf_fun, cr)
                probs = probs + tmp_probs

            lr = probs.max() # Maximize LR

            # Get a small region around the maximum
            sel = np.logical_and(probs > 0, probs > 0.99 * lr)
            if np.any(sel):
                ranges = get_ranges(a[sel], b[sel], max_ins)
                support_range = ranges[0]
                sv_len = support_range[1][0] - support_range[0][1]
                output_types = set([])

                # Finally, get the SV type based on the final range
                for rp in type_readpairs:
                    rp_range = rp.get_break_ranges(max_ins)
                    does_overlap = overlaps(rp_range[0], support_range[0], 0) and overlaps(rp_range[1], support_range[1], 0)
                    if does_overlap:
                        output_types.add(rp.sv_type)
                        num_split += rp.is_split()
                        num_pairs += rp.is_pair()
                        support_readpairs.append(rp)


                if len(ranges) > 1:
                    sv_type = COMPLEX_STR
                else:
                    assert(list(output_types) == [svt])
                    sv_type = svt
                res.append((int(lr), num_split, num_pairs, sv_len, support_range, sv_type, support_readpairs))
    return res


def get_full_sv_lr(cand_ranges, readpair_groups, ins_logsf_fun, max_ins, bc_weights,
              chimera_rate = 1e-9, in_sv_type = DEL_STR):
    pos_pairs = set([])
    num_split = 0
    num_pairs = 0
    lr = 0
    sv_len = 0

    pos_pairs = []
    for g in readpair_groups:
        for rp in g:
            if rp.sv_type != NORMAL_STR:
                range1, range2 = rp.get_break_ranges(max_ins)
                tmp_pairs = [(a, b) for (a, b) in product(np.arange(range1[0], range1[1]),
                                                          np.arange(range2[0], range2[1]))]
                pos_pairs.extend(tmp_pairs)
    pos_pairs = list(set(pos_pairs))

    if len(pos_pairs) == 0:
        return (lr, num_split, num_pairs, sv_len, (None, None), np.array([]), np.array([]))

    a = np.array([p[0] for p in pos_pairs])
    b = np.array([p[1] for p in pos_pairs])

    # Compute joint probabilities in the above ranges
    probs = np.zeros((len(a),))
    group_probs = np.zeros((len(a), len(readpair_groups)))
    for i, g in enumerate(readpair_groups):
        #tmp_probs = 0
        for rp in g:
            group_probs[:, i] += rp.get_break_lr(a, b, max_ins, ins_logsf_fun, chimera_rate, in_sv_type)
        #    tmp_probs += rp.get_break_lr(a, b, max_ins, ins_logsf_fun, chimera_rate, in_sv_type)
        if not bc_weights is None:
            probs += tk_stats.logaddexp([group_probs[:, i] + bc_weights[i],
                                tk_stats.log_1minus(np.exp(bc_weights[i]))])
        #    probs += logaddexp([tmp_probs + bc_weights[i],
        #                        log_1minus(np.exp(bc_weights[i]))])
        else:
            probs += group_probs[:, i]
        #    probs += tmp_probs

    lr = probs.max() # Maximize LR
    # Get a small region around the maximum
    sel = np.logical_and(probs > 0, probs > 0.99 * lr)
    if not np.any(sel):
        return (lr, num_split, num_pairs, sv_len, (None, None), np.array([]), np.array([]))

    ranges = get_ranges(a[sel], b[sel], max_ins)
    support_range = ranges[0]
    sv_len = support_range[1][0] - support_range[0][1]

    # Get the groups that support the SV
    support_groups = np.where(np.any(group_probs[sel, :] > 0, 0))[0]
    normal_groups = np.where(np.any(group_probs[sel, :] < 0, 0))[0]
    for gidx in support_groups:
        for rp in readpair_groups[gidx]:
            rp_range = rp.get_break_ranges(max_ins)
            if rp.sv_type != NORMAL_STR and overlaps(rp_range[0], support_range[0], 0) and overlaps(rp_range[1], support_range[1], 0):
                num_split += rp.is_split()
                num_pairs += rp.is_pair()

    return (lr, num_split, num_pairs, sv_len, support_range, support_groups, normal_groups)


def overlaps(range1, range2, max_dist = 1):
    direct_ov = not(range1[0] > range2[1] or range1[1] < range2[0])
    return direct_ov or 0 <= range2[0] - range1[1] <= max_dist or 0 <= range1[0] - range2[1] <= max_dist


class ReadPair(object):
    """"Class representing a junction-supporting pair of reads (i.e. improperly mapped or split)."""
    def __init__(self, read1, read2, max_insert=1000):
        assert is_pair(read1, read2) or is_split(read1, read2)
        # self.read1 will be the read that comes earlier on the genome
        if read1.tid == read2.tid and read1.pos > read2.pos:
            self.read1, self.read2 = read2, read1
        else:
            self.read1, self.read2 = read1, read2
        self.chrom1 = self.read1.tid
        self.chrom2 = self.read2.tid
        self.qstart1 = map_qstart(self.read1)
        self.qstart2 = map_qstart(self.read2)
        self.sv_type = self.get_sv_type(max_insert)
        self.max_insert = max_insert
        bc1 = tk_io.get_read_barcode(read1)
        bc2 = tk_io.get_read_barcode(read2)
        assert (bc1 is None and bc2 is None) or bc1 == bc2
        self.barcode = bc1

    def __str__(self):
        return str(self.read1.qname) + " " + str(self.is_split())

    def is_first_half_first(self):
        return self.qstart1 < self.qstart2

    def is_split(self):
        return is_split(self.read1, self.read2)

    def is_pair(self):
        return is_pair(self.read1, self.read2)

    def get_insert_size(self):
        return self.read2.pos - self.read1.pos + len(self.read2.seq)

    def get_sv_type(self, max_ins):
        """Returns the type of SV that the pair of reads supports."""
        if is_pair(self.read1, self.read2) and self.chrom1 == self.chrom2 and \
            not self.read1.is_reverse and self.read2.is_reverse and \
            self.get_insert_size() < max_ins:
            return NORMAL_STR

        # inter-chromosomal
        if self.chrom1 != self.chrom2 or \
            (self.chrom1 == self.chrom2 and self.get_insert_size() > MAX_FRAG_SIZE):
            if not self.read1.is_reverse and not self.read2.is_reverse:
                return TRANS_FF_STR
            if not self.read1.is_reverse and self.read2.is_reverse:
                return TRANS_FR_STR
            if self.read1.is_reverse and not self.read2.is_reverse:
                return TRANS_RF_STR
            if self.read1.is_reverse and self.read2.is_reverse:
                return TRANS_RR_STR

        # same-chromosome, pair-end
        if is_pair(self.read1, self.read2):
            if self.read1.is_reverse == self.read2.is_reverse:
                return INV_STR
            if self.read1.is_reverse:
                return TDUP_STR
            return DEL_STR

        # same chromosome, split read
        if self.read1.is_reverse != self.read2.is_reverse:
            return INV_STR
        if (self.read1.is_reverse and self.is_first_half_first()) or \
            (not self.read1.is_reverse and not self.is_first_half_first()):
            return TDUP_STR
        return DEL_STR

    @staticmethod
    def get_pair_break_range(read, other_read, max_ins):
        if read.is_reverse:
            return (max(0, read.pos - max_ins), read.pos + 1)
        else:
            return (read.aend, read.aend + max_ins + 1)


    def get_break_ranges(self, max_ins):
        """Returns the range of breakpoint coordinates supported by this ReadPair.
        Note that all positions returned at 0-based.

        Return value:
        A tuple of tuples ((start1, stop1), (start2, stop2)). A range (start1, stop1)
        means that the breakpoint happened right before basepair [start1, stop1).
        """

        if is_pair(self.read1, self.read2):
            # length of insert "between" the two reads
            max_ins = max(0, max_ins - self.read1.rlen - self.read2.rlen)

            range1 = ReadPair.get_pair_break_range(self.read1, self.read2, max_ins)
            range2 = ReadPair.get_pair_break_range(self.read2, self.read1, max_ins)
        else:
            if self.sv_type == INV_STR:
                if (self.read1.is_reverse and self.is_first_half_first()) or \
                    (not self.read1.is_reverse and not self.is_first_half_first()):
                    range1 = (self.read1.pos - SPLIT_LEN, self.read1.pos + SPLIT_LEN + 1)
                    range2 = (self.read2.pos - SPLIT_LEN, self.read2.pos + SPLIT_LEN + 1)
                else:
                    range1 = (self.read1.aend - SPLIT_LEN, self.read1.aend + SPLIT_LEN + 1)
                    range2 = (self.read2.aend - SPLIT_LEN, self.read2.aend + SPLIT_LEN + 1)
            elif self.sv_type == TDUP_STR:
                range1 = (self.read1.pos - SPLIT_LEN, self.read1.pos + SPLIT_LEN + 1)
                range2 = (self.read2.aend - SPLIT_LEN, self.read2.aend + SPLIT_LEN + 1)
            else:
                if (self.read1.is_reverse and self.qstart1 == 0) or \
                    (not self.read1.is_reverse and self.qstart1 > 0):
                    pos1 = self.read1.pos
                else:
                    pos1 = self.read1.aend
                if (self.read2.is_reverse and self.qstart2 == 0) or \
                    (not self.read2.is_reverse and self.qstart2 > 0):
                    pos2 = self.read2.pos
                else:
                    pos2 = self.read2.aend
                range1 = (pos1 - SPLIT_LEN, pos1 + SPLIT_LEN + 1)
                range2 = (pos2 - SPLIT_LEN, pos2 + SPLIT_LEN + 1)
        return (range1, range2)


    def get_break_lr(self, b1, b2, max_ins, ins_logsf_fun, chimera_rate=1.0e-9, in_sv_type=None):
        """Get the log-probability of the reads assuming that there is an SV
        of the same type as the pair at positions b1 and b2.

        Args:
        - b1/b2: array-like objects with lists of breakpoint positions. For every
        pair of positions in zip(b1, b2) it will compute the probability of the
        read pair given a breakpoint at this pair of positions.
        - mean_ins, std_ins, max_ins: mean and standard deviation of the insert
        size distribution and maximum insert size.
        """

        b1 = np.array(b1)
        b2 = np.array(b2)
        assert(len(b1) == len(b2))

        range1, range2 = self.get_break_ranges(max_ins)
        in_range_ind1 = np.logical_and(b1 >= range1[0], b1 < range1[1])
        in_range_ind2 = np.logical_and(b2 >= range2[0], b2 < range2[1])
        in_range_ind = np.logical_and(in_range_ind1, in_range_ind2)
        in_range_idx = np.where(in_range_ind)[0]

        probs = np.zeros(b1.shape)

        if is_pair(self.read1, self.read2):
            q1 = max(0.1, self.read1.mapq)
            q2 = max(0.1, self.read2.mapq)
            prob1_correct = tk_stats.log_prob_correct_from_qual(q1)
            prob2_correct = tk_stats.log_prob_correct_from_qual(q2)
            prob_pair_correct = prob1_correct + prob2_correct
            prob_pair_wrong = tk_stats.log_1minus(np.exp(prob_pair_correct))

            if self.sv_type == NORMAL_STR:
                ins_len = self.get_insert_size()
                sv_lr = prob_pair_wrong - tk_stats.logaddexp([prob_pair_wrong, prob_pair_correct + ins_logsf_fun([ins_len])[0]])
                if in_sv_type == DEL_STR:
                    # Anything inside the deletion should be false
                    in_range_idx = np.logical_and(self.read1.pos < b2, self.read2.pos > b1)
                elif in_sv_type == INV_STR:
                    in_range_idx = np.logical_or(np.logical_and(self.read1.pos < b1, self.read2.pos > b1),
                                                 np.logical_and(self.read1.pos < b2, self.read2.pos > b2))
                elif in_sv_type == TDUP_STR:
                    in_range_idx = np.zeros(b1.shape, dtype = np.bool)
                elif in_sv_type == TRANS_STR:
                    in_range_idx = np.logical_or(np.logical_and(self.read1.pos < b1, self.read2.pos > b1),
                                                 np.logical_and(self.read1.pos < b2, self.read2.pos > b2))
                probs[in_range_idx] = sv_lr
                return probs

            if len(in_range_idx) > 0:
                # Insert size implied by the given breakpoints
                ins_len = get_pair_break_dist(self.read1, b1) + get_pair_break_dist(self.read2, b2)
                # Probability under no SV (and not chimeric)
                if self.chrom1 == self.chrom2:
                    prob_norm = tk_stats.logaddexp([prob_pair_wrong, prob_pair_correct + ins_logsf_fun([self.get_insert_size()])[0]])
                else:
                    prob_norm = tk_stats.logaddexp([prob_pair_wrong, prob_pair_correct + np.log(chimera_rate)])
                probs[in_range_idx] = tk_stats.logaddexp([prob_pair_wrong, prob_pair_correct + ins_logsf_fun(ins_len[in_range_idx])]) - prob_norm
        else:
            # Either both sides are wrong or both are correct. The MAPQ (I believe) reflects
            # the mapping quality of the entire read (both sides).
            q = max(0.1, min(self.read1.mapq, self.read2.mapq))
            prob_pair_wrong = tk_stats.log_prob_wrong_from_qual(q)
            prob_pair_correct = tk_stats.log_prob_correct_from_qual(q)
            probs[in_range_idx] = prob_pair_correct - tk_stats.logaddexp([prob_pair_wrong, np.log(chimera_rate)])

        return probs


    def get_break_lr_new(self, b1, b2, max_ins, ins_logsf_fun, chimera_rate=1.0e-9):
        """Get the log-probability of the reads assuming that there is an SV
        of the same type as the pair at positions b1 and b2.

        Args:
        - b1/b2: array-like objects with lists of breakpoint positions. For every
        pair of positions in zip(b1, b2) it will compute the probability of the
        read pair given a breakpoint at this pair of positions.
        - mean_ins, std_ins, max_ins: mean and standard deviation of the insert
        size distribution and maximum insert size.
        """

        b1 = np.array(b1)
        b2 = np.array(b2)
        assert(len(b1) == len(b2))

        range1, range2 = self.get_break_ranges(max_ins)
        in_range_ind1 = np.logical_and(b1 >= range1[0], b1 < range1[1])
        in_range_ind2 = np.logical_and(b2 >= range2[0], b2 < range2[1])
        in_range_ind = np.logical_and(in_range_ind1, in_range_ind2)
        in_range_idx = np.where(in_range_ind)[0]

        probs = np.zeros(b1.shape)

        if len(in_range_idx) > 0:
            prob_pair_correct = self.get_prob_correct()
            prob_pair_wrong = tk_stats.log_1minus(np.exp(prob_pair_correct))

            if is_pair(self.read1, self.read2):

                # Probability under no SV (and not chimeric)
                if self.chrom1 == self.chrom2:
                    norm_prob = np.logaddexp(prob_pair_correct + ins_logsf_fun([self.get_insert_size()])[0],
                                             prob_pair_wrong + MIN_LOG_PROB)
                else:
                    # This is (prob_correct * epsilon) + (prob_incorrect * epsilon)
                    norm_prob = MIN_LOG_PROB
                # Probabability under no SV
                norm_prob = np.logaddexp(np.log(1 - chimera_rate) + norm_prob,
                                         np.log(chimera_rate) + MIN_LOG_PROB)

                # Insert size implied by the given breakpoints
                ins_len = get_pair_break_dist(self.read1, b1)
                ins_len += get_pair_break_dist(self.read2, b2)
                ins_prob = ins_logsf_fun(ins_len[in_range_idx])

                # Probability under SV (and not chimeric)
                sv_prob = np.logaddexp(prob_pair_correct + ins_prob,
                                       prob_pair_wrong + MIN_LOG_PROB)

                # Probability of observing read pair:
                # Either read pair is a chimera
                # - then the probability of observing this is epsilon
                # or it is not a chimera
                # - in which case it's either correctly mapped
                # -- then the probability is determined by the insert distribution
                # - or incorrectly mapped
                # -- then the probability of observing this is epsilon
                #
            else:
                norm_prob = MIN_LOG_PROB

                sv_prob = np.logaddexp(prob_pair_correct, prob_pair_wrong + MIN_LOG_PROB)

            # If inter-chromosomal read-pair, this will be positive if
            # sv_prob > MIN_LOG_PROB
            # Otherwise, it will be positive if
            # ins_prob > ins_logsf_fun([self.get_insert_size()])[0]
            # If split read, it will be always be positive.
            probs[in_range_idx] = np.logaddexp(np.log(1 - chimera_rate) + sv_prob,
                                               np.log(chimera_rate) + MIN_LOG_PROB) - norm_prob

        return probs


    def get_prob_correct(self):
        """Log-probability of a ReadPair being mapped correctly"""
        q1 = max(0.1, self.read1.mapq)
        q2 = max(0.1, self.read2.mapq)
        prob1_correct = tk_stats.log_prob_correct_from_qual(q1)
        prob2_correct = tk_stats.log_prob_correct_from_qual(q2)
        return prob1_correct + prob2_correct


    def get_chimera_rate(self, chimera_rates):
        return chimera_rates.get(self.sv_type, chimera_rates['TRANS'])


    def get_ref_prob(self, ins_logsf_fun, chimera_rates):
        """Probability of observing the read under the reference model"""

        if is_pair(self.read1, self.read2):

            chimera_rate = self.get_chimera_rate(chimera_rates)

            # Probability under no SV (and not chimeric)
            if self.sv_type == NORMAL_STR or self.sv_type == DEL_STR:
                prob_pair_correct = self.get_prob_correct()
                prob_pair_wrong = tk_stats.log_1minus(np.exp(prob_pair_correct))

                norm_prob = np.logaddexp(prob_pair_correct + ins_logsf_fun([self.get_insert_size()])[0],
                                         prob_pair_wrong + MIN_LOG_PROB)
            else:
                # This is (prob_correct * epsilon) + (prob_incorrect * epsilon)
                norm_prob = MIN_LOG_PROB
            # Probabability under no SV
            return np.logaddexp(np.log(1 - chimera_rate) + norm_prob,
                                np.log(chimera_rate) + MIN_LOG_PROB)
        else:
            return MIN_LOG_PROB


    def get_sv_prob(self, norm_prob, breaks1, breaks2, sv_type, max_ins, ins_logsf_fun, chimera_rates):
        breaks1 = np.array(breaks1)
        breaks2 = np.array(breaks2)
        assert len(breaks1) == len(breaks2)
        assert sv_type != NORMAL_STR

        chimera_rate = self.get_chimera_rate(chimera_rates)

        range1, range2 = self.get_break_ranges(max_ins)
        in_range_ind1 = np.logical_and(breaks1 >= range1[0], breaks1 < range1[1])
        in_range_ind2 = np.logical_and(breaks2 >= range2[0], breaks2 < range2[1])
        in_range_ind = np.logical_and(in_range_ind1, in_range_ind2)
        in_range_idx = np.where(in_range_ind)[0]

        probs = np.ones(breaks1.shape) * norm_prob

        if sv_type == self.sv_type:
            if len(in_range_idx) > 0:
                prob_pair_correct = self.get_prob_correct()
                prob_pair_wrong = tk_stats.log_1minus(np.exp(prob_pair_correct))

                if is_pair(self.read1, self.read2):

                    # Insert size implied by the given breakpoints
                    ins_len = get_pair_break_dist(self.read1, breaks1)
                    ins_len += get_pair_break_dist(self.read2, breaks2)
                    ins_prob = ins_logsf_fun(ins_len[in_range_idx])

                    # Probobility under SV (and not chimeric)
                    sv_prob = np.logaddexp(prob_pair_correct + ins_prob,
                                           prob_pair_wrong + MIN_LOG_PROB)
                else:
                    sv_prob = np.logaddexp(prob_pair_correct, prob_pair_wrong + MIN_LOG_PROB)
                probs[in_range_idx] = np.logaddexp(np.log(1 - chimera_rate) + sv_prob,
                                                   np.log(chimera_rate) + MIN_LOG_PROB)

        elif self.sv_type == NORMAL_STR:
            bad_idx = np.zeros(breaks1.shape, dtype = np.bool)

            if sv_type == INV_STR:
                bad_idx = np.logical_or(np.logical_and(self.read1.pos < breaks1, self.read2.pos > breaks1),
                                        np.logical_and(self.read1.pos < breaks2, self.read2.pos > breaks2))
            elif sv_type in ALL_TRANS_STRS:
                in_range_idx = np.logical_or(np.logical_and(self.read1.pos < breaks1, self.read2.pos > breaks1),
                                             np.logical_and(self.read1.pos < breaks2, self.read2.pos > breaks2))
            probs[bad_idx] = MIN_LOG_PROB

        if sv_type == DEL_STR:
            del_idx = np.logical_or(np.logical_and(self.read1.aend > breaks1, self.read1.pos < breaks2),
                                    np.logical_and(self.read2.aend > breaks1, self.read2.pos < breaks2))
            probs[del_idx] = MIN_LOG_PROB

        return probs


def eval_sv_em(read_groups, breaks1, breaks2, sv_type, chimera_rates,
               phase_set1, phase_set2,
               bc_phase_set_dict1, bc_phase_set_dict2,
               max_ins, ins_logsf_fun, em_iters = 5):

    nbcs = len(read_groups)
    loci = list(product(breaks1, breaks2))
    breaks1 = [locus[0] for locus in loci]
    breaks2 = [locus[1] for locus in loci]
    nloci = len(loci)

    # Get prior haplotype probabilities
    prior_hap_probs = init_hap_probs(read_groups, bc_phase_set_dict1, bc_phase_set_dict2)
    hap_probs = np.array(prior_hap_probs)

    no_sv_probs = np.zeros((nbcs,))
    sv_probs = np.zeros((nbcs, nloci))
    # het_sv_probs[b, n, h1, h2] is the probability of observing the data
    # for barcode b under the assumption that there is an sv at locus n
    # (or the n-th pair of loci), on haplotypes h1 and h2
    het_sv_probs = np.zeros((nbcs, nloci, 2, 2))

    if phase_set1 == phase_set2:
        het_sv_probs[:, :, 0, 1] = -np.inf
        het_sv_probs[:, :, 1, 0] = -np.inf

    for i, (bc, read_pairs) in enumerate(sorted(read_groups.iteritems())):
        for rp in read_pairs:
            norm_prob = rp.get_ref_prob(ins_logsf_fun, chimera_rates)
            no_sv_probs[i] += norm_prob
            sv_probs[i, :] += rp.get_sv_prob(norm_prob, breaks1, breaks2,
                                             sv_type, max_ins, ins_logsf_fun, chimera_rates)
            
    for em_iter in range(em_iters):
        print >> sys.stderr, 'EM iteration', em_iter
        for i, (_, read_pairs) in enumerate(sorted(read_groups.iteritems())):
            for hap1, hap2 in product([0, 1], [0, 1]):
                if phase_set1 == phase_set2 and hap1 != hap2:
                    continue
                hp = min(hap_probs[i, hap1], hap_probs[i, hap2 + 2])
                hp_neg = 1 - hp
                het_sv_probs[i, :, hap1, hap2] = np.logaddexp(np.log(hp) + sv_probs[i, :],
                                                              np.log(hp_neg) + no_sv_probs[i])
        new_hap_probs = posterior_hap_probs(no_sv_probs, sv_probs, het_sv_probs, prior_hap_probs)
        max_change = np.abs(np.max(new_hap_probs - hap_probs))
        print >> sys.stderr, 'Max change in phase probabilities', max_change
        if max_change < 1e-4:
            break
        hap_probs = new_hap_probs

    # Find optimal model
    # Sum over barcodes
    total_no_sv_probs = np.sum(no_sv_probs, axis=0)
    no_sv_max = np.max(total_no_sv_probs)

    total_sv_probs = np.sum(sv_probs, axis=0)
    sv_max_idx = np.argmax(total_sv_probs)
    sv_max = total_sv_probs[sv_max_idx]

    total_het_sv_probs = np.sum(het_sv_probs, axis=0)
    het_sv_max_idx = np.argmax(total_het_sv_probs)
    res = np.unravel_index(het_sv_max_idx, (nloci, 2, 2))
    het_sv_max_idx, het_sv_max_hap1, het_sv_max_hap2 = res
    het_sv_max = total_het_sv_probs[het_sv_max_idx, het_sv_max_hap1, het_sv_max_hap2]

    if sv_max >= het_sv_max:
        support = sv_probs[:, sv_max_idx] - no_sv_probs
        zygosity = sv_call.Zygosity.hom
        max_hap = (sv_call.Haplotype.unknown, sv_call.Haplotype.unknown)
        max_locus = loci[sv_max_idx]
    else:
        support = sv_probs[:, het_sv_max_idx] - no_sv_probs
        zygosity = sv_call.Zygosity.het
        max_hap = (sv_call.Haplotype.hap1 if het_sv_max_hap1 else sv_call.Haplotype.hap0,
                   sv_call.Haplotype.hap1 if het_sv_max_hap2 else sv_call.Haplotype.hap0)
        max_locus = loci[het_sv_max_idx]

    return ((no_sv_max, sv_max, het_sv_max), max_locus, zygosity,
            max_hap, prior_hap_probs, hap_probs, support)


def posterior_hap_probs(no_sv_probs, sv_probs, het_sv_probs, hap_probs):
    # Find optimal breakpoints and phasing
    # Sum across barcodes
    bc_probs = np.sum(het_sv_probs, axis = 0)
    res = np.unravel_index(np.argmax(bc_probs), bc_probs.shape)
    best_locus_idx, best_hap_idx1, best_hap_idx2 = res

    new_hap_probs = np.array(hap_probs)
    for i in range(het_sv_probs.shape[0]):
        # Mixed hap, don't touch.
        if 1.0 - hap_probs[i, 0] - hap_probs[i, 1] > 0.1 or 1.0 - hap_probs[i, 2] - hap_probs[i, 3] > 0.1:
            continue
        hp = min(hap_probs[i, best_hap_idx1], hap_probs[i, 2 + best_hap_idx2])
        hp_neg = 1 - hp
        p1 = sv_probs[i, best_locus_idx] + np.log(hp)
        p2 = no_sv_probs[i] + np.log(hp_neg)
        norm_factor = np.logaddexp(p1, p2)
        new_hap_probs[i, best_hap_idx1] = np.exp(p1 - norm_factor)
        new_hap_probs[i, 2 + best_hap_idx2] = np.exp(p1 - norm_factor)
        new_hap_probs[i, 1 - best_hap_idx1] = np.exp(p2 - norm_factor)
        new_hap_probs[i, 3 - best_hap_idx2] = np.exp(p2 - norm_factor)
    return new_hap_probs


def init_hap_probs(read_groups, bc_phase_set_dict1, bc_phase_set_dict2):
    nbcs = len(read_groups)
    # hap_probs[bc, 0] and hap_probs[bc, 1] are going to be the
    # probabilities of barcode bc being on haplotypes 0 or 1 respectively
    # around the first breakpoint
    # hap_probs[bc, 2] and hap_probs[bc, 3] are the same for the second
    # breakpoint.
    hap_probs = np.zeros((nbcs, 4))

    #bcs_on_hap0 = np.sum([v[2][0] > v[2][1] for v in bc_phase_set_dict1.values()])
    default_hap0_prob =  0.5 #bcs_on_hap0 / float(len(bc_phase_set_dict))

    for i, (bc, _) in enumerate(sorted(read_groups.iteritems())):
        if bc in bc_phase_set_dict1:
            hap_probs[i, 0] = bc_phase_set_dict1[bc][2][0]
            hap_probs[i, 1] = bc_phase_set_dict1[bc][2][1]
        else:
            hap_probs[i, 0] = default_hap0_prob
            hap_probs[i, 1] = 1 - default_hap0_prob

        if bc in bc_phase_set_dict2:
            hap_probs[i, 2] = bc_phase_set_dict2[bc][2][0]
            hap_probs[i, 3] = bc_phase_set_dict2[bc][2][1]
        else:
            hap_probs[i, 2] = default_hap0_prob
            hap_probs[i, 3] = 1 - default_hap0_prob

    return hap_probs
