#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make sparse matrices of bcs X genomic windows
#
import cPickle
import numpy as np
import tenkit.pandas as pd
import scipy.sparse as sp
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.hdf5 as tk_hdf5
import tenkit.reference as tk_reference
from longranger.sv.utils import *
import os
import os.path
import shutil

import martian

__MRO__ = """
stage COUNT_READS_BCS(
    in  bam      input,
    in  string   sex,
    in  string   reference_path,
    in  tsv      blacklist,
    in  h5       barcodes,
    in  pickle   inv_bc_map,
    in  bed      targets,
    in  int      target_extend,
    in  string   restrict_locus,
    in  int      window_size,
    in  int      step,
    in  int      min_reads,
    in  int      max_merge_dist,
    in  int      min_mapq,
    in  bool     read1_only,
    in  bool     no_split,
    in  bool     slide,
    out pickle[] bc_pos_mats,
    out pickle   inv_bc_map,
    out pickle   bc_counts,
    out pickle   win_counts,
    out pickle   loci,
    src py       "stages/structvars/count_reads_bcs",
) split using (
    in  string   chrom,
    in  int[]    starts,
    in  int[]    stops,
)
"""

MIN_TARGET_FRAC = 0.2

def split(args):
    win = args.window_size
    input_bam = tk_bam.create_bam_infile(args.input)
    chroms = input_bam.references
    chrom_lengths = input_bam.lengths
    chrom_len_map = {}
    for i, chrom in enumerate(chroms):
        chrom_len_map[chrom] = chrom_lengths[i]
    input_bam.close()

    max_mem_in_gb = 4 # Be a little conservative
    chunk_size = get_max_chunk(win, max_mem_in_gb)
    if not args.restrict_locus is None:
        locus_chrom, locus_start, locus_stop = tk_io.get_locus_info(args.restrict_locus)
        assert(locus_chrom in chrom_len_map)
        locus_start = max(0, locus_start)
        locus_stop = min(locus_stop, chrom_len_map[locus_chrom])

    primary_contigs = tk_reference.load_primary_contigs(args.reference_path)
    genome_size = np.sum([length for chrom, length in chrom_len_map.iteritems() if chrom in primary_contigs])
    
    prev_chrom = ''
    tot_bp = 0
    starts = []
    stops = []
    chunks = []

    # Genome-wide windows
    if not args.restrict_locus is None:
        chunks.append({'chrom':locus_chrom, 'starts':[locus_start], 'stops':[locus_stop], '__mem_gb':8})
    else:
        for chrom, length in chrom_len_map.iteritems():
            if not args.sex is None and args.sex.lower() in ['f', 'female'] and chrom in ['Y', 'chrY']:
                continue
            if not chrom in primary_contigs:
                continue
            nchunks = int(np.ceil(length / float(chunk_size)))
            # Divide as evenly as possible the windows across the chunks
            # This also makes sure that all chunks except the last will
            # have sizes that are multiples of the window size.
            win_per_chunk = int(np.ceil(length / float(nchunks * win)))
            new_chunk_size = win_per_chunk * win
            for c in range(nchunks):
                chunk_start = c * new_chunk_size
                chunk_stop = min((c + 1 ) * new_chunk_size, length)
                chunks.append({'chrom':chrom, 'starts':[chunk_start], 'stops':[chunk_stop], '__mem_gb':8})

    # Target-centered windows. If the targets (plus the extent) cover too much of the
    # genome, then skip these.
    if not args.targets is None and not args.target_extend is None:
        target_regions = []
        bed_iterator = tk_io.get_bed_iterator(args.targets)
        for chrom, start, stop in bed_iterator:
            if not args.sex is None and args.sex.lower() in ['f', 'female'] and chrom in ['Y', 'chrY']:
                continue
            if not chrom in primary_contigs:
                continue
            stop = min(chrom_len_map[chrom], stop)
            if args.restrict_locus is None or (chrom == locus_chrom and overlaps(start, stop, locus_start, locus_stop)):
                target_regions.append((chrom, start, stop))

        target_regions = sort_and_merge(target_regions, args.target_extend)
        target_size = np.sum([stop - start for _, start, stop in target_regions])
        
        if target_size / float(genome_size) < MIN_TARGET_FRAC:
            for (chrom, start, stop) in target_regions:
                if (prev_chrom != chrom and prev_chrom != '') or tot_bp > 1 * 1e7:
                    chunks.append({'chrom':str(prev_chrom), 'starts':starts, 'stops':stops, '__mem_gb':8})
                    starts = []
                    stops = []
                    tot_bp = 0
                    tot_bp += (stop - start)
                prev_chrom = chrom
                starts.append(start)
                stops.append(stop)

            if prev_chrom != '':
                chunks.append({'chrom':str(prev_chrom), 'starts':starts, 'stops':stops, '__mem_gb':8})
                
    return {'chunks': chunks}


def overlaps(start, stop, locus_start, locus_stop):
    return not(start > locus_stop or stop < locus_start)


def get_max_chunk(window_size, max_mem):
    # Assuming that the BC overlap matrices use 64 bit elements,
    # the size of the ovelap matrix will be
    # 8*(chunk_size / window_size)^2
    # Given max_mem Gb of memory, we get:
    # max_mem * (2^30) = 8*(chunk_size / window_size)^2
    # I'm actually accounting for 4 such matrices.
    return int(np.sqrt(max_mem * (2**30) * (window_size**2) / 32))


def get_max_mem(nwin):
    return int(np.ceil(3 * 8 * nwin * nwin / 2**30))


def main(args, outs):
    main_make_bc_locus_matrices(args, outs)


def join(args, outs, chunk_defs, chunk_outs):
    if args.inv_bc_map is None:
        shutil.copyfile(chunk_outs[0].inv_bc_map, outs.inv_bc_map)
    else:
        shutil.copyfile(args.inv_bc_map, outs.inv_bc_map)

    with open(outs.inv_bc_map, 'rb') as f:
        inv_bc_map = cPickle.load(f)

    nbcs = len(inv_bc_map)
    # bc_counts[chrom] will be the number of BCs in each window of chrom
    bc_counts = {}
    read_counts = np.zeros((nbcs,)) # num reads with each barcode
    win_counts = np.zeros((nbcs,)) # num windows with each barcode
    loci = []
    win_on_target = 0
    tot_win = 0
    bc_pos_mat_list = []

    for cidx, chunk in enumerate(chunk_outs):
        with open(chunk.bc_pos_mats, 'rb') as f:
            locus = cPickle.load(f)

        chrom = locus[0]
        loci.append(locus)
        out_bc_file = outs.bc_pos_mats.strip('.pickle') + '_chunk' + str(cidx) + '.pickle'
        shutil.copyfile(chunk.bc_pos_mats, out_bc_file)
        bc_pos_mat_list.append(out_bc_file)

        with open(chunk.win_counts, 'rb') as f:
            read_counts += cPickle.load(f)
            win_counts += cPickle.load(f)
            win_on_target += cPickle.load(f)
            tot_win += cPickle.load(f)
            tmp_bc_counts = cPickle.load(f)

        if not chrom in bc_counts:
            bc_counts[chrom] = tmp_bc_counts
        else:
            bc_counts[chrom] = np.concatenate((bc_counts[chrom], tmp_bc_counts))

    outs.bc_pos_mats = bc_pos_mat_list

    with open(outs.win_counts, 'wb') as f:
        cPickle.dump(read_counts, f)
        cPickle.dump(win_counts, f)
        cPickle.dump(win_on_target, f)
        cPickle.dump(tot_win, f)

    with open(outs.bc_counts, 'wb') as f:
        cPickle.dump(bc_counts, f)

    with open(outs.loci, 'wb') as f:
        cPickle.dump(loci, f)


def main_make_bc_locus_matrices(args, outs):
    window_size = args.window_size
    step = args.step
    min_mapq = args.min_mapq
    min_reads = args.min_reads
    max_merge_dist = args.max_merge_dist

    if window_size != step or args.targets is None:
        # Sliding/merging windows doesn't work for non-targeted data.
        min_reads = 0
        martian.log_info('Window merging not supported. Setting min_reads to 0.')

    if args.blacklist is None or not os.path.isfile(args.blacklist):
        blacklist = set([])
    else:
        black_df = pd.read_csv(args.blacklist, sep = '\t', index_col = False)
        blacklist = set(list(black_df['bc']))

    chrom = args.chrom
    starts = args.starts
    stops = args.stops
    if args.inv_bc_map is None:
        if args.barcodes is None:
            whitelist = set([])
        else:
            whitelist = set(tk_hdf5.read_data_frame(args.barcodes).bc)
        bc_map, inv_bc_map = get_bc_map(whitelist, blacklist)
    else:
        with open(args.inv_bc_map, 'r') as f:
            inv_bc_map = cPickle.load(f)
        bc_map = {}
        for i, b in inv_bc_map.iteritems():
            bc_map[b] = i

    bc_mat, win_starts, win_stops  = create_bc_matrix_step(args.input, chrom, starts, stops,
                                                           window_size, step, bc_map, min_mapq,
                                                           read1_only=args.read1_only,
                                                           no_split=args.no_split)

    if min_reads > 0:
        bc_mat, win_starts, win_stops = merge_bc_mat(bc_mat, win_starts, win_stops,
            min_reads, max_merge_dist, slide = args.slide)

    locus = (chrom, win_starts, win_stops)
    with open(outs.bc_pos_mats, 'wb') as outfile:
        cPickle.dump(locus, outfile)
        cPickle.dump(bc_mat, outfile)

    with open(outs.inv_bc_map, 'wb') as outfile:
        cPickle.dump(inv_bc_map, outfile)

    win_for_stats = get_non_overlapping_wins(locus[1], locus[2])
    win_on_target = len(win_for_stats)
    tot_win = len(locus[1])

    filtered_bc_mat = bc_mat.tocsc()
    filtered_bc_mat = filtered_bc_mat[:, win_for_stats].tolil()
    read_counts = np.array(filtered_bc_mat.sum(axis = 1)).flatten()
    filtered_bc_mat[filtered_bc_mat > 1.0] = 1
    win_counts = np.array(filtered_bc_mat.sum(axis = 1)).flatten()
    bc_mat[bc_mat > 1.0] = 1
    bc_counts = np.array(bc_mat.sum(axis = 0)).flatten()

    with open(outs.win_counts, 'wb') as f:
        cPickle.dump(read_counts, f)
        cPickle.dump(win_counts, f)
        cPickle.dump(win_on_target, f)
        cPickle.dump(tot_win, f)
        cPickle.dump(bc_counts, f)


def get_bc_map(whitelist, blacklist = set([])):
    """
    - whitelist: set of whitelisted barcodes
    - blacklist: set of blacklisted barcodes
    Return values:
    - bc_map: bc_map[bc] is a unique index for barcode bc
    - inv_bc_map: inv_bc_map[idx] is the barcode with index idx
    """
    bc_map = {}
    inv_bc_map = {}
    for bc in sorted(list(set(whitelist))):
        if bc in blacklist:
            continue
        inv_bc_map[len(bc_map)] = bc
        bc_map[bc] = len(bc_map)
    return (bc_map, inv_bc_map)


def create_bc_matrix_step(bam_filename, chrom, starts, stops,
                          win, step, bc_map, min_mapq=30, read1_only=False,
                          no_split=False):
    """ Creates a (BCs X Windows) sparse matrix of barcodes (columns) versus windowed locations
    """

    step = min(max(1, step), win)
    in_bam = tk_bam.create_bam_infile(bam_filename)
    nsteps_per_win = int(np.ceil(win / float(step)))
    nbcs = len(bc_map.values())
    bc_idx = []
    win_idx = []
    win_starts = []
    win_stops = []

    for start, stop in zip(starts, stops):
        nwin = int(np.ceil((stop - start)/float(step)))
        for read in in_bam.fetch(str(chrom), start, stop):
            mapq = read.mapq
            if mapq < min_mapq or (read1_only and read.is_read2) or read.is_duplicate or read.is_secondary:
                continue

            if no_split and ('SA' in [t[0] for t in read.tags]):
                continue

            bc = tk_io.get_read_barcode(read)
            if bc is None or bc == '' or not bc in bc_map:
                continue

            pos = read.pos
            if pos < start or pos > stop:
                continue

            last_win_idx = int(np.floor((pos - start)/float(step)))
            first_win_idx = max(0, last_win_idx - nsteps_per_win + 1)
            last_win_idx = min(last_win_idx, nwin - 1)
            first_win_idx = min(first_win_idx, nwin - 1)
            bc_idx.extend([bc_map[bc] for i in range(first_win_idx, last_win_idx + 1)])
            win_idx.extend([i + len(win_starts) for i in range(first_win_idx, last_win_idx + 1)])
        win_starts.extend([start + i * step for i in range(nwin)])
        win_stops.extend([min(stop, start + i * step + win) for i in range(nwin)])
    bc_idx = np.reshape(np.array(bc_idx, ndmin = 2), (1, len(bc_idx)))
    win_idx = np.reshape(np.array(win_idx, ndmin = 2), (1, len(win_idx)))
    bc_mat = sp.csc_matrix((np.ones((bc_idx.shape[1],), dtype = np.float32),
                            np.concatenate((bc_idx, win_idx), axis = 0)),
                           (nbcs, len(win_starts))).tolil()
    win_starts = np.array(win_starts).flatten()
    win_stops = np.array(win_stops).flatten()
    return (bc_mat, win_starts, win_stops)


def get_non_overlapping_wins(starts, stops):
    assert(len(starts) == len(stops))
    if len(starts) == 0:
        return np.array([])

    sel_wins = [0]
    for i in range(1, len(starts)):
        if starts[i] >= stops[sel_wins[-1]]:
            sel_wins.append(i)
    return np.array(sel_wins)


def merge_bc_mat(bc_mat, starts, stops, min_count, max_merge_dist, slide = False):
    nbcs, nwin = bc_mat.shape
    assert(nwin == len(starts))
    assert(nwin == len(stops))
    assert(np.all(starts == sorted(starts)))
    assert(np.all(stops == sorted(stops)))

    max_merge_dist = max(max_merge_dist, 0)

    rows = []
    cols = []
    vals = []
    new_starts = []
    new_stops = []
    is_good = np.ones((nwin, ), dtype = np.bool)
    bc_mat = bc_mat.tocsc()
    bc_counts = np.array(bc_mat.sum(axis = 0)).flatten()
    starts = np.array(starts).flatten()
    stops = np.array(stops).flatten()

    # merge windows to achieve >= min_count BCs per window
    for i in range(nwin):
        if not is_good[i]:
            continue
        assert np.all(starts[(i+1):] >= stops[i]), 'Cannot merge if the windows are overlapping'
        cand_win = np.where(starts[i:] - stops[i] < max_merge_dist)[0]
        last_idx = np.where(np.cumsum(bc_counts[cand_win + i]) > min_count)[0]
        if len(last_idx) > 0:
            last_idx = last_idx[0] + i
        else:
            last_idx = np.min([nwin - 1, cand_win[-1] + i])
        new_counts = np.reshape(bc_mat[:, list(range(i, last_idx + 1))].sum(axis = 1), (nbcs, 1))
        # flag windows that have been merged so we don't touch them again
        if not slide:
            is_good[(i + 1):(last_idx + 1)] = False
        r, c, v = sp.find(new_counts)
        vals.extend(v)
        rows.extend(r)
        cols.extend(np.ones(c.shape) * np.sum(is_good[0:i]))
        new_starts.append(starts[i])
        new_stops.append(stops[last_idx])

    rows = np.reshape(np.array(rows, ndmin = 2), (1, len(rows)))
    cols = np.reshape(np.array(cols, ndmin = 2), (1, len(cols)))
    new_bc_mat = sp.csc_matrix((vals, np.concatenate((rows, cols), axis = 0)),
                        (nbcs, np.sum(is_good))).tolil()

    assert(len(new_starts) == np.sum(is_good))
    new_starts = np.array(new_starts).flatten()
    new_stops = np.array(new_stops).flatten()
    return (new_bc_mat, new_starts, new_stops)
