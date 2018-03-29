#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
## Code to output a file containing barcodes mapped to locations
#
import os.path
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.tabix as tk_tabix

import subprocess

import martian

__MRO__ = """
stage INFER_CONTIGS(
    in  string restrict_locus,
    in  bam    input,
    in  int    window_size,
    out tabix  contig_output,
    src py     "stages/reads/infer_fragments",
) split using (
    in  string chrom,
)
"""

# THIS SHOULDNT BE NECESSARY!!
BUFFER = 500

def split(args):
    bam_in = tk_bam.create_bam_infile(args.input)

    if args.restrict_locus is None:
        chunk_defs = [{'chrom':chrom} for chrom in bam_in.references]
    else:
        chrom, start, stop = tk_io.get_locus_info(args.restrict_locus)
        chunk_defs = [{'chrom':chrom}]

    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    contig_outs = open(outs.contig_output, 'w')
    contig_outs.write('\t'.join(['#CHROM', 'START', 'END', 'BC_SEQ', 'NUM_READS']) + '\n')

    for chunk in chunk_outs:
        with open(chunk.contig_output, 'r') as chunk_contigs:
            for l in chunk_contigs:
                contig_outs.write(l)

    contig_outs.close()


def main(args, outs):
    """ Outputs barcode file """
    args.coerce_strings()
    bam_in = tk_bam.create_bam_infile(args.input)
    unsorted_temp_name = martian.make_path(outs.contig_output + '_TEMPUNSORTED')
    sorted_temp_name = martian.make_path(outs.contig_output + '_TEMPSORTED')
    base_dir = os.path.dirname(outs.contig_output)
    unsorted_temp_file = open(unsorted_temp_name, 'w')
    contig_output_file = open(outs.contig_output, 'w')
    window_size = args.window_size

    chroms = bam_in.references

    # Output the raw poses
    unsorted_temp_file.write('\t'.join(['#CHROM', 'START', 'END', 'BC_SEQ']) + '\n')
    if args.restrict_locus is None:
        bam_iter = bam_in.fetch(args.chrom)
    else:
        restrict_chrom, restrict_start, restrict_stop = tk_io.get_locus_info(args.restrict_locus)
        assert(args.chrom == restrict_chrom)
        bam_iter = bam_in.fetch(restrict_chrom, restrict_start, restrict_stop)

    for read in bam_iter:
        chrom = chroms[read.tid]
        start = read.pos
        end = read.aend

        if end is None:
            end = start + len(read.seq)

        bc = tk_io.get_read_barcode(read)

        if not(bc is None):
            unsorted_temp_file.write('\t'.join([chrom, str(start), str(end), bc]) + '\n')

    # Sort the poses
    unsorted_temp_file.close()
    tk_tabix.sort_bc_loc_tabix(unsorted_temp_name, sorted_temp_name, temp_dir_name=base_dir)

    # Infer the contig locations
    # This header is written during join
    #contig_output_file.write('\t'.join(['#CHROM', 'START', 'END', 'BC_SEQ', 'NUM_READS']) + '\n')
    sorted_temp_file = open(sorted_temp_name, 'r')
    sorted_temp_file.readline()
    old_bc_seq = None
    bc_poses = []
    for line in sorted_temp_file:
        (chrom, start, end, bc_seq) = line.strip('\n').split('\t')
        start = int(start)
        end = int(end)

        if not(bc_seq == old_bc_seq):
            if not(old_bc_seq is None):
                frags = infer_fragments(bc_poses, window_size)
                for (frag_chrom, frag_start, frag_end, num_reads) in frags:
                    contig_output_file.write('\t'.join([frag_chrom, str(frag_start - BUFFER), str(frag_end + BUFFER), old_bc_seq, str(num_reads)]) + '\n')
            bc_poses = []
        old_bc_seq = bc_seq
        bc_poses.append((chrom, start, end))

    # Output for the last barcode
    if not(old_bc_seq is None):
        frags = infer_fragments(bc_poses, window_size)
        for (frag_chrom, frag_start, frag_end, num_reads) in frags:
            contig_output_file.write('\t'.join([frag_chrom, str(frag_start - BUFFER), str(frag_end + BUFFER), old_bc_seq, str(num_reads)]) + '\n')

    sorted_temp_file.close()
    subprocess.check_call(['rm', sorted_temp_name])
    subprocess.check_call(['rm', unsorted_temp_name])
    contig_output_file.close()

# SHOULD MOVE TO HAVING THE WINDOW LENGTH CAN BE INFERRED FROM THE DATA
def infer_fragments(poses, window_size):
    """ Clusters reads within window distance of each other. """
    poses_by_chrom = {}
    for (chrom, start, end) in poses:
        chrom_poses = poses_by_chrom.setdefault(chrom, [])
        chrom_poses.append((start, end))

    frags = []
    for chrom, chrom_poses in poses_by_chrom.iteritems():
        sorted_poses = sorted(chrom_poses)
        this_start, this_end = sorted_poses[0]
        num_reads = 1

        for (pos, end) in sorted_poses[1:]:
            if pos > this_end + window_size:
                frags.append((chrom, this_start, this_end, num_reads))
                this_start = pos
                this_end = end
                num_reads = 0
            this_end = end
            num_reads += 1

        frags.append((chrom, this_start, this_end, num_reads))

    return frags
