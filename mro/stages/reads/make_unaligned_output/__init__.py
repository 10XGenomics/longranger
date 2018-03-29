#!/usr/bin/env python
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import subprocess
import pysam
import numpy as np
import json
import pandas

import tenkit.bam
import tenkit.stats as tk_stats
import martian
from tenkit.constants import PROCESSED_BARCODE_TAG, SAMPLE_INDEX_TAG, SAMPLE_INDEX_QUAL_TAG, TRIM_TAG, TRIM_QUAL_TAG, RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG


__MRO__ = '''
 stage MAKE_UNALIGNED_OUTPUT(
     in  string     sample_id,
     in  string     output_format,
     in  string     read_group,
     in  string[]   read_groups,
     in  fastq.gz[] reads,
     out fastq.gz   barcoded,
     out bam        barcoded_unaligned,
     out csv        basic_stats,
 ) split using (
     in  fastq.gz read_chunk,
     out int      num_pairs,
     out int      correct_barcode_pairs,
 )
'''

def split(args):
    chunks = [{'read_chunk': c} for (idx, c) in enumerate(args.reads)]
    return {'chunks': chunks}


def main(args, outs):

    def read_clusters(it):
        while True:
            head = str(it.next().strip())
            seq1 = str(it.next().strip())
            trim1 = seq1[:args.trim_length]
            seq1 = seq1[args.trim_length:]
            q1 = str(it.next().strip())
            trim_q1 = q1[:args.trim_length]
            q1 = q1[args.trim_length:]
            seq2 = str(it.next().strip())
            q2 = str(it.next().strip())
            bc = str(it.next().strip())
            rx = bc.split(',')
            bc = rx[0]
            if len(rx) > 1:
                rx = rx[1]
            else:
                rx = rx[0]
                bc = None
            bcq = str(it.next().strip())

            si = str(it.next().strip())
            siq = str(it.next().strip())

            lines = [head, seq1, q1, seq2, q2, rx, bc, bcq, si, siq, trim1, trim_q1]
            yield (bc, lines)


    try:
        version = martian.get_pipelines_version()
    except NameError:
        version = 'unknown'


    def make_rg_header(packed_rg_string):
        '''Make the RG header, matching how it's done in Lariat.'''
        result = packed_rg_string.split(':')
        if len(result) != 5:
            raise Exception("RG string must have this format - sample_id:library_id:gem_group:flowcell:lane")
        sample_id, library_id, gem_group, flowcell, lane = result
        return { 'ID': packed_rg_string, 'SM': sample_id, 'LB': library_id, 'PU': gem_group, 'PL': 'ILLUMINA' }
        #return '@RG\\tID:{0}\\tSM:{1}\\tLB:{2}.{3}\\tPU:{0}\\tPL:ILLUMINA'.format(packed_rg_string, sample_id, library_id, gem_group)

    header = {
            'HD': {'VN': '1.3', 'SO': 'unknown'},
            'RG': [make_rg_header(rg_string) for rg_string in args.read_groups],
            'PG': [{'ID': 'make_unaligned_bam', 'PN': '10X longranger/make_unaligned_bam', 'VN': version}],
            'CO': ['10x_bam_to_fastq:R1(RX:QX,TR:TQ,SEQ:QUAL)', '10x_bam_to_fastq:R2(SEQ:QUAL)', '10x_bam_to_fastq:I1(BC:QT)']
    }


    if args.output_format == "bam":
        out_bam = pysam.Samfile(outs.barcoded_unaligned, mode='wb', header=header)
        out_fastq = None

    elif args.output_format == "fastq":
        out_fastq = open(outs.barcoded, 'w')
        out_bam = None

    else:
        martian.exit("MAKE_UNALIGNED_OUPUT: invalid output format: '%s'" % args.output_format)

    def wfq(head, seq, qual):
        out_fastq.write(head)
        out_fastq.write("\n")
        out_fastq.write(seq)
        out_fastq.write("\n+\n")
        out_fastq.write(qual)
        out_fastq.write("\n")

    # Open FASTQ input chunk
    proc = subprocess.Popen(["gunzip", "--stdout", args.read_chunk], stdout=subprocess.PIPE)
    reader = proc.stdout

    num_pairs = 0
    correct_bc_pairs = 0

    for (bc, fields) in read_clusters(reader):

        (head, seq1, q1, seq2, q2, rx ,bc, bcq, si, siq, trim, trim_qual) = fields
        head_parts = head.split(" ")
        qname = head_parts[0]
        rg = head_parts[-1]

        tags1 = [
            ('RG', str(rg)),
            (SAMPLE_INDEX_TAG, si),
            (SAMPLE_INDEX_QUAL_TAG, siq)]
        tags2 = [
            ('RG', str(rg)),
            (SAMPLE_INDEX_TAG, si),
            (SAMPLE_INDEX_QUAL_TAG, siq)]

        if len(trim) > 0:
            tags1.append((TRIM_TAG, str(trim)))
            tags1.append((TRIM_QUAL_TAG, str(trim_qual)))

        num_pairs += 1

        if bc:
            tags1.append((PROCESSED_BARCODE_TAG, bc))
            tags2.append((PROCESSED_BARCODE_TAG, bc))
            correct_bc_pairs += 1
        tags1.append((RAW_BARCODE_TAG, rx))
        tags1.append((RAW_BARCODE_QUAL_TAG, bcq))
        tags2.append((RAW_BARCODE_TAG, rx))
        tags2.append((RAW_BARCODE_QUAL_TAG, bcq))

        if out_bam is not None:
            # Read 1
            a = pysam.AlignedRead()
            a.qname = qname
            a.seq = seq1
            a.qual = q1

            # Unmapped R1
            a.is_unmapped = True
            a.is_read1 = True
            a.tid = -1
            a.pos = -1
            a.mapq = 0
            a.cigar = [(4, len(seq1))]

            a.mrnm = -1
            a.mpos = -1
            a.tlen = -1

            a.tags = tags1

            # Read 2
            b = pysam.AlignedRead()
            b.qname = qname
            b.seq = seq2
            b.qual = q2

            b.is_unmapped = True
            b.is_read2 = True
            b.tid = -1
            b.pos = -1
            b.mapq = 0
            b.cigar = [(4, len(seq2))]

            b.mrnm = -1
            b.mpos = -1
            b.tlen = -1

            b.tags = tags2

            out_bam.write(a)
            out_bam.write(b)

        if out_fastq is not None:
            header = qname
            if bc:
                bc_header = "%s:Z:%s" % (PROCESSED_BARCODE_TAG, bc)
                header = header + " " + bc_header

            wfq(header, seq1, q1)
            wfq(header, seq2, q2)

    if out_bam is not None:
        out_bam.close()

    if out_fastq is not None:
        out_fastq.close()

    outs.num_pairs = num_pairs
    outs.correct_bc_pairs = correct_bc_pairs

def join(args, outs, chunk_defs, chunk_outs):
    if args.output_format == 'bam':
        tenkit.bam.concatenate(outs.barcoded_unaligned, [c.barcoded_unaligned for c in chunk_outs])
        outs.barcoded = None

    elif args.output_format == 'fastq':
        fqs = [c.barcoded for c in chunk_outs]
        subprocess.check_call('cat ' + ' '.join(fqs) + ' | bgzip -c > ' + outs.barcoded, shell=True)
        outs.barcoded_unaligned = None

    # Make a basic set of metrics
    num_pairs = sum(c.num_pairs for c in chunk_outs)
    correct_bc_pairs = sum(c.correct_bc_pairs for c in chunk_outs)

    stats = {}
    stats['num_read_pairs'] = num_pairs
    stats['bc_on_whitelist'] = tk_stats.robust_divide(float(correct_bc_pairs), num_pairs)

    if args.bc_counts is not None:
        # Load the bc counts for this GEM group
        counts = json.load(open(args.bc_counts, 'r'))
        count_arrays = [ np.array(gem_group['bc_counts'], dtype=np.float) for gem_group in counts.values() ]

        # Compute effective BC diversity and n90 bc count
        bc_df = pandas.DataFrame({'bc_num_reads': np.concatenate(count_arrays)})

        # Read-based effective diversity
        reads = bc_df.bc_num_reads.values
        sum_sq = (reads**2.0).sum()
        effective_diversity = tk_stats.robust_divide((reads.sum()**2.0), float(sum_sq))
        stats['barcode_diversity'] = effective_diversity
    else:
        stats['barcode_diversity'] = None


    basic_stats = pandas.DataFrame(stats, index=[0])
    basic_stats.to_csv(outs.basic_stats, index=False)
