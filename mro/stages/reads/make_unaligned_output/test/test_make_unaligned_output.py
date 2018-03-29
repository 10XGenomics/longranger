#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
import tenkit.test as tk_test
from .. import *
import pysam

import martian

martian.test_initialize("")

# Patrick Marks
# Simple test of deduper

IN_FASTQ = tk_test.in_path('test_bc_sorted_big_fastq.fastq.gz')
OUT_BAM = tk_test.out_path('test_unaligned_out.bam')


class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_make_unaligned(self):
        args = martian.Record({'sample_id': 1234, 'output_format': "bam", 'read_group': "RG", 'read_chunk': IN_FASTQ})
        outs = martian.Record({'barcoded_unaligned': OUT_BAM})
        main(args, outs)

        out_bam = pysam.Samfile(OUT_BAM, check_sq=False)
        reads = list(out_bam)

        assert(len(reads) == 2000)
