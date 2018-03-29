#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

# Test attach bcs
import pysam
import tenkit.test as tk_test
import tenkit.bio_io as tk_io
from .. import *

import martian

IN_BAM = tk_test.in_path("attach_phasing/bam_test.bam")
IN_FRAGS = tk_test.in_path("attach_phasing/fragments_test.tsv.gz")
OUT_BAM = tk_test.out_path("test_attach_phasing.bam")

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_attach_phasing(self):

        args = martian.Record({ 'input': IN_BAM, 'fragment_phasing': IN_FRAGS, 'chunk_start': 0, 'chunk_end': 1<<32 })
        outs = martian.Record({ 'phased_possorted_bam': OUT_BAM, 'phased_possorted_bam_index': OUT_BAM + ".bai" })

        main(args, outs)

        # Ensure each read has a barcode
        out_bam = pysam.Samfile(OUT_BAM)
        bam_reads = list(out_bam)

        '''
        chr1    628490  701466  10565419        565419  711255  GTACACAGAGTGTT-1        0.9996837673    0.000316232700235       5.00029071137e-61
        chr1    628789  678258  10565419        565419  711255  CGAACTCACTCCAA-1        0.999800468729  0.000199531270515       5.00029083957e-61
        chr1    628958  726129  10565419        565419  711255  AGGCTTCATCAGAA-1        3.01287901923e-08       0.999999969871  2.50113486814e-61
        chr1    630911  726153  10565419        565419  711255  CTAAGCAGGTTTAG-1        0.998004731897  0.00199526810283        5.00029096742e-61
        '''

        def check_reads(bc, start, end, haplotype):
            for r in bam_reads:
                if tk_io.get_read_barcode(r) != bc:
                    continue

                tags = { t:v for (t,v) in r.tags }
                if r.pos >= start and r.pos < end:
                    self.assertTrue(tags.has_key('HP'))
                    self.assertEqual(tags['HP'] == haplotype)

        check_reads("AGGCTTCATCAGAA-1", 565419, 711255, 1)
        check_reads("CTAAGCAGGTTTAG-1", 565419, 711255, 0)
        check_reads("CTAAGCAGGTTTAG-1", 565419, 711255, 1)


        self.assertEqual(len(bam_reads) > 0, True)
