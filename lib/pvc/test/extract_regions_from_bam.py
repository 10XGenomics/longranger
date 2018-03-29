#!/usr/bin/env python

import pysam
import sys

bam = pysam.Samfile(sys.argv[1])
bam_out = pysam.Samfile(sys.argv[2], 'wb', template=bam)

for l in open(sys.argv[3]):
    chrom, start, end = l.split('\t')
    for r in bam.fetch(chrom, int(start), int(end)):
        bam_out.write(r)

bam.close()
bam_out.close()
