# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.regions as tk_regions
import numpy as n
import math
import os
import tenkit.bam as tk_bam
import tenkit.constants
import tenkit.bio_io as tk_io
import pysam 

#from subprocess import Popen, PIPE


__MRO__ = """
stage GET_COV_BEDGRAPH(
    in  bam              bam_infile,
    out bedGraph.gz      hp1,
    out bedGraph.gz      hp2,
    out bedGraph.gz      hp0,
    out bedGraph.gz.tbi  hp1_idx,
    out bedGraph.gz.tbi  hp2_idx,
    out bedGraph.gz.tbi  hp0_idx,
) split using (
    in  string locus,
)
"""

MaxRead = 100000
Amp = 1.5

def Harmonic(x):
    return 2.0 / 3.0 * math.pow(x, 1.5) + math.sqrt(x) / 2.0 - 0.2


def InvHarmonic(y):
    return math.pow(1.5 * (y + 0.2), 2.0 / 3.0)


def getBins(n):
    numBins = int(round(InvHarmonic(n * 1.0 / Amp))) + 1
    bins = [int(round(Harmonic(x) * Amp)) for x in range(1, numBins)]
    if bins[0]==1:
        bins=[0]+bins
    else:
        bins=[0,1]+bins
    return bins


Bins = getBins(MaxRead)
Values = [0] + [int(round((Bins[i-1]+Bins[i])/2.0)) for i in range(1,len(Bins))]


def discretize(num):
    idx = int(round(InvHarmonic(num * 1.0 / Amp))) + 2
    if num > Bins[idx]:
        print "Wrong! idx in discreteize function is too low", num, Bins[idx]
        return -1
    while idx > 0 and num <= Bins[idx - 1]:
        idx = idx - 1
    return Values[idx]

# discretizing
Disc=n.vectorize(discretize) # vectorize the function

#bedtools = "/mnt/home/wei/longranger/dev/install/devpipes/bedtools/v2.21.0/bedtools"


def split(args):
    input_bam = tk_bam.create_bam_infile(args.bam_infile)

    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    loci = []
    for (chrom, length) in zip(chroms, chrom_lengths):
        bad_chrom = ('random' in chrom or 'U' in chrom or 'hap' in chrom)
        if bad_chrom or chrom[:3] != 'chr' or chrom == 'chrM'or chrom == 'chrY': continue
        start = 0
        while start + tenkit.constants.PARALLEL_LOCUS_SIZE < length:
            stop = start + tenkit.constants.PARALLEL_LOCUS_SIZE
            loci.append({'locus': tk_io.create_locus_info( chrom, start, stop)})
            start += tenkit.constants.PARALLEL_LOCUS_SIZE
        loci.append({'locus': tk_io.create_locus_info( chrom, start, length)})

    return {'chunks': loci, 'join': {'__mem_gb': 12.0}}


def main(args, outs):

    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    chrom = str(chrom)
    
    # further split each chunk into regions of 100kB
    regionSize = 10000
    regionStarts = range(start, stop, regionSize)
    regionEnds = [x for x in regionStarts[1:]]
    regionEnds.append(stop)

    # read the bam file
    samfile = pysam.Samfile(args.bam_infile, "rb")
    fouts = [open(outs.hp0.strip(".gz"), 'w'), open(outs.hp1.strip(".gz"), 'w'), open(outs.hp2.strip(".gz"), 'w')]

    for rStart, rEnd in zip(regionStarts, regionEnds):
        print "\n", rStart, rEnd
        barcode_reads = {}
        barcode_hp = {}
        # initialize the coverage track
        coverage = n.zeros((3,rEnd-rStart))

        for r in samfile.fetch(chrom, rStart, rEnd):
            if not r.is_proper_pair: continue
            if r.is_duplicate or not r.is_paired or r.is_qcfail or \
                r.is_secondary or r.is_unmapped or r.mapq<20 or r.tid != r.rnext or \
                abs(r.pos - r.pnext)>5000:
                continue
            
            tags = dict(r.tags)
            if not "MI" in tags: continue
            mid = tags["MI"] 
            hp = 0
            if "HP" in tags: hp = tags["HP"]
            if not mid in barcode_hp:
                barcode_hp[mid] = hp
                barcode_reads[mid] = tk_regions.Regions()
            barcode_reads[mid].add_region((r.pos, r.aend))

        # add the unique sequenced mol coverage
        for bc in barcode_reads.keys():
            regions = barcode_reads[bc]
            hp = barcode_hp[bc]
            for rgs in regions:
                rgs_start = max(rStart, rgs[0])
                rgs_end = min(rEnd, rgs[1])
                coverage[hp][rgs_start-rStart:(rgs_end-rStart)]+=1

        for hp in range(3):
            disc_cov = Disc(coverage[hp])
            #disc_cov = [ discretize(x) for x in coverage[hp] ]
            sel = (disc_cov[:-1]!=disc_cov[1:])
            print sel.sum()
            pos=n.arange(len(disc_cov))
            boundaries=n.append(0, pos[sel], len(disc_cov))
            print disc_cov.size, sel.sum(), boundaries.size
            #print coverage[hp][:10]
            #print disc_cov[:min(10,len(disc_cov))]
            for i in range(len(boundaries)-1):
                fouts[hp].write( "%s\t%d\t%d\t%d\n" % (chrom, boundaries[i]+rStart, boundaries[i+1]+rStart, disc_cov[boundaries[i]]))
    #disc_cov = Disc(coverage) # discereitzed coverage

    for hp in range(3):
        fouts[hp].close()
        #bedGraph[hp] = [(chrom, boundaries[i]+start, boundaries[i+1]+start, b[boundaries[i]]) for i in range(len(boundaries)-1)]
    # write out bedGraph


def join(args, outs, chunk_defs, chunk_outs):
    outs_hp = [ outs.hp0.strip(".gz"), outs.hp1.strip(".gz"), outs.hp2.strip(".gz") ]
    chunks_hp = [ [out.hp0.strip(".gz")  for out in chunk_outs],
        [out.hp1.strip(".gz")  for out in chunk_outs],
        [out.hp2.strip(".gz")  for out in chunk_outs]]
    
    for hp in range(3):
        cmd1 = "cat "+ " ".join(chunks_hp[hp]) + " > " + outs_hp[hp]#+"_tmp"
        os.system(cmd1)
        cmd2 = "bgzip "+outs_hp[hp]
        cmd3 = "tabix -p bed "+outs_hp[hp]+".gz"
        os.system(cmd2)
        os.system(cmd3)
        #cmd2 = "sort -k1,1 -k2,2n "+ outs_hp[hp]+"_tmp" + " > " + outs_hp[hp]
        #os.system(cmd2)
        #os.remove(outs_hp[hp]+"_tmp")
