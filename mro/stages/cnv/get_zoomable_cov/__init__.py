# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import numpy as np
import math
import tenkit.bam as tk_bam
import tenkit.constants
import tenkit.bio_io as tk_io
import pysam
import pyBigWig
import pandas as pd

############################################################
# Viewing WGS bam file with IGV at large scale is          #
# challenging. This stage will generate zoomable coverage  #
# track for the two haploids and the unphased one.         #
# Meanwhile, discretization scheme is used to compress the #
# file size. The discretization scheme generates bins of   #
# size sqrt for the coverage level x. It would reduce the  #
# Poisson Process noise while compress the data amount to  #
# ~ 4%.                                                    #
#                                                          #
# In calculating the coverage, the total number of DNAs    #
# with overlapping reads is computed. Compared to the      #
# number of covering reads, this measure avoids            #
# overcounting if the covering DNA is covered by multiple  #
# overlapping reads at a given location.                   #
############################################################

### This stage does not work, unless we find a good zoomable
### coverage foramt without a license issue

__MRO__ = """
stage GET_COV_BEDGRAPH(
    in  bam         possorted_bam,
    in  tsv.gz      fragment_phasing,
    out bw          hp_read_1,      # hp 1 read count
    out bw          hp_read_2,      # hp 2 read count
    out bw          hp_read_0,      # no phasing read count
    out bw          hp_read_t,      # total read count
    out bw          hp_bc_1,        # hp 1 bc count
    out bw          hp_bc_2,        # hp 2 bc count
    out bw          hp_bc_0,        # no phasing bc count
    out bw          hp_bc_t,        # total bc count
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
    numBins = int(round(InvHarmonic(n * 1.0 / Amp))) + 3
    bins = [int(round(Harmonic(x) * Amp)) for x in range(1, numBins)]
    if bins[0]==1:
        bins=[0]+bins
    else:
        bins=[0,1]+bins
    return bins


Bins = getBins(MaxRead)
Values = [0] + [int(round((Bins[i-1]+Bins[i])/2.0)) for i in range(1,len(Bins))]


def discretize(num):
    if num > MaxRead:
        num = MaxRead

    idx = int(round(InvHarmonic(num * 1.0 / Amp))) + 2
    if num > Bins[idx]:
        print "Wrong! idx in discreteize function is too low", num, Bins[idx]
        return -1
    while idx > 0 and num <= Bins[idx - 1]:
        idx = idx - 1
    return Values[idx]

# discretizing
Disc=np.vectorize(discretize) # vectorize the function

#bedtools = "/mnt/home/wei/longranger/dev/install/devpipes/bedtools/v2.21.0/bedtools"


def split(args):
    input_bam = tk_bam.create_bam_infile(args.possorted_bam)

    # chromosomes should be ordered as chr1, chr11, chr12, ..., chr19, chr2
    chroms0 = input_bam.references
    chrom_lengths0 = input_bam.lengths

    chroms, chrom_lengths =(list(t) for t in zip(*sorted(zip(chroms0, chrom_lengths0))))

    loci = []
    for (chrom, length) in zip(chroms, chrom_lengths):
        start = 0
        while start + tenkit.constants.PARALLEL_LOCUS_SIZE < length:
            stop = start + tenkit.constants.PARALLEL_LOCUS_SIZE
            loci.append({'locus': tk_io.create_locus_info( chrom, start, stop), '__mem_gb': 4})
            start += tenkit.constants.PARALLEL_LOCUS_SIZE
        loci.append({'locus': tk_io.create_locus_info( chrom, start, length), '__mem_gb': 4})

    return {'chunks': loci, 'join': {'__mem_gb': 8.0}}


def main(args, outs):

    (chrom, start, stop) = tk_io.get_locus_info(args.locus)
    chrom = str(chrom)

    dtypes = {'#chrom': "str", "frag_start": "int64", "frag_end": "int64", "h0": "float64", "h1": "float64"}
    frags = pd.read_csv(args.fragment_phasing, sep="\t", compression='gzip',\
       usecols=["#chrom","frag_start","frag_end","h0","h1"], dtype=dtypes)
    frags = frags[(frags["#chrom"]==chrom) & (frags["frag_end"]>=start) & (frags["frag_start"]<=stop)]
    frags["hp"]=0
    frags.loc[frags["h0"] >= 0.95, 'hp'] = 1
    frags.loc[frags["h1"] >= 0.95, 'hp'] = 2
    del frags["h0"]
    del frags["h1"]

    # further split each chunk into regions of 100kB
    regionSize = 10000
    regionStarts = range(start, stop, regionSize)
    regionEnds = [x for x in regionStarts[1:]]
    regionEnds.append(stop)


    # read the bam file
    samfile = pysam.Samfile(args.possorted_bam, "rb")
    fouts = [ [ open(outs.hp_read_0.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_read_1.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_read_2.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_read_t.strip(".bw")+".bedGraph", 'w')],
              [ open(outs.hp_bc_0.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_bc_1.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_bc_2.strip(".bw")+".bedGraph", 'w'),
                open(outs.hp_bc_t.strip(".bw")+".bedGraph", 'w')]
            ]

    # convert h5 file to filtered csv file
    def filter_func(df):
        return (df["chrom"] == chrom) & (df["start_pos"] <= stop) & (df["end_pos"] >= start)


    # work with small region at a time to avoid large memory
    for rStart, rEnd in zip(regionStarts, regionEnds):
        frags2=frags[(frags["frag_end"]>=rStart) & (frags["frag_start"]<=rEnd)]
        coverage   = [np.zeros((4,rEnd-rStart)), np.zeros((4,rEnd-rStart))]
        #bc_2_phase = {}
        print "\n", rStart, rEnd
        # initialize the coverage track

        ## read count
        for r in samfile.fetch(chrom, rStart, rEnd):
            if not r.is_proper_pair: continue
            if r.is_duplicate or (not r.is_paired) or r.is_qcfail or \
                r.is_secondary or r.is_unmapped or r.mapq<30 or r.tid != r.rnext or \
                abs(r.pos - r.pnext)>5000:
                continue

            tags = dict(r.tags)
            if not "MI" in tags: continue
            hp = 0
            if "HP" in tags:
                hp = tags["HP"]

            s = max(rStart, r.pos)
            e = min(rEnd,r.aend)
            if s >= e: continue
            #print max(rStart, r.pos), min(rEnd,r.aend)
            coverage[0][hp][(s-rStart):(e-rStart)]+=1
            coverage[0][3][(s-rStart):(e-rStart)]+=1
            #bc = tags["BX"]
            #if not bc in bc_2_phase:
            #    bc_2_phase[bc] = hp

        ## bc count
        for _, row in frags2.iterrows():
            s = max(rStart, row["frag_start"])
            e = min(rEnd, row["frag_end"])

            if s >= e: continue;
            coverage[1][3][(s-rStart):(e-rStart)]+=1
            coverage[1][row["hp"]][(s-rStart):(e-rStart)]+=1

        # discretization and then print out in the bedGraph format
        for kind in range(2): ## read and then bc counts
            for hp in range(4):
                disc_cov = Disc(coverage[kind][hp])
                #disc_cov = [ discretize(x) for x in coverage[hp] ]
                sel = (disc_cov[:-1]!=disc_cov[1:])
                print sel.sum()
                pos=np.arange(len(disc_cov))
                boundaries=np.append(0, [x+1 for x in pos[sel]])
                boundaries=np.append(boundaries, len(disc_cov))
                print disc_cov.size, sel.sum(), boundaries.size
                #print coverage[hp][:10]
                #print disc_cov[:min(10,len(disc_cov))]
                for i in range(len(boundaries)-1):
                    fouts[kind][hp].write( "%s\t%d\t%d\t%d\n" % (chrom, boundaries[i]+rStart, boundaries[i+1]+rStart, disc_cov[boundaries[i]]))
    #disc_cov = Disc(coverage) # discereitzed coverage

    for kind in range(2): ## read and then bc counts
        for hp in range(4):
            fouts[kind][hp].close()
        #bedGraph[hp] = [(chrom, boundaries[i]+start, boundaries[i+1]+start, b[boundaries[i]]) for i in range(len(boundaries)-1)]
    # write out bedGraph


def join(args, outs, chunk_defs, chunk_outs):

    input_bam = tk_bam.create_bam_infile(args.possorted_bam)
    chroms0 = input_bam.references
    chrom_lengths0 = input_bam.lengths
    header_info =[t for t in sorted(zip(chroms0, chrom_lengths0))]

    chunks_hp = [           [out.hp_read_0.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_read_1.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_read_2.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_read_t.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_bc_0.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_bc_1.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_bc_2.strip(".bw")+".bedGraph"  for out in chunk_outs],
                            [out.hp_bc_t.strip(".bw")+".bedGraph"  for out in chunk_outs]]
    outs_hp         = [     outs.hp_read_0,
                            outs.hp_read_1,
                            outs.hp_read_2,
                            outs.hp_read_t,
                            outs.hp_bc_0,
                            outs.hp_bc_1,
                            outs.hp_bc_2,
                            outs.hp_bc_t]


    for hp in range(len(chunks_hp)):
        bw = pyBigWig.open(outs_hp[hp], "w")
        bw.addHeader(header_info)
        # chromosomes should be ordered as chr1, chr11, chr12, ..., chr19, chr2
        for chk in chunks_hp[hp]:
            cov_bed = pd.read_csv(chk, sep='\t', header=None)
            bw.addEntries(chroms=cov_bed[[0]].values.astype('str'), starts=cov_bed[[1]].values,
                ends=cov_bed[[2]].values, values=cov_bed[[3]].values.astype(np.float64))

        bw.close()
        del bw
        del cov_bed

