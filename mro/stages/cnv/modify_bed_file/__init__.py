#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# modify the queryregions to add near-target regions and off target regions (WES)
#
# if two neighboring target regions have exonID (targetID) 15 and 25
# then all near-target regions have exonID 20
# the off-target regions have exonDI 21.
#

#import os
#import json

__MRO__ = """
stage MODIFY_BED_FILE(
    in  bed      in_bed,
    in  int      bin,
    in  int[]    nearTargetBin,
    in  bool     addOffTarget,
    out bed      modified,
    src py       "stages/modify_bed_file",
)
"""

def split(args):
    return {'chunks': [{'__mem_gb': 6, '__threads':1}]}

def join(args, outs, chunk_defs, chunk_outs):
    outs.modified = chunk_outs[0].modified

# bed file may overlap and be adjacent to each other
def newRemains(rms, regionID, binID):
    return [regionID, binID] + rms[2:]


def main(args, outs):
    #args.coerce_strings()
    #outs.coerce_strings()

    MIN_SZ = args.bin          # minimal bin size
    MAX_SZ = MIN_SZ + 40       # maximal bin size
    NearBins = args.nearTargetBin  # nearBins
    numNearBins = len(NearBins)
    BoundaryForward = [0]*(numNearBins+1)
    BoundaryBackward = [0]*(numNearBins+1)
    for i in range(numNearBins):
        BoundaryForward[i+1]=BoundaryForward[i]+NearBins[i]
        BoundaryBackward[numNearBins-i-1]=BoundaryBackward[numNearBins-i]-NearBins[i]
    maxSizes = [-x for x in BoundaryBackward]
    numB = len(BoundaryBackward)

    def regionSplit(start, end):
        size = end-start
        if size <= MAX_SZ:
            return [(start, end)]
        else:
            fold = size / MIN_SZ +1
            if size % MIN_SZ == 0: fold -= 1
            chunk = size/fold
            rst=[]
            for i in range(fold-1):
                rst.append((start+i*chunk, start+(i+1)*chunk))
            rst.append((start+(fold-1)*chunk, end))
            return rst


    def getLeftNearTargetBins(preEnd, curStart):
        # NearBins list of near target bin sizes
        # curStart=-1 if this is the first target (region)
        #
        # return a list of off-/near-target regions, and the index corresponding to off-target bin or None

        if preEnd==-1:
            return ([(curStart+BoundaryBackward[i], curStart+BoundaryBackward[i+1]) for i in range(numB-1)], None)
        else:
            if curStart-preEnd >= 2*maxSizes[0]:
                OffTarget = (preEnd+BoundaryForward[numB-1],curStart+BoundaryBackward[0])
                return ([(preEnd+BoundaryForward[i], preEnd+BoundaryForward[i+1]) for i in range(numB-1)] + \
                    [OffTarget] + \
                    [(curStart+BoundaryBackward[i], curStart+BoundaryBackward[i+1]) for i in range(numB-1)], \
                        numB-1)
            else:
                gap = curStart-preEnd
                halfGap = gap/2.0
                for max2nd in range(len(maxSizes)):
                    if halfGap >= maxSizes[max2nd]:
                        break

                size2nd = maxSizes[max2nd]

                sizeInMid = 2*(halfGap-size2nd)
                #print gap, halfGap, max2nd, size2nd, sizeInMid
                if sizeInMid > maxSizes[max2nd-1]:
                    MidBins = [(preEnd+size2nd, preEnd+int(halfGap)),(preEnd+int(halfGap), curStart-size2nd)]
                elif curStart-size2nd > preEnd+size2nd:
                    MidBins = [(preEnd+size2nd, curStart-size2nd)]
                else:
                    MidBins = []

                return ([(preEnd+BoundaryForward[i], preEnd+BoundaryForward[i+1]) for i in range(0,numB-1-max2nd)] +\
                    MidBins +\
                    [(curStart+BoundaryBackward[i], curStart+BoundaryBackward[i+1]) for i in range(max2nd,numB-1)], \
                        None)


    def getRightNearTargetBins(curEnd):
        return [(curEnd+BoundaryForward[i], curEnd+BoundaryForward[i+1]) for i in range(numB-1)]


    # read in original bed file info
    # remains stores some ID information: exon ID and bin ID
    regions = {}
    with open(args.in_bed) as f:
        for l in f:
            if l.startswith('browser') or l.startswith('track') or l.startswith('-browser') or l.startswith('-track') or l.startswith('#'):
                continue
            info = l.split()
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            remains = info[3:]
            if chrom not in regions:
                regions[chrom]=[]
            regions[chrom].append((start, end, remains))


    new_regions = {}
    # new modified bed files, including possible near-target regions
    # here is how exon ID got transformed:
    #     id : old ID
    #     new id = id*10+5
    #     near target bin id = new id -5 and new id +5

    binID=-1
    fout = open(outs.modified, "w")

    ordered_chrom = sorted(regions.keys())
    for chrom in ordered_chrom:
        regionID=-5
        rgs = regions[chrom]
        new_regions[chrom]=[]
        new_regions_by_chr = []
        ln = len(rgs)

        for idx in range(ln):
            regionID+=10

            if len(NearBins)>0:
                # generate left nearTarget Bins
                if idx==0:
                    middleBins, idxOffTarget = getLeftNearTargetBins(-1, rgs[idx][0])
                else:
                    middleBins, idxOffTarget = getLeftNearTargetBins(rgs[idx-1][1], rgs[idx][0])
                for ii in range(len(middleBins)):
                    r=middleBins[ii]
                    if idxOffTarget and ii==idxOffTarget:
                        if args.addOffTarget:
                            binID+=1
                            new_regions_by_chr.append((r[0], r[1], newRemains(rgs[idx][2], regionID-4, binID)))
                    else:
                        binID+=1
                        new_regions_by_chr.append((r[0], r[1], newRemains(rgs[idx][2], regionID-5, binID)))
            elif len(NearBins)==0 and args.addOffTarget and idx > 0:
                binID+=1
                new_regions_by_chr.append((rgs[idx-1][1], rgs[idx][0], newRemains(rgs[idx][2], regionID-4, binID)))


            # add the split bins from the original target (region)
            for r in regionSplit(rgs[idx][0], rgs[idx][1]):
                binID+=1
                new_regions_by_chr.append((r[0], r[1], newRemains(rgs[idx][2], regionID, binID)))

            # add the right near target bins for the last target
            if idx==ln-1 and len(NearBins)>0:
                for r in getRightNearTargetBins(rgs[idx][1]):
                    binID+=1
                    new_regions_by_chr.append((r[0], r[1], newRemains(rgs[idx][2], regionID+5, binID)))

        for r in new_regions_by_chr:
            fout.write("%s\t%d\t%d\t%s\n" % (chrom, r[0], r[1], "\t".join([str(ee) for ee in r[2]])))

    fout.flush()
    fout.close()
