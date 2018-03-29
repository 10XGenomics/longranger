# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
from scipy.stats import poisson
import math
import longranger.cnv.rpm_utils as ru
from  shutil import copyfileobj
import os

__MRO__ = """
stage HOMO_DEL(
    in  string[] rPMFiles,
    in  int   minNumExons,
    in  txt   covStat,
    out bed   homo_del,
    src py    "stages/homo_del",
) split using ()
"""

LENGTH_CAP=2000000
MAX_NREAD=2
MAX_PVAL=1.0e-5


class ExonRecord(object):
    def __init__(self, records, MolPerSZ):
        self.exonID = records[0].exonID
        self.chrom = records[0].chrom
        isHomoDel = True
        minStart = None
        maxEnd  = None
        self.Nread = 0.0

        for r in records:
            self.Nread += r.Nread
            if not r.isHomoDel:
                isHomoDel = False
            if minStart is None:
                minStart = r.start
                maxEnd = r.end
            else:
                minStart = min(minStart, r.start)
                maxEnd = max(maxEnd, r.end)

        self.isHomoDel = isHomoDel
        self.start = minStart
        self.end = maxEnd
        self.pval = 1
        if self.isHomoDel: self.pval=1e-99

        if (not self.isHomoDel) and (self.exonID%10==5) and (self.Nread <=MAX_NREAD):
            ExpectedReads = int(math.ceil((self.end-self.start)*MolPerSZ*0.5*0.5))
            pval = poisson.cdf(int(math.ceil(self.Nread)), ExpectedReads)

            if pval < MAX_PVAL:
                self.isHomoDel = True
                self.pval=pval
                print self.Nread, self.isHomoDel, self.end-self.start, self.start, self.end
                print ExpectedReads, pval, self.exonID
                print


def split(args):
    # split by rPMFiles as the result of the STAGE GET_RPMS
    mem_gb = 18
    if args.rPMFiles:
        files = [{"rPM":rpm, "__mem_gb": mem_gb} for rpm in args.rPMFiles]
        return {'chunks': files}
    else:
        return {'chunks':[{"rPM":None}]}



def main(args, outs):

    # return None if no input rPMFiles
    if not args.rPMFiles:
        outs.homo_del = None
        return

    # dict to carry all rPM information
    all_rPM = {}

    # read the coverage statistics from args.covStat
    MolPerSZ=0
    with open(args.covStat) as fin:
        fin.readline()
        line = fin.readline()
        infos = line.strip().split()
        MolPerSZ = float(infos[3])

    # load rPM info to all_rPM
    with open(args.rPM) as fin:
        for l in fin:
            record = ru.rPMRecord(l)
            chrom = record.chrom
            exonID = record.exonID
            if chrom not in all_rPM:
                all_rPM[chrom]={}
                print "reading chromosome "+chrom
            if exonID in all_rPM[chrom]:
                all_rPM[chrom][exonID].append(record)
            else:
                all_rPM[chrom][exonID]=[record]

    # get all exon_ids by chromosome
    all_exon_ids = {}
    for chrom in all_rPM:
        all_exon_ids[chrom] = all_rPM[chrom].keys()
        all_exon_ids[chrom].sort()

    all_exons = {} # event detection at the exon level
    numExon = 0
    numExonHomoDel = 0
    numExonBigHomoDel = 0

    for chrom in all_exon_ids:
        all_exons[chrom]={}
        numExon += len(all_exon_ids[chrom])
        for id in all_exon_ids[chrom]:
            # merge all bins coming from the same exon id into a single bin
            er = ru.MergedBin(all_rPM[chrom][id], "homo_del", MolPerSZ,
                              MAX_NREAD, MAX_PVAL)
            all_exons[chrom][id]=er
            if er.isHomoDel: numExonHomoDel +=1

    # merge nerighbin bins with homo deletions into a single merged one
    # neighboring exons with homo del signal will be merged together, regardless how far it is apart
    with open(outs.homo_del, "w") as fout:
        for chrom in all_exons:
            if chrom in ["chrY","chrX", "Y", "X"]:
                continue
            Nexons = len(all_exon_ids[chrom])

            print "Nexons", Nexons

            idx =0
            while idx < Nexons:
                if all_exon_ids[chrom][idx]%10 == 5 and \
                    all_exons[chrom][all_exon_ids[chrom][idx]].isHomoDel:


                    numHetExon = 1
                    pval = all_exons[chrom][all_exon_ids[chrom][idx]].HomoDelPVal
                    idx0 = idx
                    Nread = all_exons[chrom][all_exon_ids[chrom][idx]].Nread
                    idx +=1

                    while idx < Nexons:
                        # search stops when current exon is not a homo del
                        if not all_exons[chrom][all_exon_ids[chrom][idx]].isHomoDel: break
                        if all_exon_ids[chrom][idx]%10 == 5:
                            numHetExon += 1
                            pval *= all_exons[chrom][all_exon_ids[chrom][idx]].HomoDelPVal
                            Nread += all_exons[chrom][all_exon_ids[chrom][idx]].Nread
                        idx += 1

                    #print pval
                    exons_this_call = []
                    for ii in range(idx0, idx):
                        exons_this_call.append(all_exons[chrom][all_exon_ids[chrom][ii]])
                    merged_bin_this_call = ru.MergedBin(exons_this_call, "homo_del", MolPerSZ)

                    event_length = all_exons[chrom][all_exon_ids[chrom][idx-1]].end - all_exons[chrom][all_exon_ids[chrom][idx0]].start


                    if numHetExon >= args.minNumExons and event_length <= LENGTH_CAP:
                        numExonBigHomoDel += 1
                        fout.write("%s\t%d\t%d\tPASS\t%.4g\t%.1f\t%d\n" % (chrom, all_exons[chrom][all_exon_ids[chrom][idx0]].start, \
                            all_exons[chrom][all_exon_ids[chrom][idx-1]].end, pval, Nread, merged_bin_this_call.ttlMol))

                idx+=1



def join(args, outs, chunk_defs, chunk_outs):
    homoDelFiles = [out.homo_del for out in chunk_outs if out.homo_del and os.path.exists(out.homo_del) ]
    with open(outs.homo_del,"w") as fout:
        for infile in homoDelFiles:
            with open(infile, "rb") as fin:
                copyfileobj(fin, fout)
