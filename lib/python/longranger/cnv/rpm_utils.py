# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
from sets import Set
import math
import gzip
import bisect
from scipy.stats import binom_test
from scipy.stats import poisson

def GetMID(token):
    if token=="": return []
    mids = token.split(":")
    return mids[:(len(mids)-1)]


class rPMRecord(object):
    def __init__(self, line=None):
        if line is None:
            self.chrom = None
            self.start = None
            self.end   = None
            self.exonID = None
            self.binID = None
            self.ttlMol = None
            self.ttlMolMs = None
            self.hp1Mol = None
            self.hp1MolMs = None
            self.hp2Mol = None
            self.hp2MolMs = None
            self.NakeReads = None  # unbarcoded reads
            self.MIDAll = []
            self.MIDAllSeen = []
            self.MIDHP1 = []
            self.MIDHP1Seen = []
            self.MIDHP2 = []
            self.MIDHP2Seen = []
            self.MIDs = [self.MIDAll, self.MIDAllSeen, self.MIDHP1, self.MIDHP1Seen, self.MIDHP2, self.MIDHP2Seen]
            self.DM = None
            self.AS = None
            self.size = None
            self.Nread = None

        else:
            infos = line.strip().split(",")
            self.chrom = infos[0]
            self.start = int(infos[1])
            self.end   = int(infos[2])
            self.exonID = int(infos[3])
            self.binID = int(infos[4])
            self.ttlMol = int(infos[5])
            self.ttlMolMs = int(infos[6])
            self.hp1Mol = int(infos[7])
            self.hp1MolMs = int(infos[8])
            self.hp2Mol = int(infos[9])
            self.hp2MolMs = int(infos[10])
            self.NakeReads = int(infos[11])  # unbarcoded reads
            self.MIDAll = GetMID(infos[12])
            self.MIDAllSeen = GetMID(infos[13])
            self.MIDHP1 = GetMID(infos[14])
            self.MIDHP1Seen = GetMID(infos[15])
            self.MIDHP2 = GetMID(infos[16])
            self.MIDHP2Seen = GetMID(infos[17])
            self.MIDs = [self.MIDAll, self.MIDAllSeen, self.MIDHP1, self.MIDHP1Seen, self.MIDHP2, self.MIDHP2Seen]
            self.size = self.end - self.start
            self.Nread = self.ttlMol - self.ttlMolMs

            self.DM = None
            self.AS = None

            DM = float(infos[18])
            AS = float(infos[19])
            if DM > - 990 and AS > -990:
                self.DM = DM
                self.AS = AS

            self.isHomoDel = (self.hp1Mol == self.hp1MolMs) and \
                (self.hp2Mol == self.hp2MolMs)


    def testMissingHP(self):
        # return two values:
        # +1 yes, -1 no, 0, no coverage, -2 wrong
        # 1 haploid 1 is missing
        # 2 haploid 2 is missing
        # 0 both are missing
        if self.hp1Mol == self.hp1MolMs and self.hp2Mol == self.hp2MolMs:
            return [0,0]
        elif self.hp1Mol == self.hp1MolMs and self.hp2Mol > self.hp2MolMs:
            return [1,1]
        elif self.hp1Mol > self.hp1MolMs and self.hp2Mol == self.hp2MolMs:
            return [1,2]
        elif self.hp1Mol > self.hp1MolMs and self.hp2Mol > self.hp2MolMs:
            return [-1, 0]
        else:
            return [-2,-2]



def ReadRPM(RPMFile):
    rpm_by_chrom = {}
    with open(RPMFile,"rb") as f:
        for line in f:
            r = rPMRecord(line)
            if r.chrom not in rpm_by_chrom:
                rpm_by_chrom[r.chrom]={}
            rpm_by_chrom[r.chrom][r.binID]=r

    return rpm_by_chrom


class MergedBin (rPMRecord):
    def __init__(self, binList, mode="het_del", MolPerSZ=None,
                 max_num_contm_reads=2, homo_pval_thrs = 1e-5):
        sz = len(binList)
        rPMRecord.__init__(self)

        if sz==0:
            raise Exception("No bins to be merged")

        self.chrom = binList[0].chrom
        self.exonID = binList[0].exonID
        for i in range(1,sz):
            if self.chrom != binList[i].chrom:
                raise Exception("bins not come from the same chromosome")


        self.start = min([b.start for b in binList])
        self.end = max([b.end for b in binList])
        self.size=0
        self.NakeReads = 0
        FinalMIDs = [Set(), Set(), Set(), Set(), Set(), Set()]
        DMSum, DMCnt, ASSum, ASCnt = 0.0, 0.0, 0.0, 0.0
        for b in binList:
            for store, info in zip(FinalMIDs, b.MIDs):
                store.update(info)
            self.NakeReads = max(self.NakeReads, b.NakeReads)
            if b.DM is not None and b.DM > -990:
                DMSum += b.DM
                DMCnt += 1
            if b.AS is not None and b.AS > -990:
                ASSum += b.AS
                ASCnt += 1
            self.size += b.size

        if DMCnt ==0:
            self.DM = -999
        else:
            self.DM = DMSum / DMCnt
        if ASCnt ==0:
            self.AS = -999
        else:
            self.AS = ASSum / ASCnt


        self.MIDAll, self.MIDAllSeen, \
            self.MIDHP1, self.MIDHP1Seen, \
            self.MIDHP2, self.MIDHP2Seen = \
            [sorted(list(s)) for s in FinalMIDs]
        self.ttlMol, self.ttlMolMs = len(self.MIDAll), len(self.MIDAll)-len(self.MIDAllSeen)
        self.hp1Mol, self.hp1MolMs = len(self.MIDHP1), len(self.MIDHP1)-len(self.MIDHP1Seen)
        self.hp2Mol, self.hp2MolMs = len(self.MIDHP2), len(self.MIDHP2)-len(self.MIDHP2Seen)
        self.Nread = len(self.MIDAllSeen)

        self.MissingHP=None
        self.HetDelPVal = None
        self.WildCov = None
        self.WildMol, self.WildMolMs, self.DelMol, self.DelMolMs = None, None, None, None

        self.isHomoDel = False
        self.HomoDelPVal = None

        if mode=="het_del":
            self.WildRead = max(self.hp1Mol-self.hp1MolMs, self.hp2Mol -self.hp2MolMs)
            if self.hp1Mol == self.hp1MolMs:
                p = 1.0
                if self.hp2Mol >0: p = self.hp2MolMs*1.0/self.hp2Mol
                self.HetDelPVal = math.pow(p,self.hp1MolMs)
                self.MissingHP = 1
                self.WildRead = (self.hp2Mol-self.hp2MolMs)
                self.WildMol, self.WildMolMs, self.DelMol, self.DelMolMs = \
                    self.hp2Mol, self.hp2MolMs, self.hp1Mol, self.hp1MolMs
            elif self.hp2Mol == self.hp2MolMs:
                p = 1.0
                if self.hp1Mol >0: p = self.hp1MolMs*1.0/self.hp1Mol
                self.HetDelPVal = math.pow(p,self.hp2MolMs)
                self.MissingHP = 2
                self.WildRead = (self.hp1Mol-self.hp1MolMs)
                self.WildMol, self.WildMolMs, self.DelMol, self.DelMolMs = \
                    self.hp1Mol, self.hp1MolMs, self.hp2Mol, self.hp2MolMs

        elif mode=="homo_del":
            if self.Nread == 0:
                self.isHomoDel= True
                self.HomoDelPVal = 1e-99
            elif self.exonID%10==5 and self.Nread <= max_num_contm_reads:
                ExpectedReads = int(math.ceil((self.end-self.start)*MolPerSZ*0.5*0.5))
                self.HomoDelPVal = poisson.cdf(int(math.ceil(self.Nread)), ExpectedReads)
                if self.HomoDelPVal <= homo_pval_thrs: self.isHomoDel= True


class VarDensity(object):
    def __init__(self, vcfgz):
        self.allVars= {}
        self.ttlVars = 0
        with gzip.open(vcfgz, "rb") as f:
            for l in f:
                if l[0]=="#": continue
                infos = l.split("\t",7)
                if infos[6] != "PASS": continue
                chrom = infos[0]
                pos = int(infos[1])
                if not chrom in self.allVars:
                    self.allVars[chrom] = []
                else:
                    self.allVars[chrom].append(pos)
                self.ttlVars += 1
        print "Have read %d variants" % self.ttlVars


    def getVarDen(self, chrom, start, end):
        if not chrom in self.allVars or end<=start: return None
        numVarInRegion = bisect.bisect_left(self.allVars[chrom], end) \
            - bisect.bisect_left(self.allVars[chrom], start)
        return numVarInRegion*1000.0/(end-start)



def findFlankingIDRange (minID, maxID):
    # if no flanking bins found, extend to farther regions
    if minID%10!=0: minID -= 5
    if maxID%10!=0: maxID +=5
    return (minID, maxID)

def unparityPval(r1, r2):
    bn_test = binom_test([r1,r2], p=0.5)
    #mu = (r1+r2)/2.0
    #ps=poisson(mu=mu)
    #delta = abs(r1-r2)/2.0
    #ps_test = ps.cdf(x=(mu-delta))*2
    return bn_test#, ps_test


def calBackgroundStats (bins, ids):
    # what to test?
    # is it as het del?
    # if not, does it look like a het del with contamination?
    # if not, return
    #     coverages: # reads *100/bin size
    #         total, hp1, hp2,
    #     min/max coverages:
    #         total, hp1, hp2
    #     fraction of phasing
    #         total, ph1, ph2
    #     min/max fraction of phasing
    #         total, ph1, ph2
    raw_size = sum([bins[id].end - bins[id].start for id in ids])
    ids = [id for id in ids if bins[id].ttlMol > bins[id].ttlMolMs] # only use positive ids
    if len(ids)==0: return "0\t0\t0\t0\t0\t0\t0"

    sizes = [bins[id].end - bins[id].start for id in ids]
    totalRead = [bins[id].ttlMol - bins[id].ttlMolMs for id in ids]
    hp1Read = [bins[id].hp1Mol - bins[id].hp1MolMs for id in ids]
    hp2Read = [bins[id].hp2Mol - bins[id].hp2MolMs for id in ids]
    hp0Read = [ttl-hp1-hp2 for ttl, hp1, hp2 in zip(totalRead, hp1Read, hp2Read)]
    cov_ttl, cov_hp1, cov_hp2  = [sum(read)*100.0/sum(sizes) for read in [totalRead, hp1Read, hp2Read]]
    covs_ttl, covs_hp1, covs_hp2 = [[read[i]*100.0/sizes[i] for i in range(len(ids))] for read in [totalRead, hp1Read, hp2Read]]
    phaseFrac = 1-sum(hp0Read)*1.0/sum(totalRead)

    parity=unparityPval(sum(hp1Read), sum(hp2Read))

    return "%6d\t%8.2g\t%6d\t%6d\t%6d\t%8.2f\t%8.4f" % \
        (raw_size, parity, sum(totalRead), sum(hp1Read), sum(hp2Read), cov_ttl, phaseFrac)
    #return "%d\t%d/%d/%d\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.4f/%.4f/%.4f" % \
    #    (raw_size, \
    #     sum(totalRead), sum(hp1Read), sum(hp2Read), \
    #     cov_ttl, cov_hp1, cov_hp2, \
    #     min(covs_ttl), min(covs_hp1), min(covs_hp2), \
    #     max(covs_ttl), max(covs_hp1), max(covs_hp2), \
    #     phaseFrac, min(phaseFracs), max(phaseFracs)
    #    )


def calMDAS(bins, background_id_range):
    ttlDM = 0.0
    ttlAS = 0.0
    DMCnt = 0
    ASCnt = 0
    for id in range(background_id_range[0], background_id_range[1]+1):
        if not id in bins: continue
        r = bins[id]
        numReads = r.ttlMol - r.ttlMolMs
        if numReads > 0:
            if r.DM:
                DMCnt += numReads
                ttlDM += (numReads * r.DM)
            if r.AS:
                ASCnt += numReads
                ttlAS += (numReads * r.AS)

    outString = ""
    if DMCnt > 0:
        outString += ("%10.4f" % (ttlDM/DMCnt))
    else:
        outString += ("%10.4f" % (-999.0))

    if ASCnt > 0:
        outString += ("\t%10.4f" % (ttlAS/ASCnt))
    else:
        outString += ("\t%10.4f" % (-999.0))

    return outString
