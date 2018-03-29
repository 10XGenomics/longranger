# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
from collections import namedtuple
from longranger.cnv.event_bed import SimpleEvent
import tenkit.bio_io as tk_io
import tenkit.regions as tk_regions
from tenkit.stats import robust_divide

OverlapInfo = namedtuple("Info", "isOverlap size queryFraction databaseFraction sizeDB numEntries name")


class TruthStatus(object):
    def __init__(self, iT = None, iF =None ):
        self.isTrue=iT
        self.isFalse=iF


class Status(object):
    def __init__(self, td=[], trackInfo=[], iP=None):
        self.TDStatus = td
        self.trackInfo = trackInfo
        self.isPass=iP


def read_bed_file(bed):
    # output chrom_idx_truth is a chr-indexed dictionary with tk_regions.Regions value
    with open(bed) as f:
        bed_dict = tk_io.get_target_regions(f)
    return bed_dict


def find_overlap(chrom, start, end, database_bed, trackName):
    if start > end: raise Exception("errRegionForamt","the stop position is smaller than the start position "+" ".join([start,end]))
    # database_bed is the object returned by read_database_bed
    # return
    #    1. True or False for finding or not finding
    #    2. Total overlapping base pairs
    #    3. The percentage of overlap of the query
    #    4. The percentage of overlap of the largest region in database
    #    5. Number of regions overlapped

    if chrom not in database_bed:
        return OverlapInfo(False, 0, 0, 0, 0, 0, trackName)
    else:
        overlapping_regions = database_bed[chrom].overlapping_regions(start,end)
        if len(overlapping_regions) == 0:
            return OverlapInfo(False, 0, 0, 0, 0, 0, trackName)
        region_sizes = [r[1]-r[0] for r in overlapping_regions]
        overlapping_sizes = [min(end,r[1])-max(start,r[0]) for r in overlapping_regions]
        overlapping_fractions=[robust_divide(o*1.0, s) for s, o in zip(region_sizes, overlapping_sizes)]
        total_overlap_size=sum(overlapping_sizes)
        fraction_as_query=robust_divide(total_overlap_size*1.0, (end-start))
        #print region_sizes, overlapping_fractions
        return OverlapInfo(total_overlap_size>0, total_overlap_size, fraction_as_query, max(overlapping_fractions), max(region_sizes), len(overlapping_sizes), trackName)



class TruthSet(object):
    def __init__(self, bedfile, name="Truth", regionsAvoided = {}):
        self.content = read_bed_file(bedfile)
        self.found = {}
        self.name=name


        if len(regionsAvoided) > 0:
            tracks = {}
            for name, t  in regionsAvoided.iteritems():
                tracks[name] = read_bed_file(t)
            old_content = self.content
            self.content = {}
            for chrom, rgs in old_content.iteritems():
                if not chrom in self.content:
                    self.content[chrom] = tk_regions.Regions()
                for st, en in zip(rgs.starts, rgs.ends):
                    isPass= True
                    for name, t in tracks.iteritems():
                        if find_overlap(chrom, st, en, t, name)[0]:
                            isPass = False
                            break
                    if isPass:
                        self.content[chrom].starts.append(st)
                        self.content[chrom].ends.append(en)

        for chrom, regions in self.content.iteritems():
            for st, en in regions:
                key="_".join([chrom, str(st), str(en)])
                self.found[key]=0

    def size(self):
        return len(self.found)

    def getNumConfirmed(self):
        cnt = 0
        for k, v in self.found.iteritems():
            if v >0: cnt+=1
        return cnt

    def getSensitivity(self):
        return robust_divide(self.getNumConfirmed()*1.0, len(self.found))

    def resetFound(self):
        for k in self.found:
            self.found[k]=0

    def checkOverlap(self, chrom, start, end):
        if start > end: raise Exception("errRegionForamt","the stop position is smaller than the start position "+" ".join([start,end]))
    # database_bed is the object returned by read_database_bed
    # return
    #    1. True or False for finding or not finding
    #    2. Total overlapping base pairs
    #    3. The percentage of overlap of the query
    #    4. The percentage of overlap of the largest region in database
    #    5. Number of regions overlapped

        if chrom not in self.content:
            return OverlapInfo(False, 0, 0, 0, 0, 0, self.name)
        else:
            overlapping_regions = self.content[chrom].overlapping_regions(start,end)
            if len(overlapping_regions) == 0:
                return OverlapInfo(False, 0, 0, 0, 0, 0, self.name)
            for r in overlapping_regions:
                key = "_".join([chrom, str(r[0]), str(r[1])])
                self.found[key]=1
            region_sizes = [r[1]-r[0] for r in overlapping_regions]
            overlapping_sizes = [min(end,r[1])-max(start,r[0]) for r in overlapping_regions]
            overlapping_fractions=[robust_divide(o*1.0, s) for s, o in zip(region_sizes, overlapping_sizes)]
            total_overlap_size=sum(overlapping_sizes)
            fraction_as_query=robust_divide(total_overlap_size*1.0, (end-start))
        #print region_sizes, overlapping_fractions
            return OverlapInfo(total_overlap_size>0, total_overlap_size, fraction_as_query, \
                max(overlapping_fractions), max(region_sizes), len(overlapping_sizes), self.name)

    def getMissed(self):
        return [k for k in self.found if self.found[k]==0]




class Calls(object):
    def __init__(self, het_del, ts, genomeTracks={}, overlapThr={}):
        # genomeTracks is a map of various genomeTracks
        self.genomeTracks = {}
        self.genomeTracksOverlapThr = {}
        if len(overlapThr)==0:
            for name in genomeTracks:
                self.genomeTracksOverlapThr[name] = 0.0
        else:
            self.genomeTracksOverlapThr = overlapThr

        #print self.genomeTracksOverlapThr

        for name, f in genomeTracks.iteritems():
            self.genomeTracks[name]=read_bed_file(f)

        self.AllEvents = []
        self.AllRawData = []
        self.TruthData = ts
        self.NumTD = len(self.TruthData) if self.TruthData is not None else 0
        self.NPos = [0]*self.NumTD
        self.Status = []
        with open(het_del) as f:
            for l in f:
                evt = SimpleEvent(l)
                if evt.isPass:
                    self.AllEvents.append(evt)
                    self.AllRawData.append(l.strip())
                    TDStatus = [TruthStatus() for i in range(self.NumTD)]
                    trackOverlapMap = {}
                    for name, track in self.genomeTracks.iteritems():
                        trackOverlapMap[name]=find_overlap(evt.chrom, evt.start, evt.end, track, name)
                    self.Status.append(Status(TDStatus, trackOverlapMap))

        self.TotalEvent = len(self.AllEvents)
        self.AllPos = 0
        self.HasPerformance = False

        # variables not initialized here



    def clearStatus(self):
        for t in self.TruthData:
            t.resetFound()
        for i in range(self.TotalEvent):
            for j in range(self.NumTD):
                self.Status[i].TDStatus[j] = TruthStatus(None,None)
        self.NPos = [0]*self.NumTD
        self.AllPos = 0
        self.HasPerformance = False


    def checkCallPerformance(self, SENIDX=0, PPVIDX=2, trackAvoided=[], overlapThr={}):
        self.clearStatus()
        overlapThrsGood = True
        for n in self.genomeTracks:
            if not n in overlapThr:
                overlapThrsGood = False
                break
            if overlapThr[n] <0.0 or overlapThr[n] > 1.0:
                overlapThrsGood = False
                break

        if not overlapThrsGood:
            overlapThr = self.genomeTracksOverlapThr


        validTrackAvoided=[]
        for n in trackAvoided:
            if n in self.genomeTracks:
                validTrackAvoided.append(n)


        for i in range(self.TotalEvent):
            evt = self.AllEvents[i]
            isPassTrackFilter = True
            for n in self.genomeTracks:
                if self.Status[i].trackInfo[n].queryFraction > overlapThr[n] :
                    isPassTrackFilter = False
                    break

            if isPassTrackFilter:
                self.Status[i].isPass = True

                self.AllPos+=1
                for j in range(self.NumTD):
                    if self.TruthData[j].checkOverlap(evt.chrom, evt.start, evt.end)[0]:
                        self.NPos[j]+=1
                        self.Status[i].TDStatus[j].isTrue = True
                        self.Status[i].TDStatus[j].isFalse = False
                    else:
                        self.Status[i].TDStatus[j].isTrue = False
                        self.Status[i].TDStatus[j].isFalse = True
            else:
                self.Status[i].isPass = False
        self.HasPerformance = True

        #print len(self.AllEvents), numFailDMAS, numFailTrack, numFail

        if self.AllPos <= 0:
            return self.TruthData[SENIDX].getSensitivity(), 0.0, self.AllPos
        else:
            return self.TruthData[SENIDX].getSensitivity(), robust_divide(self.NPos[PPVIDX]*1.0, self.AllPos), self.AllPos



    def getFalsePos(self, whichTD=1, num=20, outfile=None):
        if not outfile is None:
            fout = open(outfile,"w")
            num=9999999

        if not self.HasPerformance:
            raise Exception("Err","has not beeen compared to a truth set yet")

        cnt = 0
        for i in range(self.TotalEvent):
            if self.Status[i].isPass and self.Status[i].TDStatus[whichTD].isFalse:
                cnt+=1
                if not outfile is None:
                    fout.write(self.AllRawData[i]+"\n")
                else:
                    print self.AllRawData[i]
                if cnt >= num:
                    break

        if not outfile is None:
            fout.flush()
            fout.close()


    def getFalseNegative(self, whichTD=0, num=20, outfile=None):
        if not outfile is None:
            fout = open(outfile,"w")
            num=9999999

        if not self.HasPerformance:
            raise Exception("Err","has not beeen compared to a truth set yet")

        td = self.TruthData[whichTD]
        allKeys = td.found.keys()
        allKeys.sort()
        cnt = 0
        for k in allKeys:
            v=td.found[k]
            if v==0:
                cnt+=1
                rg = k.replace("_","\t")
                if not outfile is None:
                    fout.write(rg+"\n")
                else:
                    print rg
                if cnt >= num:
                    break

        if not outfile is None:
            fout.flush()
            fout.close()

