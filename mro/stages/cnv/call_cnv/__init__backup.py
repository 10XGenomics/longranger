# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
""" this is an old version without calculating metrics for background sanity check """


import math
import csv

__MRO__ = """
stage CALL_CNV(
    in  csv       rPM,
    in  int       maxUnknownGap,
    in  int       maxNormalGap,
    in  int       maxHomoGap,
    in  string    wgsmode,
    out bed       het_del,
    src py       "stages/call_cnv",
) split using ()
"""

def testMissingHP(info):
    # return two values:
    # +1 yes, -1 no, 0, no coverage
    # 1 haploid 1 is missing
    # 2 haploid 2 is missing
    # 0 both are missing
    if info[2]==info[3] and info[4]==info[5]:
        return [0,0]
    elif info[2]==info[3] and info[4]!=info[5]:
        return [1,1]
    elif info[2]!=info[3] and info[4]==info[5]:
        return [1,2]
    elif info[2]!=info[3] and info[4]!=info[5]:
        return [-1,0]

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    with open(args.rPM,"rb") as f:
        reader = csv.reader(f)
        data_list = list(reader)

    data_by_chrom = {}
    for r in data_list:
        chrom, start, end, exonID, ID = r[:5]
        start = int(start)
        end = int(end)
        ID = int(ID)
        info = [int(a) for a in r[5:]]
        if chrom not in data_by_chrom:
            data_by_chrom[chrom]={}
        data_by_chrom[chrom][ID]=(start, end, exonID, info)

    het_del=open(outs.het_del,'w')
    print outs.het_del
    #het_del_data=open(outs.het_del_data,'w')

    #numCNV = 0

    sorted_chrom = sorted(data_by_chrom.keys())
    for chrom in sorted_chrom:
        print chrom
        r = data_by_chrom[chrom]
        IDs = sorted(r.keys())
        ln = len(IDs)

        idx = 0
        while idx < ln:
            # find stretch of potential het_del
            info = r[IDs[idx]][3]
            hasMissing, whichMissing = testMissingHP(info)
            #print hasMissing, whichMissing
            
            if hasMissing == 1:
                # initialize a het_del test
                id_list = [IDs[idx]]
                hp_list = [whichMissing]
                lastHet = r[IDs[idx]][1] # the end of the last bin showing sign of het del
                lastMissing = whichMissing
                MergedStart = r[IDs[idx]][0]
                MergedEnd = r[IDs[idx]][1]

                idx += 1
                # check whether there is homo del bins in front of it.
                # DO it later. No impoact on final calculation
                #idx_old = idx-1
                #while idx_old >=0:
                #    hasM, whichM = testMissingHP(r[IDs[idx_old]][3])
                #    if hasM !=0 or whichM != 0 and r[IDs[idx]][1]:
                #        break
                

                while idx < ln:
                    info = r[IDs[idx]][3]
                    hasMissing, whichMissing = testMissingHP(info)
                    #print hasMissing, whichMissing, lastMissing, r[IDs[idx]][1], r[IDs[idx]][0] - r[IDs[idx-1]][1]
                    
                    # when the region is extended:
                    #    1. next bin is the same type and the gap is smaller than 200 bp
                    #    2. reads on both haps for the next bin, next bin is adjacenct to the previous one, the strech is no larger than 200 bp
                    #    3. "homo" type bin, gap is no larger than 200 bp, the strech can not be largert than 1kb??
                    #    4. NO extension when the missing hap got swtiched, say from hap1 to hap2. 
                    #
                    #    a strech starts from the end pos of the last het  lastHet
                    #
                    #    if normal and homo bins alternate, whichever breaks the condition break the region continuation
                    
                    gapToPre = r[IDs[idx]][0] - r[IDs[idx-1]][1]
                    if gapToPre > args.maxUnknownGap: break
                    
                    if args.wgsmode == "wgs":
                        if hasMissing == 1: # het del
                            if whichMissing != lastMissing:
                                break
                        elif hasMissing == 0: # homo
                            if r[IDs[idx]][1] - lastHet > args.maxHomoGap:
                                break
                        elif hasMissing == -1: # none
                            if r[IDs[idx]][1] - lastHet > args.maxNormalGap:
                                break
                    else:
                        if abs(int(r[IDs[idx]][2])-int(r[IDs[idx-1]][2])) > 10:
                            break

                    if hasMissing ==1:
                        id_list.append(IDs[idx])
                        hp_list.append(whichMissing)
                        lastHet = r[IDs[idx]][1]
                    if hasMissing >=0:
                        MergedEnd = r[IDs[idx]][1]

                    idx += 1

                # trim off if the last several bins with reads on both
                backward = len(hp_list) -1
                while hp_list[backward]==-1:
                    backward -= 1
                id_list=id_list[:backward+1]
                hp_list=hp_list[:backward+1]

                # calculate the merged info
                MergedEnd=r[id_list[len(id_list)-1]][1]
                totalMol, totalMolNoRead, MissingMol, WildMol, WildMolNoRead, ttlCoveredRegion = 0,0,0,0,0,0
                for i in range(len(id_list)):
                    ID = id_list[i]
                    whichMissing = hp_list[i]
                    info = r[ID][3]
                    if whichMissing == 0: continue
                    else:
                        totalMol += info[0]
                        totalMolNoRead += info[1]
                        MissingMol += info[whichMissing*2]
                        WildMol += info[(3-whichMissing)*2]
                        WildMolNoRead += info[(3-whichMissing)*2+1]
                        ttlCoveredRegion += (r[ID][1]-r[ID][0])

                if totalMol ==0: continue
                p = WildMolNoRead*1.0/WildMol
                binom = math.pow(p,MissingMol)
                if binom > 1e-2: continue
                #het_del.write("%s\t%d\t%d\t1.0\t%.3e\t%.2e\t%d\t%d\t%d\t%d\n" %
                #        (chrom, MergedStart, MergedEnd, binom, (WildMol-WildMolNoRead)*100.0/ttlCoveredRegion,totalMol, totalMolNoRead, MissingMol, WildMol))
                het_del.write("%s\t%d\t%d\t1.0\t%.3e\t%.2e\t%d\t%d\t%d\t%d\t%d\n" %
                        (chrom, MergedStart, MergedEnd, binom, (WildMol-WildMolNoRead)*100.0/ttlCoveredRegion,ttlCoveredRegion,totalMol, totalMolNoRead, MissingMol, WildMol))
                #for ID in id_list:
                #    het_del.write("%s\t%d\t%d\t1.0\t%.3e\t%d\t%d\t%d\t%d\t%d\n" %
                #            (chrom, r[ID][0], r[ID][1], binom, WildMolNoRead, totalMol, totalMolNoRead, MissingMol, WildMol))
            idx+=1

    het_del.close()
