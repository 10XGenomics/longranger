# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import csv
import os
##########################################
## v6 core final 
##########################################

badRegions = {}
MaxRegions = {}
BadRegionsByExtID = {}
fullData= {}
MinExtIDRequired = range(0,11)
PER_THR=0.20
MinMax = 0.33

excel="zero_cov_target/summaries_2016_06_06__EX_1.1_VV.csv"
sample_info={} # external_id -> lena ID
with open(excel, 'rb') as file:
    sampleInfo = csv.reader(file)
    sampleInfo.next()
    for row in sampleInfo:
        lenaID=row[1]
        externalID = row[19]
        if externalID not in sample_info:
            sample_info[externalID]=[]
        sample_info[externalID].append(lenaID)

sample_info_test ={'D270216': ['20262', '20261'], 'HCC1143': ['20257', '20258'], 'HCC1143 BL': ['20259', '20260'], 'HCC38': ['20433', '20432']}
SI = sample_info
numExIDs = len(SI)

i = 0
exIDSorted = SI.keys()
exIDSorted.sort()
print exIDSorted

for exID in exIDSorted:
    lenaID = SI[exID]
    i+=1
    print exID
    allData={}
    allCnt={}
    for s in lenaID:
        rpm = "/mnt/yard/wei/data_CNV/WES_probe_assessment/CNV_"+s+"/outs/rPM.csv"
        print rpm
        if not os.path.isfile(rpm): continue
        ttlRead = 0
        ttlSize = 0
    
        localData = {}
        with open(rpm) as f:
            for l in f:
                infos = l.split(",")
                sz = int(infos[2])-int(infos[1])
                numRead = int(infos[5]) - int(infos[6])
                ttlRead+=numRead
                ttlSize += sz
                key="_".join(infos[:3])
                localData[key]=numRead*1.0/sz
        
            mn = ttlRead*1.0/ttlSize
            for k, v in localData.iteritems():
                if not k in allData:
                    allData[k]=0
                    allCnt[k]=0
                allData[k]+=(v/mn)
                allCnt[k]+=1
                
                if k not in fullData:
                    fullData[k] = {}
                if exID not in fullData[k]:
                    fullData[k][exID]={}
                fullData[k][exID][s]=(v/mn)
           

    BadRegionsByExtID[exID]={}
    for k, v in allData.iteritems():
        perc = v / allCnt[k]
        if k not in MaxRegions:
            MaxRegions[k]=0
        if perc > MaxRegions[k]:
            MaxRegions[k] = perc
        if perc < PER_THR:
            if k not in badRegions:
                badRegions[k] = 1
            else:
                badRegions[k] += 1
        BadRegionsByExtID[exID][k]=perc

        
for thr in MinExtIDRequired:
    maxCov = 0
    fout=open("/mnt/home/wei/ipython_dir/data_result/v6_core_final_black_list/v6_core_final_black_list_"+str(thr)+"_"+str(MinMax)+"_extIDs.bed", "w")
    numBad = 0
    numHighMax = 0
    for k, v in badRegions.iteritems():
        if v>=thr and MaxRegions[k]<=MinMax: 
            numBad += 1
            maxCov += MaxRegions[k]
            if MaxRegions[k] > 0.25: 
                numHighMax += 1
            fout.write(k.replace("_","\t")+"\n")
    
    fout.close()
    print thr, numBad, maxCov/numBad, numHighMax
                
