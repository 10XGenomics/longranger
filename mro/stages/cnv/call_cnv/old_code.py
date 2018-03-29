#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# convert fragments.h5 file into csv with selected columns
#
import pandas as pd
from scipy import stats
from scipy.stats import poisson
from scipy.special import gammaln
import numpy as np
import math

__MRO__ = """
stage CALL_CNV(
    in  csv       rPM,
    out bed       het_del,
    out txt       het_del_data,
    out bed       tandem_dup,
    out bed       reg_dup,
    src py       "stages/call_cnv",
) split using () 
"""

def postP(n,m,l):
    return math.exp(gammaln(m+l+1)+gammaln(n+2)-gammaln(l+1)-gammaln(m+n+2))


def getFreq(numzero, rPM):
    m = max(rPM)
    cnt = [0]*(m+1)
    cnt[0]=numzero
    for i in rPM:
        cnt[i]+=1

    sz = (numzero+len(rPM))*1.0
    return [x/sz for x in cnt]

def getGoodFreq(cnt, n, h):
    # cnt   the original frequncy
    # n     the desized size of the dist.
    # h     the parameter for the Poisson kernel
    p = [0]*n
    for i in range(len(cnt)):
        for j in range(n):
            p[j] += cnt[i]*poisson.pmf(j, i+h)

    return p


def probConv(p1, p2): # numpy.convolve
    n1 = len(p1)-1
    n2 = len(p2)-1
    n = n1+ n2
    p = [0.0]*(n+1)
    for i in xrange(n+1):
        for j in xrange(i+1):
            k = i-j
            if j>n1 or k>n2: continue
            p[i]+= p1[j]*p2[k]

    return p


def testProb(p1):
    p2 = np.convolve(p1, p1)
    p3 = np.convolve(p1, p2)
    p4 = np.convolve(p1, p3)
    p5 = np.convolve(p1, p4)
    return [p1,p2,p3,p4,p5]


def llhTest(pt, ps, nTest):
    # log likelihood test
    m = len(ps)

    maxL = -1e6
    idx = -1

    allL = []

    for t in range(m):
        p = ps[t]
        n=len(p)
        ll = 0.0
        for i in range(len(pt)):
            fi = 1.0/n
            if i<n and p[i]>0: fi=p[i]
            ll += (math.log(fi)*pt[i])
        ll *= (2*nTest)
        allL.append(ll)
        if ll > maxL:
            maxL = ll
            idx = t+1
       # print t+1, ll
    
    #allL.sort(reverse=True)
    #print "%.4f\t%.4f\t%.4f" % (maxL, allL[0], maxL-allL[0])
    # assume 6 or more
    return idx, maxL-allL[0]


def numLines(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    a = pd.read_json(args.rPM)
    het_del=open(outs.het_del,'w')
    het_del_data=open(outs.het_del_data,'w')
    tandem_dup=open(outs.tandem_dup,'w')
    reg_dup=open(outs.reg_dup,'w')
    sz = a.shape[0]
    # tandem duplication test

    numCNV = 0

    for i in xrange(sz):
        b = a.iloc[i]

        if i%10000==0: print i

        if b[3]==0 or b[4]<0 or b[6]==0 or b[7]<0 or b[4]>b[3] or b[7]>b[6]: continue
        l1, l2 = len(b[5]), len(b[8])

        if l1==0 and l2>0 or l1>0 and l2==0:
            # test het deletion
            p = 0.0
            binom = 0.0
            num_DNA_control = 0.0 # umber of DNA covering the non zero haploid
            #postPV = 0.0
            #size_ref = 10000
            if l1==0:
                p = b[7]*1.0/b[6]
                binom = math.pow(p,b[4])
                num_DNA_control = b[6]
                #size_ref = b[6]
                #postPV = postP(b[6],b[3],b[7])
            elif l2==0:
                p = b[4]*1.0/b[3]
                binom = math.pow(p,b[7])
                num_DNA_control = b[3]
                #size_ref = b[3]
                #postPV = postP(b[3],b[6],b[4])
        
            pvalue = binom
            #if binom > 1e-2:# or size_ref <= 20 : continue
            if pvalue > 1e-2:# or size_ref <= 20 : continue
                continue

            d, pv = stats.ks_2samp([0]*b[4]+a.iloc[i][5], [0]*b[6]+a.iloc[i][8])
            #print b[0], b[1], b[2]
            #print b[3], b[4], b[5]
            #print b[6], b[7], b[8]
            
            # write out the output file with varous threshold

            het_del.write("%s\t%d\t%d\t%.3e\t%.3e\t%d\n" % (b[9],b[10],b[11], pv, pvalue, num_DNA_control))
            het_del_data.write("%s\t%d\t%d\t%.3e\t%.3e\t" % (b[9],b[10],b[11], pv, pvalue))
            het_del_data.write("%s\t%d\t%d\t" % (b[9],b[10],b[11]))
            het_del_data.write("%d\t" % (b[0]))
            het_del_data.write("%d\t" % (b[1]))
            for ii in b[2]:
                het_del_data.write("%d," % (ii))
            het_del_data.write("\t")
            het_del_data.write("%d\t" % (b[3]))
            het_del_data.write("%d\t" % (b[4]))
            for ii in b[5]:
                het_del_data.write("%d," % (ii))
            het_del_data.write("\t")
            het_del_data.write("%d\t" % (b[6]))
            het_del_data.write("%d\t" % (b[7]))
            for ii in b[8]:
                het_del_data.write("%d," % (ii))
            het_del_data.write("\n")
            #het_del.write(b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],"\n")
            #print i, d, pv, p, binom
            #print
        else:
            # test dup

            rPM1 = [0]*b[4]+b[5]
            rPM2 = [0]*b[7]+b[8]
            #mu1, mu2 = np.mean(rPM1), np.mean(rPM2)
            #z0, z1, z2 = b[1]*1.0/b[0], b[4]*1.0/b[3], b[7]*1.0/b[6]
            z1, z2 = b[4]*1.0/b[3], b[7]*1.0/b[6]

            d, pv = stats.ks_2samp(rPM1, rPM2)
            if pv < 0.20:
                #print mu1, mu2
                #print "%.2g" % z0,
                #print b[0], b[1], b[2]
                #print "%.2g" % z1,
                #print b[3], b[4], b[5]
                #print "%.2g" % z2,
                #print b[6], b[7], b[8]
                #print i, d, pv
                #nz1, nz2 = b[4], b[7]
                pt, p1 = getFreq(b[4],b[5]), getFreq(b[7],b[8])
                nTest = b[3]
                if z2 < z1:  # the one with lower zero freq is assumed to be the test freq
                    pt, p1 = p1, pt
                    nTest = b[6]
                cp, ll = llhTest(pt, testProb(getGoodFreq(p1,len(pt), 0.1)), nTest)
                #print cp, ll
                if ll > 6:
                    numCNV+=1
                    tandem_dup.write("%s\t%d\t%d\t%.2g\t%d\n" % (b[9], b[10], b[11], ll, cp))
                    #print "one more tandem dup", numCNV, (i+1)
                #print

    het_del.close()
    het_del_data.close()
    tandem_dup.close()
    reg_dup.close()




