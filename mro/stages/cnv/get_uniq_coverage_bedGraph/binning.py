# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import math


def Harmonic(x):
    return 2.0 / 3.0 * math.pow(x, 1.5) + math.sqrt(x) / 2.0 - 0.2


def InvHarmonic(y):
    return math.pow(1.5 * (y + 0.2), 2.0 / 3.0)


def getBins(n, amp):
    numBins = int(round(InvHarmonic(n * 1.0 / amp))) + 1
    bins = [int(round(Harmonic(x) * amp)) for x in range(numBins)]
    return bins


def discretize(num, bins, amp):
    idx = int(round(InvHarmonic(num * 1.0 / amp)))
    if num > bins[idx]:
        print "Wrong! idx in discreteize function is too low", num, bins[idx]
        return -1
    while idx > 0 and num <= bins[idx - 1]:
        idx = idx - 1
    return idx
