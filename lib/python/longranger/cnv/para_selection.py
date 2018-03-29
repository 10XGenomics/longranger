# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import math

def cal_homo_del_parameters():
    return (1.0e-10, 4.0)

def cal_het_del_parameters(REPORT_LENGTH_MASS):
    if REPORT_LENGTH_MASS:
        wildcov = REPORT_LENGTH_MASS*2.073 + 14.43
        pval = math(10.0, -2.03*REPORT_LENGTH_MASS-5.48)
        return (pval, wildcov)
    else:
        return (1.0e-08, 15.0)
