# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#import json
#import math

___MRO___ = """
stage COMPUTE_CUTOFFS(
    in  map   length_mass,
    in  float inferred_total_mass_ng,
    out float het_del_pval_flt,
    out float het_del_wildcov_flt,
    out float homo_del_pval_flt,
    out float homo_del_contamination_flt,
    src py "stages/cnv/compute_cutoffs",
)
"""

def main(args, outs):
    #inferred_total_mass_ng = None
    #if args.inferred_total_mass_ng: inferred_total_mass_ng = args.inferred_total_mass_ng
    #else:
    #    if args.length_mass:
    #        with open(args.length_mass) as data_file:
    #            length_mass_data = json.loads(data_file.read())
    #            if "inferred_total_mass_ng" in length_mass_data:
    #                inferred_total_mass_ng = length_mass_data["inferred_total_mass_ng"]

    outs.homo_del_pval_flt = 1.0e-9
    outs.homo_del_contamination_flt = 3.0
    #if inferred_total_mass_ng:
    #    ### the linear formula are estimated using weight liner regressions
    #    ### on the cutoffs values passing the market specs (0.5 sen and 0.8 ppv)
    #    outs.het_del_wildcov_flt = 30.922 - inferred_total_mass_ng*6.883
    #    outs.het_del_pval_flt = math.pow(10, -9.2108-0.1244*inferred_total_mass_ng)
    #else:
    outs.het_del_wildcov_flt = 18.0
    outs.het_del_pval_flt =1e-8
