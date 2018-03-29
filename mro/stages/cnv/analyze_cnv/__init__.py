# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import longranger.cnv.analysis as cnv_tool
from tenkit import safe_json

__MRO__ = """
stage CNV_ANALYSIS_WGS(
    in  bool        wgsmode,
    in  bed         het_del,
    in  bed         het_del_truth_sensitivity,
    in  bed         het_del_truth_ppv,
    in  map         het_del_stratification,
    in  map         blacklist_map,
    in  bed         homo_del,
    in  bed         homo_del_truth_sensitivity,
    in  bed         homo_del_truth_ppv,
    out json        summary,
    out bed         het_del_falsePositives,
    out bed         het_del_falseNegatives,
    out bed         homo_del_falsePositives,
    out bed         homo_del_falseNegatives,
)
"""


def main(args, outs):

    if not args.het_del:
        outs.summary = None
        return


    output_dict = {}

    # No ground truth -- just report number of calls
    if (not args.het_del_truth_sensitivity)  or (not args.het_del_truth_ppv):
        myHets = cnv_tool.Calls(args.het_del, None)
        output_dict["het_del_numPosCalls"] = myHets.TotalEvent

        myHoms = cnv_tool.Calls(args.homo_del, None)
        output_dict["homo_del_numPosCalls"] = myHoms.TotalEvent

        output_dict["total_del_numPosCalls"] = myHets.TotalEvent + myHoms.TotalEvent

        fsummary = open(outs.summary, "w")
        fsummary.write(safe_json.safe_jsonify(output_dict))
        fsummary.close()

        return

   ############################################################
   ## het del                                                 #
   ############################################################

   # set up the true data variable truthData
    truthNames = ["TDsens", "TDPPV"]
    truths = [args.het_del_truth_sensitivity, args.het_del_truth_ppv]
    if args.het_del_stratification:
        truthData = [cnv_tool.TruthSet(t,n, regionsAvoided=args.het_del_stratification) for t,n in zip(truths, truthNames)]
    else:
        truthData = [cnv_tool.TruthSet(t,n) for t,n in zip(truths, truthNames)]


    # set up stratification tracks
    tracksAvoiding = {}
    if args.het_del_stratification:
        for k, v in args.het_del_stratification.iteritems():
            tracksAvoiding[k]=v


    mycalls = cnv_tool.Calls(args.het_del, truthData)

    sen, ppv, numPos = mycalls.checkCallPerformance(SENIDX=0, PPVIDX=1)

    output_dict["het_del_sen"]=sen
    output_dict["het_del_ppv"]=ppv
    output_dict["het_del_numPosCalls"]=numPos

    mycalls.getFalseNegative(outfile=outs.het_del_falseNegatives, whichTD=0)
    mycalls.getFalsePos(outfile=outs.het_del_falsePositives, whichTD=1)

   ############################################################
   ## homo del                                                 #
   ############################################################

    if not args.wgsmode:
        truthNames = ["TDsens", "TDPPV"]
        truths = [args.homo_del_truth_sensitivity, args.homo_del_truth_ppv]
        truthData = [cnv_tool.TruthSet(t,n) for t,n in zip(truths, truthNames)]

        mycalls = cnv_tool.Calls(args.homo_del, truthData)

        sen, ppv, numPos = mycalls.checkCallPerformance(SENIDX=0, PPVIDX=1)

        output_dict["homo_del_sen"]=sen
        output_dict["homo_del_ppv"]=ppv
        output_dict["homo_del_numPosCalls"]=numPos

        mycalls.getFalseNegative(outfile=outs.homo_del_falseNegatives, whichTD=0)
        mycalls.getFalsePos(outfile=outs.homo_del_falsePositives, whichTD=1)


    output_dict["total_del_numPosCalls"] = output_dict.get("homo_del_numPosCalls", 0) + output_dict.get("het_del_numPosCalls", 0)
    fsummary = open(outs.summary, "w")
    fsummary.write(safe_json.safe_jsonify(output_dict))
    fsummary.close()
