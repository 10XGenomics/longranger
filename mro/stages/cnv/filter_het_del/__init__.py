# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.

import longranger.cnv.event_bed as cnv_event
#import longranger.cnv.para_selection as para_selection

__MRO__ = """
stage FILTER_HET_DEL(
    in  bool        wgsmode,
    in  bed         het_del,
    in  float       pvalFlt,
    in  float       WildCovFlt,
    in  float       BGImparityPvalFlt,
    in  float       BGBaseCovFlt,
    in  float       BGPhaseFracFlt,
    in  map         blacklist_map,
    out bed         het_del_filtered,
    out bedpe       het_del_filtered_bedpe,
    src py          "stages/cnv/filter_het_del",
)
"""

def main(args, outs):
    if not args.het_del:
        outs.het_del_filtered = None
        outs.het_del_filtered_bedpe = None
        return

    #if not args.pvalFlt or not args.WildCovFlt:
    #    pval_flt, wild_cov_flt = para_selection.cal_het_del_parameters()
    #else:
    pval_flt, wild_cov_flt = args.pvalFlt, args.WildCovFlt

    with open(args.het_del) as f:
        events = [cnv_event.Event.from_native(l) for l in f]

    for (event_num, e) in enumerate(events):
        e.name = "del%d" % event_num
        e.set_filters(pval_flt, wild_cov_flt, args.BGImparityPvalFlt, args.BGBaseCovFlt, args.BGPhaseFracFlt)

    # currently no blacklist filtering for wgs het deletions
    if not args.wgsmode and args.blacklist_map:
        extra_blacklist = {}
        for k, v in args.blacklist_map.iteritems():
            if k not in ["SEGDUP","LowCov", "whitelist"]:
                extra_blacklist[k] = v
        if extra_blacklist: cnv_event.do_blacklist_filtering(events, extra_blacklist)
        #if "whitelist" in args.blacklist_map and "het" in args.blacklist_map["whitelist"]:
        #    cnv_event.do_whiltelist_filtering(events, args.blacklist_map["whitelist"]["het"], "1000G_ACCS")

    with open(outs.het_del_filtered, "w") as fout:
        for e in events:
            if e.filter_rst[0]:
                fout.write(e.get_bed_output()+"\n")

    with open(outs.het_del_filtered_bedpe, "w") as fout:
        for e in events:
            fout.write(e.get_bedpe_output()+"\n")
