# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import longranger.cnv.event_bed as cnv_event
import longranger.sv.io as sv_io
import shutil
import os

__MRO__ = """
stage FILTER_HOMO_DEL(
    in  bed        homo_del,
    in  float      homo_pval_flt,
    in  float      homo_ctm_flt, # number of contamination reads
    in  map        blacklist_map,
    in  bedpe      het_del_filtered_bedpe,
    in  string     reference_path,
    in  string     sample_name,
    out bed        homo_del_filtered,
    out bedpe      homo_del_filtered_bedpe,
    out bedpe      deletions,
    out vcf.gz     deletions_vcf,
    out vcf.gz.tbi deletions_vcf_index,
    src py         "stages/filter_homo_del",
) split using ()
"""

def main(args, outs):
    if not args.homo_del:
        outs.homo_del_filtered = None
        outs.homo_del_filtered_bedpe = None
        if args.het_del_filtered_bedpe:
            shutil.copyfile(args.het_del_filtered_bedpe, outs.deletions)
        else:
            outs.deletions = None
        return

    #if not args.homo_pval_flt or not args.homo_ctm_flt:
    #    homo_pval_flt, homo_ctm_flt = para_selection.cal_homo_del_parameters()
    #else:
    homo_pval_flt, homo_ctm_flt = args.homo_pval_flt, args.homo_ctm_flt

    with open(args.homo_del) as f:
        events = [cnv_event.HomoEvent.from_native(line) for line in f]

    for e in events:
        print e.ttlMolMs


    #if "whitelist" in args.blacklist_map and "homo" in args.blacklist_map["whitelist"]:
    #    cnv_event.do_homo_whiltelist_filtering(events, args.blacklist_map["whitelist"]["homo"], "1000G_ACCS")




    with open(outs.deletions, "w") as fout:
        event_num = 0
        with open(args.het_del_filtered_bedpe) as fhet:
            for line in fhet:
                fout.write(line)
                event_num += 1

        for e in events:
            e.name = "del%d" % event_num
            event_num += 1

            #e.check_blacklist(blacklist, "LowProbeCoverage")
            if e.pval > homo_pval_flt:
                if e.filters == "PASS":
                    e.filters = "Large_Pval"
                else:
                    e.filters += ";Large_Pval"
            if e.Nread > homo_ctm_flt:
                if e.filters == "PASS":
                    e.filters = "HighContamination"
                else:
                    e.filters += ";HighContamination"

            fout.write(e.get_bedpe_output()+"\n")

    with open(outs.homo_del_filtered_bedpe, "w") as fout:
        for e in events:
            fout.write(e.get_bedpe_output()+"\n")

    with open(outs.homo_del_filtered, "w") as fout:
        for e in events:
            fout.write(e.get_bed_output()+"\n")

    sv_io.bedpe_to_vcf(outs.deletions, outs.deletions_vcf.strip('.gz'),
                       args.sample_name, "10X_HAP_CNV", args.reference_path)
    sv_io.index_sv_vcf(outs.deletions_vcf.strip(".gz"))
    outs.deletions_vcf_index = outs.deletions_vcf + '.tbi'

    # delete the non-gzipped file
    os.remove(outs.deletions_vcf.strip('.gz'))

