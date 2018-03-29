# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.bed_utils as bedtools
import tenkit.bio_io as tk_io
#import tenkit.regions as tk_regions
import shutil
import tenkit.constants
import tenkit.reference as tk_reference

__MRO__ = """
stage GENERATE_EXOME_GROUND_TRUTH (
    in  path    reference_path,
    in  bedpe   wgs_deletions_gt,
    in  bool    wgsmode,
    in  bed     target_regions,
    in  int     min_overlap,
    in  bed     low_coverage_blacklist,
    in  map     blacklist_map,
    out bed     het_gt_sen,
    out bed     het_gt_ppv,
    out bed     hom_gt_sen,
    out bed     hom_gt_ppv,
    out bed     het_del_query_region,
    out bed     hom_del_query_region,
    out map     blacklist_map,
    src py      "stages/cnv/generate_exome_ground_truth",
)
"""

def main(args, outs):

    if not args.wgsmode:
        # blacklist curation
        if args.blacklist_map:
            # If we got an explicit blacklist, use it
            blacklist_map = args.blacklist_map
        else:
            blacklist_map = {}
        
            
        # We did not get an explicit blacklist -- in this case, combine the built-in segdup file and the homo_del_blacklist
        # Combine the homo_del_blacklist with the internal segdup file to get the filter set
        if not "LowCov" in blacklist_map:
            blacklist_map["LowCov"] = args.low_coverage_blacklist

        genome = tk_reference.get_genome(args.reference_path)
        if not "SEGDUP" in blacklist_map:
            cnv_segdups = tenkit.constants.find_sv_file(genome, "cnv_segdup_filter.bed")
            blacklist_map["SEGDUP"] = cnv_segdups
            
        if not "whitelist" in blacklist_map:  blacklist_map["whitelist"] = {}
        if not "homo" in blacklist_map["whitelist"]:
            pilot_accs = tenkit.constants.find_sv_file(genome, "20141020.pilot_mask.whole_genome.bed")
            blacklist_map["whitelist"]["homo"] = pilot_accs
        if not "het" in blacklist_map["whitelist"]:
            strict_accs = tenkit.constants.find_sv_file(genome, "20141020.strict_mask.whole_genome.bed")
            blacklist_map["whitelist"]["het"] = strict_accs
        
        outs.blacklist_map = blacklist_map
            

        # generate filtered target regions for het del and homo del callers
        ## homo
        if "LowCov" in blacklist_map:
            bedtools.no_overlap(args.target_regions, blacklist_map["LowCov"], outs.hom_del_query_region+"_tmp.bed")
        else:
            shutil.copyfile(args.target_regions, outs.hom_del_query_region+"_tmp.bed")
        bedtools.overlap(outs.hom_del_query_region+"_tmp.bed", blacklist_map["whitelist"]["homo"], outs.hom_del_query_region)
        #shutil.copyfile(outs.hom_del_query_region+"_tmp.bed", outs.hom_del_query_region)
            
        ## het
        bedtools.no_overlap(outs.hom_del_query_region+"_tmp.bed", blacklist_map["SEGDUP"], outs.het_del_query_region+"_tmp.bed")
        bedtools.overlap(outs.het_del_query_region+"_tmp.bed", blacklist_map["whitelist"]["het"], outs.het_del_query_region)
        
    
#         ## 
#         with open(blacklist_map["whitelist"]["homo"]) as f:
#             accs_1000g_pilot = tk_io.get_target_regions(f)
#         with open(blacklist_map["whitelist"]["het"]) as f:
#             accs_1000g_strict = tk_io.get_target_regions(f)
    
    
        if args.wgs_deletions_gt:
            ## het evetns
            with open(outs.het_del_query_region) as f:
                bed_target = tk_io.get_target_regions(f)

            fhet_sen = open(outs.het_gt_sen, "w")
            fhet_ppv = open(outs.het_gt_ppv, "w")

            with open(args.wgs_deletions_gt) as f:
                for line in f:
                    if line[0]=="#": continue
                    infos = line.strip().split()
                    chrom, start, end = infos[0], int(infos[1]), int(infos[5])
                    if chrom in bed_target:
                        overlappings = bed_target[chrom].overlapping_regions(start,end)
                        overlap_size = 0
                        for s, e in overlappings:
                            overlap_size += (min(e, end)-max(s,start))
                        if overlap_size >= args.min_overlap: #and \
                            #chrom in accs_1000g_strict and\
                            #accs_1000g_strict[chrom].overlaps_region(start, end):
                            record = "\t".join((infos[0], infos[1], infos[5]))+"\n"
                            if ("HET" in line) and ("TIER=1" in line):
                                fhet_sen.write(record)
                            fhet_ppv.write(record)

            
            ## homo events
            with open(outs.hom_del_query_region) as f:
                bed_target = tk_io.get_target_regions(f)

            fhom_sen = open(outs.hom_gt_sen, "w")
            fhom_ppv = open(outs.hom_gt_ppv, "w")

            with open(args.wgs_deletions_gt) as f:
                for line in f:
                    if line[0]=="#": continue
                    infos = line.strip().split()
                    chrom, start, end = infos[0], int(infos[1]), int(infos[5])
                    if chrom in bed_target:
                        overlappings = bed_target[chrom].overlapping_regions(start,end)
                        has_full_exon = False
                        for s, e in overlappings:
                            if start <= s+1 and end >= e-1:
                                print start, end, s, e
                                has_full_exon = True
                                break

                        if has_full_exon: #and \
                            #chrom in accs_1000g_pilot and\
                            #accs_1000g_pilot[chrom].overlaps_region(start, end):
                            record = "\t".join((infos[0], infos[1], infos[5]))+"\n"
                            if ("HOM" in line) and ("TIER=1" in line):
                                fhom_sen.write(record)
                            fhom_ppv.write(record)

            fhet_sen.flush(); fhet_sen.close()
            fhet_ppv.flush(); fhet_ppv.close()
            fhom_sen.flush(); fhom_sen.close()
            fhom_ppv.flush(); fhom_ppv.close()

        else:
            outs.het_gt_sen = None
            outs.het_gt_ppv = None
            outs.hom_gt_sen = None
            outs.hom_gt_ppv = None

    else:
        outs.hom_gt_sen = None
        outs.hom_gt_ppv = None
        outs.het_del_query_region = None
        outs.hom_del_query_region = None
        outs.blacklist_map = None
        fhet_sen = open(outs.het_gt_sen, "w")
        fhet_ppv = open(outs.het_gt_ppv, "w")

        if args.wgs_deletions_gt:
            with open(args.wgs_deletions_gt) as f:
                for line in f:
                    if line[0]=="#": continue
                    infos = line.strip().split()
                    record = "\t".join((infos[0], infos[1], infos[5]))+"\n"
                    if "TIER=1" in line:
                        fhet_sen.write(record)
                    fhet_ppv.write(record)

            fhet_sen.flush(); fhet_sen.close()
            fhet_ppv.flush(); fhet_ppv.close()
        else:
            outs.het_gt_sen = None
            outs.het_gt_ppv = None
