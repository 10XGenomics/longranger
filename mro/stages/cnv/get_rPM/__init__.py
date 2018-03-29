#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# convert fragments.h5 file into csv with selected columns
#

import os
import subprocess

import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.chunk_utils as tk_chunks
from tenkit.constants import PARALLEL_LOCUS_SIZE

__MRO__ = """
stage GET_RPM(
    in  csv      fragments,
    in  tsv      fragment_phasing,
    in  bam      phased_possorted_bam,
    in  bed      queryRegions,
    in  bool     wgsmode,
    in  int      mapq,
    in  int      overlap,
    out bed      rPM,
    out string[] rPMFiles,
    out txt      stat,
    out txt      covStat,
    src py       "stages/get_rPM",
)
"""

def findChroms(bedfile):
    chroms = {}
    with open(bedfile) as f:
        for l in f:
            items = l.split()
            ch = items[0]
            if not ch in chroms:
                chroms[ch]=1
    return chroms



def split(args):
    input_bam = tk_bam.create_bam_infile(args.phased_possorted_bam)
    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    loci = tk_chunks.chunk_by_locus(chroms, chrom_lengths, PARALLEL_LOCUS_SIZE, extra_args = {'__mem_gb': 12, '__threads':4})
    return {'chunks': loci}


def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    if not args.fragments:
        outs.rPM = None
        outs.rPMFiles = None
        outs.stat = None
        outs.fracPhased = None
        outs.covStat = None
        return


    tool = "molecular_count"
    cmd0 = tool
    if args.wgsmode:
        cmd0 += " --wgs"

    tmp_bed_file = outs.stat+"_tmp.bed"
    (chrom, start, end) = tk_io.get_locus_info(args.locus)

    nrows = 0
    with open(tmp_bed_file, "w") as fout:
        with open(args.queryRegions) as f:
            for l in f:
                items=l.split()
                reg_chrom = items[0]
                reg_start = int(items[1])
                reg_end = int(items[2])

                if chrom == reg_chrom and reg_start >= start and reg_end < end:
                    nrows += 1
                    fout.write(l)

    if nrows == 0:
        outs.rPM = None
        return

    # specify the target file
    cmd0 += " -b "+tmp_bed_file

    cmd = cmd0 + \
          " -f " + args.fragments + \
          " -p " + args.fragment_phasing + \
          " -m " + args.phased_possorted_bam + \
          " --rPM " + outs.rPM + " --stat " + outs.stat + \
          " --fracPhased " + outs.fracPhased + " --mapq "+str(args.mapq) + \
          " --overlap " + str(args.overlap) + " --covstat " + str(outs.covStat)

    print "Running cmd:"
    print cmd

    subprocess.check_call(cmd, shell=True)

    #try:
    #    outtext=subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    #    print outtext
    #except subprocess.CalledProcessError as err:
    #    print err.output
    #    raise Exception ("error in running molecular_count")


def join(args, outs, chunk_defs, chunk_outs):

    if not args.fragments:
        outs.rPM = None
        outs.rPMFiles = None
        outs.stat = None
        outs.fracPhased = None
        outs.covStat = None
        return

    rPMFiles = [out.rPM for out in chunk_outs if out.rPM and os.path.exists(out.rPM) ]
    outs.rPMFiles = rPMFiles

    with open(outs.rPM, "w") as fout:
        fout.write("chrom\tstart\tend\tnumBC\tbc_hap1\tbc_hap2\n")
        for fl in rPMFiles:
            with open(fl) as fin:
                #fin.readline()
                for l in fin:
                    infos = l.split(',')
                    out_infos = [infos[0], infos[1], infos[2], infos[5], infos[7], infos[9]]
                    fout.write("\t".join(out_infos)+"\n")

    frac_files = [out.fracPhased for out in chunk_outs if out.fracPhased and os.path.exists(out.fracPhased)]
    with open(outs.fracPhased, "w") as fout:
        fout.write("exonName,fracPhased,phasedRead,totalRead\n")
        for fl in frac_files:
            with open(fl) as fin:
                fin.readline()
                for l in fin:
                    fout.write(l)

    stat_files = [out.stat for out in chunk_outs if out.stat and os.path.exists(out.stat)]
    with open(outs.stat, "w") as fout:
        for fl in stat_files:
            with open(fl) as fin:
                for l in fin:
                    fout.write(l)

    ttlTargetSz = 0.0
    ttlMol = 0.0
    ttlMolSeen = 0.0
    for f in  [out.covStat for out in chunk_outs if out.covStat and os.path.exists(out.covStat) ]:
        with open(f) as fin:
            fin.readline()
            line = fin.readline()
            infos = line.strip().split()
            molSeen, mol, sz = float(infos[0]), int(infos[1]), int(infos[2])
            ttlTargetSz += sz
            ttlMol += mol
            ttlMolSeen += molSeen
    MolSeenPerSZ = ttlMolSeen/ttlTargetSz

    with open(outs.covStat, "w") as fout:
        fout.write("OverallMolecularSeen\tOverallMolecular\tOverallTargetSize\tMolSeenPerSZ\n")
        fout.write("%.1f\t%d\t%d\t%.6f\n" % (ttlMolSeen, ttlMol, ttlTargetSz, MolSeenPerSZ))
