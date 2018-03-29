#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# This is a shim to call the loupe preprocessor form python

import json
import martian
import tenkit.reference
import subprocess
import os
import shutil

__MRO__ = """
stage LOUPE_PREPROCESS(
    in  string reference_path,
    in  vcf.gz input_vcf,
    in  bam    sorted_deduplicated_bam,
    in  bedpe  structvar_data,
    in  bedpe  shortstructvar_data
    in  tsv    bkpt_details,
    in  json   fragment_histogram,
    in  json   summarize_output,
    in  json   single_partition_data,
    in  json   coverage_data,
    in  json   phasing_quality_data,
    in  json   alarms,
    in  bed    targets,
    in  tsv.gz fragments,
    in  string sample_desc,
    in  bool   noloupe,
    out loupe  output_for_loupe,
    src py     "stages/loupe",
 ) split using
 )
"""

# Use split to request more memory
def split(args):
    mem_gb = 16
    chunks = [{'__mem_gb': mem_gb}]
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    if chunk_outs[0].output_for_loupe is None:
        # Set output to null if noloupe is set.
        outs.output_for_loupe = None
    else:
        # Copy loupe file to output directory.
        shutil.copy(chunk_outs[0].output_for_loupe, outs.output_for_loupe)

def main(args, outs):
    # Skip stage if noloupe is set.
    if args.noloupe:
        outs.output_for_loupe = None
        return

    # Path to store temporary goloupe data. (We'll delete this at the end of
    # this function.)
    temp_path = martian.make_path("temporary")
    os.mkdir(temp_path);

    uncompressed_vcf_path = martian.make_path("temporary/phased_snps.vcf");
    uncompressed_vcf_file = open(uncompressed_vcf_path, "w");

    # Uncompress VCF file for goloupe.
    ret = subprocess.call(["gunzip", "-c", args.input_vcf], stdout=uncompressed_vcf_file)
    uncompressed_vcf_file.close()

    if ret != 0:
        raise Exception("Cannot uncompress VCF file for goloupe")

    pipeline_params_path = martian.make_path("temporary/pipeline_params.json")
    pipeline_params_file = open(pipeline_params_path, "w")

    pipeline_params = martian.get_invocation_args()
    pipeline_params['pipelines_version'] = martian.get_pipelines_version()

    pipeline_params_file.write(json.dumps(pipeline_params))

    pipeline_params_file.close()


    # These need to be our files for hg19
    genes_path = tenkit.reference.get_loupe_genes(args.reference_path)
    contigs_path = tenkit.reference.get_fasta_contig_index_path(args.reference_path)

    summary_files = [pipeline_params_path + ":pipeline_params"]
    if args.summarize_output is not None:
        summary_files.append(args.summarize_output + ":summary_data")
    if args.single_partition_data is not None:
        summary_files.append(args.single_partition_data + ":partition_data")
    if args.coverage_data is not None:
        summary_files.append(args.coverage_data + ":coverage_data")
    if args.phasing_quality_data is not None:
        summary_files.append(args.phasing_quality_data + ":phasing_data")
    if args.fragment_histogram is not None:
        summary_files.append(args.fragment_histogram + ":fragment_histogram")
    if args.length_mass_data is not None:
        summary_files.append(args.length_mass_data + ":length_mass_data")
    if args.alarms is not None:
        summary_files.append(args.alarms + ":alarms")

    call = ["goloupe",
        "-vcf=" + uncompressed_vcf_path,
        "-output=" + outs.output_for_loupe,
        "-bam=" + args.sorted_deduplicated_bam,
        "-tempdir=" + temp_path,
        "-contigs=" + contigs_path,
        "-summary=" + str.join(",", summary_files)]

    if os.path.isfile(genes_path):
        call.append("-bed=" + genes_path)

    if args.structvar_data is not None:
        call.append("-sv=" + args.structvar_data)

    if args.shortstructvar_data is not None:
        call.append("-ssv=" + args.shortstructvar_data)

    if args.bkpt_details is not None:
        call.append("-bkpt=" + args.bkpt_details)

    if args.targets is not None:
        call.append("-targets=" + args.targets)

    # Call goloupe which will produce a .loupe file from a bunch of inputs
    ret = subprocess.call(call);

    # If goloupe fails, crash only if this is an internal pipeline. 
    if ret != 0:
        try:
            import kitten
            kitten = kitten
            raise Exception("goloupe failed. Sorry")
        except ImportError:
            outs.output_for_loupe = None

    # Delete all of our temporary stuff
    shutil.rmtree(temp_path);
