# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import subprocess

__MRO__='''
stage COPY_FRAGMENT_INFO(
    in  h5   fragments,
    in  h5   barcodes,
    in  json fragment_histogram,
    in  json single_partition_results,
    in  json barcode_histogram,
    out h5   fragments,
    out h5   barcodes,
    out json fragment_histogram,
    out json single_partition_results,
    out json barcode_histogram,
    src py "stages/reporter/copy_fragment_info",
)
'''

def main(args, outs):

    if args.barcodes is not None and args.fragments is not None:
        subprocess.check_call(['cp', args.barcodes, outs.barcodes])
        subprocess.check_call(['cp', args.fragments, outs.fragments])
        subprocess.check_call(['cp', args.fragment_histogram, outs.fragment_histogram])
        subprocess.check_call(['cp', args.single_partition_results, outs.single_partition_results])
        subprocess.check_call(['cp', args.barcode_histogram, outs.barcode_histogram])
    else:
        outs.barcodes = None
        outs.fragments = None
        outs.fragment_histogram = None
        outs.single_partition_results = None
        outs.barcode_histogram = None
