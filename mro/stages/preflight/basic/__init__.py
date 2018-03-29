#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import re
import socket
import martian
import tenkit.preflight as tk_preflight

__MRO__ = """
stage BASIC_PREFLIGHT(
    in string sample_id,
    in string read_group,
    in string output_format,
    in map[]  sample_def,
    in bool   check_executables,
    src py     "stages/preflight/basic",
)
"""

def main(args, outs):
    hostname = socket.gethostname()

    if args.output_format == 'bam' and args.read_group is None:
        martian.exit("Please specify a read_group to populate the @RG field of the BAM file")

    if args.sample_id is not None:
        if not re.match("^[\w-]+$", args.sample_id):
            martian.exit("Sample name may only contain letters, numbers, underscores, and dashes: " + args.sample_id)

    for sample_def in args.sample_def:
        read_path = sample_def["read_path"]
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, longranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        library_id = sample_def.get("library_id")
        if library_id is not None:
            if not re.match("^[\w-]+$", library_id):
                martian.exit("Library name may only contain letters, numbers, underscores, and dashes: " + library_id)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not tk_preflight.is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

        ok, msg = tk_preflight.check_sample_indices(sample_def)
        if not ok:
            martian.exit(msg)

    # Check open file handles limit
    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

    martian.log_info(tk_preflight.record_package_versions())


