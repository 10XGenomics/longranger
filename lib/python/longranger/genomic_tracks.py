# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import os
import tenkit.reference as tk_reference

GENOMIC_TRACKS_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "genomic_tracks"))
print "genomic_tracks_dir", GENOMIC_TRACKS_DIR

## sv_blacklist terminal_cnv_blacklist etc.


# find the right blacklist according to the user input information, mode, reference_path, default_file_name
# user input information can either be specific file path or none
# mode takes effect only when user_input is noen. It can specify the exact modes
# genome_version spficifies the specific refernece genomic version
# default_file_name is the file name for the default genomic track
def get_genomic_track(user_input, mode, reference_path, default_file_name):
    ## return a file path pointing to a valid blacklist file
    if user_input is not None:
        return user_input
    
    # use one of the shipped blacklist files
    genome_version = tk_reference.get_genome(reference_path)
    if genome_version is None:
        return None
    print GENOMIC_TRACKS_DIR, genome_version, mode, default_file_name
    file_path = os.path.join(GENOMIC_TRACKS_DIR, genome_version, mode, default_file_name)
    return file_path if os.path.exists(file_path) else None
