#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import subprocess
import os
import os.path
import tenkit.test as tk_test
import tenkit.bio_io as tk_io
from tenkit.constants import TEST_FILE_IN_DIR, TEST_FILE_OUT_DIR
from .. import *

import martian

def infile(path):
    return os.path.join(TEST_FILE_IN_DIR, "populate_info_fields", path)

def outfile(path):
    return os.path.join(TEST_FILE_OUT_DIR, path)

mrs_job_name = "populate_info_fields"
base_out_dir = outfile("")
job_dir = os.path.join(base_out_dir, mrs_job_name,  "POPULATE_INFO_FIELDS", "fork0", "chnk0", "files")

class TestFunctions(tk_test.UnitTestBase):

    def setUp(self):
        self.clear_directory()
        martian.test_initialize(outfile(""))

        self.args = {
                'bam': infile("chr10.bam"),
                'variants': "null",
                'input': infile("chr10.vcf")
                }
    def write_mro(self, args):
        tpl = """
            @include "_snpindel_caller_stages.mro"
            call POPULATE_INFO_FIELDS(
                bam = "%(bam)s",
                input = "%(input)s",
                variant_mode = "freebayes",
                vc_precalled = null,
                reference_path = "/mnt/opt/refdata_new/hg19",
                min_mapq_attach_bc = 60,
                default_indel_qual = 43,
            )
        """
        fn = os.path.join(base_out_dir, "test.mro")
        with open(fn, 'w') as f:
            f.write(tpl % args)

        return fn


    def run_stage(self, args):
        mro_file = self.write_mro(args)
        proc = subprocess.Popen(['mrs', mro_file, mrs_job_name], cwd=base_out_dir, stdout=subprocess.PIPE)
        out, err = proc.communicate()

        if proc.returncode != 0:
            print out
            raise Exception("mrs failed on during test")


    def test_basic(self):
        self.run_stage(self.args)

        # Load the output file
        vfr = tk_io.VariantFileReader(os.path.join(job_dir,"default.vcf"))
        for r in vfr.record_getter():
            pos = tk_io.get_record_pos(r)
            barcodes = tk_io.get_record_barcodes(r)
            if pos == 26357747 or pos == 26357748:

                print barcodes
                assert(barcodes[1][0] =='1-ATAGGAGTTCAGGG_63')
                print tk_io.get_record_alt_allele_counts(r)
                assert(int(tk_io.get_record_alt_allele_counts(r)[0]) in [31,32])
                assert(int(tk_io.get_record_ref_allele_count(r)) == 0)
            if pos == 26501280:
                print barcodes
                assert(barcodes[0][0] == '1-TGAAGACATAACCC_61_61')
                assert(int(r.INFO['POSTHPC'][0]) == 10)


