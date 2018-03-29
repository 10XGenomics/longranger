// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * Provide a simple "main" wrapper to implement the actual executable
 */

package main

import (
	"flag"
	"io/ioutil"
	"log"
	"loupe/formats"
	"loupe/preprocessor"
	"os"
	"runtime/pprof"
)

/* Command line arguments */
var input_vcf_path = flag.String("vcf", "", "input VCF file [required]")
var input_bed_path = flag.String("bed", "", "input BED file [required]")
var output_loupe_path = flag.String("output", "", "output LOUPE file [required]")
var profile_path = flag.String("profile", "", "Write gprof profile data to file")
var extra_summary_paths = flag.String("summary", "", "optional list of comma separated summary files")
var input_refseq_path = flag.String("refseq", "", "input refseq file")
var input_bam_path = flag.String("bam", "", "sorted, deduplicated bam file")
var temp_dir_path = flag.String("tempdir", "", "Temporary directory to do work in")
var mem_profile_path = flag.String("mprofile", "", "Write memprof data to this file")
var input_structvar_path = flag.String("sv", "", "Structural variant data")
var input_shortstructvar_path = flag.String("ssv", "", "Short structural variant data")
var input_breakpoints_path = flag.String("bkpt", "", "Structural variant breakpoint details")
var input_targets_path = flag.String("targets", "", "Targets File")
var input_fragments_path = flag.String("fragments", "", "Fragments File")
var input_contig_path = flag.String("contigs", "", "Contigs File");

func main() {
	flag.Parse()

	/* Check to make sure everything is filled out */
	if *output_loupe_path == "" {
		log.Printf("Bad arguments")
		os.Exit(1)
	}

	if *profile_path != "" {
		profile_file, err := os.Create(*profile_path)
		if err != nil {
			log.Printf("Error: %v", err)
			os.Exit(1)
		}
		pprof.StartCPUProfile(profile_file)
	}

	var extra_summary_files map[string]string

	if extra_summary_paths != nil && *extra_summary_paths != "" {
		extra_summary_files = formats.ParseSummarySpecification(*extra_summary_paths)
	} else {
		extra_summary_files = map[string]string{}
	}

	/*
	 * Make a new temporary directory (or use the one specified on the command line)
	 */
	var tmpdir string
	if *temp_dir_path == "" {
		var err error
		tmpdir, err = ioutil.TempDir(".", "loupe-tmp")
		if err != nil {
			log.Printf("cannot get temporary directory: %v", err)
			os.Exit(1)
		}
	} else {
		tmpdir = *temp_dir_path
	}

	/* Process the data */
	/* TODO: We should report some basic statistics WRT how much
	 * data was processed and how long it took.
	 */

	log.Printf("Preprocessor converting %v (with %v and %v) to %v",
		*input_vcf_path,
		*input_bed_path,
		*input_refseq_path,
		*output_loupe_path)

	err := preprocessor.PreprocessRunData(tmpdir,
		*input_vcf_path,
		*input_bed_path,
		*input_refseq_path,
		*input_bam_path,
		*input_structvar_path,
		*input_shortstructvar_path,
		*input_breakpoints_path,
		*input_targets_path,
		*input_fragments_path,
                *input_contig_path,
		*output_loupe_path,
		extra_summary_files)

	if *profile_path != "" {
		pprof.StopCPUProfile()
	}
	if *mem_profile_path != "" {
		mprofile_file, err := os.Create(*mem_profile_path)
		if err != nil {
			log.Printf("Error: %v", err)
			os.Exit(1)
		}

		pprof.WriteHeapProfile(mprofile_file)
		mprofile_file.Close()

	}

	if err == nil {
		log.Printf("Success!")
		os.Exit(0)
	} else {
		log.Printf("Preprocessing failed: %v", err)
		os.Exit(1)
	}

}
