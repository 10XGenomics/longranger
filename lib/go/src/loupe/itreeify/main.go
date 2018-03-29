// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This is a helper program that computes "track" data for exons
 * from the genes/exons file and stores that as a separate file.
 */
package main

import (
	"flag"
	"log"
	"loupe/formats"
	"os"
)

var input_path = flag.String("bed", "", "input BED file")
var output_path = flag.String("out", "", "output itree file")

func main() {
	flag.Parse()

	if *input_path == "" || *output_path == "" {
		log.Printf("Bad arguments")
		os.Exit(1)
	}

	log.Printf("Loading exon data from: %v", *input_path)
	gene_data := formats.ReadBedFile(*input_path)

	log.Printf("Computing itree")
	ExonTrackData := formats.GenerateTrackIndex("exon-track", gene_data.RawExonData)

	log.Printf("Writing output to: %v", output_path)
	formats.WriteTrackDataToFile(*output_path, ExonTrackData)
}
