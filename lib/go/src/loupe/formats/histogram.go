// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This implements a mechanism to report the histogram of
 * block lengths.
 */

package formats

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
)

type PhaseBlockHistogramData struct {
	PhaseBlockHistogram map[string]int
}

/*
 * Add a single datum to a histogram. We quantize everything by the quantization
 * factor and add the "end-start" to said histogram.
 */
func WriteOut(m map[string]int, quantization int, id int, start int, end int) {
	length := end - start
	/* Ignore phase blocks with only one SNP in them... those dont count */
	if length > 0 {
		bin := fmt.Sprintf("%v", length-(length%quantization))
		m[bin]++
	}
}

/*
 * Return a map of phase-block-size to freequency.
 */
func BuildHistogramFromVCF(vcf_array []*SimpleVCFRow, quantization int) map[string]int {

	histogram := make(map[string]int)
	var lastchr string
	var last_pb_id int
	var last_pb_start int
	var last_pb_end int

	/* Loop over every variant*/
	for _, variant := range vcf_array {
		/* Are we at the end of a phase block */
		if variant.WellPhased && (variant.Sequence[0].Sequence != variant.Sequence[1].Sequence) {
			if lastchr != variant.Chromosome || last_pb_id != variant.PhaseId {
				lastchr = variant.Chromosome
				WriteOut(histogram,
					quantization,
					last_pb_id,
					last_pb_start,
					last_pb_end)
				last_pb_start = variant.Position
				last_pb_id = variant.PhaseId
				last_pb_end = variant.Position
			} else {
				last_pb_end = variant.Position
			}
		}
	}
	WriteOut(histogram,
		quantization,
		last_pb_id,
		last_pb_start,
		last_pb_end)

	return histogram
}

func WritePBHistogram(path string, vcf_array []*SimpleVCFRow) {
	m := BuildHistogramFromVCF(vcf_array, 1000)
	var pbh PhaseBlockHistogramData

	pbh.PhaseBlockHistogram = m
	js, err := json.Marshal(pbh)
	if err != nil {
		panic(err)
	}

	fp, err := os.Create(path)
	if err != nil {
		log.Printf("Cannot open file")
	}
	defer fp.Close()

	fp.Write(js)
}
