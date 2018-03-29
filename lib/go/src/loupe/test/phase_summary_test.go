// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

func npb(chr string, start int, stop int, a int, b int, n int, id int) GenericTrackData {
	extra := &ExtraPhaseBlockData{a, b, n}
	return GenericTrackData{chr, start, stop, fmt.Sprintf("%v", id), extra}
}

func TestPB1(t *testing.T) {
	bf := ReadBedFile("inputs/pb_bed_1.bed")

	bedtrack := GenerateTrackIndex("silly", bf.RawGeneData)

	vcf, err := ReadVCFToArray("inputs/pb_vcf1.vcf")
	if err != nil {
		t.Error("VCF??")
	}

	a, _ := GeneratePhaseBlockSummary(bedtrack, vcf)
	answer := []GenericTrackData{
		npb("chr1", 1200, 1531, 0, 5, 0, 1),
		npb("chr1", 2200, 3111, 3, 1, 0, 3),
		npb("chr1", 3123, 3134, 0, 2, 0, 4),
		npb("chr1", 4143, 4154, 2, 0, 0, 8),
		npb("chr1", 4163, 5154, 3, 0, 0, 9),
		npb("chr2", 500, 901, 3, 0, 0, 20),
		npb("chr3", 500, 551, 2, 0, 0, 21)}
	/*
			for i := 0; i < len(answer); i++ {
		                log.Printf("%v", i);
				if a.TrackData[i] != answer[i] {
					log.Printf("***")
				}

				log.Printf("%v (%v) vs %v", a.TrackData[i], a.TrackData[i].Info, answer[i])
			}
	*/
	if !reflect.DeepEqual(answer, a.TrackData) {
		t.Error(fmt.Printf("Unequal: %v versus %v", answer, a.TrackData))
	}
}
