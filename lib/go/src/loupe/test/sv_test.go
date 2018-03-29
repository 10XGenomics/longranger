// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	. "loupe/formats"
	"testing"
)

func TestSV1(t *testing.T) {
	svs, err := LoadBreakpoints("inputs/sv_2.bed")

	check(t, err, nil, "uhoh")

	check(t, len(svs), 7, "uhoh2")

	sv2 := StructuralVariant{"C",
		"chr1", 1500, 2500,
		"chr3", 1000, 2000, 10,
		[]string{"T1", "T2"},
		map[string]string{"NBCS1": "5", "NBCS2": "4", "NOV": "1"},
		nil}

	check(t, svs[2], &sv2, "uhoh3")

	err = WriteSVDataToFile("sv-data.json", svs)

	check(t, err, nil, "uhoh4")

	svt := Variants2Track(svs)

	check(t, svt.TrackData[0].Start, 1000, "uhoh5")
	check(t, svt.TrackData[1].Start, 1010, "uhoh6")
	check(t, svt.TrackData[4].Info, &sv2, "uhoh8")

	check(t, SearchTrackForRange(svt, "chr4", 1001)[0].Name, "D", "uhoh9")

}

func TestSV2(t *testing.T) {
	svs, _ := LoadBreakpoints("inputs/nasty.bedpe")

	Variants2Track(svs)
}

func TestSV3(t *testing.T) {
	details, err := LoadDetails("inputs/bkpt_details_1.tsv")

	check(t, err, nil, "uhoh")

	check(t, len(details), 3, "uhoh2")

	check(t, details[0].Name, "A", "uhoh3")
	check(t, details[2].Name, "D", "uhoh4")

	check(t, details[0].Details["bkpt_number"], "1", "uhoh5")
	check(t, details[0].Details["phase_set"], ".", "uhoh6")

	check(t, details[2].Details["called_hap"], "2", "uhoh7")
}

func TestSV4(t *testing.T) {
	variants, err := LoadBreakpointsWithDetails("inputs/sv_2.bed", "inputs/bkpt_details_1.tsv")

	check(t, err, nil, "uhoh")
	check(t, len(variants), 7, "uhoh2")

	names := []string{"A", "B", "D"}
	for _, variant := range variants {
		for _, name := range names {
			if variant.Name == name {
				check(t, len(variant.Details), 1, "uhoh3")
				if len(variant.Details) > 0 {
					check(t, len(variant.Details[0]), 9, "uhoh4")
				}
			}
		}
	}

	if len(variants[0].Details) > 0 {
		details := variants[0].Details[0]
		check(t, details["bkpt_number"], "1", "uhoh5")
		check(t, details["phase_set_start"], ".", "uhoh6")
	}

	if len(variants[3].Details) > 0 {
		details := variants[3].Details[0]
		check(t, details["called_hap"], "2", "uhoh7")
	}
}
