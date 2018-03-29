// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This provides facilities to encode a list of candidate structural variants
 * in a bed-like format into a JSON object that loupe can understand.
 *
 * The input is a bed file of breakpoint regions. The "name" field in said bed file
 * will be shared by exactly the two regions that participate in the structural
 * variant. Thus our job is to parse the file and collate rows based on the
 * value of the name field.
 */

package formats

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"strings"
)

/*
 * Holds data related to structural variant phasing
 */
type BreakpointDetails struct {
	Name    string
	Details map[string]string
}

/*
 * Holds data for a single structural variant.
 */
type StructuralVariant struct {
	/* An (arbitrary) name that uniquely identifies this variant */
	Name string
	/* First breakpoint */
	Chromosome1 string
	Start1      int
	Stop1       int

	/* Second breakpoint */
	Chromosome2 string
	Start2      int
	Stop2       int

	/* Other data */
	Quality int

	Filters []string
	Info    map[string]string
	Details []map[string]string
}

/*
 * Load a file of structural variants breakpoints and return a
 * list of SVs.
 */
func LoadSVs(svPath string) (*BundledTrackData, error) {
	/* Grab the contents of the file */
	variants, err := LoadBreakpoints(svPath)
	if err != nil {
		return nil, err
	}
	track := Variants2Track(variants)
	return track, nil
}

func LoadSVsWithDetails(svPath string, detailsPath string) (*BundledTrackData, error) {
	variants, err := LoadBreakpointsWithDetails(svPath, detailsPath)
	if err != nil {
		return nil, err
	}
	track := Variants2Track(variants)
	return track, nil
}

func LoadBreakpointsWithDetails(svPath string, detailsPath string) ([]*StructuralVariant, error) {
	variants, err := LoadBreakpoints(svPath)
	if err != nil {
		return nil, err
	}
	if details, err := LoadDetails(detailsPath); err == nil {
		/* Use a naive O(n^2) algorithm to match SVs with their breakpoint details.
		 * The number of SVs and details is small (< 1000) so this should be ok. */
		for _, variant := range variants {
			variant.Details = []map[string]string{}
			found := false
			for _, detail := range details {
				if variant.Name == detail.Name {
					variant.Details = append(variant.Details, detail.Details)
					found = true
				}
			}
			if !found {
				log.Printf("Warning: Cannot find structural variant %s in breakpoint details file",
					variant.Name)
			}
		}
	} else {
		log.Printf("Cannot load structural variants with breakpoint details: %v", err)
		return nil, err
	}
	return variants, nil
}

/*
 * Serialize a StructuralVariant array to JSON
 */
func WriteSVDataToFile(path string, svs []*StructuralVariant) error {
	fp, err := os.Create(path)
	if err != nil {
		log.Printf("cannot create file: %v", err)
		return err
	}

	defer fp.Close()

	js, err := json.Marshal(svs)
	if err != nil {
		panic(err)
	}
	fp.Write(js)
	return nil
}

func Variants2Track(svs []*StructuralVariant) *BundledTrackData {
	gtd := make([]GenericTrackData, 0, 2*len(svs))

	for index := 0; index < len(svs); index++ {
		sv := svs[index]
		gtd = append(gtd, GenericTrackData{NormalizeChromosomeName(sv.Chromosome1),
			sv.Start1,
			sv.Stop1,
			sv.Name,
			sv})
		gtd = append(gtd, GenericTrackData{NormalizeChromosomeName(sv.Chromosome2),
			sv.Start2,
			sv.Stop2,
			sv.Name,
			sv})
	}
	return GenerateTrackIndex("structvars", gtd)
}

/*
 * Load a breakpoint details file from disk.
 */
func LoadDetails(path string) ([]*BreakpointDetails, error) {
	details := []*BreakpointDetails{}
	fieldnames := []string{"chrom1", "start1", "end1", "chrom2", "start2", "end2", "sv", "bkpt_number",
		"phase_set", "phase_set_start", "phase_set_end", "called_hap", "qual", "sv_bcs", "hap1_bcs",
		"hap2_bcs"}
	err := ReadGenericTsvFile(path, fieldnames, func(line_number int, info map[string]string) error {
		if len(info) < len(fieldnames) {
			return &TsvParseError{path, line_number}
		}

		/* Filter out header line */
		header := true
		for _, fieldname := range fieldnames {
			if info[fieldname] != fieldname {
				header = false
				break
			}
		}
		if !header {
			/* If this line contains data, add it to breakpoint details. */
			name := info["sv"]

			/* Ignore fields "chrom1" through "sv" since those are already
			 * contained in StructuralVariant struct. */
			detail := map[string]string{}
			for i := 7; i < len(fieldnames); i++ {
				detail[fieldnames[i]] = info[fieldnames[i]]
			}
			details = append(details, &BreakpointDetails{name, detail})
		}

		return nil
	})
	return details, err
}

/*
 * Load a breakpoints file from disk.
 */
func LoadBreakpoints(path string) ([]*StructuralVariant, error) {

	fp, err := os.Open(path)
	if err != nil {
		log.Printf("Cannot open: %v: %v", path, err)
		return nil, err
	}

	defer fp.Close()

	variants := make([]*StructuralVariant, 0, 0)
	input_buffer := bufio.NewReader(fp)
	var line_number = 0
	for {
		line_number++
		line, err := input_buffer.ReadString('\n')

		if err != nil {
			break
		}

		if line[0] == '#' {
			continue
		}

		bp := new(StructuralVariant)
		var junk1, info, filters string

		_, err = fmt.Sscanf(line, "%s%d%d%s%d%d%s%d%s%s%s%s",
			&(bp.Chromosome1),
			&(bp.Start1),
			&(bp.Stop1),
			&(bp.Chromosome2),
			&(bp.Start2),
			&(bp.Stop2),
			&(bp.Name),
			&(bp.Quality),
			&(junk1),
			&(junk1),
			&(filters),
			&(info))

		if err != nil {
			log.Printf("Error: %v. Cannot understand SV file (%v) at line %v: %v",
				err, path, line_number, line)
			continue
		}

		bp.Filters = ParseFilterBlock(filters)
		bp.Info = ParseInfoBlock(info)

		variants = append(variants, bp)
	}

	if len(variants) == 0 {
		log.Printf("Warning: Did not successfully parse any SV breakpoints.")
	}
	return variants, nil
}

func ParseFilterBlock(filters string) []string {
	filter_array := strings.Split(filters, ";")

	if filters == "." || filters == "" || len(filter_array) == 0 {
		return []string{}
	}

	return filter_array
}

func ParseInfoBlock(info string) map[string]string {

	result := make(map[string]string)
	fields := strings.Split(info, ";")

	if len(fields) == 0 || info == "" || info == "." {
		return result
	}

	for _, kv := range fields {
		kva := strings.Split(kv, "=")
		if (len(kva) > 1) {
			result[kva[0]] = kva[1]
		}
	}
	return result
}
