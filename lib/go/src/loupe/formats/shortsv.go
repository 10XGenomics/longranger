// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This provides facilities to encode a list of short structural variants
 * in a bed-like format into a JSON object that loupe can understand.
 *
 * The input is a bed file of breakpoint regions. The "name" field in said bed file
 * will be shared by exactly the two regions that participate in the structural
 * variant. Thus our job is to parse the file and collate rows based on the
 * value of the name field.
 */

package formats

/*
 * Load a file of structural variants breakpoints and return a
 * list of SVs. Short SVs are encoded like SVs in the bedpe file.
 */
func LoadShortSVs(svPath string) (*BundledTrackData, error) {
	/* Grab the contents of the file */
	variants, err := LoadBreakpoints(svPath)
	if err != nil {
		return nil, err
	}
	track := ShortVariants2Track(variants)
	return track, nil
}

func LoadShortSVsWithDetails(svPath string, detailsPath string) (*BundledTrackData, error) {
	variants, err := LoadBreakpointsWithDetails(svPath, detailsPath)
	if err != nil {
		return nil, err
	}
	track := ShortVariants2Track(variants)
	return track, nil
}

func ShortVariants2Track(shortSVs []*StructuralVariant) *BundledTrackData {
	gtd := make([]GenericTrackData, 0, 2*len(shortSVs))

	for index := 0; index < len(shortSVs); index++ {
		shortSV := shortSVs[index]
		gtd = append(gtd, GenericTrackData{NormalizeChromosomeName(shortSV.Chromosome1),
			shortSV.Start1,
			shortSV.Stop2,
			shortSV.Name,
			shortSV})
	}
	return GenerateTrackIndex("shortstructvars", gtd)
}
