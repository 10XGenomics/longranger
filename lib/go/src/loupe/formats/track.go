// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
package formats

import (
	"fmt"
)

/*
 * This module generates loupe-client ready "track" data from an array
 * of ranges.  Each range is specified by an instance of GenericTrackData
 * which at a minimum associates a range with a name, and possibly some
 * other arbitrary data.
 *
 * A bundled-track data object is a collection of GenericTrackData ranges
 * that is indexed for storage in a .loupe file.  A BundledTrackData object
 * can be jsonified and interpreted by the loupe client.
 *
 * A BundledTrackData object is created from an array of GenericTrackData ranges
 * with the GenerateTrackIndex function.  It can be written to a file with the
 * WritrTrackDataToFile function.
 */

/*
 * A bundle of everything loupe needs to render a track. This can
 * be output (in JSON) as a section of the loupe file.
 */
type BundledTrackData struct {
	/* What do we call this track */
	Name string
	/* Associate from names to indices into TrackData */
	Dictionary map[string]int
	/* This array holds the "raw" track data */
	TrackData []GenericTrackData
	/* Associate from an arbitary position or range to a set of
	 * indexes into track data.
	 */
	Itree map[string]*IntervalTree
}

/*
 * This structure holds the basic data we need for each range in a loupe
 * track
 */
type GenericTrackData struct {
	/* What chromosome are we on */
	Chromosome string
	/* What range within the chromosome */
	Start int
	Stop  int
	/* What do we call it? */
	Name string
	/* What other data is along for the ride? This is not interpreted by
	 * the track mechanism but may be used to pass data through to the UI
	 */
	Info interface{}
}

/*
 * Search a BundledTrackData instance for all of the ranges that overlap
 * a given position on a chromosome
 */
func SearchTrackForRange(track *BundledTrackData,
	chromosome string,
	position int) []*GenericTrackData {

	/* Make an array to store the results in */
	results := make([]*GenericTrackData, 0, 0)

	/* Find the itree for the chromosome in question */
	tree := track.Itree[chromosome]
	if tree == nil {
		return results
	}

	/* Find all of the IDs that overlap */
	indexes := SearchIntervalTreeToArray(tree, position)

	/* Convert those IDs to GenericTrackData objects by looking them
	 * up int he TrackData field.
	 */
	for i := range indexes {
		var id int
		fmt.Sscanf(indexes[i].Name, "%d", &id)
		results = append(results, &track.TrackData[id])
	}
	return results
}

/*
 * Index and generate a BdundledTrackData from an array of GenericTrackData
 * objects.
 */
func GenerateTrackIndex(name string,
	track_data []GenericTrackData) *BundledTrackData {

	/* Step 1: Build a name index */

	name_index := make(map[string]int)

	for index, data := range track_data {
		name_index[data.Name] = index
	}

	/*
	 * Step 2: Split up the track_data by chromosome and build
	 * interval structures for each range.
	 */
	track_data_itree := make(map[string][]*Interval)

	for index, data := range track_data {
		chr := data.Chromosome

		/* Instantiate a new interval */
		interval := new(Interval)
		interval.Left = data.Start
		interval.Right = data.Stop
		/* For JS compatability reasons, we store the numerical
		 * ID as a string. For Historical reasons, the field is called
		 * Name
		 */
		interval.Name = fmt.Sprintf("%v", index)

		/* Find or use an array for this chromosome */
		iarray := track_data_itree[chr]
		if iarray == nil {
			iarray = make([]*Interval, 0, 0)
		}

		track_data_itree[chr] = append(iarray, interval)
	}

	/* Step 4: Build a separate interval tree for the data
	 * on each chromosome
	 */
	itree_by_chr := make(map[string]*IntervalTree)

	for chr, iarray := range track_data_itree {
		itree_by_chr[chr] = BuildIntervalTree(iarray)
	}

	/*
	 * Step 4: Build a map with all of the data we'll need
	 */

	bundled := &BundledTrackData{
		name,
		name_index,
		track_data,
		itree_by_chr}
	return bundled

}

/*
 * Write a BundledTrackData instance to a file.
 */
func WriteTrackDataToFile(writer *CompressedWriter, track_data *BundledTrackData) LoupeSection {

	section, err := writer.WriteJSON(track_data)
	if err != nil {
		panic(err)
	}
	return section

}
