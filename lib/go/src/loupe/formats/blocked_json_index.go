// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/*
 The Blocked JSON Index implements an index very similar to the Interval Tree.
 It is a mechanism for associating data with intervals on the genome but with
 different performance characteristics. In particular, this has a much lower
 memory overhead but can yield explosive behavior if there is a lot of overlap
 between intervals.

 The blocked JSON index is structured similarly to the Indexed VCF. We write
 intervals as discrete chunks with about 100 intervals stored in each chunk.
 Then we write an index that associates larger regions with a single chunk.
 However, because a interval might cross multiple chunks, we may have to write
 it down multiple times.

 This means that this index can have parasitic behavior in the event that there
 are too many overlapping regions.
*/

package formats

import (
	"log"
	"sort"
)

/*
 This associates a particular region with an array of GenericTrackData objects
 (for teh actual intervals within that region) located at a particular offset
 in the file
*/
type BlockedJSONIndexEntry struct {
	Chromosome string
	Start      int
	End        int
	Where      LoupeSection
}

/*
 This holds all of the data we require to build (and store) a blocked index.
 TODO: Right now this holds both the "internal" state as well as any external
 state that we publish to the .loupe file.
*/
type BlockedJSONIndex struct {
	/* Name of this index */
	KindOfIndex string

	/* Per chromosome arrays of indexes */
	IndexPerChromosome map[string]*[]BlockedJSONIndexEntry

	/* This olds the GenericTrackData for a specific chunk, as it is
	   being built. */
	Temporary []*GenericTrackData
	/* What chromosome are we currently working on */
	CurrentChromosome string
}

/*
 Build a type and associated functions to sort a TrackDataArray
 by chromosome and position
*/
type GenericTrackDataArray []*GenericTrackData

func (gta GenericTrackDataArray) Len() int {
	return len(gta)
}

func (gta GenericTrackDataArray) Swap(i, j int) {
	gta[i], gta[j] = gta[j], gta[i]
}

func (gta GenericTrackDataArray) Less(i, j int) bool {
	d := ChromosomeCMP(gta[i].Chromosome, gta[j].Chromosome)
	if d == 0 {
		return (gta[i].Start < gta[j].Start)
	} else {
		return d < 0
	}

}

/*
 Sort an array of GenericTrackData objects. This returns an array of pointers,
 sorted by the chromosome and start.
*/
func SortGenericTrackData(d []GenericTrackData) []*GenericTrackData {

	var gta GenericTrackDataArray
	gta = make([]*GenericTrackData, len(d))

	for i := 0; i < len(d); i++ {
		gta[i] = &(d[i])
	}

	sort.Sort(gta)
	return gta
}

/*
 * Write the final chunk of data for a chromosome (and update the index)
 */
func FinishChromosome(writer *CompressedWriter, index *BlockedJSONIndex) {

	if len(index.Temporary) > 0 {
		/* Write the chunk to the data file */
		section, err := writer.WriteJSON(index.Temporary)
		if err != nil {
			panic(err)
		}

		/* Calculate the range this chunk overlaps */
		end := index.Temporary[0].Stop
		start := index.Temporary[0].Start
		for i := 0; i < len(index.Temporary); i++ {
			if index.Temporary[i].Stop > end {
				end = index.Temporary[i].Stop
			}
			if index.Temporary[i].Start < start {
				start = index.Temporary[i].Start
			}
		}
		/* Update the index */
		ic := index.IndexPerChromosome[index.CurrentChromosome]

		*ic = append(*ic, BlockedJSONIndexEntry{index.CurrentChromosome,
			start,
			end,
			section})
		/* Truncate scratch space */
		index.Temporary = index.Temporary[0:0]
	}
}

/*
 * Add a single GenericTrackData object to the Blocked JSON index.
 * This works by adding data to a temporary holding array. When the array
 * is large enough (or we switch chromosomes), we write the array out
 * to the data file as a single "chunk" and update the index.
 */
func AddToBlockedJSONIndex(writer *CompressedWriter, index *BlockedJSONIndex, entry *GenericTrackData, max int) {

	/* Flush holding section if we switched chromosomes */
	if index.CurrentChromosome != entry.Chromosome {
		FinishChromosome(writer, index)
	}

	/* Get an index array for this chromosome, creating it if it does not
	 * exist.
	 */

	index_for_chr, exists := index.IndexPerChromosome[entry.Chromosome]
	if !exists {
		log.Printf("Chromosome: %v", entry.Chromosome)
		q := make([]BlockedJSONIndexEntry, 0, 0)
		index_for_chr = &q
		index.IndexPerChromosome[entry.Chromosome] = index_for_chr
	}

	/* Do we need to flush index.Temporary to disk? */
	if len(index.Temporary) >= max {
		/* Yes. Copy it to disk and update the index */
		section, err := writer.WriteJSON(index.Temporary)
		if err != nil {
			panic(err)
		}
		*index_for_chr = append(*index_for_chr, BlockedJSONIndexEntry{entry.Chromosome,
			index.Temporary[0].Start,
			entry.Start,
			section})

		/*
		 * Yuck! Some intervals might extend past the end of this section!!!
		 * take those intervals and copy them into the new section.
		 */
		var j = 0
		for i := 0; i < len(index.Temporary); i++ {
			if index.Temporary[i].Start > entry.Start {
				index.Temporary[j] = index.Temporary[i]
				j++
			}
		}
		/* Truncate 'Temporary' to only include the carryover data */
		index.Temporary = index.Temporary[0:j]
	}

	/* Finally..... Append the new data to temporary holding area and update the
	 * current chromosome as it may have changed.
	 */
	index.Temporary = append(index.Temporary, entry)
	index.CurrentChromosome = entry.Chromosome
}

/* Build a block index for an array of GenericTrackData object.
 * data_path must be the path of the "data file" that goes at the end of the loupe file.
 * index_path must be the path to the index file that will be written with the contents of the
 * index calculated above.
 */
func BuildBlockedJSONIndex(writer *CompressedWriter, track []GenericTrackData) LoupeSection {

	log.Printf("Sorting data...")
	/* Need track sorted by chromosome and position */
	sorted_track := SortGenericTrackData(track)

	index := NewBlockedIndex()

	log.Printf("adding data")
	/* Iterate over GenericTrackData, copying it to data_path and our index */
	for i := 0; i < len(sorted_track); i++ {
		tptr := sorted_track[i]
		AddToBlockedJSONIndex(writer, index, tptr, 10000)
	}
	/* Don't forget to write out the last chunk */
	FinishChromosome(writer, index)

	section, err := writer.WriteJSON(index)
	if err != nil {
		panic(err)
	}
	return section

}

func NewBlockedIndex() *BlockedJSONIndex {

	index := new(BlockedJSONIndex)
	index.IndexPerChromosome = make(map[string]*[]BlockedJSONIndexEntry)

	return index

}
