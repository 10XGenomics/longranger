// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This implements a simple VCF to JSON converted and indexer. The idea is to
 * create a large JSON representation of a VCF file. Because this would be
 * unwieldy to parse as a single JSON blob, we cut the file into multiple
 * sections with a few thousand entries each.  Each section, by its self is a
 * complete JSON object.  We compute an index into that file (also JSON) which
 * associates a logical (chr+offset) range of VCF data to a physical offset
 * into the file.
 *
 * TODO: There are a few idioms in use in this file that are pretty mean
 *       to the GC.
 */

package formats

import (
	"fmt"
)

/*
 This struct holds entries that go in the VCF index
*/
type VCFIndexEntry struct {
	/* What is the first coordinate of this chunk */
	FirstPosition string
	/* What is the one-after-last coordinate of this chunk */
	LastPosition string
	/* Where does this chunk start in the JSONified VCF file */
	Where LoupeSection
}

/*
 * This structure contains all of the data that we include on a per-SNP
 * basis in the VCF summary.  This is a copy of the data in the SimpleVCFRow
 * structure with additional fields that are computed by other indexing work.
 * (for example, Genes)
 *
 * TODO: This is gross. Is there a better way to augment the VCF structure
 * that to have to copy it?
 */
type SNPData struct {
	Chromosome   string               // Same as VCF
	Position     int                  // Same as VCF
	Id           string               // Same as VCF
	Reference    string               // Same as VCF
	Sequence     [2]PhasedAlternative // Same as VCF
	WellPhased   bool                 // Same as VCF
	PhaseId      int                  // Same as VCF
	Genes        []string             // What genes overlap this SNP
	PhaseQuality float64              // Same as VCF
}

/*
 * Compute the data that goes into a single SNP record in a loupe
 * file.
 */
func BuildSNPData(v *SimpleVCFRow, s *SNPData, gene_track *BundledTrackData, refseq map[string]string) {

	/* Copy fields from VCF structure that we care about */
	s.Chromosome = v.Chromosome
	s.Position = v.Position
	s.Id = v.Id
	s.Reference = v.Reference
	s.Sequence = v.Sequence
	s.WellPhased = v.WellPhased
	s.PhaseId = v.PhaseId
	s.PhaseQuality = v.PhaseQuality

	/*
	 * If the input VCF file does not specify an ID for this variant,
	 * try to associate it with a rsID.
	 *
	 * TODO: The way that we overwrite ID here is.... a nasty kludge.
	 */
	if s.Id == "" || s.Id == "." {
		s.Id = GetRefIdForVCFRow(refseq, v)
	}
	/* Augment it with an array of genes that overlap this SNP */
	if gene_track != nil {
		s.Genes = make([]string, 0, 0)
		genes := SearchTrackForRange(gene_track, s.Chromosome, s.Position)

		for _, g := range genes {
			/* XXX This is gross! We mutate the original Gene data
			 * here to keep track of the per-gene SNP count!!!
			 */
			g.Info.(*GeneExtraData).SNPCount++
			s.Genes = append(s.Genes, g.Name)
		}
	}
}

/*
 * Write single chink of VCF entries as a complete JSON blob to a file. Additionally,
 * update "index" to include an index record for those entries.
 * This also replaces "this_segment" with a new, blank array suitable for future
 * data.
 */
func WriteOneChunk(
	index *[]VCFIndexEntry,
	this_segment *[]SimpleVCFRow,
	writer *CompressedWriter,
	gene_track *BundledTrackData,
	refseq map[string]string) {

	// Step 0. Annotate data with gene information
	annotated_data := make([]SNPData, len(*this_segment))
	for index, vcf_data := range *this_segment {
		BuildSNPData(&vcf_data, &annotated_data[index], gene_track, refseq)
	}

	section, err := writer.WriteJSON(annotated_data)

	if err != nil {
		panic(err)
	}

	// Step 2: Update the index array with the logical and
	// physical positions of the new data.
	var index_entry VCFIndexEntry
	lastdata := len(*this_segment) - 1

	index_entry.FirstPosition = fmt.Sprintf("%s+%d",
		(*this_segment)[0].Chromosome,
		(*this_segment)[0].Position)

	index_entry.LastPosition = fmt.Sprintf("%s+%d",
		(*this_segment)[lastdata].Chromosome,
		(*this_segment)[lastdata].Position)

	index_entry.Where = section

	// More pain to the GC
	*index = append(*index, index_entry)

	// Step 3: reset this_segment
	// Screw you too, garbage collector
	*this_segment = make([]SimpleVCFRow, 0, 0)
}

func VCFtoJSON(vcf_array []*SimpleVCFRow,
	writer *CompressedWriter,
	chunk_size int,
	gene_track *BundledTrackData,
	refseq map[string]string) LoupeSection {

	/* This index data; built as we parse the VCF file*/
	index := make([]VCFIndexEntry, 0, 0)
	/* The data for a single output JSON segment */
	this_segment := make([]SimpleVCFRow, 0, 0)

	/* Read the VCF file one line at a time. Accumulate data in
	 * this this_segment array. Once it gets big enough, write out
	 * that data as a single JSON chunk, update the index, and reset
	 * this_segment.  This effectively chunks the data into small
	 * easily digestible JSON segments.
	 */
	for _, data := range vcf_array {
		// Do we have enough data?
		if len(this_segment) > chunk_size {
			/*
			 * Yes we do. Write this chunk of data out and update
			 * the index
			 */
			WriteOneChunk(&index,
				&this_segment,
				writer,
				gene_track,
				refseq)
		}

		// No seriously. I'm just doing to this make the GC's
		// life more painful.
		this_segment = append(this_segment, *data)
	}

	/* Write the final, partial chunk */
	if len(this_segment) > 0 {
		WriteOneChunk(&index, &this_segment, writer, gene_track, refseq)
	}

	section, err := writer.WriteJSON(index)
	if err != nil {
		panic(err)
	}

	return section
}
