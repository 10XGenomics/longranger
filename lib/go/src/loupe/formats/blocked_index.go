// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * Block indexes are a particular (and annoyingly wierd) format for
 * storing large amounts of genomic annotation data. A blocked index has
 * two parts: An "index" part and a bunch of "data" segments.
 * Each data segments is an array of 64-bit structures in the following
 * format:
 *  6   5         4         3         2         1         0
 *  3210987654321098765432109876543210987654321098765432109876543210]
 * [0XLRIIIIIIIIIIIIIIIIIIIIIIII00TTTTTTEEEEEEEESSSSSSSSSSSSSSSSSSSS]
 * 0: ALways 0
 * Q: Set of read quality is < 60
 * L: Read from haplotype 0
 * R: Read from haplotyepe 1
 * X: Read has quality below 60
 * I: An identifier for this annotation.
 * E: The 'end' position of this annotation, relative to the segment
 * S: The 'start' position of this annotation, relative to the segment
 * T: exTended flags (see XFLAG defs below)
 *
 * The data within a segment is sorted by the 'Start' field.
 *
 * Each segment has either as most 16K of the above datums.
 *
 * The "index" part contains two subparts: a 'database' which is a map
 * from an indent number to more information about it. We keep that
 * as an interface{}, allowing any data to be added there.
 *
 * The index also contains a per-chromosome array of Chromosome-Offset -->
 * File-Offset associations.  Those associations allow lookup of the segment
 * containing data for a specific range.
 *
 * TODO/BUG: It is possible for a range to overlap multiple segments.  Which will cause
 * results to be missed.  We should ranges by 'end' and put them into all segments that
 * they overlap. The L and R flags can be used to help keep track of this.
 */

package formats

const (
	MAX_RECORDS_PER_SEGMENT = 16*1024 - 1
	MAX_OFFSET_PER_SEGMENT  = 1024*1024 - 1
	BYTES_PER_RECORD        = 8
	FLAG_HAPLOTYPE0         = 1
	FLAG_HAPLOTYPE1         = 2
	FLAG_LOW_QUALITY        = 4 // Set if MAPQ < 60
	FLAG_ALLWAYS_ZERO       = 8
	XFLAG_DIRECTION_RC      = 1
	XFLAG_ALIGN_UNIQUE      = 0
	XFLAG_ALIGN_LQFIXED     = 2
	XFLAG_ALIGN_ZQFIXED     = 4
	XFLAG_ALIGN_REMAPPED    = 6
	XFLAG_ALIGN_LOCAL       = 8
	XFLAG_REALLY_LOW_QUAL   = 16 // Set if MAPQ < 30

)

/*
 * The IndexEntry object is used to make an (sorta) array that associates logical
 * "chromosome" offsets with a physical file offset. Loupe can bsearch a (small)
 * array of these to find the range within the .loupe file that will have the data
 * it is interested in.
 */
type IndexEntry struct {
	/* What is the logical offset */
	ChromosomeOffset int
	/* Where does data for this offset start */
	Where LoupeSection
	//FileOffset int64
	//Length     int
}

/*
 * This forms the entirety of the "index" section of a blocked annotation
 * data structure.
 */
type BlockedIndex struct {
	/* This maps from "idents" to opaque information about this. For right now
	 * the interface is always a string corresponding to the actual barcode.
	 */
	Database []string
	/* An array of IndexEntries for each chromosome.*/
	Index map[string][]IndexEntry
}

/*
 * This object contains all of the intermediate data we need to build up the index
 * and data segments of the blocked index.
 */
type BlockIndexBuilder struct {
	/* The index (that will eventually end up in the .loupe file) */
	IndexDataSoFar BlockedIndex
	/* The file that we are appending compress data to */
	Writer *CompressedWriter
	//DataFile *os.File
	/* The current chromosome we're working on */
	CurChromosome string
	/* How many ranges have been added to this segment */
	OutstandingUnIndexed int
	/* What is the genome offset of the current segment */
	IndexStartsAt int
	/* Keep a map of barcodes to idents */
	BarcodeToIdent map[string]int
	/* What is the next ident to use for a new barcode */
	NextBarcodeIdent int
	/* A pointer to the last completed index record. */
	Last *IndexEntry
	/* A buffer holds data bound for the data section to
	 * reduce calls to Write()
	 */
	Buffer [1024 * 1024 * 8]byte
	/* Where should the next byte written to buffer go */
	UsedBuffer int64
}

/*
 * This assembles data into a single 64-bit integer that we can store
 */
func Assemble(start int, end int, ident int, flags int, xflags int) int64 {
	dist := end - start
	if start < 0 || dist < 0 || end < 0 || ident < 0 || flags < 0 {
		panic("You know what? You don't even deserve an error message for this kind of crap.")
	}
	if start >= 1<<20 || dist >= 1<<8 || ident >= 1<<24 || flags >= 1<<4 || xflags >= 1<<8 {
		panic("offset too big")
	}
	return (int64)(start) | (int64)(dist)<<20 | (int64)(ident)<<36 | (int64)(flags)<<60 | (int64)(xflags)<<28
}

/*
 * Find a ident number for this barcode.  If we've seen the barcode before
 * recycle the number that we've used. Otherwise allocate a new one and write
 * down our choice.
 */
func (bib *BlockIndexBuilder) IdentForBarcode(barcode string) int {
	i, present := bib.BarcodeToIdent[barcode]

	if !present {
		bib.IndexDataSoFar.Database = append(bib.IndexDataSoFar.Database, barcode)
		i = len(bib.IndexDataSoFar.Database) - 1
		bib.BarcodeToIdent[barcode] = i
	}

	if i >= (1 << 24) {
		/* We'll run out of bits in Assemble if we go past this */
		panic("OMG!!!!! You have more than a 16 million barcodes. Thats a lot of barcodes!")
	}
	return i
}

/*
 * Open a file to store the data section of a blocked index and prepare a BlockIndexBuilder
 * data structure.
 * The generate approach here is to call OpenBlockIndex(), then call Add() a bunch of times
 * before finally calling Close(). You can use WriteIndex() to get the index section after
 * you've finished with the block section.
 */
func OpenBlockIndex(writer *CompressedWriter) (*BlockIndexBuilder, error) {

	bf := new(BlockIndexBuilder)
	bf.Writer = writer

	bf.IndexDataSoFar.Database = append(make([]string, 0, 1024*1024), "")
	bf.IndexDataSoFar.Index = make(map[string][]IndexEntry)
	bf.BarcodeToIdent = make(map[string]int, 1024*1024)

	return bf, nil
}

/*
 * Add a new barcode annotation
 */
func (bib *BlockIndexBuilder) Add(chromosome string,
	start int,
	end int,
	barcode string,
	flags int,
	xflags int) error {

	/* If we switch chromosomes, write too many things or go to far, add a new index entry */

	if chromosome != bib.CurChromosome ||
		(end-bib.IndexStartsAt) > MAX_OFFSET_PER_SEGMENT ||
		bib.OutstandingUnIndexed > MAX_RECORDS_PER_SEGMENT {

		bib.FlushBuffer()
		/* Reset some bean counters */
		bib.IndexStartsAt = start
		bib.OutstandingUnIndexed = 0
		bib.CurChromosome = chromosome
	}

	/* Now that we've adjusted the index (if needed) add the packed data */
	bib.OutstandingUnIndexed++

	if start < bib.IndexStartsAt || end < start {
		panic("CANT GO BACKWARDS NO NO NO NO NO")
	}

	/*
	 * Compute the packed representation of this record
	 */
	packed_data := Assemble(start-bib.IndexStartsAt,
		end-bib.IndexStartsAt,
		bib.IdentForBarcode(barcode),
		flags,
		xflags)

	/* Flush bib.Buffer if it is out of space */
	/* Convert to a byte array and write it to the file */
	for j := uint64(0); j < 8; j++ {
		bib.Buffer[bib.UsedBuffer] = byte((packed_data >> (j * 8)) & 0xff)
		bib.UsedBuffer++
	}

	return nil
}

/*
 * Flush the contents of bib.Buffer to the data file.  If force==false,
 * we only flush data if we're out of space.
 */
func (bib *BlockIndexBuilder) FlushBuffer() {

	section, err := bib.Writer.WriteChunk(bib.Buffer[0:bib.UsedBuffer])
	if err != nil {
		panic(err)
	}

	/* Get (or create) a map for this chromosome */
	map_for_chr := bib.IndexDataSoFar.Index[bib.CurChromosome]
	if map_for_chr == nil {
		map_for_chr = make([]IndexEntry, 0)
	}

	/* Add this chr_offset-->file_offset to the index */
	bib.IndexDataSoFar.Index[bib.CurChromosome] = append(map_for_chr,
		IndexEntry{bib.IndexStartsAt, section})

	//bib.DataFile.Write(bib.Buffer[0:bib.UsedBuffer])
	bib.UsedBuffer = 0
}

/*
 * Finalize a BlockBuilder.
 */
func (bib *BlockIndexBuilder) Close() {
	/* Write any outstanding data to the data file */
	bib.FlushBuffer()
}

/*
 * Write the index part to a file
 */
func (bib *BlockIndexBuilder) WriteIndex() LoupeSection {
	section, err := bib.Writer.WriteJSON(bib.IndexDataSoFar)
	if err != nil {
		panic(err)
	}

	return section
}
