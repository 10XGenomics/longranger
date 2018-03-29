// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This file implements a module that parses a BAM file (using the biogo parser)
 * and creates a .loupe-style index and data chunk for said BAM file.
 */

package formats

import (
	"bufio"
	"bytes"
	"code.google.com/p/biogo.bam"
	"io"
	"log"
	"os/exec"
)

var (
	/* What field of the BAM record has the 10X barcode */
	BARCODE_FIELD_NAME   = []byte{'B', 'X'}
	HAPLOTYPE_FIELD_NAME = []byte{'H', 'P'}
)

/*
 * We need a "ReReader" reader object. A ReReader contains a source reader
 * but implements a read function that retries if the source reader returned
 * fewer bytes than requested.
 */
type ReReader struct {
	/* The exact reader that was passed to us */
	Source io.Reader
	/* A buffered reader that we use as our source */
	Buffered io.Reader
}

/*
 * Implement a "ReRead" function. It reads from a Reader and it retries
 * if it gets fewer bytes than it wanted.  It will always completely fill
 * the "data" array or return an error.
 */
func (r ReReader) Read(data []byte) (int, error) {

	/* How much data do we wish to read */
	want := len(data)
	/* Did something go wrong? */
	var err error
	/* Where are we currently reading to in data */
	var offset int
	/* How much data did we just read*/
	var read_len int
	/* Real programmers never use curley braces. They just write all their code as
	 * increment and check conditions in a body-less for loop.
	 */
	for offset = 0; offset < want && err == nil; read_len, err = r.Buffered.Read(data[offset:]) {
		offset += read_len
	}

	return offset, err
}

/*
 * Make a new ReReader from any Reader object
 */
func MakeReReader(src io.Reader) *ReReader {
	return &ReReader{src, bufio.NewReader(src)}
}

/*
 * Search the auxilary tags for the barcode tag, BX.
 */
func FindBarcodeAux(record *bam.Record) string {
	/* Iterate over each AUx tag */
	for _, aux_tag := range record.AuxTags {
		/* The first two bytes of the aux tag
		 * are the tag name. See if they match
		 */
		if bytes.Compare(BARCODE_FIELD_NAME, aux_tag[0:2]) == 0 {
			/* The 3rd byte is the tag type: we don't care
			 * about that.  Convert the data we care about
			 * into a string and return that.
			 * XXX: This conversion is going to cause GC problems
			 */
			return string(aux_tag[3:])
		}
	}
	return ""
}

func ReadTags(record *bam.Record) (string, int, int) {
	var barcode string

	indexes := make([][]byte, 4)
	searches := []string{"AS", "XS", "HP", "XT"}

	for _, aux_tag := range record.AuxTags {
		if bytes.Compare(BARCODE_FIELD_NAME, aux_tag[0:2]) == 0 {
			barcode = string(aux_tag[3:])
		} else {
			for i := 0; i < len(searches); i++ {
				searchfor := searches[i]
				if searchfor[0] == aux_tag[0] && searchfor[1] == aux_tag[1] {
					indexes[i] = aux_tag[3:]
				}
			}
		}
	}
	as := BytesToFloat(indexes[0])
	xs := BytesToFloat(indexes[1])
	haplotype := BytesToInt(indexes[2])
	xt := BytesToInt(indexes[3])

	var flag int
	if haplotype == 1 {
		flag = FLAG_HAPLOTYPE0
	} else if haplotype == 2 {
		flag = FLAG_HAPLOTYPE1
	}

	var xflag int
	if as-xs > 4 {
		xflag |= XFLAG_ALIGN_UNIQUE
		xt = 0
	} else if as-xs > 0 {
		xflag |= XFLAG_ALIGN_LQFIXED
	} else if as-xs == 0 {
		xflag |= XFLAG_ALIGN_ZQFIXED
	} else {
		xflag |= XFLAG_ALIGN_REMAPPED
	}

	if xt > 0 {
		xflag |= XFLAG_ALIGN_LOCAL
	}

	if record.Flags&bam.Reverse != 0 {
		xflag |= XFLAG_DIRECTION_RC
	}
	if record.MapQ < 60 {
		flag |= FLAG_LOW_QUALITY
	}

	if record.MapQ < 30 {
		xflag |= XFLAG_REALLY_LOW_QUAL
	}

	return barcode, flag, xflag
}

/*
 * Implement an "uncompressing" reader that uses the system's
 * gunzip program to uncompress the input file.  We do this
 * because the system gzip is much much faster than the go library
 */
func FastZipReader(path string) (io.Reader, *exec.Cmd, error) {
	/* Setup a "Cmd" structure to describe the command we
	 * want to execute
	 * */
	cmd := exec.Command("gunzip", "-c", path)
	stdout, err := cmd.StdoutPipe()

	if err != nil {
		return nil, nil, err
	}
	/* Start the command */
	err = cmd.Start()

	if err != nil {
		return nil, nil, err
	}

	/* Make a reader that will return data with no short reads from
	 * gzip's stdout */
	return MakeReReader(stdout), cmd, nil
}

/*
 * .Loupe-ify a BAM file.  This involves copying the "interesting" data from the bam
 * file into the .loupe data section using the BlockedIndex mechanism.  Then, we
 * write out the index to a separate file.
 */
func BuildBAMIndexAndData(bam_path string,
	writer *CompressedWriter,
	fragments *[]GenericTrackData) (LoupeSection, *BlockedIndex) {

	//var fragment_haplotype_data map[string]*IntervalTree

	sofar := 0
	/* Get a blocked_index writer object that we'll use to build
	 * up our data, one datum at a time.
	 */

	bib, err := OpenBlockIndex(writer)

	if err != nil {
		log.Printf("ERR:  %v", err)
		panic(err)
	}

	bam_file_reader, cmd, err := FastZipReader(bam_path)
	if err != nil {
		log.Printf("ERR:  %v", err)
		panic(err)
	}
	_ = cmd
	/*
	 * XXX This is troublesome XXX
	 * We would like to Wait() for the gzip subprocess when we're done to
	 * avoid zombie processies. If the process is still trying to generate
	 * data we'll deadlock.  Closing the Reader should make the gunzip give
	 * up but it doesn't work.  To prevent deadlock's, we don't wait for
	 * gzip at all and leak the PID.
	 */
	/*
		        defer cmd.Wait();
			defer bam_file_reader.(*ReReader).Source.(io.ReadCloser).Close();
	*/

	/* Setup a reader for said BAM file */
	bam_reader, err := bam.NewUncompressedReader(bam_file_reader)

	if err != nil {
		log.Printf("ERR:  %v", err)
		panic(err)
	}

	/* Grab the first record from our BAM reader */
	var next_record bam.Record
	_, err = bam_reader.Read(&next_record)
	if err != nil {
		panic("Cant read even a single valid record from your BAM file")
	}

	/* Iterate over the bam file, processing every record (and not ignoring
	 * the record we just got).
	 */
	for ; err == nil; _, err = bam_reader.Read(&next_record) {

		fixed_name := NormalizeChromosomeName(next_record.Ref.Name())

		barcode, flags, xflag := ReadTags(&next_record)
		if barcode != "" {
			/* Add it to the blocked index */
			bib.Add(fixed_name,
				next_record.Pos+1,
				next_record.End()+1,
				barcode,
				flags,
				xflag)
		}

		/* Produce occational progress reports */
		sofar++
		if (sofar % 10000000) == 0 {
			log.Printf("Processed %v reads", sofar)
		}

	}
	/* Did we exit the above loupe because something went wrong? */
	if err != nil && err != io.EOF {
		log.Printf("BAM parsing failed: %v after %v records", err, sofar)
		panic("OMGOMGOMG")
	}

	log.Printf("Parsed %v BAM records", sofar)
	/* Close teh blocked_index writer. This will ensure that all data
	 * has been written and that the final index record is correct.
	 */
	bib.Close()

	/* Write the index */
	section := bib.WriteIndex()
	return section, &bib.IndexDataSoFar
}
