// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
package formats

/*
 * This implements functions that load a refseq database that
 * associates coordinates with IDs from a VCF file and stores them
 * in a map so that the inverse association can be stored with
 * each SNP in a loupe file.
 */

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

func WriteRSIndexToFile(writer *CompressedWriter, rsindex map[string]string) LoupeSection {
	section, err := writer.WriteJSON(rsindex)

	if err != nil {
		panic(err)
	}
	return section
}

/*
 * This function iterates through a "vcf-like" file that defines rsIDs for
 * genomic coordinates and associates those with the SNPs in the run VCF.
 * It returns a map with the "names" of every SNP that has a name in the
 * format {"chr4+54321": "rs12345", ...}
 */
func ComputeVCFAnnotation(refseq_path string,
	vcf_array []*SimpleVCFRow) (map[string]string, error) {
	output := make(map[string]string)

	/* Give up now if we don't even have any VCF data. BsearchVCFFile
	 * doesn't work of there are no elements in vcf_array.
	 */
	if len(vcf_array) == 0 {
		return output, nil
	}

	/* Iterate over every line in the reference file */
	err := RefseqIterate(refseq_path, func(line int,
		chr string,
		position int,
		name string) {

		offset := BsearchVCFFile(vcf_array, chr, position)
		vcf_entry := vcf_array[offset]
		/* Does this rsID correspond to something in our VCF? */
		if vcf_entry.Position == position && vcf_entry.Chromosome == chr {
			/* Yes! Add it to the map */
			coordinate := fmt.Sprintf("%v+%v", chr, position)
			output[coordinate] = name
		}
	})

	if err != nil {
		return nil, err
	} else {
		return output, nil
	}
}

/*
 * This is a rediculously simple VCF parser that only knows how to parse the
 * first three elements of each line to produce a chr+offset-->rsid mapping
 * from a VCF file waith said data.
 */
func RefseqIterate(path string,
	callback func(line int, chr string, position int, name string)) error {

	/* Open the file */
	fp, err := os.Open(path)
	if err != nil {
		return err
	}

	defer fp.Close()
	line_number := 0
	input_buffer := bufio.NewReader(fp)

	/* Process the file one line at a time */
	for {
		line_number++
		/* Grab the next line from the file */
		line, err := input_buffer.ReadString('\n')

		/* This is probably just an EOF but could be worse. */
		if err != nil {
			break
		}

		/* Ignore comment lines */
		if line[0] == '#' {
			continue
		}

		var chromosome string
		var offset int
		var refid string

		/* Parse the data we care about from this line */
		_, err = fmt.Sscanf(line, "%s%d%s",
			&chromosome,
			&offset,
			&refid)
		if err != nil {
			log.Printf("Trouble parsing %v at %v: %v",
				path, line_number, err)
			continue
		}

		/* Fix chromosome names */
		if !strings.HasPrefix(chromosome, "chr") {
			chromosome = "chr" + chromosome
		}

		if line_number != 0 && (line_number%5000000 == 0) {
			log.Printf("Processed %v lines so far", line_number)
		}

		/* TODO XXX BUG: It is unclear to me if there is an off-by-one
		 * bug here! Do these formats use the same indexing scheme?
		 */
		callback(line_number, chromosome, offset, refid)

	}
	return nil
}

/*
 * This is just a wrapper around looking up something in a map
 */
func GetRefIdForVCFRow(ref_ids map[string]string, row *SimpleVCFRow) string {

	coordinate := fmt.Sprintf("%s+%d", row.Chromosome, row.Position)
	id := ref_ids[coordinate]

	if id != "" {
		return id
	} else {
		return row.Id
	}
}
