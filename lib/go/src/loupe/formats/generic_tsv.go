// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/*
 * This implements a module for parsing generic BED files. A generic BED file
 * has the format:
 *  Chromosome \t Start \t End \t field1 \t field2 \t..... \n
 * The Chromosome, Start, and End are interpreted and genomic coordinates and
 * the fields are interpreted as raw strings.  We use an array of strings
 * to assign names to the different fields.
 */

package formats

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type TsvCallback func(int, map[string]string) error

type TsvParseError struct {
	Path       string
	LineNumber int
}

func (self *TsvParseError) Error() string {
	return fmt.Sprintf("Line %d in bed file %s does not contain proper fields", self.LineNumber, self.Path)
}

/*
 * Does a character separate two fields in a BED file?
 */
func ShouldSplit(r rune) bool {
	return r == ' ' || r == '\t'
}

/*
 * Extract all of the optional fields in a bed file into a map.
 * filenmaes is an array of strings that specifies what to call them
 */
func ExtractFields(fieldnames []string, line string) map[string]string {

	parsed_fields := make(map[string]string)

	trimmed_line := strings.Trim(line, "\n")
	splitted := strings.FieldsFunc(trimmed_line, ShouldSplit)

	for i := 0; i < len(splitted) && i < len(fieldnames); i++ {
		parsed_fields[fieldnames[i]] = splitted[i]
	}
	return parsed_fields
}

func ReadEvenMoreGenericTsvFile(path string, callback func(int, *bufio.Reader) error) error {

	fp, err := os.Open(path)
	if err != nil {
		log.Printf("Error: %v", err)
		return err
	}

	defer fp.Close()

	buffer := bufio.NewReader(fp)
	log.Printf("%v  %v", *fp, path)
	line_number := 0
	for {
		chr, _, err := buffer.ReadRune()

		if err != nil {
			break
		}
		line_number++

		/* Ignore lines starting with # */
		if chr == '#' {
			continue
		}

		buffer.UnreadRune()

		if err := callback(line_number, buffer); err != nil {
			log.Printf("Warning (gen): Error when parsing line %d in tsv file: %v", line_number, err)
		}
		if line_number%1000000 == 0 {
			log.Printf("Read %v lines", line_number)
		}
	}
	return nil
}

func ReadGenericTsvFile(path string, fieldnames []string, callback TsvCallback) error {

	fp, err := os.Open(path)
	if err != nil {
		log.Printf("Error: %v", err)
		return err
	}

	defer fp.Close()

	buffer := bufio.NewReader(fp)
	line_number := 0
	for {
		line, err := buffer.ReadString('\n')
		if err != nil {
			break
		}
		line_number++

		/* Ignore lines starting with # */
		if line[0] == '#' {
			continue
		}

		/* Extract fields */
		info := ExtractFields(fieldnames, line)

		if err := callback(line_number, info); err != nil {
			log.Printf("Warning: Error when parsing line %d in tsv file: %v", line_number, err)
		}
	}
	return nil
}

/*
 * Read a generic bed file into a GenericTrackData array.
 * fieldnames is an array of strings that specify the names of the fields
 * after the first three mandatory fields.  They will be copied into the info
 * section of the GenericTrackData objects as a string-->string map.
 */
func ReadGenericBedFile(path string, fieldnames []string) ([]GenericTrackData, error) {
	fieldnames = append([]string{"chrom1", "start1", "end1"}, fieldnames...)
	gtd := make([]GenericTrackData, 0, 0)
	err := ReadGenericTsvFile(path, fieldnames, func(line_number int, info map[string]string) error {
		if len(info) < 3 {
			return &TsvParseError{path, line_number}
		}
		chromosome := info["chrom1"]
		start, err := strconv.Atoi(info["start1"])
		if err != nil {
			return err
		}
		end, err := strconv.Atoi(info["end1"])
		if err != nil {
			return err
		}

		/* Remove chromosome, start, end from info dictionary */
		bedInfo := map[string]string{}
		for i := 3; i < len(fieldnames); i++ {
			bedInfo[fieldnames[i]] = info[fieldnames[i]]
		}

		gtd = append(gtd, GenericTrackData{chromosome, start, end, fmt.Sprintf("%v", line_number), bedInfo})
		return nil
	})
	return gtd, err
}
