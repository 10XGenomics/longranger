// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package formats

import (
	"bufio"
	"fmt"
)

type Fragment struct {
	Barcode    TenXBarcode
	PhaseBlock int
	Haplotype  byte
}

/*
This function loads all of the fragments from a fragments.tsv file into RAM */
func LoadFragmentsFile(path string) ([]GenericTrackData, error) {

	/* Allocate objects out of these giant (and expandable arrays) to keep
	 * allocation overhead down
	 */
	gtd := make([]GenericTrackData, 0, 0)
	gtd_i := 0
	lastchr := ""
	var chromosome string
	var start int
	var end int
	var barcode string
	var phaseblock int
	var prob_0 float32
	var prob_1 float32
	var prob_mixed float32
	var pb_start int
	var pb_end int
	var count int
	err := ReadEvenMoreGenericTsvFile(path, func(lineno int, b *bufio.Reader) error {

		_, err := fmt.Fscanln(b,
			&chromosome,
			&start,
			&end,
			&phaseblock,
			&pb_start,
			&pb_end,
			&barcode,
			&prob_0,
			&prob_1,
			&prob_mixed,
			&count)
		if err != nil {
			return err
		}

		/* This is some GC black magic. When we have two copies of the same string,
		 * adjust the second to be the first (and share memory with it) so that the
		 * latter can be GC'ed
		 */
		chromosome = NormalizeChromosomeName(chromosome)
		if chromosome == lastchr {
			chromosome = lastchr
		} else {
			lastchr = chromosome
		}

		/*
		 * Decide what haplotype this fragment is on
		 */
		var hap byte
		if prob_0 > 0.995 {
			hap = 0
		} else if prob_1 > 0.995 {
			hap = 1
		} else {
			hap = 2
		}

		ftd := Fragment{TenXBarcodeInit(barcode), phaseblock, hap}
		gtd = append(gtd, GenericTrackData{chromosome, start, end, fmt.Sprintf("%s", lineno), &ftd})

		/* Append data to the gtd array */
		/*
			                gtd = gtd[0 : gtd_i+1]
					ftd = ftd[0 : gtd_i+1]
					gtd[gtd_i].Chromosome = chromosome
					gtd[gtd_i].Start = start
					gtd[gtd_i].Stop = end
					gtd[gtd_i].Name = fmt.Sprintf("%d", lineno)
					gtd[gtd_i].Info = &ftd[gtd_i]
					ftd[gtd_i].Haplotype = hap
					ftd[gtd_i].PhaseBlock = phaseblock
					ftd[gtd_i].Barcode.SetBarcode(barcode)
		*/
		gtd_i++
		return nil
	})

	return gtd[0:gtd_i], err

}
