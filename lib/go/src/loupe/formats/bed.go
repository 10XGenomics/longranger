// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This is a parser for a bed-like file format that contains all of the
 * data we need to construct gene and exon maps.  We expect to parse a file
 * in the following format:
 *
 * <refseq gene ID> chr<chromosome-number> <direction> <start> <end>\
 *     <exonstart>,<exonstart>,... <exonend>,<exonend>,... <name>,\
 *     <read-frame>,<read-frame>,...
 *
 * Lines starting with # should be ignored but currently arent.
 *
 * The main entry points are ReadBedFile() which parses a BED file and
 * SearchForEntry() which searches the parsed representation for genes
 * overlapping a given coordinate.
 *
 */

package formats

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"sort"
	"strings"
)

/*
 * This structure represents the entire contents of a BED file.
 */
type BEDFile struct {
	/* An array of every exon in the file */
	RawExonData []GenericTrackData
	/* An array of every gene in the file */
	RawGeneData []GenericTrackData
}

/*
 * Data we stick in the GenericTrackData.Info field for genes
 */
type GeneExtraData struct {
	Direction     string
	AlternateName string
	SNPCount      int
	Row           int
	GeneIndex     int
}

type ExonExtraData struct {
	//ExonNumber int
	//Id         string
	GeneIndex int
}

type RawGene struct {
	AlternateName                     string
	Chromosome                        string
	Start                             int
	End                               int
	Direction                         string
	Name                              string
	exon_starts, exon_ends, exon_orfs string
}

type RawGeneArray []RawGene

func (g RawGeneArray) Len() int {
	return len(g)
}

func (g RawGeneArray) Swap(i, j int) {
	g[i], g[j] = g[j], g[i]
}

func (g RawGeneArray) Less(i, j int) bool {
	gi := g[i]
	gj := g[j]

	if gi.Name < gj.Name {
		return true
	} else if gi.Name > gj.Name {
		return false
	} else {
		if gi.Chromosome < gj.Chromosome {
			return true
		} else if gi.Chromosome > gj.Chromosome {
			return false
		} else {
			return gi.Start < gj.Start
		}
	}
}

func ReadBedFile(path string) *BEDFile {

	r := LoadBEDRecords(path)
	if r == nil {
		panic("ERROR")
	}

	return InterpretBEDRecords(r)
}

func LoadBEDRecords(path string) []RawGene {
	rawdata := make([]RawGene, 0, 0)

	var line_number int

	/* Open the file. */
	fp, err := os.Open(path)
	if err != nil {
		log.Printf("Error: %v", err)
		return nil
	}
	defer fp.Close()

	buffer := bufio.NewReader(fp)

	/* Iterate until we run out of file. */
	for {
		line, err := buffer.ReadString('\n')
		if err != nil {
			break
		}

		line_number++

		var datum RawGene
		/*
		 * Parse data from a line of the file
		 */
		_, err = fmt.Sscanf(line, "%s%s%s%d%d%s%s%s%s",
			&datum.AlternateName,
			&datum.Chromosome,
			&datum.Direction,
			&datum.Start,
			&datum.End,
			&datum.exon_starts,
			&datum.exon_ends,
			&datum.Name,
			&datum.exon_orfs)

		if err != nil {
			log.Printf("Trouble parsing %v at %v: %v",
				path,
				line_number,
				err)
			continue
		}

		/* Normalize gene names to upper case */
		datum.Name = strings.ToUpper(datum.Name)

		datum.Chromosome = NormalizeChromosomeName(datum.Chromosome)

		rawdata = append(rawdata, datum)
	}

	return rawdata
}

func InterpretBEDRecords(rawdata []RawGene) *BEDFile {

	raw_gene_data := make([]GenericTrackData, 0, 0)

	raw_exon_data := make([]GenericTrackData, 0, 0)

	/* A map of gene names that will need to be munged */
	genes_to_munge := make(map[string]bool)

	var rgd_index = 0

	/* Step 1 sort things */

	sort.Sort(RawGeneArray(rawdata))

	/* step 2:Handle the first gene */
	raw_gene_data = append(raw_gene_data, NewGene(&rawdata[0], 0))
	AppendExonsToTrack(&raw_exon_data, &rawdata[0], 0)

	/* Step 3: Iterate over every record */
	for i := 1; i < len(rawdata); i++ {
		/* rgd_index tracks the index into rawdata of the first datum that's part of
		 * this gene.*/
		newgene := rawdata[i]
		lastgene := rawdata[rgd_index]

		/* Compare this record to the previous one to decide if we're supposed to
		 * merge with that record or add a new one.*/
		if newgene.Name != lastgene.Name || newgene.Chromosome != lastgene.Chromosome || newgene.Start > lastgene.End {
			raw_gene_data = append(raw_gene_data, NewGene(&rawdata[i], i))

			/* If this gene name shows up on multiple chromosomes, note that for
			 * future name munging.
			 */
			if newgene.Name == lastgene.Name && newgene.Chromosome != lastgene.Chromosome {
				genes_to_munge[newgene.Name] = true
			}
			rgd_index = i
		}

		/* Extend the entry that we're working on */
		last := len(raw_gene_data) - 1
		if newgene.End > raw_gene_data[last].Stop {
			raw_gene_data[last].Stop = newgene.End
		}

		AppendExonsToTrack(&raw_exon_data, &rawdata[rgd_index], rgd_index)

	}

	/* Munge names for genes that show up on multiple chromosomes */
	for i := range raw_gene_data {
		gene := &raw_gene_data[i]
		if genes_to_munge[gene.Name] {
			gene.Name = gene.Name + "_" + strings.ToUpper(gene.Chromosome)
		}

	}
	AssignGeneRow(raw_gene_data)

	bf := new(BEDFile)

	bf.RawGeneData = raw_gene_data
	bf.RawExonData = raw_exon_data
	return bf
}

/* Convert a raw gene entry into a GenericTrackData object */
func NewGene(rawdata *RawGene, index int) GenericTrackData {
	var d GenericTrackData

	d.Chromosome = rawdata.Chromosome
	d.Start = rawdata.Start
	d.Stop = rawdata.End
	d.Name = rawdata.Name
	d.Info = &GeneExtraData{rawdata.Direction, rawdata.AlternateName, 0, 0, index}

	return d

}

func AppendExonsToTrack(ed *[]GenericTrackData, raw *RawGene, index int) {
	parseexons(ed,
		raw.Chromosome,
		raw.exon_starts,
		raw.exon_ends,
		raw.exon_orfs,
		raw.Name,
		raw.Direction,
		raw.AlternateName,
		index)

}

/*
 * Add a single exon to the exon list
 */
func newexon(exon_list *[]GenericTrackData,
	chromosome string,
	start_s string,
	stop_s string,
	orf_s string,
	gene_name string,
	exd *ExonExtraData) {

	var start, stop, orf int

	/* Try to parse the start, stop, and open-read-frame-offset data */
	_, err1 := fmt.Sscanf(start_s, "%d", &start)
	_, err2 := fmt.Sscanf(stop_s, "%d", &stop)
	_, err3 := fmt.Sscanf(orf_s, "%d", &orf)

	/* If everything parsed OK, add the exon to the exon list */
	if err1 == nil && err2 == nil && err3 == nil {
		*exon_list = append(*exon_list, GenericTrackData{chromosome,
			start,
			stop,
			gene_name,
			exd})
	}
}

/*
 * Add parse the exon fields for a single gene and add them to the exon list
 */
func parseexons(exon_list *[]GenericTrackData,
	chromosome string,
	starts string,
	stops string,
	orfs string,
	gene_name string,
	direction string,
	alternative_name string,
	GeneIndex int) {

	/*
	 * Exon data is stored as three comma separated lists for the start,
	 * end, and ORF of each exon.  We split those lists up and then iterate
	 * over the corresponding entries
	 */
	starts_array := strings.Split(starts, ",")
	stops_array := strings.Split(stops, ",")
	orfs_array := strings.Split(orfs, ",")

	for i := 0; i < len(starts_array); i++ {
		extra_data := new(ExonExtraData)
		if direction == "+" {
			//extra_data.ExonNumber = i + 1
		} else {
			/* WHy this -1? I have no idea! */
			//extra_data.ExonNumber = len(starts_array) - i - 1
		}
		//extra_data.Id = alternative_name
		extra_data.GeneIndex = GeneIndex
		newexon(exon_list,
			chromosome,
			starts_array[i],
			stops_array[i],
			orfs_array[i],
			gene_name,
			extra_data)
	}
}

/*
 * Add a new entry in the genes array for this gene */
func newgene(gene_list *[]GenericTrackData,
	chromosome string,
	direction string,
	start int,
	stop int,
	name string,
	alternate_name string) *GenericTrackData {

	info := GeneExtraData{direction, alternate_name, 0, 0, 0}
	g := GenericTrackData{chromosome,
		start,
		stop,
		name,
		&info}
	*gene_list = append(*gene_list, g)
	return &((*gene_list)[len(*gene_list)-1])

}

type TrackDataSorter []*GenericTrackData

func (t TrackDataSorter) Len() int {
	return len(t)
}

func (t TrackDataSorter) Less(i, j int) bool {
	if t[i].Chromosome < t[j].Chromosome {
		return true
	} else if t[i].Chromosome > t[j].Chromosome {
		return false
	} else {
		return t[i].Start < t[j].Start
	}
}

func (t TrackDataSorter) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

/*
 * This computes a global track positions for all genes. The results are filled in
 * as the "row" field in the GeneExtraData structure. This attempts to lay out genes
 * in a way that minimizes confusing visual overlap between overlapping genes.
 */
func AssignGeneRow(genes []GenericTrackData) {
	sorted_genes := TrackDataSorter(make([]*GenericTrackData, len(genes)))

	/* Step 1: Compute a list of all genes sorted by chromosome and starting
	 * position.
	 */
	for i := 0; i < len(genes); i++ {
		sorted_genes[i] = &genes[i]
	}
	sort.Sort(sorted_genes)

	/* There are three gene tracks. lastend stores the right-most extent used by
	 * each track.
	 */
	var lastend [3]int
	var lastchr string

	/* Iterate over every gene*/
	for i := 0; i < len(sorted_genes); i++ {
		gene := sorted_genes[i]

		/* If we switch chromosomes, reset everything */
		if gene.Chromosome != lastchr {
			lastend[0] = 0
			lastend[1] = 0
			lastend[2] = 0
			lastchr = gene.Chromosome
		}

		goodness := math.MinInt32
		row := 0
		/* Find the 'best' track on which to place this gene.*/
		for j := 0; j < 3; j++ {
			var g = gene.Start - lastend[j]
			if g > goodness {
				goodness = g
				row = j
			}
		}
		/* Place it up and update lastend */
		gene.Info.(*GeneExtraData).Row = row
		lastend[row] = gene.Stop
	}
}
