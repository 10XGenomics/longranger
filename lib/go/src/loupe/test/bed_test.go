// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"fmt"
	"io/ioutil"
	"log"
	. "loupe/formats"
	"loupe/preprocessor"
	"testing"
)

type Lookup struct {
	Location   int
	Chromosome string
	Name       string
}

func TestBestX(t *testing.T) {

	bf := ReadBedFile("inputs/bed_test1.bed")
	stuff1 := GenerateTrackIndex("awesome", bf.RawExonData)
	stuff2 := GenerateTrackIndex("awesome", bf.RawGeneData)
	_ = stuff1
	_ = stuff2

	/*
	   TODO: Add some actual um.... tests....
	*/
}

func TestExon1(t *testing.T) {
	bf := ReadBedFile("inputs/bed_exons1.bed")

	stuff1 := GenerateTrackIndex("awesome", bf.RawExonData)

	s := SearchTrackForRange(stuff1, "chr1", 16775588)

	if s[0].Name != "NECAP2" {
		t.Error("0")
	}

	s2 := SearchTrackForRange(stuff1, "chr1", 92145898)

	if len(s2) != 0 {
		t.Error(1)
	}

	s3 := SearchTrackForRange(stuff1, "chr1", 92145901)

	if len(s3) != 4 {
		t.Error(2)
	}

	s4 := SearchTrackForRange(stuff1, "chr1", 33560313)
	if s4[0].Name != "AZIN2" {
		t.Error("1")
	}
}

/*
 * Test our ability to read a BED file and index its intervals
 * correctly
 */
func TestBED1(t *testing.T) {
	bf := ReadBedFile("inputs/bed_test1.bed")

	//log.Printf("%v", bf)

	/* This maps a chromosome and offset to a name that's in
	 * the inputs/bed_test1.bed file.  We'll search for the
	 * chr+offset pair and see if we get back the name
	 */
	questions := []Lookup{
		{84864224, "chr1", "CHR1:84864215-84864372:DNASE2B"},
		{43545, "chr99", ""},
		{20978630, "chr14", "CHR14:20978631-20979280:RNASE10"},
		{20978650, "chr14", "CHR14:20978631-20979280:RNASE10"},
		{20978630, "chr14", "CHR14:20978631-20979280:RNASE10"},
		{20979279, "chr14", "CHR14:20978631-20979280:RNASE10"},
		{20979280, "chr14", ""},
		{20979280, "chr14", ""},
		{9999999999, "chr14", ""},
		{153640389, "chrX", "CHRX:153640290-153640422:DNASE1L1"},
	}
	stuff1 := GenerateTrackIndex("awesome", bf.RawGeneData)

	for _, q := range questions {
		so := SearchTrackForRange(stuff1, q.Chromosome, q.Location)

		/* TODO: There has got to be a way to simplify this
		 * logic */
		if len(so) == 0 {
			if q.Name != "" {
				t.Error(fmt.Sprintf("Fail! Search for %v+%v should return no results. got %v",
					q.Chromosome, q.Location, so))
			}
		} else {
			if len(so) > 1 {
				t.Error(fmt.Sprintf("Got too many results for %v+%v: %v",
					q.Chromosome, q.Location, so))
			}

			if so[0].Name != q.Name {
				t.Error(fmt.Sprintf("Got wrong answer: %v vs %v", so[0], q))
			}
		}
	}

}

func TestBED2(t *testing.T) {
	bf := ReadBedFile("inputs/bed_test7.bed")

	check(t, bf.RawGeneData[0], GenericTrackData{"chr1", 10, 500, "GENEA_CHR1", &GeneExtraData{"-", "x", 0, 0, 0}}, "uhoh0")
	check(t, bf.RawGeneData[1], GenericTrackData{"chr2", 1, 600, "GENEA_CHR2", &GeneExtraData{"-", "x", 0, 0, 1}}, "uhoh1")
	check(t, bf.RawGeneData[3], GenericTrackData{"chr2", 1, 600, "GENEB", &GeneExtraData{"-", "x", 0, 1, 3}}, "uhoh1")
	check(t, bf.RawGeneData[2], GenericTrackData{"chr4_zing", 1, 600, "GENEA_CHR4_ZING", &GeneExtraData{"-", "x", 0, 0, 2}}, "uhoh1")
}

/*
 * Test our ability to read a BED and VCF file and spit ot a gene index */
func TestIntegration1(t *testing.T) {

	temp_dir_name, err := ioutil.TempDir(".", "loupe-tmp")
	if err != nil {
		log.Printf("cannot get temporary directory: %v", err)
	}

	preprocessor.PreprocessRunData(temp_dir_name, "inputs/bigfile.vcf", "inputs/integration1.bed", "", "", "", "", "", "", "", "awesome-output.json", map[string]string{})
	/* XXX: The output here needs to be tested... somehow */
}

/*
 * Build a gene index from some sample files and check that it is correct
 */
func TestIntegration2(t *testing.T) {

	x := ReadBedFile("inputs/bed_test2.bed")

	indexed_bed := GenerateTrackIndex("awesome", x.RawGeneData)

	vcf_data, _:= ReadVCFToArray("inputs/vcf_test2.vcf")

	VCFtoJSON(vcf_data, GetTestWriter("vcftest2.dat"), 5, indexed_bed, nil)


	wanted_counts := map[string]int{
		"C10.700":    5,
		"B600.800":   1,
		"B2000.4000": 4,
		"B3000.3200": 1,
		"B10.3100":   4,
		"C5000.6000": 0,
	}

	for k := range wanted_counts {
		id := indexed_bed.Dictionary[k]
		record := indexed_bed.TrackData[id]
                if record.Info == nil {
                        t.Error("Missing data");
                }
	}

}

func TestIntegration3(t *testing.T) {
	temp_dir_name, err := ioutil.TempDir(".", "loupe-tmp")
	if err != nil {
		log.Printf("cannot get temporary directory: %v", err)
	}

	preprocessor.PreprocessRunData(temp_dir_name, "inputs/vcf_test4.vcf", "inputs/bed_test4.bed", "inputs/refseq1.vcf", "inputs/bam1.bam", "inputs/sv_2.bed", "inputs/bkpt_details_1.tsv", "inputs/targets_test1.bed", "inputs/fragment-1.bed", "", "integration3.json", map[string]string{})

	/* XXX: Some way to test this output? */

}

func TestIntegration4(t *testing.T) {
	temp_dir_name, err := ioutil.TempDir(".", "loupe-tmp")
	if err != nil {
		log.Printf("cannot get temporary directory: %v", err)
	}

	preprocessor.PreprocessRunData(temp_dir_name, "inputs/empty.vcf", "inputs/bed_test4.bed", "inputs/refseq1.vcf", "inputs/bam1.bam", "inputs/sv_2.bed", "inputs/bkpt_details_1.tsv", "inputs/targets_test1.bed", "", "", "integration4.json", map[string]string{})

	/* XXX: Some way to test this output? */

}

