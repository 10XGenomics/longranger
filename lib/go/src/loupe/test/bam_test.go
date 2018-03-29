// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"bufio"
	"code.google.com/p/biogo.bam"
	. "loupe/formats"
	"os"
	"testing"
)

func TestBam1(t *testing.T) {

	fp, err := os.Open("inputs/bam1.bam")

	if err != nil {
		t.Error("LEARN TO SPELL")
	}

	reader := bufio.NewReader(fp)

	bam_reader, err := bam.NewReader(reader)

	if err != nil {
		t.Error("bam.NewReader -- it is crying: %v", err)
	}

	records := make([]bam.Record, 10)

	for i := 0; i < 10; i++ {
		_, err := bam_reader.Read(&records[i])

		if err != nil {
			t.Error("bam_reader.Read no workie")
		}

		//fmt.Printf("%v\n", (records[i]))
	}
	check(t, records[0].Pos, 1811055, "1")
	check(t, string(records[0].Seq.Expand()), "ATGGCGTGAACCCAGGAGGCGGAGCTTGCAGTGAGCCAAGATCGCACCACTGCACTCCAGCCTGGGTGACAGAGTGAG", "2")
	check(t, FindBarcodeAux(&records[0]), "1-GCAAGGTCTTTACC", "3")
	check(t, records[9].Pos, 8300365, "4")
	check(t, records[9].Ref.Name(), "chr1", "5")
	check(t, string(records[9].Seq.Expand()), "GCCCCACGTGGTGGCTCGCACCTGCATTTCAGCACTTTGGAAGCTGGAGGTGGGAGCATTGCTTGGGCCCAGCGATGG", "6")
	check(t, FindBarcodeAux(&records[9]), "1-TTCGCCTCAGGTTC", "7")

}

func GetTestWriter(file string) *CompressedWriter {
        x, err := NewCompressedWriter(file, 4096);
        if (err != nil) {
                panic(err);
        }

        return x;
}

func TestBam2(t *testing.T) {
	_, bl := BuildBAMIndexAndData("inputs/bam1.bam", GetTestWriter("bam1.dat"), nil)


	check(t, bl.Database[1], "1-GCAAGGTCTTTACC", "a")
	check(t, bl.Database[2], "1-TTGGCAGTGCGTAG", "b")
	check(t, len(bl.Index["chr1"]) > 1, true, "c")

}

func TestEmptyBam(t *testing.T) {
	// Can we read a bam file without crashing?
	_, _ = BuildBAMIndexAndData("inputs/bam1.bam", GetTestWriter("bam2.dat"), nil)

}

func TestFragments1(t *testing.T) {

	d := []GenericTrackData{
		{"chr1", 500, 600, "A", &Fragment{TenXBarcodeInit("ACCT"), 12345, 0}},
		{"chr2", 500, 600, "A", &Fragment{TenXBarcodeInit("TCCT"), 12346, 1}},
		{"chr1", 700, 900, "A", &Fragment{TenXBarcodeInit("ACCT"),12347, 1}},
		{"chr1", 550, 750, "A", &Fragment{TenXBarcodeInit("GCCT"), 12348, 1}}}

	dt := BuildFragmentIndex(d)

	check(t, LookupFragmentHaplotype(dt, "ACCT", "chr1", 575), 1, "a")
	check(t, LookupFragmentHaplotype(dt, "TCCT", "chr1", 575), 0, "b")
	check(t, LookupFragmentHaplotype(dt, "ACCT", "chr1", 275), 0, "c")
	check(t, LookupFragmentHaplotype(dt, "TCCT", "chr2", 575), 2, "d")
	check(t, LookupFragmentHaplotype(dt, "ACCT", "chr1", 775), 2, "e")

}
