// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"encoding/json"
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

/* Test that ReadVCFWithCallback dies if given a bad file path */
func Test1(t *testing.T) {

	err := ReadVCFWithCallback("badpath", func(a int, b int, d *SimpleVCFRow) bool { return false })

	if err == nil {
		t.Error("incorrect behavior for invalid path")
	}

}

/* Confirm that rendering an array of SimpleVCFRow's to JSON produces
 * plausable output
 */
func Test3(t *testing.T) {
	stuff := []SimpleVCFRow{
		{
			"chr1",
			10501,
			".",
			"T",
			[2]PhasedAlternative{{Sequence: "T"}, {Sequence: "T"}},
			false,
			0,
			0,
			0,
		},
		{
			"chr1",
			10583,
			".",
			"G",
			[2]PhasedAlternative{{Sequence: "A"}, {Sequence: "G"}},
			true,
			0,
			0,
			0,
		},
	}

	str := RenderVCFSliceToJSON(stuff)
	if len(str) < 20 {
		t.Error("JSON text '%v' cant possibly be right", str)
	}
	//log.Printf("JSON: %v", str)
}

/*
 * Open a small VCF file and parse it. Make sure the data meets
 * this hard-coded expectation
 */
func Test2(t *testing.T) {

	rows := make([]SimpleVCFRow, 0, 0)

	ReadVCFWithCallback("inputs/vcf_test1.vcf",
		func(l int, p int, data *SimpleVCFRow) bool {
			rows = append(rows, *data)
			return true
		})

	expected_result :=
		[]SimpleVCFRow{
			{
				"chr1",
				10492,
				".",
				"C",
				[2]PhasedAlternative{{Sequence: "C"}, {Sequence: "C"}},
				false,
				0,
				1778.2,
				0,
			},
			{
				"chr1",
				10501,
				".",
				"T",
				[2]PhasedAlternative{{Sequence: "T"}, {Sequence: "T"}},
				false,
				0,
				1778.2,
				0,
			},
			{
				"chr1",
				10583,
				".",
				"G",
				[2]PhasedAlternative{{Sequence: "A"}, {Sequence: "G"}},
				false,
				0,
				1238.9,
				0,
			},
			{
				"chr1",
				10621,
				".",
				"GTTGCAAAGGCGCGCCGCGCCG",
				[2]PhasedAlternative{{Sequence: "G"}, {Sequence: "G"}},
				false,
				0,
				20330.4,
				0,
			},
		}

	if reflect.DeepEqual(expected_result, rows) {
		// yay!
	} else {
		t.Error(fmt.Printf("Uh oh!!!! \n %v \n IS NOT \n %v\n", rows, expected_result))
	}

}

/*
 * Index a large VCF file.
 * TODO... someone conform that the output is reasonable
 */
func Test4(t *testing.T) {

	vcf, err := ReadVCFToArray("inputs/bigfile.vcf")
	if err != nil {
		t.Error("VCF???")
	}

	VCFtoJSON(vcf, GetTestWriter("temp-vcf.dat"), 100, nil, nil)

}

/*
 * Test code that extracts barcode data
 */
func Test5(t *testing.T) {

	b, ok := ParseBarcode("AAAAA_5;AAAAC_99;AAAAG_50,AAAAT_2")
	if !reflect.DeepEqual(b[0], []TenXBarcode{TenXBarcodeInit("AAAAA"),
		TenXBarcodeInit("AAAAC"),
		TenXBarcodeInit("AAAAG")}) {
		t.Error(fmt.Sprintf("Uh oh! %v is wrong", b[0]))
	}

	if !reflect.DeepEqual(b[1], []TenXBarcode{TenXBarcodeInit("AAAAT")}) {
		t.Error(fmt.Sprintf("Uh oh %v is wrong", b[1]))
	}

	c, ok := ParseBarcode(",")

	if ok != true || len(c[0]) != 0 || len(c[1]) != 0 {
		t.Error("Blank barcode shuuld succeed: %v %v %v", c[0], c[1], ok)
	}

	d, ok := ParseBarcode(",AAAAT_2")
	if !reflect.DeepEqual(d[1], []TenXBarcode{TenXBarcodeInit("AAAAT")}) {
		t.Error(fmt.Sprintf("Uh oh %v is wrong", d[1]))
	}

	if len(d[0]) != 0 {
		t.Error("Bad answer")
	}
}

func Test6(t *testing.T) {

	rows := make([]SimpleVCFRow, 0, 0)

	ReadVCFWithCallback("inputs/vcf_test4.vcf",
		func(l int, p int, data *SimpleVCFRow) bool {
			rows = append(rows, *data)
			return true
		})

	if rows[0].Sequence[0].BarcodeEvidence[0].GetBarcode() != "AAA" {
		t.Error("Bad barcode 1")
	}
	if rows[4].Sequence[1].BarcodeEvidence[1].GetBarcode() != "GCTA" {
		t.Error("Bad barcode 2")
	}
	if rows[4].Sequence[0].BarcodeEvidence[0].GetBarcode() != "AATA" {
		t.Error("Bad barcode 3")
	}
	if rows[4].Sequence[0].BarcodeEvidence[2].GetBarcode() != "ACC" {
		t.Error("Bad barcode 4: %v", rows[4].Sequence[0].BarcodeEvidence[1].Barcode)
	}

}

func Test7(t *testing.T) {
	a, _ := ReadVCFToArray("inputs/vcf_test5.vcf")
	b, _ := ReadVCFToArray("inputs/vcf_test4.vcf")

	/*
		for i := 0; i < len(a); i++ {
			log.Printf("%v vs %v", *a[i], *b[i])
		}
	*/
	if !reflect.DeepEqual(a, b) {
		t.Error(fmt.Sprintf("Not sorting right %v vs %v", a, b))
	}
}

func Test8(t *testing.T) {
	a, _ := ReadVCFToArray("inputs/vcf_test6.vcf")

	b := []SimpleVCFRow{
		{
			"chr1",
			20,
			".",
			"C",
			[2]PhasedAlternative{{Sequence: "T"}, {Sequence: "C"}},
			false,
			0,
			10.123,
			7.4,
		},
	}

	if !reflect.DeepEqual(a[0], &b[0]) {
		t.Error(fmt.Sprintf("Not right right %v vs %v", a[0], b[0]))
	}

}

func Test9(t *testing.T) {
	a, _ := ReadVCFToArray("inputs/vcf_test7.vcf")

	check(t, a[0].Sequence[0].Frequency, 22, "a")
	check(t, a[0].Sequence[1].Frequency, 16, "b")
	check(t, a[1].Sequence[0].Sequence, "A", "c")
	check(t, a[1].Sequence[1].Sequence, "G", "d")
	check(t, a[1].Sequence[0].Frequency, 5, "e")
	check(t, a[1].Sequence[1].Frequency, 25, "f")
	check(t, a[2].Sequence[0].Frequency, 70, "g")
	check(t, a[2].Sequence[1].Frequency, 1, "h")
}

func Test10(t *testing.T) {
	seqs := []string{"ACGGGCTT-12", "9-CCCGT", "ACGT"}
	for _, seq := range seqs {
		barcode := TenXBarcodeInit(seq)

		if barcode.GetBarcode() != seq {
			t.Error(fmt.Sprintf("Barcode sequence %s is wrong", barcode.GetBarcode()))
		}

		barcodeJson, err := json.Marshal(barcode)
		if err == nil {
			var decodedBarcode map[string]string
			if err = json.Unmarshal(barcodeJson, &decodedBarcode); err == nil {
				if _, ok := decodedBarcode["Barcode"]; ok {
					if decodedBarcode["Barcode"] != seq {
						t.Error(fmt.Sprintf("Barcode sequence JSON %v has wrong 'Barcode' field", decodedBarcode))
					}
				} else {
					t.Error(fmt.Sprintf("Barcode sequence JSON %v does not contain 'Barcode' field", decodedBarcode))
				}
			}
		}

		if err != nil {
			t.Error(fmt.Sprintf("Barcode sequence %s could not be encoded in JSON: %s", seq, err.Error()))
		}
	}
}

func testit_vcf(t *testing.T, vcf []*SimpleVCFRow, chr string, position int, chrwant string, positionwant int) {

	i := BsearchVCFFile(vcf, chr, position)
	found := vcf[i].Chromosome == chrwant && vcf[i].Position == positionwant

	if !found {
		t.Error(fmt.Sprintf("search: %v:%v got: %v:%v wnat: %v:%v",
			chr, position,
			vcf[i].Chromosome,
			vcf[i].Position,
			chrwant,
			positionwant))
	}
}

func Test11(t *testing.T) {
	a, _ := ReadVCFToArray("inputs/vcf_test8.vcf")

	testit_vcf(t, a, "chr2_y", 100, "chr2_y", 100)
	testit_vcf(t, a, "chr2_y", 99, "chr2_x", 500)
	testit_vcf(t, a, "chr2_y", 101, "chr2_y", 100)
	testit_vcf(t, a, "chr2_y", 300, "chr2_y", 300)
	testit_vcf(t, a, "chr2", 300, "chr2", 300)
	testit_vcf(t, a, "chr2", 500, "chr2", 500)
	testit_vcf(t, a, "chr2", 100, "chr2", 100)

}
