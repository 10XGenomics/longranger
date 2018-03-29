// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

func testeq(t *testing.T, got int, want int, what string) {
	if got != want {
		t.Error("%v: %v != %v", what, got, want)
	}
}

func TestRS1(t *testing.T) {
	vcf, err := ReadVCFToArray("inputs/vcf_test2.vcf")
	if err != nil {
		t.Error("VCF???")
	}

	testeq(t, BsearchVCFFile(vcf, "chr1", 600), 3, "1")
	testeq(t, BsearchVCFFile(vcf, "chr1", 604), 3, "2")
	testeq(t, BsearchVCFFile(vcf, "chr1", 599), 2, "3")
	testeq(t, BsearchVCFFile(vcf, "chr7", 6), 11, "4")
	testeq(t, BsearchVCFFile(vcf, "chr2", 2), 5, "5")
	testeq(t, BsearchVCFFile(vcf, "chr2", 6), 6, "5")
	testeq(t, BsearchVCFFile(vcf, "chr2", 2499), 7, "6")
	testeq(t, BsearchVCFFile(vcf, "chr2", 2500), 8, "7")
	testeq(t, BsearchVCFFile(vcf, "chr2", 2501), 9, "8")
	testeq(t, BsearchVCFFile(vcf, "chr2", 2502), 9, "9")

}

func TestAll(t *testing.T) {
	vcf, err := ReadVCFToArray("inputs/vcf_test2.vcf")
	if err != nil {
		t.Error("VCF???")
	}

	annotation, _ := ComputeVCFAnnotation("inputs/refseq1.vcf", vcf)
	answer := map[string]string{
		"chr1+30": "RS1",
		"chr1+40": "RS2"}

	if !reflect.DeepEqual(annotation, answer) {
		t.Error(fmt.Sprintf("WRONG ANSWER: %v %v", answer, annotation))
	}

	//log.Printf("AAA: %v", annotation)
}
