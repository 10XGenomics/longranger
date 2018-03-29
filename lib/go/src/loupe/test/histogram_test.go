package test

import (
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

func TestHist(t *testing.T) {
	vcf, err := ReadVCFToArray("inputs/vcfhist.vcf")

	if err != nil {
		t.Error("VCF???")
	}

	answer := BuildHistogramFromVCF(vcf, 100)

	expected := map[string]int{
		"0":    1,
		"600":  1,
		"1000": 3,
	}

	if !reflect.DeepEqual(answer, expected) {
		t.Error(fmt.Sprintf("UHOH %v != %v", answer, expected))
	}
}
