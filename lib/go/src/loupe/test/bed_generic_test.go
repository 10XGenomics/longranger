package test

import (
	"log"
	. "loupe/formats"
	"testing"
)

func tassert(t *testing.T, ok bool, str string) {
	if !ok {
		t.Error(str)
	}
}

func TestTargets1a(t *testing.T) {
	i1, err := ReadGenericBedFile("inputs/generic-bed1.bed", []string{"phase_set", "ps_start", "ps_end", "bc"})

	if err != nil {
		t.Error("huh??")
	}

	tassert(t, i1[3].Chromosome == "chr1", "a")
	tassert(t, i1[3].Start == 725026, "b")
	m := (i1[3].Info).(map[string]string)
	log.Printf("%v", m)
	log.Printf("%v", i1)
	tassert(t, m["phase_set"] == "10725156", "c")
	tassert(t, m["bc"] == "GGCAAGAGCCAGGT-1", "d")
}

func TestBJI1(t *testing.T) {
	var gtt = []GenericTrackData{
		{"chr1", 500, 600, "a", nil},
		{"chr2", 800, 900, "a", nil},
		{"chr2", 801, 2000, "a", nil},
		{"chr2", 900, 990, "a", nil},
		{"chr2", 1000, 1200, "a", nil},
		{"chr2", 1500, 1800, "a", nil},
		{"chr3", 800, 900, "a", nil},
		{"chr4", 800, 900, "a", nil}}

	f := GetTestWriter("bj1-file.json")
	idx := NewBlockedIndex()
	AddToBlockedJSONIndex(f, idx, &gtt[0], 2)
	log.Printf("%v", idx)

	tassert(t, idx.CurrentChromosome == "chr1", "1")
	tassert(t, len(idx.Temporary) == 1, "2")
	tassert(t, len(*idx.IndexPerChromosome["chr1"]) == 0, "3")

	AddToBlockedJSONIndex(f, idx, &gtt[1], 2)
	log.Printf("%v", idx)
	tassert(t, idx.CurrentChromosome == "chr2", "4")
	tassert(t, len(idx.Temporary) == 1, "5")
	tassert(t, len(*idx.IndexPerChromosome["chr1"]) == 1, "6")
	tassert(t, len(*idx.IndexPerChromosome["chr2"]) == 0, "7")

	AddToBlockedJSONIndex(f, idx, &gtt[2], 2)
	log.Printf("%v", idx)
	AddToBlockedJSONIndex(f, idx, &gtt[3], 2)
	log.Printf("%v", idx)
	log.Printf("%v", *idx.IndexPerChromosome["chr2"])
	tassert(t, len(*idx.IndexPerChromosome["chr2"]) == 1, "8")

	AddToBlockedJSONIndex(f, idx, &gtt[4], 2)
	log.Printf("%v", *idx.IndexPerChromosome["chr2"])
	tassert(t, len(*idx.IndexPerChromosome["chr2"]) == 1, "9")
	log.Printf("%v", idx)
	AddToBlockedJSONIndex(f, idx, &gtt[5], 2)
	AddToBlockedJSONIndex(f, idx, &gtt[6], 2)
	AddToBlockedJSONIndex(f, idx, &gtt[7], 2)
	log.Printf("%v", idx)

}

func TestBJI2(t *testing.T) {

	i2, _ := ReadGenericBedFile("inputs/targets.bed", []string{})

        f := GetTestWriter("bj2-data.dat");
	BuildBlockedJSONIndex(f, i2)

}
