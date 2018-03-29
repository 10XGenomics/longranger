// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

func check(t *testing.T, got interface{}, want interface{}, txt string) {
	if !reflect.DeepEqual(got, want) {
		s := fmt.Sprintf("%v: got: %v wanted: %v", txt, got, want)
		t.Error(s)
	}
}

func TestPacking(t *testing.T) {

	check(t, Assemble(0, 0, 0, 0), int64(0), "a")
	//                                     FIIIIIIEEEESSSSS
	check(t, Assemble(1, 3, 7, 0), int64(0x0000007000200001), "b")
	check(t, Assemble(0, 0, 0, 1), int64(0x1000000000000000), "c")
}

func TestBarcode(t *testing.T) {

	bib, _ := OpenBlockIndex(GetTestWriter("silly"))

	a1 := bib.IdentForBarcode("AAAA")
	a2 := bib.IdentForBarcode("GGGG")
	a3 := bib.IdentForBarcode("AAAA")
	a4 := bib.IdentForBarcode("CCCC")

	if a1 != a3 {
		t.Error("HUH???")
	}
	if a1 == a2 || a2 == a4 || a1 == a4 {
		t.Error("nope!")
	}

	check(t, bib.IndexDataSoFar.Database[a1], "AAAA", "1")
	check(t, bib.IndexDataSoFar.Database[a2], "GGGG", "2")
}

func TestBlock1(t *testing.T) {

	bib, _ := OpenBlockIndex(GetTestWriter("my-block-index"))


	bib.Add("c1", 0x500, 0x600, "ACC", 4)
	bib.Add("c1", 0x600, 0x700, "ACC", 0)
	bib.Add("c1", 0x5000000, 0x5005000, "ACT", 0)
	bib.Add("c2", 0x100, 0x200, "ACC", 0)
	bib.Add("c2", 0x500, 0x600, "ACCT", 0)

	bib.Close()

        /*

	contents, err := ioutil.ReadFile("my-block-index")

	if err != nil {
		t.Error("2")
	}
	result1 := []int64{
		//0xFIIIIIIEEEESSSSS
		0x4000001010000000,
		0x0000001010000100,
		0x0000002500000000,
		0x0000001010000000,
		0x0000003010000400,
	}
	var result [40]byte

	for i := 0; i < len(result1); i++ {
		for j := 0; j < 8; j++ {
			result[i*8+j] = byte((result1[i] >> (uint64(8 * j))) & 0xff)
		}
	}

	check(t, contents[10:], result[0:40], "OMG")

	check(t, bib.IndexDataSoFar.Database[1], "ACC", "blah1")

	check(t, bib.IndexDataSoFar.Index["c1"][0], IndexEntry{0x500, 10, 2}, "t1")
	check(t, bib.IndexDataSoFar.Index["c1"][1], IndexEntry{0x5000000, 26, 1}, "t2")
	check(t, bib.IndexDataSoFar.Index["c2"][0], IndexEntry{0x100, 34, 2}, "t3")
        */
}
