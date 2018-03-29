// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"reflect"
	"testing"
)

var bs BedSet = BedSet{
	Starts: map[string][]int{"chr1": []int{1000, 1100, 1200, 1300, 1400, 1500}},
	Ends:   map[string][]int{"chr1": []int{1120, 1220, 1320, 1420, 1520, 1620}},
}

var fs FragmentSet = FragmentSet{
	Chroms: []string{"chr1", "chr1", "chr1"},
	Starts: []int{0, 700, 100},
	Ends:   []int{3000, 2000, 200},
}

func init() {
	OverlapThr = 50
	FragEndBuffer = 500
}

//func TestInnerBaits(t *testing.T) {
//	bs.ResetIdx()
//	allInnerBeds := make([][]int, 0)
//	_, inner := bs.InnerBaits("chr1", 0, 3000)
//	allInnerBeds = append(allInnerBeds, inner)
//	_, inner = bs.InnerBaits("chr1", 700, 2000)
//	allInnerBeds = append(allInnerBeds, inner)
//	_, inner = bs.InnerBaits("chr1", 1000, 2000)
//	allInnerBeds = append(allInnerBeds, inner)
//
//	allExpected := [][]int{
//		[]int{0, 1, 2, 3, 4, 5},
//		[]int{2, 3},
//		[]int{},
//	}
//
//	for i, innerBed := range allInnerBeds {
//		expected := allExpected[i]
//		if !reflect.DeepEqual(innerBed, expected) {
//			t.Error("(bs *BedSet) InnerBaits miss some target regions")
//			t.Error("expected", expected)
//			t.Error("returned", innerBed)
//		}
//	}
//}
//
//func TestInnerFragCov(t *testing.T) {
//	_, observed := bs.InnerFragCov(&fs)
//	expected := map[string][]int{
//		"chr1": []int{1, 1, 2, 2, 1, 1},
//	}
//	if !reflect.DeepEqual(observed, expected) {
//		t.Error("(bs *BedSet) InnerFragCov has wrong results")
//	}
//
//}

func TestReadOverlap(t *testing.T) {
	allTests := []string{
		"chr1\t800\t1100\t0\t3000",
		"chr1\t1100\t1600\t0\t3000",
		"chr1\t1400\t1600\t0\t3000",
		"chr1\t1100\t1600\t1100\t1600",
		"chr1\t500\t800\t0\t3000",
		"chr2\t1100\t1600\t0\t3000",
	}
	allExpected := [][]int{
		[]int{0},
		[]int{0}, // 800, 1100, 0, 3000
		[]int{1, 2, 3, 4, 5},
		[]int{1, 2, 3, 4, 5}, // 1100, 1600, 0, 3000
		[]int{4, 5},
		[]int{4, 5}, // 1400, 1600, 0, 3000
		[]int{1, 2, 3, 4, 5},
		[]int{}, // 1100, 1600, 1100, 1600
		[]int{},
		[]int{}, // 500, 800, 0, 3000
		[]int{},
		[]int{}, // chr2, 1100, 1600, 0, 3000
	}

	allObserved := make([][]int, 0)
	bs.ResetIdx()
	all, inner := bs.ReadOverlap("chr1", 800, 1100, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlap("chr1", 1100, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlap("chr1", 1400, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlap("chr1", 1100, 1600, 1100, 1600)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlap("chr1", 500, 800, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlap("chr2", 1100, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)

	for i, test := range allTests {
		observed1, observed2 := allObserved[2*i], allObserved[2*i+1]
		expected1, expected2 := allExpected[2*i], allExpected[2*i+1]
		if !reflect.DeepEqual(observed1, expected1) ||
			!reflect.DeepEqual(observed2, expected2) {
			t.Error("(bs *BedSet) ReadOverlap has wrong results")
			t.Error("for test case ", test)
			t.Error("expected all", expected1)
			t.Error("expected inner", expected2)
			t.Error("observed all", observed1)
			t.Error("observed inner", observed2)
			t.Error()
		}

	}
	//Silence = true

}

func TestReadOverlapUniq(t *testing.T) {
	allTests := []string{
		"chr1\t800\t1100\t0\t3000",
		"chr1\t1100\t1600\t0\t3000",
		"chr1\t1400\t1600\t0\t3000",
		"chr1\t1100\t1600\t1100\t1600",
		"chr1\t500\t800\t0\t3000",
		"chr2\t1100\t1600\t0\t3000",
		"chr1\t1150\t1350\t0\t3000",
	}
	allExpected := [][]int{
		[]int{0},
		[]int{0}, // 800, 1100, 0, 3000
		[]int{1},
		[]int{1}, // 1100, 1600, 0, 3000
		[]int{4},
		[]int{4}, // 1400, 1600, 0, 3000
		[]int{1},
		[]int{}, // 1100, 1600, 1100, 1600
		[]int{},
		[]int{}, // 500, 800, 0, 3000
		[]int{},
		[]int{}, // chr2, 1100, 1600, 0, 3000
		[]int{2},
		[]int{2}, // 1150, 1350, 0, 3000
	}

	allObserved := make([][]int, 0)
	bs.ResetIdx()
	all, inner := bs.ReadOverlapUniq("chr1", 800, 1100, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr1", 1100, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr1", 1400, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr1", 1100, 1600, 1100, 1600)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr1", 500, 800, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr2", 1100, 1600, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()
	all, inner = bs.ReadOverlapUniq("chr1", 1150, 1350, 0, 3000)
	allObserved = append(allObserved, all, inner)
	bs.ResetIdx()

	for i, test := range allTests {
		observed1, observed2 := allObserved[2*i], allObserved[2*i+1]
		expected1, expected2 := allExpected[2*i], allExpected[2*i+1]
		if !reflect.DeepEqual(observed1, expected1) ||
			!reflect.DeepEqual(observed2, expected2) {
			t.Error("(bs *BedSet) ReadOverlapUniq has wrong results")
			t.Error("for test case ", test)
			t.Error("expected all", expected1)
			t.Error("expected inner", expected2)
			t.Error("observed all", observed1)
			t.Error("observed inner", observed2)
			t.Error()
		}

	}
	//Silence = true

}

func TestCheckFragPos(t *testing.T) {
	allExpected := [][]int{
		[]int{1000, 2000},
		[]int{1000, 2000},
	}
	allObserved := make([][]int, 0)
	fS, fE := checkFragPos(1000, 2000)
	allObserved = append(allObserved, []int{fS, fE})
	fS, fE = checkFragPos(2000, 1000)
	allObserved = append(allObserved, []int{fS, fE})

	for i, expected := range allExpected {
		observed := allObserved[i]
		if !reflect.DeepEqual(observed, expected) {
			t.Error("checkFragPos has wrong results")
		}
	}
}
