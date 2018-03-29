// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strconv"
	"strings"

	"cnv/commonFunc"
)

var (
	PhasingThr float64 = 0.995
)

type FragmentSet struct {
	Chroms []string
	Starts []int
	Ends   []int
}

type PhasedFragmentSet struct {
	Chroms []string
	Starts []int
	Ends   []int
	HPs    []int
}

func NewFragmentSet() *FragmentSet {
	fragS := FragmentSet{}
	fragS.Starts = make([]int, 0)
	fragS.Ends = make([]int, 0)
	fragS.Chroms = make([]string, 0)
	return &fragS
}

func ReadFrags(csvFile string) *FragmentSet {
	fragS := FragmentSet{}

	ffrag, err1 := os.Open(csvFile)
	//fmt.Println(csvFile)
	defer ffrag.Close()
	commonFunc.CheckErr(err1)

	scanFrag := bufio.NewScanner(ffrag)
	fragS.Starts = make([]int, 0)
	fragS.Ends = make([]int, 0)
	fragS.Chroms = make([]string, 0)

	cnt := 0
	scanFrag.Scan()
	for scanFrag.Scan() {
		its := strings.Split(scanFrag.Text(), ",")
		chrom := its[1]
		mid, _ := strconv.Atoi(its[3])
		s, _ := strconv.Atoi(its[4])
		e, _ := strconv.Atoi(its[2])

		for mid > cnt {
			fragS.Chroms = append(fragS.Chroms, "")
			fragS.Starts = append(fragS.Starts, -1)
			fragS.Ends = append(fragS.Ends, -1)
			cnt++
		}
		fragS.Chroms = append(fragS.Chroms, chrom)
		fragS.Starts = append(fragS.Starts, s)
		fragS.Ends = append(fragS.Ends, e)
		cnt++
		if cnt%1000000 == 0 {
			fmt.Println(cnt)
		}
		if TESTING && cnt == testNumFrag {
			break
		}
	}

	fmt.Println("read in ", len(fragS.Starts), "fragments")
	return &fragS
}

func ReadPhasedFrags(csvFile string) *PhasedFragmentSet {
	fragS := PhasedFragmentSet{}

	ffrag, err1 := os.Open(csvFile)
	defer ffrag.Close()
	commonFunc.CheckErr(err1)

	gfrag, gzerr := gzip.NewReader(ffrag)
	defer gfrag.Close()
	commonFunc.CheckErr(gzerr)

	scanFrag := bufio.NewScanner(gfrag)
	fragS.Starts = make([]int, 0)
	fragS.Ends = make([]int, 0)
	fragS.Chroms = make([]string, 0)
	fragS.HPs = make([]int, 0)

	cnt := 0
	scanFrag.Scan()
	for scanFrag.Scan() {
		its := strings.Split(scanFrag.Text(), "\t")
		chrom := its[0]
		s, _ := strconv.Atoi(its[1])
		e, _ := strconv.Atoi(its[2])
		var probs [2]float64
		probs[0], _ = strconv.ParseFloat(its[7], 64)
		probs[1], _ = strconv.ParseFloat(its[8], 64)

		hp := -1
		if probs[0] >= PhasingThr {
			hp = 1
		}
		if probs[1] >= PhasingThr {
			hp = 2
		}

		fragS.Chroms = append(fragS.Chroms, chrom)
		fragS.Starts = append(fragS.Starts, s)
		fragS.Ends = append(fragS.Ends, e)
		fragS.HPs = append(fragS.HPs, hp)
		cnt++
		if cnt%1000000 == 0 {
			fmt.Println(cnt)
		}
		if TESTING && cnt == testNumFrag {
			break
		}
	}

	fmt.Println("read in ", len(fragS.Starts), "phased fragments")
	return &fragS
}
