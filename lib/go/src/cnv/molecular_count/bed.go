// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"cnv/commonFunc"
)

var (
	FragEndBuffer int  = 500
	TESTING       bool = false
	testNumBed    int  = 1000
	testNumFrag   int  = 5000
	testNumRead   int  = 50000
)

type BedSet struct {
	Starts   map[string][]int
	Ends     map[string][]int
	chromIdx string
	baitIdx  int
	Exon     map[string][]string
	ID       map[string][]string
}

func ReadBaits(file string) *BedSet {
	fbait, err1 := os.Open(file)
	defer fbait.Close()
	commonFunc.CheckErr(err1)
	scanBait := bufio.NewScanner(fbait)

	bedSet := BedSet{}

	bedSet.Starts = make(map[string][]int)
	bedSet.Ends = make(map[string][]int)
	bedSet.Exon = make(map[string][]string)
	bedSet.ID = make(map[string][]string)

	cnt := 0

	scanBait.Text()
	for scanBait.Scan() {
		its := strings.Split(scanBait.Text(), "\t")
		sz_its := len(its)
		if BEDTYPE == -1 {
			BEDTYPE = sz_its
		} else if BEDTYPE != sz_its {
			log.Fatal("Rows in the bedfile do not have the same number of fileds")
		}
		chrom := its[0]
		s, _ := strconv.Atoi(its[1])
		e, _ := strconv.Atoi(its[2])

		if _, ok := bedSet.Starts[chrom]; !ok {
			bedSet.Starts[chrom] = make([]int, 0)
			bedSet.Ends[chrom] = make([]int, 0)
			if BEDTYPE >= 4 {
				bedSet.Exon[chrom] = make([]string, 0)
			}
			if BEDTYPE >= 5 {
				bedSet.ID[chrom] = make([]string, 0)
			}
		}

		bedSet.Starts[chrom] = append(bedSet.Starts[chrom], s)
		bedSet.Ends[chrom] = append(bedSet.Ends[chrom], e)
		if sz_its >= 4 {
			bedSet.Exon[chrom] = append(bedSet.Exon[chrom], its[3])
		}
		if sz_its >= 5 {
			bedSet.ID[chrom] = append(bedSet.ID[chrom], its[4])
		}
		cnt++
		if TESTING && cnt == testNumBed {
			break
		}
	}
	fmt.Println("read in ", cnt, "bait")
	return &bedSet
}

func (bs *BedSet) InnerBaits(chrom string, fragS, fragE int) ([]int, []int) {
	fS0, fE0 := checkFragPos(fragS, fragE)
	fS, fE := fS0+FragEndBuffer, fE0-FragEndBuffer

	innerBaits := make([]int, 0)
	allBaits := make([]int, 0)

	if chrom != bs.chromIdx {
		bs.chromIdx = chrom
		bs.baitIdx = 0
	}

	ss := bs.Starts[bs.chromIdx]
	ee := bs.Ends[bs.chromIdx]
	sz := len(ee)

	for bs.baitIdx < sz && fS0 > ee[bs.baitIdx] {
		bs.baitIdx++
	}
	if bs.baitIdx >= sz {
		return allBaits, innerBaits
	}

	for i := bs.baitIdx; i < sz; i++ {
		if ss[i] > fE0 {
			return allBaits, innerBaits
		}

		ol := commonFunc.IntMin(fE0, ee[i]) -
			commonFunc.IntMax(fS0, ss[i]+1)
		if ol >= OverlapThr {
			allBaits = append(allBaits, i)
		}

		if ss[i]+1 >= fS && ee[i] <= fE {
			innerBaits = append(innerBaits, i)
			//if !Silence {
			//	fmt.Println("bed index", i, "fragment", fS, fE, "bed", ss[i], ee[i])
			//}
		}
	}
	return allBaits, innerBaits
}

func (bs *BedSet) GetAllMolInfo(MidFile string, PhasingFile string) (*FullCoverage, *FragmentSet) {
	fc := NewFullCoverage(bs)
	fragS := NewFragmentSet()

	/////////////////////////////////////
	// read in molecular with MID
	/////////////////////////////////////
	ffrag, err1 := os.Open(MidFile)
	defer ffrag.Close()
	commonFunc.CheckErr(err1)
	scanFrag := bufio.NewScanner(ffrag)

	cnt := 0
	scanFrag.Scan()
	for scanFrag.Scan() {
		bs.ResetIdx()
		its := strings.Split(scanFrag.Text(), ",")
		chrom := its[2]
		mid, _ := strconv.Atoi(its[4])
		s, _ := strconv.Atoi(its[5])
		e, _ := strconv.Atoi(its[3])
		fragS.Chroms = append(fragS.Chroms, chrom)
		fragS.Starts = append(fragS.Starts, s)
		fragS.Ends = append(fragS.Ends, e)

		allBaits, _ := bs.InnerBaits(chrom, s, e)
		for _, b := range allBaits {
			fc.AllCnt[chrom][b][0][mid] = 1
		}

		cnt++
		if cnt%1000000 == 0 {
			fmt.Println(cnt)
		}
	}

	fmt.Println("read in ", cnt, "fragments with mid")

	/////////////////////////////////////
	// read in molecular with phasing info
	/////////////////////////////////////
	ffrag2, err2 := os.Open(PhasingFile)
	defer ffrag2.Close()
	commonFunc.CheckErr(err2)
	gfrag, gzerr := gzip.NewReader(ffrag2)
	defer gfrag.Close()
	commonFunc.CheckErr(gzerr)
	scanFrag2 := bufio.NewScanner(gfrag)

	cnt = 0
	scanFrag2.Scan()
	for scanFrag2.Scan() {
		bs.ResetIdx()
		its := strings.Split(scanFrag2.Text(), "\t")
		chrom := its[0]
		fragS, _ := strconv.Atoi(its[1])
		fragE, _ := strconv.Atoi(its[2])
		phaseS, _ := strconv.Atoi(its[4])
		phaseE, _ := strconv.Atoi(its[5])

		var probs [2]float64
		probs[0], _ = strconv.ParseFloat(its[7], 64)
		probs[1], _ = strconv.ParseFloat(its[8], 64)
		mid, _ := strconv.Atoi(its[11])
		hp := -1
		if probs[0] >= PhasingThr {
			hp = 1
		} else if probs[1] >= PhasingThr {
			hp = 2
		} else {
			continue
		}

		if phaseS < fragS {
			phaseS = fragS
		}
		if phaseE > fragE {
			phaseE = fragE
		}
		allBaits, _ := bs.InnerBaits(chrom, phaseS, phaseE)
		for _, b := range allBaits {
			fc.AllCnt[chrom][b][hp][mid] = 1
		}
		cnt++
		if cnt%1000000 == 0 {
			fmt.Println(cnt)
		}
	}
	fmt.Println("read in ", cnt, "fragments with phasing information")
	return fc, fragS
}

func (bs *BedSet) InnerFragCov(fragS *FragmentSet, phasedFragS *PhasedFragmentSet) *SimpleCoverage {
	bs.ResetIdx()

	sc := &SimpleCoverage{}
	sc.AllCnt = make(map[string][][]int)
	sc.InnerCnt = make(map[string][][]int)
	sc.bedset = bs

	for k, v := range bs.Starts {
		sc.AllCnt[k] = make([][]int, len(v))
		sc.InnerCnt[k] = make([][]int, len(v))
		for i := 0; i < len(v); i++ {
			sc.AllCnt[k][i] = make([]int, 3)
			sc.InnerCnt[k][i] = make([]int, 3)
		}
	}

	cnt := 0
	for i, chrom := range fragS.Chroms {
		allBaits, innerBaits := bs.InnerBaits(chrom, fragS.Starts[i], fragS.Ends[i])
		if !Silence {
			ln := fragS.Ends[i] - fragS.Starts[i] + 1
			if ln > 2*FragEndBuffer {
				fmt.Printf("read # %6d\tchrom %s  fragment [%8d, %8d: %8d]\t%v\n",
					i, chrom, fragS.Starts[i], fragS.Ends[i],
					fragS.Ends[i]-fragS.Starts[i]+1, innerBaits)
			}
		}
		for _, b := range innerBaits {
			sc.InnerCnt[chrom][b][0]++
		}
		for _, b := range allBaits {
			sc.AllCnt[chrom][b][0]++
		}
		cnt++
		if TESTING && cnt == testNumFrag {
			break
		}
	}

	bs.ResetIdx()
	cnt = 0
	for i, chrom := range phasedFragS.Chroms {
		allBaits, innerBaits := bs.InnerBaits(chrom, phasedFragS.Starts[i], phasedFragS.Ends[i])
		hp := phasedFragS.HPs[i]
		if hp <= 0 {
			continue
		}
		for _, b := range innerBaits {
			sc.InnerCnt[chrom][b][hp]++
		}
		for _, b := range allBaits {
			sc.AllCnt[chrom][b][hp]++
		}
		cnt++
		if TESTING && cnt == testNumFrag {
			break
		}
	}

	return sc
}

func (bs *BedSet) ReadOverlap(chrom string, insertS, insertE, fragS, fragE int) ([]int, []int) {
	//fmt.Println(chrom, insertS, insertE, mid)

	// this is useful when bin size is smaller than 100 bp.
	if bs.baitIdx <= 3 {
		bs.baitIdx = 0
	} else {
		bs.baitIdx -= 3
	}

	overlappingBaits := make([]int, 0)
	innerBaits := make([]int, 0)

	var fS0, fE0, fS, fE int
	if fragS == -1 || fragE == -1 {
		fS0 = -1
		fE0 = -1
		fS = -1
		fE = -1
	} else {
		fS0, fE0 = checkFragPos(fragS, fragE)
		fS, fE = fS0+FragEndBuffer, fE0-FragEndBuffer
	}

	//fmt.Println(fS, fE)

	if chrom != bs.chromIdx {
		bs.chromIdx = chrom
		bs.baitIdx = 0
	}
	//fmt.Println(bs.chromIdx, baitIdx)

	ss := bs.Starts[bs.chromIdx]
	ee := bs.Ends[bs.chromIdx]
	sz := len(ee)

	for bs.baitIdx < sz && insertS > ee[bs.baitIdx] {
		bs.baitIdx++
		if !Silence {
			fmt.Println(sz, bs.baitIdx)
		}
	}

	//fmt.Println(bs.baitIdx, sz)

	if bs.baitIdx >= sz {
		return overlappingBaits, innerBaits
	}

	if !Silence {
		fmt.Println("bed index", bs.baitIdx, "insert", insertS, insertE,
			"bed", ss[bs.baitIdx], ee[bs.baitIdx])
	}

	for i := bs.baitIdx; i < sz; i++ {
		if ss[i]+1 > insertE {
			return overlappingBaits, innerBaits
		}
		ol := commonFunc.IntMin(insertE, ee[i]) -
			commonFunc.IntMax(insertS, ss[i]+1)
		if !Silence {
			fmt.Println("bed index", i, "insert", insertS, insertE,
				"bed", ss[i], ee[i], "frag", fS, fE, ol)
		}
		if ol >= OverlapThr {
			overlappingBaits = append(overlappingBaits, i)
			if fS != -1 && fE != -1 && ss[i]+1 >= fS && ee[i] <= fE {
				innerBaits = append(innerBaits, i)
			}
		}
	}
	return overlappingBaits, innerBaits
}

func (bs *BedSet) ReadOverlapUniq(chrom string, insertS, insertE, fragS, fragE int) ([]int, []int) {
	//fmt.Println(chrom, insertS, insertE, mid)
	overlappingBaits := make([]int, 0)
	innerBaits := make([]int, 0)

	fS0, fE0 := checkFragPos(fragS, fragE)
	fS, fE := fS0+FragEndBuffer, fE0-FragEndBuffer

	//fmt.Println(fS, fE)

	if chrom != bs.chromIdx {
		bs.chromIdx = chrom
		bs.baitIdx = 0
	}
	//fmt.Println(bs.chromIdx, baitIdx)

	ss := bs.Starts[bs.chromIdx]
	ee := bs.Ends[bs.chromIdx]
	sz := len(ee)

	for bs.baitIdx < sz && insertS > ee[bs.baitIdx] {
		bs.baitIdx++
		if !Silence {
			fmt.Println(sz, bs.baitIdx)
		}
	}

	//fmt.Println(bs.baitIdx, sz)

	if bs.baitIdx >= sz {
		return overlappingBaits, innerBaits
	}

	if !Silence {
		fmt.Println("bed index", bs.baitIdx, "insert", insertS, insertE,
			"bed", ss[bs.baitIdx], ee[bs.baitIdx])
	}

	maxOL := -1
	bestIdx := -1
	//fmt.Println()
	for i := bs.baitIdx; i < sz; i++ {
		if ss[i]+1 > insertE {
			break
		}
		ol := commonFunc.IntMin(insertE, ee[i]) -
			commonFunc.IntMax(insertS, ss[i]+1)
		if !Silence {
			fmt.Println("bed index", i, "insert", insertS, insertE,
				"bed", ss[i], ee[i], "frag", fS, fE, ol)
		}
		//if ol > 0 {
		//	fmt.Println(insertS, insertE, ss[i], ee[i], ol, maxOL, bestIdx, OverlapThr)
		//}
		if ol >= OverlapThr {
			if ol > maxOL {
				maxOL = ol
				bestIdx = i
			}
		}
	}

	//fmt.Println(maxOL, bestIdx)
	if maxOL >= OverlapThr {
		overlappingBaits = append(overlappingBaits, bestIdx)
		if ss[bestIdx]+1 >= fS && ee[bestIdx] <= fE {
			innerBaits = append(innerBaits, bestIdx)
		}
	}
	//fmt.Println(overlappingBaits, innerBaits)
	return overlappingBaits, innerBaits
}

func checkFragPos(fS, fE int) (int, int) {
	if fS > fE {
		fS, fE = fE, fS
		if !Silence {
			fmt.Println("swapping")
		}
	}
	return fS, fE
}

func (bs *BedSet) ResetIdx() {
	bs.chromIdx = ""
	bs.baitIdx = -1
}
