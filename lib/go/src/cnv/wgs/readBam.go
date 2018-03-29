// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

// read in bam, output 3 bedGraph for hap1, hap2, unphased
// read chr by chr
// parallelization? here or stage?
// get bins first
// insert level counting when overlap
// unique number of bar code coverage

// overall there are 3 phaasing tracks
// data structure for molecuar
//	mol id, chrom, start, end, phasing
//	slice of intervals
//	slice of int (bool) for covering status

// at each moment, there are 3 phasing tracks and

// when a molecular has been fully processed, the coverage info will be updated
//	to the phasing tracks

import (
	"encoding/binary"
	"fmt"
	"os"

	"commonFunc"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

var (
	//MaxMolLen int = 1000000
	MinGap int = 1000 // minimum gap required to break a molecular into clusters of reads
)

func ReadBamByChrom(file, chrom string) {
	f, err0 := os.Open(file)
	commonFunc.CheckErr(err0)
	defer f.Close()

	bf, err := bam.NewReader(f, 1)
	commonFunc.CheckErr(err)
	defer bf.Close()

	refs := bf.Header().Refs()
	numChr := len(refs)
	chroms := make([]string, numChr)
	chrSizes := make([]int, numChr)

	for i, r := range refs {
		chroms[i] = r.Name()
		chrSizes[i] = r.Len()
	}

	for {
		r, err := bf.Read()
		if err != nil {
			break
		}

		mi, ok := r.Tag([]byte("MI"))
		if !ok || mi == nil {
			continue
		}

		var mid int
		switch len(mi) {
		case 4:
			mid = int(uint8(mi[3]))
		case 5:
			mid = int(binary.LittleEndian.Uint16(mi[3:]))
		case 7:
			mid = int(binary.LittleEndian.Uint32(mi[3:]))
		default:
			mid = -1
			fmt.Println("wrong size", len(mi))
		}

		hp0, hasHP := r.Tag([]byte("HP"))
		hp := -1
		idxTrack := []int{0}
		if hasHP {
			hp = int(hp0.Value().(uint8))
			idxTrack = append(idxTrack, hp)
		}

		if TESTING && mid > testNumFrag {
			fmt.Println("read ", cnt, "reads")
			return fullCoverage
		}

		rid := r.Ref.ID()
		mrid := r.MateRef.ID()
		if rid < 0 || rid >= numChr {
			continue
		}
		if mrid != rid { //|| r.MatePos > r.Pos {
			continue
		}

		start_diff := r.Pos - r.MatePos
		if start_diff < 0 {
			start_diff = -start_diff
		}

}

func ReadBam(file string, bs *BedSet, fragS *FragmentSet) *FullCoverage {
	cnt := 0
	for {
		r, err := bf.Read()
		if err != nil {
			break
		}

		mi, ok := r.Tag([]byte("MI"))
		if !ok || mi == nil {
			continue
		}

		var mid int
		switch len(mi) {
		case 4:
			mid = int(uint8(mi[3]))
		case 5:
			mid = int(binary.LittleEndian.Uint16(mi[3:]))
		case 7:
			mid = int(binary.LittleEndian.Uint32(mi[3:]))
		default:
			mid = -1
			fmt.Println("wrong size", len(mi))
		}

		hp0, hasHP := r.Tag([]byte("HP"))
		hp := -1
		idxTrack := []int{0}
		if hasHP {
			hp = int(hp0.Value().(uint8))
			idxTrack = append(idxTrack, hp)
		}

		if TESTING && mid > testNumFrag {
			fmt.Println("read ", cnt, "reads")
			return fullCoverage
		}

		rid := r.Ref.ID()
		mrid := r.MateRef.ID()
		if rid < 0 || rid >= numChr {
			continue
		}
		if mrid != rid { //|| r.MatePos > r.Pos {
			continue
		}

		start_diff := r.Pos - r.MatePos
		if start_diff < 0 {
			start_diff = -start_diff
		}

		//isProperPairBasic := (r.Flags&sam.Paired > 0) &&
		//	(r.Flags^sam.Unmapped > 0) &&
		//	(r.Flags^sam.MateUnmapped > 0) &&
		//	(r.Flags&sam.Reverse != r.Flags&sam.MateReverse)
		//isProperPair := isProperPairBasic && (start_diff <= 500)

		//fmt.Printf("%016b\t%02b\t%v\t%v\t%v\t%v\n", r.Flags, r.Flags&sam.ProperPair, isProperPair, start_diff, r.Ref.Name(), r.Pos)

		//isProperPair := r.Flags&sam.ProperPair > 0

		//if r.MapQ < 30 {
		//	continue
		//}

		//if !isProperPair || r.Flags >= sam.Secondary {
		if r.Flags >= sam.Secondary {
			continue
		}
		// a major correction
		chrom, insertS, insertE := r.Ref.Name(), r.Pos, r.End()

		if fragS.Chroms[mid] != chrom {
			fmt.Println("Wrong")
			continue
		}
		if _, ok := fullCoverage.AllCnt[chrom]; !ok {
			continue
		}
		cnt++
		if cnt%100000 == 0 {
			fmt.Println("Get read ", cnt)
		}

		if TESTING && cnt == testNumRead {
			break
		}

		allOLBaits, innerOLBaits := bs.ReadOverlap(chrom, insertS, insertE,
			fragS.Starts[mid], fragS.Ends[mid])

		//fmt.Printf("read # %6d\tchrom %s  interval [%8d, %8d]\t[%8d, %8d: %8d] fragment \tmid %d\n",
		//	cnt, chrom, insertS, insertE, fragS.Starts[mid], fragS.Ends[mid],
		//	fragS.Ends[mid]-fragS.Starts[mid]+1, mid)
		//fmt.Println(allOLBaits)
		//fmt.Println(innerOLBaits, "\n")

		if !Silence {
			//fmt.Println(r)
			fmt.Println(chrom, insertS, insertE, mid, allOLBaits, innerOLBaits, "\n")
		}
		for _, i := range allOLBaits {
			for _, track := range idxTrack {
				if _, ok := fullCoverage.AllCnt[chrom][i][track][mid]; !ok {
					fullCoverage.AllCnt[chrom][i][track][mid] = 1
				} else {
					fullCoverage.AllCnt[chrom][i][track][mid]++
				}
			}
		}

		for _, i := range innerOLBaits {
			for _, track := range idxTrack {
				if _, ok := fullCoverage.InnerCnt[chrom][i][track][mid]; !ok {
					fullCoverage.InnerCnt[chrom][i][track][mid] = 1
				} else {
					fullCoverage.InnerCnt[chrom][i][track][mid]++
				}
			}
		}
	}
	fmt.Println("read ", cnt, "reads")
	return fullCoverage
}
