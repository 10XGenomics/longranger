// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"

	"cnv/commonFunc"
	"code.google.com/p/biogo.bam"
	//"github.com/biogo/hts/bam"
	//"github.com/biogo/hts/sam"
)

var (
	BamFile  string = "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/PHASER_SVCALLER_PD/20486/1014.0.13-0/PHASER_SVCALLER_PD/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam"
	BamIndex string = BamFile + ".bai"
)

func main() {
	Chrom := os.Args[1]
	Start, _ := strconv.Atoi(os.Args[2])
	End, _ := strconv.Atoi(os.Args[3])

	if _, err := os.Stat(BamIndex); os.IsNotExist(err) {
		log.Fatal("bam index file " + BamIndex + " does not exist!")
	}

	f, err0 := os.Open(BamFile)
	commonFunc.CheckErr(err0)
	defer f.Close()

	fIdx, err1 := os.Open(BamIndex)
	commonFunc.CheckErr(err1)
	defer fIdx.Close()

	bf, err := bam.NewReader(f)
	commonFunc.CheckErr(err)
	bif, err2 := bam.ReadIndex(fIdx)
	commonFunc.CheckErr(err2)
	//defer bf.Close()

	refs := bf.Header().Refs()
	numChr := len(refs)
	chroms := make([]string, numChr)
	chrSizes := make([]int, numChr)
	chr2ID := make(map[string]int)
	for i, r := range refs {
		chroms[i] = r.Name()
		chrSizes[i] = r.Len()
		chr2ID[r.Name()] = i
	}

	chID := chr2ID[Chrom]
	readIter, err_readIter := bf.Fetch(bif, chID, Start, End)
	commonFunc.CheckErr(err_readIter)
	for readIter.Next() {
		r := readIter.Get()
		if r == nil {
			continue
		}

		// process each read
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

		as := -9999.0
		xs := -9999.0
		dm := -9999.0
		as0, hasAS0 := r.Tag([]byte("AS"))
		xs0, hasXS0 := r.Tag([]byte("XS"))
		dm0, hasDM0 := r.Tag([]byte("DM"))
		if hasAS0 {
			as = float64(as0.Value().(float32))
			//as, err = strconv.ParseFloat(as0.Value().(string), 64)
			//commonFunc.CheckErr(err)
		}
		if hasXS0 {
			xs = float64(xs0.Value().(float32))
			//xs, err = strconv.ParseFloat(xs0.Value().(string), 64)
			//commonFunc.CheckErr(err)
		}
		if hasDM0 {
			dm, err = strconv.ParseFloat(dm0.Value().(string), 64)
			commonFunc.CheckErr(err)
		}

		fmt.Println(r.Name, mi, mid, as, dm, r.Pos, as0, dm0, mi, xs0, xs, string(as0.Type()))

		if r.Pos < Start || r.End() > End {
			continue
		}

		fmt.Println(mid, hp, r.Pos, r.End(), r.Ref.Name(), r.MapQ, r.Flags >= bam.Secondary, math.Min(float64(End), float64(r.End()))-math.Max(float64(Start), float64(r.Pos)))
	}

}
