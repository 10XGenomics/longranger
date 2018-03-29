// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	//"reflect"
	"fmt"
	"testing"
)

func init() {
	OverlapThr = 50
	FragEndBuffer = 500
}

func TestReadBam(t *testing.T) {
	//Silence = false

	bs := &BedSet{
		Starts: map[string][]int{"chr1": []int{12080, 12595, 13346, 13425, 17359}},
		Ends:   map[string][]int{"chr1": []int{12200, 12715, 13466, 13545, 17479}},
	}
	fs := ReadFrags("data/test_fragments.csv")

	allFragCnt, innerFragCnt := bs.InnerFragCov(fs)
	innerInfo, allInfo := ReadBam("data/test2.bam", bs, fs)

	allCnt := make(map[string][]int)
	innerCnt := make(map[string][]int)

	for k, v := range innerInfo {
		innerCnt[k] = make([]int, len(v))
		for i, m := range v {
			innerCnt[k][i] = len(m)
		}
	}
	for k, v := range allInfo {
		allCnt[k] = make([]int, len(v))
		for i, m := range v {
			allCnt[k][i] = len(m)
		}
	}

	for i, _ := range allCnt["chr1"] {
		fmt.Printf("%d\t%d\t%d\t%d\n", allCnt["chr1"][i], innerCnt["chr1"][i],
			allFragCnt["chr1"][i], innerFragCnt["chr1"][i])
	}
}
