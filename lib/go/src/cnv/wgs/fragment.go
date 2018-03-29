// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

type Fragment struct {
	FirstLoc   int
	HP         int
	ReadRanges []IntervalBasic
}

func NewFrag(firstLoc, hp int) *Fragment {
	f := new(Fragment)
	f.FirstLoc = firstLoc
	f.HP = hp
	f.ReadRanges = make([]IntervalBasic, 10)
}
