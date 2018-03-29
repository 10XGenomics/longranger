// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"fmt"
	"os"

	"commonFunc"
	"github.com/biogo/hts/bam"
)

func main() {
	f, err0 := os.Open(os.Args[1])
	commonFunc.CheckErr(err0)
	defer f.Close()

	bf, err := bam.NewReader(f, 1)
	commonFunc.CheckErr(err)
	defer bf.Close()
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

		fmt.Println(mi.String())
		//fmt.Println(r)
		cnt++

		if cnt == 5000 {
			break
		}
	}
}
