// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"flag"
	"fmt"
	"os"

	"commonFunc"

	"github.com/biogo/boom"
)

func main() {
	bamFile := flag.String("bam", "", "Input bam file")
	help := flag.Bool("help", false, "Print this usage message")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: %s -target <target file> \n", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	bf, err := boom.OpenBAM(*bamFile)
	commonFunc.CheckErr(err)
	fmt.Println(bf.RefNames())

	chroms := bf.RefNames()
	chrSize := len(chroms)

	var lastMI string = ""

	for cnt := 0; ; {
		r, _, err := bf.Read()
		if err != nil {
			break
		}

		mi, ok := r.Tag([]byte("MI"))
		if !ok || mi == nil {
			continue
		}
		miS := string(mi)
		if r.RefID() < 0 || r.RefID() >= chrSize {
			continue
		}
		if r.RefID() != r.NextRefID() || r.NextStart() > r.Start() {
			continue
		}
		//chrom, insertS, insertE := chroms[r.RefID()], r.NextStart(), r.End()
		cnt++
		if cnt%100000 == 0 {
			fmt.Println("Get read ", cnt)
			fmt.Println(r)
			fmt.Println()
		}
		if miS != lastMI {
			lastMI = miS
		}

		//fmt.Println(chrom, insertS, insertE, bx)
	}
}
