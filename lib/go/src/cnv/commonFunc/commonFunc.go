// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package commonFunc

import (
	"fmt"
	"os"
)

func IntMin(a, b int) int {
	if a > b {
		return b
	} else {
		return a
	}
}

func IntMax(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

func CheckErr(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}
}
