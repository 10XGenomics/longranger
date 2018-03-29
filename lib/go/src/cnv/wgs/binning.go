// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

// recommended amp 1.5

import (
	"fmt"
	"math"
)

func GetBins(n int, amp float64) (bins []int) {
	// n is the largest possible read count
	numBins := int(math.Floor(InvHarmonic(float64(n)/amp)+0.5)) + 1
	bins = make([]int, numBins)
	for i := 0; i < numBins; i++ {
		bins[i] = int(math.Floor(Harmonic(float64(i))*amp + 0.5))
		fmt.Println(i, Harmonic(float64(i))*amp)
	}
	return bins
}

func Harmonic(x float64) float64 {
	return 0.66666667*math.Pow(x, 1.5) + math.Sqrt(x)/2.0 - 0.2
}

func RHarmonic(x float64) int {
	return int(math.Floor(Harmonic(x) + 0.5))
}

func InvHarmonic(x float64) float64 {
	return math.Pow(1.5*x, 2.0/3.0)
}
