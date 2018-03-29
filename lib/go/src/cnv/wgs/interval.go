// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

// basic interval only has start and end locations
// leaving out chromosome information
type IntervalBasic struct {
	Start int
	End   int
}

// will add function to *sort* and *merge* IntervalBasic slice
