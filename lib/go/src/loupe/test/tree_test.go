// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	"fmt"
	. "loupe/formats"
	"reflect"
	"testing"
)

/*
 * This is a trivially simple interval search. We try to confirm
 * that this has the same behavior as the faster-better-probably-more-buggy
 * version.
 */
func DumbSimpleSearch(iva []*Interval, location int) []*Interval {
	output := make([]*Interval, 0, 0)

	for i := 0; i < len(iva); i++ {
		if iva[i].Left <= location && iva[i].Right > location {
			output = append(output, iva[i])
		}
	}
	SortIntervalArray(output)
	return output
}

func TestTree1(t *testing.T) {
	//log.Printf("STARTTREETEST1")
	data := []*Interval{
		&Interval{0, 20, "0.20"},
		&Interval{80, 100, "80.100"},
		&Interval{4, 12, "4.12"},
		&Interval{10, 30, "10.30"},
		&Interval{40, 55, "40.55"},
		&Interval{40, 80, "40.80"},
		&Interval{28, 35, "38.35"},
		&Interval{32, 40, "32.40"},
		&Interval{100, 400, "100.400"},
		&Interval{500, 600, "500.600"},
		&Interval{440, 550, "440.550"},
		&Interval{520, 560, "520.560"},
		&Interval{530, 700, "530.700"},
		&Interval{-10, -1, "-10.-1"},
		&Interval{530, 531, "530.531"},
	}

	r := BuildIntervalTree(data)
	for i := -20; i < 1000; i++ {
		/* Look for i using two interval searches.  Get made
		 * if they don't return the same thing
		 */
		a := DumbSimpleSearch(data, i)
		b := SearchIntervalTreeToArray(r, i)
		if !reflect.DeepEqual(a, b) {
			t.Error(fmt.Sprintf("Error for value %v: %v vs %v", i, a, b))
		}
	}

}
