// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This file implements a simple interval tree that we use for matching SNP
 * other variations against a gene database.  An interval tree is an efficient
 * mechanism for storing a set of [L...R] intervals and then computing which
 * intervals overlap a given position  The BuildIntervalTree function
 * constructs the interval tree from an array of Intervals and the
 * SearchIntervalTree function finds the intervals that overlap a given
 * position.
 *
 * TODO: The mechanism this uses to generate the interval tree might be O(N^2).
 * TODO: This does not correctly sort the intervals within a given node. It
 * should maintain two sorted arrays, one by 'left' and one by 'right' to
 * enable fast exclusion of intervals that are in the same node.
 * TODO: The Interval struct should be extended with an interface{}-typed field
 * so that we may attach arbitrary data to an interval.
 */

package formats

import (
	"encoding/json"
	"log"
	"os"
	"sort"
)

/*
 * Interval type: associate a left (inclusive) and right (exclusive) boundary
 * with a name
 */
type Interval struct {
	/* 'Left'/'Minimum' value for this interval */
	Left int
	/* 'Right'/Maximum' value for this interval */
	Right int
	/* Name this interval is associated with */
	Name string
}

/*
 * The 'IntervalArray' type and Len/Swap/Less functions define the primitives
 * so that we can call sort.Sort over intervals to sort arrays by their
 * midpoint.
 */
type IntervalArray []*Interval

func (ia IntervalArray) Len() int {
	return len(ia)
}

func (ia IntervalArray) Swap(i, j int) {
	ia[i], ia[j] = ia[j], ia[i]
}

func (ia IntervalArray) Less(i, j int) bool {
	return (ia[i].Left + ia[i].Right) < (ia[j].Left + ia[j].Right)
}

/*
 * This function defines a node in an interval tree
 */
type IntervalTree struct {
	/* Midpoint represented by this node */
	MidPoint int
	/* Array of all intervals that overlap midpoint */
	/* TODO: We should store two arrays here. One sorted by 'left' and
	 * the other sorted by 'right'. This will allow more efficient searching
	 */
	DataHere []*Interval
	/* Tree of all intervals that are fully to the left of midpoint */
	LeftTree *IntervalTree
	/* Tree of all intervals that are fully to the right of midpoint */
	RightTree *IntervalTree
}

/*
 * Sort (destructively) an array of pointers to intervals.
 */
func SortIntervalArray(input []*Interval) {
	var ia IntervalArray
	ia = input
	sort.Sort(ia)
}

/*
 * Build an interval tree from an (sorted) array of intervals.
 */
func DoBuildIntervalTree(input []*Interval) *IntervalTree {
	/* Allocate this node */
	it := new(IntervalTree)
	it.DataHere = make([]*Interval, 0, 0)
	/* Allocate sub-arrays to be passed to children */
	LeftArray := make([]*Interval, 0, 0)
	RightArray := make([]*Interval, 0, 0)

	mid_index := len(input) / 2

	/* Head of this tree is the midpoint of the data we want */
	it.MidPoint = (input[mid_index].Left + input[mid_index].Right) / 2

	/* Split the contents of interval into three arrays:
	 * -> An array of intervals that overlap it.MidPoint
	 * -> An array that is fully to the left
	 * -> An array that is fully to the right
	 *
	 * The copies that are preformed preserve the order of the input array
	 * in the output arrays such that data does not need to be resorted
	 * TODO: This is very expensive! Perhaps there is a more efficient way
	 * to do this?
	 */
	for i := 0; i < len(input); i++ {
		if input[i].Left == input[i].Right {
			log.Printf("Dropping zero-lengthed data")
		} else if input[i].Left <= it.MidPoint && input[i].Right > it.MidPoint {
			it.DataHere = append(it.DataHere, input[i])
		} else if input[i].Left > it.MidPoint {
			RightArray = append(RightArray, input[i])
		} else {
			LeftArray = append(LeftArray, input[i])
		}
	}

	/* Recursively call DoBuildIntervalTree on the left and right arrays.*/
	if len(LeftArray) > 0 {
		it.LeftTree = DoBuildIntervalTree(LeftArray)
	}
	if len(RightArray) > 0 {
		it.RightTree = DoBuildIntervalTree(RightArray)
	}
	return it
}

/*
 * Search an interval tree. Calls |callback| for each interval found that
 * overlaps with location.
 */
func SearchIntervalTree(it *IntervalTree, location int, callback func(i *Interval)) {

	if it == nil {
		return
	}

	/* Check all intervals that cross this node's midpoint. */
	for i := 0; i < len(it.DataHere); i++ {
		interval := it.DataHere[i]
		if interval.Left <= location && interval.Right > location {
			callback(interval)
		}
	}
	/* Recurse into either the left or right subtree. */
	if location > it.MidPoint && it.RightTree != nil {
		SearchIntervalTree(it.RightTree, location, callback)
	} else {
		SearchIntervalTree(it.LeftTree, location, callback)
	}
}

/*
 * Use SearchIntervalTree and an internal callback to compute an array of all
 * intervals that overlap |location| and return that array.
 */
func SearchIntervalTreeToArray(it *IntervalTree, location int) []*Interval {

	output := make([]*Interval, 0, 0)

	SearchIntervalTree(it, location, func(i *Interval) {
		output = append(output, i)
	})

	/* Sort the output by the midpoint of the intervals. This feature
	 * mostly exists to give a deterministic output for unit tests.
	 */
	SortIntervalArray(output)
	return output
}

/*
 * Build an interval tree by (destructively) sorting the input and then using
 * DoBuildIntervalTree to build the tree its self.
 */
func BuildIntervalTree(input []*Interval) *IntervalTree {
	SortIntervalArray(input)
	return DoBuildIntervalTree(input)
}

func WriteItreeToFile(path string, it map[string]*IntervalTree) error {

	fp, err := os.Create(path)

	if err != nil {
		log.Printf("Cannot create itree file: %v", err)
		return err
	}
	defer fp.Close()

	js, err := json.Marshal(it)

	if err != nil {
		panic(err)
	}

	fp.Write(js)

	return nil
}
