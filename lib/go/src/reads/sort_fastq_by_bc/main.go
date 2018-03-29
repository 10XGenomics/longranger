// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"container/heap"
	"encoding/json"
	"flag"
	"io/ioutil"
	"os"
	"reads/fastq_util"
	"sort"
	"strings"
)

var prefix = flag.String("prefix", "", "prefix of barcode that this bucket of barcodes [required].")
var barcode_bucket_json = flag.String("barcode_bucket_json", "", "filename of json describing the locations of all of the files for this barcode bucket [required].")
var output_dir = flag.String("output_dir", "", "output directory [required].")

func main() {
	flag.Parse()

	// first load each file from last stage that had this prefix and sort them and write them out again
	file_map := loadFileMap()
	sorted_filename_maps := []map[string]string{}
	for key, chunk_map := range file_map {
		bucket_chunk := loadBucketChunk(chunk_map)
		sort.Sort(BucketChunk(bucket_chunk))
		sorted_filename_map := writeBucketChunk(key, bucket_chunk, *output_dir)
		sorted_filename_maps = append(sorted_filename_maps, sorted_filename_map)
	}

	// then streamingly merge those files TODO
	output_filename := strings.Join([]string{*output_dir, "/", "reads.fastq.gz"}, "")
	read_pair_writer := fastq_util.NewReadPairWriter(output_filename)
	defer read_pair_writer.Close()

	readers := []*fastq_util.InterleavedReadPairReader{}
	for i := 0; i < len(sorted_filename_maps); i++ {
		file_map := sorted_filename_maps[i]
		read_pair_reader := fastq_util.NewInterleavedReadPairReader(file_map["R"])
		defer read_pair_reader.Close()
		readers = append(readers, read_pair_reader)
	}

	min_heap := MinHeap{}
	heap.Init(&min_heap)
	//initialize the heap with 1 read from each chunk
	for i := 0; i < len(readers); i++ {
		read_candidate, err := readers[i].GetReadPair()
		if err == nil && read_candidate != nil {
			heap.Push(&min_heap, &ReadPairHolder{read_pair: read_candidate, read_pair_reader: readers[i]})
		}
	}

	// while the min heap is not empty, grab the min, write it out and add the next read from dat reader
	for min_heap.Len() > 0 {
		min_read := heap.Pop(&min_heap).(*ReadPairHolder)
		read_pair_writer.WriteReadPair(min_read.read_pair)
		next_read_from_that_reader, err := min_read.read_pair_reader.GetReadPair()
		if err == nil {
			heap.Push(&min_heap, &ReadPairHolder{read_pair: next_read_from_that_reader, read_pair_reader: min_read.read_pair_reader})
		}
	}

	//remove temp files
	for i := 0; i < len(sorted_filename_maps); i++ {
		file_map := sorted_filename_maps[i]
		os.Remove(file_map["R"])
	}
}

// this type is for going into the heap so I keep track of which reader a read pair came from
// so that I can put the next read pair from that reader afer popping the min
type ReadPairHolder struct {
	read_pair        *fastq_util.ReadPair
	read_pair_reader *fastq_util.InterleavedReadPairReader
}

// heap backing
type MinHeap []*ReadPairHolder

// next several functions for heap interface
func (mh MinHeap) Len() int { return len(mh) }

func (mh MinHeap) Less(i, j int) bool {
	if string(*mh[i].read_pair.Barcode.Seq) == string(*mh[j].read_pair.Barcode.Seq) {

		return string(*mh[i].read_pair.Name) < string(*mh[j].read_pair.Name)
	} else {
		return string(*mh[i].read_pair.Barcode.Seq) < string(*mh[j].read_pair.Barcode.Seq)
	}
}

func (mh MinHeap) Swap(i, j int) {
	mh[i], mh[j] = mh[j], mh[i]
}

func (mh *MinHeap) Push(x interface{}) {
	item := x.(*ReadPairHolder)
	*mh = append(*mh, item)
}

func (mh *MinHeap) Pop() interface{} {
	old := *mh
	n := len(old)
	item := old[n-1]
	*mh = old[0 : n-1]
	return item
}

func loadBucketChunk(chunk_map map[string]string) []*fastq_util.ReadPair {
	read_filename := chunk_map["R"]
	reader := fastq_util.NewInterleavedReadPairReader(read_filename)
	defer reader.Close()

	read_pairs := []*fastq_util.ReadPair{}
	for {
		read_pair, err := reader.GetReadPair()
		if err != nil {
			break
		}
		read_pairs = append(read_pairs, read_pair)
	}
	return read_pairs
}

// write a bucket chunk after sorting
func writeBucketChunk(key string, bucket_chunk []*fastq_util.ReadPair, output_dir string) map[string]string {
	read_fname := strings.Join([]string{output_dir, "/", key, "R.fastq.gz"}, "")

	file_map := map[string]string{"R": read_fname}

	writer := fastq_util.NewReadPairWriter(read_fname)

	// and be a good boy and close my files
	defer writer.Close()

	for i := 0; i < len(bucket_chunk); i++ {
		r := bucket_chunk[i]
		writer.WriteReadPair(r)
	}
	return file_map
}

// filetype for sorting
type BucketChunk []*fastq_util.ReadPair

// and some sorting interface functions
func (a BucketChunk) Len() int      { return len(a) }
func (a BucketChunk) Swap(i, j int) { a[i], a[j] = a[j], a[i] }
func (a BucketChunk) Less(i, j int) bool {
	if string(*(a[i].Barcode.Seq)) == string(*(a[j].Barcode.Seq)) {
		return string(*(a[i].Name)) < string(*(a[j].Name))
	} else {
		return string(*(a[i].Barcode.Seq)) < string(*(a[j].Barcode.Seq))
	}
}

// get input from previous stage, maps from index -> map of ("R1","R2","R3","R4") -> filename
func loadFileMap() map[string]map[string]string {
	file, _ := ioutil.ReadFile(*barcode_bucket_json)
	var map_of_file_maps = map[string]map[string]string{}
	json.Unmarshal(file, &map_of_file_maps)
	return map_of_file_maps
}
