// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"bufio"
	cgzip "code.google.com/vitess/go/cgzip"
	"crypto/md5"
	"encoding/binary"
	"flag"
	"os"
	"reads/fastq_util"
	"strconv"
	"strings"
	"tenkit/barcode"
)

/* Command line arguments */
var reads = flag.String("reads", "", "fastq.gz input file containing reads [required]")
var barcodes = flag.String("barcodes", "", "fastq.gz input containing barcodes [required]")
var sample_indices = flag.String("sample_index_reads", "", "fastq.gz file containing sample index reads [required]")
var gem_group = flag.String("gem_group", "1", "lane of GemCode instrument or null/1 if only 1 lane used [required]")
var reads_interleaved = flag.Bool("interleaved", true, "are reads interleaved? [required]")
var barcodeCounts = flag.String("barcodeCounts", "", "json out of COUNT_BCS stage[required]")
var bcConfidenceThreshold = flag.Float64("bcConfidenceThreshold", 0.975, "probability threshold of posterior of barcode")
var num_buckets = flag.Int("buckets", 256, "number of buckets to create")
var output_dir = flag.String("output_directory", "", "output director [required]")
var barcode_whitelist = flag.String("barcode_whitelist", "", "list of expected barcodes [required]")
var max_expected_barcode_errors = flag.Float64("max_expected_barcode_errors", 0.75, "maximum expected errors in the barcode (from qual string) to be considered a valid barcode")
var read_group_string = flag.String("read_group_string", "", "string to append to each FASTQ record to track read group")

func main() {
	flag.Parse()
	// Get readers
	if *sample_indices == "none" {
		nonestr := ""
		sample_indices = &nonestr
	}

	read_reader := fastq_util.NewReadPairReader(*reads, *barcodes, *sample_indices)
	defer read_reader.Close()
	barcoded_sample := true
	if *barcodes == "none" {
		barcoded_sample = false
	}

	// get writers (maps of prefix to writer)
	writer_files := make(map[int]map[string]string)
	read_writers := getReadWriters(reads, output_dir, writer_files, *num_buckets)

	barcode_validator := barcode.NewBarcodeValidator(*barcode_whitelist, *max_expected_barcode_errors, *barcodeCounts, *bcConfidenceThreshold, *gem_group)
	// read barcode, choose writers, write reads
	for {
		read, err1 := read_reader.GetReadPair()
		if err1 != nil {
			break
		}
		bc, isValid := barcode_validator.ValidateBarcode(string(*read.Barcode.Seq), *read.Barcode.Qual)

		// append read group string to read name
		new_name := (*read.Name)[0 : len(*read.Name)-1]
		new_name = append(new_name, byte(' '))
		new_name = append(new_name, *read_group_string...)
		new_name = append(new_name, byte('\n'))
		read.Name = &new_name

		writer := read_writers[getWriterIndex(&bc, read.Name, barcoded_sample)]
		if isValid {
			barcode_only := bc[0 : len(bc)-1]
			new_seq := append(barcode_only, byte('-'))
			gemGroup := []byte(*gem_group)
			for i := range gemGroup {
				new_seq = append(new_seq, gemGroup[i])
			}
			new_seq = append(new_seq, byte(','))
			for _, b := range *read.Barcode.Seq {
				if b != byte('\n') {
					new_seq = append(new_seq, b)
				}
			}
			new_seq = append(new_seq, byte('\n'))
			read.Barcode.Seq = &new_seq
		}
		writer.WriteReadPair(read)
	}

	// close all of my dirty, dirty files
	for _, writer := range read_writers {
		writer.Close()
	}
}

func getWriterIndex(bc, qname *[]byte, barcoded_sample bool) int {
	if barcoded_sample {
		md5sum := md5.Sum(*bc)
		toRet := int(binary.LittleEndian.Uint32(md5sum[0:8])) % (*num_buckets)
		if toRet < 0 {
			return -toRet
		}
		return toRet
	}
	md5sum := md5.Sum(*qname)
	toRet := int(binary.LittleEndian.Uint32(md5sum[0:8])) % int(*num_buckets)
	if toRet < 0 {
		return -toRet
	}
	return toRet
}

func writeRead(reader *bufio.Reader, writer *GzWriterCloser) error {
	for i := 0; i < 4; i++ {
		line_slice, err := reader.ReadSlice('\n')
		// error check somewhere maybe?
		if err != nil {
			return err
		}
		writer.writer.Write(line_slice)
	}
	return nil
}

func writeBarcodeRead(barcode_reader *bufio.Reader, barcode_writers map[string]*GzWriterCloser, nbases int,
	barcode_name_buffer *[]byte, barcode_seq_buffer *[]byte, barcode_validator *barcode.BarcodeValidator) (string, error) {
	// read barcode name and save in buffer until we know which bucket we are writing to
	line_slice, err := barcode_reader.ReadSlice('\n')
	if err != nil {
		return "", nil
	}
	copy(*barcode_name_buffer, line_slice)
	barcode_name := (*barcode_name_buffer)[:len(line_slice)]

	// read barcode sequence and write read name and save for later
	line_slice, err = barcode_reader.ReadSlice('\n')
	copy(*barcode_seq_buffer, line_slice)
	barcode_seq := (*barcode_seq_buffer)[:len(line_slice)]

	// read and toss +
	line_slice, err = barcode_reader.ReadSlice('\n')

	// read qual string and decide if barcode is valid
	line_slice, err = barcode_reader.ReadSlice('\n')
	barcode_seq, is_valid := barcode_validator.ValidateBarcode(string(barcode_seq), line_slice)
	var prefix string
	if is_valid {
		prefix = string(barcode_seq[:nbases])
	} else {
		count := nbases
		data := make([]byte, count)
		for i := 0; i < count; i++ {
			data[i] = byte('N')
		}
		prefix = string(data) //non valid barcodes go in the N chunk
	}
	barcode_writer := barcode_writers[prefix]
	barcode_writer.writer.Write(barcode_name)
	barcode_writer.writer.Write(barcode_seq)
	barcode_writer.writer.Write([]byte("+\n"))
	barcode_writer.writer.Write(line_slice)

	return prefix, err
}

type GzWriterCloser struct {
	writer *cgzip.Writer
	closer *os.File
}

func getReadWriters(name *string, outputdir *string, writer_names map[int]map[string]string, num_buckets int) map[int]*fastq_util.ReadPairWriter {
	writers := make(map[int]*fastq_util.ReadPairWriter)
	full_name := strings.Split(*name, "/")
	base_name := strings.Split(full_name[len(full_name)-1], ".fastq")[0]
	for i := 0; i < int(num_buckets); i++ {
		filename := strings.Join([]string{*outputdir, "/", base_name, "_", strconv.FormatInt(int64(i), 10), "_", *gem_group, ".fastq.gz"}, "")
		writers[i] = fastq_util.NewReadPairWriter(filename)
		if bucket, has := writer_names[i]; has {
			bucket["R"] = filename
			bucket["gem_group"] = *gem_group
		} else {
			writer_names[i] = make(map[string]string)
			writer_names[i]["R"] = filename
			writer_names[i]["gem_group"] = *gem_group
		}
	}
	return writers
}

// func getSeqs(nbases int) []string {
// 	nucs := []string{"A", "C", "G", "T", "N"}
// 	if nbases == 1 {
// 		return nucs
// 	}
// 	old_seqs := getSeqs(nbases - 1)
// 	new_seqs := []string{}
// 	for i := 0; i < len(old_seqs); i++ {
// 		for j := 0; j < len(nucs); j++ {
// 			new_seqs = append(new_seqs, strings.Join([]string{old_seqs[i], nucs[j]}, ""))
// 		}
// 	}
// 	return new_seqs
// }
