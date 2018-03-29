// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastq_util

import (
	"bufio"
	cgzip "code.google.com/vitess/go/cgzip"
	"fmt"
	"os"
)

// wrapper for reading from 4 files (or 3 if the sample index file is missing)
type ReadPairReader struct {
	read1_reader        *bufio.Reader
	read2_reader        *bufio.Reader
	barcode_reader      *bufio.Reader
	sample_index_reader *bufio.Reader // might be null
	read1_file          *os.File
	read2_file          *os.File
	barcode_file        *os.File
	sample_index_file   *os.File // might be null
}

type InterleavedReadPairReader struct {
	reader    *bufio.Reader
	read_file *os.File
}

// to make reading easy
func (r *ReadPairReader) GetReadPair() (*ReadPair, error) {
	return getReadPair(r.read1_reader, r.read2_reader, r.barcode_reader, r.sample_index_reader)
}

func (r *InterleavedReadPairReader) GetReadPair() (*ReadPair, error) {
	return getInterleavedReadPair(r.reader)
}

// easy to deal with closing all of those files
func (r *ReadPairReader) Close() {
	r.read1_file.Close()
	r.read2_file.Close()
	if r.barcode_file != nil {
		r.barcode_file.Close()
	}
	r.sample_index_file.Close()
}

func (r *InterleavedReadPairReader) Close() {
	r.read_file.Close()
}

// wrapper for writing files, note that the constructor for this
// only makes 1 file now and interleaves everything
type ReadPairWriter struct {
	read_writer *cgzip.Writer
	read_file   *os.File
}

// to make writing reads easy
func (r *ReadPairWriter) WriteReadPair(read *ReadPair) {
	writeReadPair(r.read_writer, read)
}

// remember guys, its just one file now, many reads, one file
func (r *ReadPairWriter) Close() {
	r.read_writer.Flush()
	r.read_writer.Close()
	r.read_file.Close()
}

// all the same file...
func NewReadPairWriter(filename string) *ReadPairWriter {
	r_file, _ := os.Create(filename)

	r_writer := cgzip.NewWriter(r_file)
	return &ReadPairWriter{
		read_writer: r_writer,
		read_file:   r_file,
	}
}

// write your 4 reads
func writeReadPair(r_writer *cgzip.Writer, read_pair *ReadPair) {
	writeReadName(r_writer, read_pair.Name)
	writeRead(r_writer, read_pair.Read1)
	writeRead(r_writer, read_pair.Read2)
	writeRead(r_writer, read_pair.Barcode)
	writeRead(r_writer, read_pair.Sample_index)
}

func writeReadName(writer *cgzip.Writer, name *[]byte) {
	_, err := writer.Write(*name)
	if err != nil {
		panic(fmt.Sprintf("cant write file"))
	}
}

// write a single read
func writeRead(writer *cgzip.Writer, read *Read) {
	writer.Write(*read.Seq)
	writer.Write(*read.Qual)
}

// ReadPairReader constructor
func NewReadPairReader(read_filename, barcode_filename, sample_index_filename string) *ReadPairReader {
	r_file, _ := os.Open(read_filename)
	r_reader := bufio.NewReader(r_file)
	var bc_file *os.File
	var bc_reader *bufio.Reader
	if barcode_filename != "none" {
		bc_file, _ = os.Open(barcode_filename)
		bc_reader = bufio.NewReader(bc_file)
	}
	var si_file *os.File
	var si_reader *bufio.Reader
	if sample_index_filename != "" {
		si_file, _ := os.Open(sample_index_filename)
		si_reader = bufio.NewReader(si_file)
	}
	toReturn := ReadPairReader{
		read1_reader:        r_reader,
		read2_reader:        r_reader,
		barcode_reader:      bc_reader,
		sample_index_reader: si_reader,
		read1_file:          r_file,
		read2_file:          r_file,
		barcode_file:        bc_file,
		sample_index_file:   si_file,
	}
	return &toReturn
}

func NewInterleavedReadPairReader(read_filename string) *InterleavedReadPairReader {
	r_file, _ := os.Open(read_filename)
	r_z_reader, _ := cgzip.NewReader(r_file)
	r_reader := bufio.NewReader(r_z_reader)
	toReturn := InterleavedReadPairReader{
		reader:    r_reader,
		read_file: r_file,
	}
	return &toReturn
}

// read all 4 reads (or 3 if there is no sample index)
func getReadPair(read1_reader, read2_reader, bc_reader, si_reader *bufio.Reader) (*ReadPair, error) {
	read_name, read_1, err := getRead(read1_reader)
	if err != nil {
		return nil, err
	}
	_, read_2, err := getRead(read2_reader)
	_, bc, err := getRead(bc_reader)

	var si_read *Read
	if si_reader != nil {
		_, si_read, err = getRead(si_reader)
	} else {
		seq := []byte("\n")
		qual := []byte("\n")
		si_read = &Read{Seq: &seq, Qual: &qual}
	}

	read_pair := ReadPair{Name: read_name, Read1: read_1, Read2: read_2, Barcode: bc, Sample_index: si_read}
	return &read_pair, nil
}

func getInterleavedReadPair(reader *bufio.Reader) (*ReadPair, error) {
	read_name, err := reader.ReadBytes('\n')
	if err != nil {
		return nil, err
	}
	read1, _ := getInterleavedRead(reader)
	read2, _ := getInterleavedRead(reader)
	bc, _ := getInterleavedRead(reader)
	si, _ := getInterleavedRead(reader)
	read_pair := ReadPair{Name: &read_name, Read1: read1, Read2: read2, Barcode: bc, Sample_index: si}
	return &read_pair, nil
}

// read a single read
func getRead(reader *bufio.Reader) (*[]byte, *Read, error) {
	if reader == nil {
		seq := []byte("\n")
		qual := []byte("\n")
		return nil, &Read{Seq: &seq, Qual: &qual}, nil
	} //for the case of no barcode read
	read_name, err := reader.ReadBytes('\n')
	if err != nil {
		seq := []byte("\n")
		qual := []byte("\n")
		return nil, &Read{Seq: &seq, Qual: &qual}, err
	}
	sequence, _ := reader.ReadBytes('\n')
	//throw away line
	reader.ReadSlice('\n')
	quality, _ := reader.ReadBytes('\n')
	return &read_name, &Read{Seq: &sequence, Qual: &quality}, nil
}

func getInterleavedRead(reader *bufio.Reader) (*Read, error) {
	read, err := reader.ReadBytes('\n')
	if err != nil {
		return nil, err
	}
	qual, err := reader.ReadBytes('\n')
	return &Read{Seq: &read, Qual: &qual}, nil
}

// keep your reads together and store only 1 name
type ReadPair struct {
	Name         *[]byte
	Read1        *Read
	Read2        *Read
	Barcode      *Read
	Sample_index *Read
}

// simple read data structure
type Read struct {
	Seq  *[]byte
	Qual *[]byte
}
