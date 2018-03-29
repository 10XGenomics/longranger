// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

/*
 This file implements a compressed file writer that writes compressed chunks
 to the loupe file. It also makes space for a preamble at the begining and
 knows how to write to that*/

package formats

import (
	"bytes"
	"compress/gzip"
	"encoding/json"
	"io"
	"log"
	"os"
)

/*
 * This is our standard representation for a "section" of a loupe file
 */
type LoupeSection struct {
	Name  string
	Start int
	End   int
}

/*
 A Compressed writer object */
type CompressedWriter struct {
	CurrentOffset int
	HeaderSpace   int
	Writer        *os.File
}

/*
 Compress and write a chunk of data to a loupe file.
 This returns a LoupeSection object for the new data.
*/
func (c *CompressedWriter) WriteChunk(data []byte) (LoupeSection, error) {
	var b bytes.Buffer

	g := gzip.NewWriter(&b)

	_, err := g.Write(data)

	if err != nil {
		return LoupeSection{"", 0, 0}, err
	}

	g.Flush()
	g.Close()

	_, err = b.WriteTo(c.Writer)

	if err != nil {
		return LoupeSection{"", 0, 0}, err
	}

	new_offset, err := c.Writer.Seek(0, 1)
	if err != nil {
		return LoupeSection{"", 0, 0}, err
	}

	ret := LoupeSection{"", c.CurrentOffset, int(new_offset)}
	c.CurrentOffset = int(new_offset)
	return ret, nil

}

/*
 * This copies a file into the compressed Loupe file and returns
 * loupesection object for that file. No transformations (compression or
 * otherwise) are done on the input.
 */
func (c *CompressedWriter) WriteFromFile(s_path string) (LoupeSection, error) {

	source, err := os.Open(s_path)

	if err != nil {

		log.Printf("Can't open your file! %v", err)
		return LoupeSection{"", 0, 0}, err

	}

	defer source.Close()

	_, err = io.Copy(c.Writer, source)

	if err != nil {
		log.Printf("Problem copying data %v", err)
		return LoupeSection{"", 0, 0}, err
	}

	new_offset, err := c.Writer.Seek(0, 1)
	if err != nil {
		return LoupeSection{"", 0, 0}, err
	}

	ret := LoupeSection{"", c.CurrentOffset, int(new_offset)}
	c.CurrentOffset = int(new_offset)

	return ret, nil
}

/*
 * Serialize, compress, and write an object as JSON data to a loupe file.
 */
func (c *CompressedWriter) WriteJSON(jdata interface{}) (LoupeSection, error) {
	data, err := json.Marshal(jdata)
	if err != nil {
		return LoupeSection{"", 0, 0}, err
	}

	return (c.WriteChunk(data))
}

func (c *CompressedWriter) Close() {
	c.Writer.Close()
	c.Writer = nil
}

/* Write a header to a loupe file.  Many '\n' will be appended to the
   header to make things more intelligable.
*/
func (c *CompressedWriter) WriteHeader(data []byte) error {

	padding := []byte("\n\n\n\n\n\n\n\n\n\n\n")
	if len(data)+len(padding) > c.HeaderSpace {
		panic("OOPS")
	}

	f := c.Writer
	_, err := f.Seek(0, 0)

	if err != nil {
		return err
	}

	_, err = f.Write(data)
	if err != nil {
		return err
	}
	f.Write(padding)

	f.Seek(0, 2)

	return nil

}

/*
 Create a new loupe file. |headerspace| is the number of bytes to reserve at the
 begining of the file for  the preamble
*/
func NewCompressedWriter(path string, headerspace int) (*CompressedWriter, error) {

	f, err := os.Create(path)

	if err != nil {
		return nil, err
	}

	_, err = f.Seek(int64(headerspace), 0)

	if err != nil {
		return nil, err
	}

	cw := &CompressedWriter{headerspace, headerspace, f}
	return cw, nil
}
