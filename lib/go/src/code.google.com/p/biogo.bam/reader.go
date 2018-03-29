// Copyright ©2012 The bíogo.bam Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package bam

import (
	"code.google.com/p/biogo.bam/bgzf"
	"encoding/binary"
	"errors"
	"io"
	"unsafe"
)

type Reader struct {
	r io.Reader
	h *Header
}

func NewReader(r io.Reader) (*Reader, error) {
	bg, err := bgzf.NewReader(r)
	if err != nil {
		return nil, err
	}
	br := &Reader{
		r: bg,
		h: &Header{
			seenRefs:   set{},
			seenGroups: set{},
			seenProgs:  set{},
		},
	}
	err = br.h.read(br.r)
	if err != nil {
		return nil, err
	}
	return br, nil
}

func NewUncompressedReader(r io.Reader) (*Reader, error) {
	br := &Reader{
		r: r,
		h: &Header{
			seenRefs:   set{},
			seenGroups: set{},
			seenProgs:  set{},
		},
	}
	err := br.h.read(br.r)
	if err != nil {
		return nil, err
	}
	return br, nil
}

func (br *Reader) Header() *Header {
	return br.h
}

// BAM record layout.
type bamRecordFixed struct {
	blockSize int32
	refID     int32
	pos       int32
	nLen      uint8
	mapQ      uint8
	bin       uint16
	nCigar    uint16
	flags     Flags
	lSeq      int32
	nextRefID int32
	nextPos   int32
	tLen      int32
}

var (
	lenFieldSize      = binary.Size(bamRecordFixed{}.blockSize)
	bamFixedRemainder = binary.Size(bamRecordFixed{}) - lenFieldSize
)

/*
 * Read data into a record.  Return a pointer to the same record
 */
func (br *Reader) Read(rec *Record) (*Record, error) {
	r := errReader{r: br.r}
	bin := binaryReader{r: &r}

	// Read record header data.
	blockSize := int(bin.readInt32())
	r.n = 0 // The blocksize field is not included in the blocksize.

	rec.buf_reset()

	refID := bin.readInt32()
	rec.Pos = int(bin.readInt32())
	nLen := bin.readUint8()
	rec.MapQ = bin.readUint8()
	_ = bin.readUint16()
	nCigar := bin.readUint16()
	rec.Flags = Flags(bin.readUint16())
	lSeq := bin.readInt32()
	nextRefID := bin.readInt32()
	rec.MatePos = int(bin.readInt32())
	rec.TempLen = int(bin.readInt32())
	if r.err != nil {
		return nil, r.err
	}

	// Read variable length data.
	name := rec.buf_alloc(int(nLen))
	if nf, _ := r.Read(name); nf != int(nLen) {
		return nil, errors.New("bam: truncated record name")
	}
	rec.Name = string(name[:len(name)-1]) // The BAM spec indicates name is null terminated.

	rec.Cigar = readCigarOps(&bin, nCigar)
	if r.err != nil {
		return nil, r.err
	}

	seq := rec.buf_alloc(int((lSeq + 1) >> 1))
	if nf, _ := r.Read(seq); nf != int((lSeq+1)>>1) {
		return nil, errors.New("bam: truncated sequence")
	}
	rec.Seq = NybbleSeq{Length: int(lSeq), Seq: *(*nybblePairs)(unsafe.Pointer(&seq))}

	rec.Qual = rec.buf_alloc(int(lSeq))

	if nf, _ := r.Read(rec.Qual); nf != int(lSeq) {
		return nil, errors.New("bam: truncated quality")
	}

	auxTags := rec.buf_alloc(blockSize - r.n)
	r.Read(auxTags)
	if r.n != blockSize {
		return nil, errors.New("bam: truncated auxilliary data")
	}
	rec.AuxTags = parseAux(auxTags, rec.aux_buffer[0:0])

	if r.err != nil {
		return nil, r.err
	}

	refs := int32(len(br.h.Refs()))
	if refID != -1 {
		if refID < -1 || refID >= refs {
			return nil, errors.New("bam: reference id out of range")
		}
		rec.Ref = br.h.Refs()[refID]
	} else {
		rec.Ref = nil
	}
	if nextRefID != -1 {
		if nextRefID < -1 || nextRefID >= refs {
			return nil, errors.New("bam: mate reference id out of range")
		}
		rec.MateRef = br.h.Refs()[nextRefID]
	} else {
		rec.MateRef = nil
	}

	return rec, nil
}

type FetchIter struct {
	leftMostOffset bgzf.Offset
	br             *Reader
	rid, beg, end  int
	hasStarted     bool
	rec            *Record
}

func (br *Reader) Fetch(idx *Index, rid, beg, end int) (FetchIter, error) {
	// Index is specified as an input, better to be included in the Reader class
	overlappingChunks := idx.Chunks(rid, beg, end)
	//for i, v := range overlappingChunks {
	//	fmt.Println(i, v.Begin, v.End)
	//}
	if len(overlappingChunks) <= 0 {
		return FetchIter{}, nil
	}
	leftMostOffset := overlappingChunks[0].Begin
	//fmt.Println(leftMostOffset)
	bgzfR, ok := br.r.(*bgzf.Reader)
	//bgzfR, err := br.r.(io.ReadSeeker)
	//bgzfR, err := (&br.r).(*bgzf.Reader)
	//bgzfR, err := br.r.(bgzf.Reader)
	if !ok {
		return FetchIter{}, errors.New("the input bam file is not a proper bgzf file")
	}
	bgzfR.Seek(leftMostOffset, 0)

	iter := FetchIter{leftMostOffset: leftMostOffset,
		br:  br,
		rid: rid,
		beg: beg,
		end: end,
		rec: &Record{},
	}
	return iter, nil
}

func (iter FetchIter) Next() bool {
    // Bail out in case of empty iterator
    if iter.br == nil {
        return false
    }

	_, err := iter.br.Read(iter.rec)
	if err != nil || iter.rec.Ref == nil {
		//fmt.Println("find a nil")
		return false
	}

	refID := iter.rec.ReferenceID()
	for refID < iter.rid || (refID == iter.rid && iter.rec.End() < iter.beg) {
		_, err = iter.br.Read(iter.rec)
		if err != nil {
			//fmt.Println("find a nil in getting to the right chromosome")
			return false
		}
		refID = iter.rec.ReferenceID()
	}
	if refID > iter.rid || refID == -1 || (refID == iter.rid && iter.rec.Pos > iter.end) {
		//fmt.Println("reach the end of query region", iter.rid, refID, iter.end, iter.rec.Pos)
		return false
	}
	if refID == iter.rid && iter.rec.End() >= iter.beg {
		return true
	} else {
		//fmt.Println("Explain why it is not true", iter.rid, refID, iter.end, iter.rec.End())
		return false
	}
}

func (iter FetchIter) Get() *Record { return iter.rec }

func readCigarOps(br *binaryReader, n uint16) []CigarOp {
	co := make([]CigarOp, n)
	for i := range co {
		co[i] = CigarOp(br.readUint32())
		if br.r.err != nil {
			return nil
		}
	}
	return co
}

type errReader struct {
	r   io.Reader
	n   int
	err error
}

func (r *errReader) Read(p []byte) (int, error) {
	if r.err != nil {
		return 0, r.err
	}
	var n int
	n, r.err = r.r.Read(p)
	r.n += n
	return n, r.err
}

type binaryReader struct {
	r   *errReader
	buf [4]byte
}

func (r *binaryReader) readUint8() uint8 {
	r.r.Read(r.buf[:1])
	return r.buf[0]
}

func (r *binaryReader) readUint16() uint16 {
	r.r.Read(r.buf[:2])
	return binary.LittleEndian.Uint16(r.buf[:2])
}

func (r *binaryReader) readInt32() int32 {
	r.r.Read(r.buf[:4])
	return int32(binary.LittleEndian.Uint32(r.buf[:4]))
}

func (r *binaryReader) readUint32() uint32 {
	r.r.Read(r.buf[:4])
	return binary.LittleEndian.Uint32(r.buf[:4])
}
