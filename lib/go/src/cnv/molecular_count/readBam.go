// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"strconv"

	"cnv/commonFunc"
	"code.google.com/p/biogo.bam"
	//"github.com/biogo/hts/bam"
	//"github.com/biogo/hts/sam"
)

func ReadBam(file string, bs *BedSet, fragS *FragmentSet) *FullCoverage {
	bs.ResetIdx()

	fullCoverage := NewFullCoverage(bs)

	f, err0 := os.Open(file)
	commonFunc.CheckErr(err0)
	defer f.Close()

	bf, err := bam.NewReader(f)
	commonFunc.CheckErr(err)
	//defer bf.Close()

	refs := bf.Header().Refs()
	numChr := len(refs)
	chroms := make([]string, numChr)
	chrSizes := make([]int, numChr)

	for i, r := range refs {
		chroms[i] = r.Name()
		chrSizes[i] = r.Len()
	}

	cnt := 0
	r := &bam.Record{}
	for {
		_, err := bf.Read(r)
		if err != nil {
			break
		}

		var mid int
		mi, ok := r.Tag([]byte("MI"))
		if !ok || mi == nil {
			mid = -1
		} else {
			switch len(mi) {
			case 4:
				mid = int(uint8(mi[3]))
			case 5:
				mid = int(binary.LittleEndian.Uint16(mi[3:]))
			case 7:
				mid = int(binary.LittleEndian.Uint32(mi[3:]))
			default:
				fmt.Println("wrong size", len(mi))
				continue
			}
		}

		hp0, hasHP := r.Tag([]byte("HP"))
		hp := -1
		idxTrack := []int{0} //0 for total, 1 for hap1, 2 for hap2
		if hasHP {
			hp = int(hp0.Value().(uint8))
			idxTrack = append(idxTrack, hp)
		}

		if TESTING && mid > testNumFrag {
			fmt.Println("read ", cnt, "reads")
			return fullCoverage
		}

		rid := r.Ref.ID()
		mrid := r.MateRef.ID()
		if rid < 0 || rid >= numChr {
			continue
		}
		if mrid != rid { //|| r.MatePos > r.Pos {
			continue
		}

		start_diff := r.Pos - r.MatePos
		if start_diff < 0 {
			start_diff = -start_diff
		}

		//isProperPairBasic := (r.Flags&sam.Paired > 0) &&
		//	(r.Flags^sam.Unmapped > 0) &&
		//	(r.Flags^sam.MateUnmapped > 0) &&
		//	(r.Flags&sam.Reverse != r.Flags&sam.MateReverse)
		//isProperPair := isProperPairBasic && (start_diff <= 500)

		//fmt.Printf("%016b\t%02b\t%v\t%v\t%v\t%v\n", r.Flags, r.Flags&sam.ProperPair, isProperPair, start_diff, r.Ref.Name(), r.Pos)

		//isProperPair := r.Flags&sam.ProperPair > 0

		if r.MapQ < MAPQ {
			continue
		}

		//if !isProperPair || r.Flags >= sam.Secondary {
		if r.Flags >= bam.Secondary {
			continue
		}
		// a major correction
		chrom, insertS, insertE := r.Ref.Name(), r.Pos, r.End()

		if mid != -1 && fragS.Chroms[mid] != chrom {
			fmt.Println("Wrong mid", mid)
			continue
		}
		if _, ok := fullCoverage.AllCnt[chrom]; !ok {
			continue
		}
		cnt++
		if cnt%100000 == 0 {
			fmt.Println("Get read ", cnt)
		}

		if TESTING && cnt == testNumRead {
			break
		}

		// grap other information:
		// AS alingment score;
		// DM the mean NM among all reads in a molecular
		// AM, XM, XT : primary and seconary alignments are in molecular; XT whether it is tandem region

		// check Overlapping
		var allOLBaits, innerOLBaits []int
		if mid == -1 {
			allOLBaits, innerOLBaits = bs.ReadOverlapUniq(chrom, insertS, insertE,
				-1, -1)
		} else {
			allOLBaits, innerOLBaits = bs.ReadOverlapUniq(chrom, insertS, insertE,
				fragS.Starts[mid], fragS.Ends[mid])
		}

		//fmt.Printf("read # %6d\tchrom %s  interval [%8d, %8d]\t[%8d, %8d: %8d] fragment \tmid %d\n",
		//	cnt, chrom, insertS, insertE, fragS.Starts[mid], fragS.Ends[mid],
		//	fragS.Ends[mid]-fragS.Starts[mid]+1, mid)
		//fmt.Println(allOLBaits)
		//fmt.Println(innerOLBaits, "\n")

		if !Silence {
			//fmt.Println(r)
			fmt.Println(chrom, insertS, insertE, mid, allOLBaits, innerOLBaits, "\n")
		}
		for _, i := range allOLBaits {
			if mid == -1 {
				fullCoverage.NumNadeRead[chrom][i]++
			} else {
				for _, track := range idxTrack {
					if _, ok := fullCoverage.AllCnt[chrom][i][track][mid]; !ok {
						fullCoverage.AllCnt[chrom][i][track][mid] = 1
					} else {
						fullCoverage.AllCnt[chrom][i][track][mid]++
					}
				}
			}
		}

		for _, i := range innerOLBaits {
			for _, track := range idxTrack {
				if _, ok := fullCoverage.InnerCnt[chrom][i][track][mid]; !ok {
					fullCoverage.InnerCnt[chrom][i][track][mid] = 1
				} else {
					fullCoverage.InnerCnt[chrom][i][track][mid]++
				}
			}
		}

		as0, hasAS0 := r.Tag([]byte("AS"))
		dm0, hasDM0 := r.Tag([]byte("DM"))
		if hasAS0 {
			as := float64(as0.Value().(float32))
			for _, i := range allOLBaits {
				fullCoverage.AllAS[chrom][i] += as
				fullCoverage.AllAScnt[chrom][i]++
			}
		}
		if hasDM0 {
			dm, errDM := strconv.ParseFloat(dm0.Value().(string), 64)
			commonFunc.CheckErr(errDM)
			for _, i := range allOLBaits {
				fullCoverage.AllDM[chrom][i] += dm
				fullCoverage.AllDMcnt[chrom][i]++
			}
		}
	}
	fmt.Println("read ", cnt, "reads")
	return fullCoverage
}

func ReadBamViaFetch(bamFile string, bs *BedSet, fragS *FragmentSet, ncpu int) *FullCoverage {
	fmt.Println("OverlapThr", OverlapThr)
	// current implementation does not use parallelization
	idxFile := bamFile + ".bai"
	if _, err := os.Stat(idxFile); os.IsNotExist(err) {
		log.Fatal("bam index file " + idxFile + " does not exist!")
	}

	bs.ResetIdx()
	fullCoverage := NewFullCoverage(bs)

	f, err0 := os.Open(bamFile)
	commonFunc.CheckErr(err0)
	defer f.Close()

	fIdx, err1 := os.Open(idxFile)
	commonFunc.CheckErr(err1)
	defer fIdx.Close()

	bf, err := bam.NewReader(f)
	commonFunc.CheckErr(err)
	bif, err2 := bam.ReadIndex(fIdx)
	commonFunc.CheckErr(err2)
	//defer bf.Close()

	refs := bf.Header().Refs()
	numChr := len(refs)
	chroms := make([]string, numChr)
	chrSizes := make([]int, numChr)
	chr2ID := make(map[string]int)
	for i, r := range refs {
		chroms[i] = r.Name()
		chrSizes[i] = r.Len()
		chr2ID[r.Name()] = i
	}

	cnt := 0
	for ch, v := range bs.Starts {
		chID := chr2ID[ch]
		for idx, beg := range v {
			end := bs.Ends[ch][idx]
			//fmt.Println(ch, chID, beg, end)
			readIter, err := bf.Fetch(bif, chID, beg, end)
			if err != nil {
				continue
			}
			for readIter.Next() {
				r := readIter.Get()
				if r == nil {
					continue
				}

				var mid int
				// process each read
				mi, ok := r.Tag([]byte("MI"))
				if !ok || mi == nil {
					mid = -1
				} else {
					switch len(mi) {
					case 4:
						mid = int(uint8(mi[3]))
					case 5:
						mid = int(binary.LittleEndian.Uint16(mi[3:]))
					case 7:
						mid = int(binary.LittleEndian.Uint32(mi[3:]))
					default:
						fmt.Println("wrong size", len(mi))
						continue
					}
				}

				hp0, hasHP := r.Tag([]byte("HP"))
				hp := -1
				idxTrack := []int{0}
				if hasHP {
					hp = int(hp0.Value().(uint8))
					idxTrack = append(idxTrack, hp)
				}

				if TESTING && mid > testNumFrag {
					fmt.Println("read ", cnt, "reads")
					return fullCoverage
				}

				rid := r.Ref.ID()
				mrid := r.MateRef.ID()
				if rid < 0 || rid >= numChr {
					continue
				}
				if mrid != rid || rid != chID { //|| r.MatePos > r.Pos {
					continue
				}

				start_diff := r.Pos - r.MatePos
				if start_diff < 0 {
					start_diff = -start_diff
				}

				//isProperPairBasic := (r.Flags&sam.Paired > 0) &&
				//	(r.Flags^sam.Unmapped > 0) &&
				//	(r.Flags^sam.MateUnmapped > 0) &&
				//	(r.Flags&sam.Reverse != r.Flags&sam.MateReverse)
				//isProperPair := isProperPairBasic && (start_diff <= 500)

				//fmt.Printf("%016b\t%02b\t%v\t%v\t%v\t%v\n", r.Flags, r.Flags&sam.ProperPair, isProperPair, start_diff, r.Ref.Name(), r.Pos)

				//isProperPair := r.Flags&sam.ProperPair > 0

				if r.MapQ < MAPQ {
					continue
				}

				//if !isProperPair || r.Flags >= sam.Secondary {
				if r.Flags >= bam.Secondary {
					continue
				}
				// a major correction
				chrom, insertS, insertE := r.Ref.Name(), r.Pos, r.End()

				//fmt.Println(chrom, insertS, insertE, mid, hp)
				if mid != -1 && fragS.Chroms[mid] != chrom {
					fmt.Println("Wrong")
					continue
				}
				if _, ok := fullCoverage.AllCnt[chrom]; !ok {
					continue
				}
				cnt++
				if cnt%100000 == 0 {
					fmt.Println("Get read ", cnt)
				}

				if TESTING && cnt == testNumRead {
					break
				}

				// check overlapping
				var allOLBaits, innerOLBaits []int
				if mid == -1 {
					allOLBaits, innerOLBaits = bs.ReadOverlap(chrom, insertS, insertE,
						-1, -1)
				} else {
					allOLBaits, innerOLBaits = bs.ReadOverlap(chrom, insertS, insertE,
						fragS.Starts[mid], fragS.Ends[mid])
				}

				if !Silence {
					//fmt.Println(r)
					fmt.Println(chrom, insertS, insertE, mid, allOLBaits, innerOLBaits, "\n")
				}
				//fmt.Println(chrom, insertS, insertE, mid, hp, r.MapQ, allOLBaits, innerOLBaits, "\n")
				for _, i := range allOLBaits {
					if mid == -1 {
						fullCoverage.NumNadeRead[chrom][i]++
					} else {
						for _, track := range idxTrack {
							if _, ok := fullCoverage.AllCnt[chrom][i][track][mid]; !ok {
								fullCoverage.AllCnt[chrom][i][track][mid] = 1
							} else {
								fullCoverage.AllCnt[chrom][i][track][mid]++
							}
						}
					}
				}

				for _, i := range innerOLBaits {
					for _, track := range idxTrack {
						if _, ok := fullCoverage.InnerCnt[chrom][i][track][mid]; !ok {
							fullCoverage.InnerCnt[chrom][i][track][mid] = 1
						} else {
							fullCoverage.InnerCnt[chrom][i][track][mid]++
						}
					}
				}

				as0, hasAS0 := r.Tag([]byte("AS"))
				dm0, hasDM0 := r.Tag([]byte("DM"))
				if hasAS0 {
					as := float64(0.0)
					switch as0.Type() {
					case 'f':
						as = float64(as0.Value().(float32))
					case 'c':
						as = float64(as0.Value().(int8))
					case 'C':
						as = float64(as0.Value().(uint8))
					case 's':
						as = float64(as0.Value().(int16))
					case 'S':
						as = float64(as0.Value().(uint16))
					case 'i':
						as = float64(as0.Value().(int32))
					case 'I':
						as = float64(as0.Value().(uint32))
					}
					for _, i := range allOLBaits {
						fullCoverage.AllAS[chrom][i] += as
						fullCoverage.AllAScnt[chrom][i]++
					}
				}
				if hasDM0 {
					dm, errDM := strconv.ParseFloat(dm0.Value().(string), 64)
					commonFunc.CheckErr(errDM)
					for _, i := range allOLBaits {
						fullCoverage.AllDM[chrom][i] += dm
						fullCoverage.AllDMcnt[chrom][i]++
					}
				}

			} // end fetch for loop
		} // end indexing in bs
	} // end chrom for loop
	fmt.Println("read ", cnt, "reads")
	return fullCoverage
}
