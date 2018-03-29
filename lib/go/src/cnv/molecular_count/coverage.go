// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"

	"cnv/commonFunc"
)

type FullCoverage struct {
	// CHANGE Made on June 8 2016
	// Adding the number of unbarcoded reads
	//
	// chrom -> probe/bin -> haploid: 0, 1, 2 -> read histogram
	AllCnt map[string][][]map[int]int
	// number of unbarcoded reads
	// chrom -> probe/bin
	NumNadeRead map[string][]int
	InnerCnt    map[string][][]map[int]int
	AllDM       map[string][]float64
	AllDMcnt    map[string][]int
	AllAS       map[string][]float64
	AllAScnt    map[string][]int
	bedset      *BedSet
}

type SimpleCoverage struct {
	// CHANGE Made on June 8 2016
	// Adding haploid "3" for unbarcoded reads
	//
	// chrom -> probe/bin -> read counts on each haploid: 0(total), 1, 2
	AllCnt map[string][][]int
	// number of unbarcoded reads
	// chrom -> probe/bin
	NumNadeRead map[string][]int
	InnerCnt    map[string][][]int
	AllDM       map[string][]float64
	AllDMcnt    map[string][]int
	AllAS       map[string][]float64
	AllAScnt    map[string][]int
	bedset      *BedSet
}

func (fc *FullCoverage) WriteOutPhaseFrac(file string) {
	f, err := os.Create(file)
	commonFunc.CheckErr(err)
	defer f.Close()
	fWriter := bufio.NewWriter(f)
	defer fWriter.Flush()

	// exon id -> [total reads, unphased reads]
	exonReadCnt := make(map[string][]int)
	exonList := make([]string, 0)
	for k, v := range fc.AllCnt {
		//fmt.Println(k, len(v), len(fc.bedset.Exon[k]))
		for i, _ := range v {
			exonID := fc.bedset.Exon[k][i]
			if exonID == "" {
				continue
			}
			phased := len(v[i][1]) + len(v[i][2])
			total := len(v[i][0])
			if _, ok := exonReadCnt[exonID]; ok {
				exonReadCnt[exonID][0] += total
				exonReadCnt[exonID][1] += phased
			} else {
				exonList = append(exonList, exonID)
				exonReadCnt[exonID] = []int{total, phased}
			}
		}
	}

	fmt.Fprint(fWriter, "exonName,fracPhased,phasedRead,totalRead\n")
	for _, id := range exonList {
		total, phased := exonReadCnt[id][0], exonReadCnt[id][1]
		if total == 0 {
			continue
		}
		frac := float64(phased) / float64(total)
		fmt.Fprintf(fWriter, "%s,%f,%d,%d\n", id, frac, phased, total)
	}

}

func NewFullCoverage(bs *BedSet) *FullCoverage {
	fc := &FullCoverage{}

	fc.AllCnt = make(map[string][][]map[int]int)
	fc.InnerCnt = make(map[string][][]map[int]int)
	fc.AllDM = make(map[string][]float64)
	fc.AllDMcnt = make(map[string][]int)
	fc.AllAS = make(map[string][]float64)
	fc.AllAScnt = make(map[string][]int)
	fc.NumNadeRead = make(map[string][]int)

	for k, v := range bs.Starts {
		l := len(v)
		fc.NumNadeRead[k] = make([]int, l)
		fc.AllCnt[k] = make([][]map[int]int, l)
		fc.InnerCnt[k] = make([][]map[int]int, l)
		fc.AllDM[k] = make([]float64, l)
		fc.AllDMcnt[k] = make([]int, l)
		fc.AllAS[k] = make([]float64, l)
		fc.AllAScnt[k] = make([]int, l)
		for i := 0; i < len(v); i++ {
			fc.AllCnt[k][i] = []map[int]int{
				make(map[int]int), make(map[int]int), make(map[int]int)}
			fc.InnerCnt[k][i] = []map[int]int{
				make(map[int]int), make(map[int]int), make(map[int]int)}
		}
	}

	fc.bedset = bs
	return fc
}

func (fc *FullCoverage) ToSimpleCount() *SimpleCoverage {
	sc := &SimpleCoverage{}
	sc.bedset = fc.bedset
	sc.NumNadeRead = fc.NumNadeRead

	sc.AllCnt = make(map[string][][]int)
	sc.InnerCnt = make(map[string][][]int)
	for k, v := range sc.bedset.Starts {
		sc.AllCnt[k] = make([][]int, len(v))
		sc.InnerCnt[k] = make([][]int, len(v))
		for i := 0; i < len(v); i++ {
			fcA := fc.AllCnt[k][i]
			fcI := fc.InnerCnt[k][i]
			sc.AllCnt[k][i] = []int{len(fcA[0]), len(fcA[1]), len(fcA[2])}
			sc.InnerCnt[k][i] = []int{len(fcI[0]), len(fcI[1]), len(fcI[2])}
		}
	}
	sc.AllDM = fc.AllDM
	sc.AllDMcnt = fc.AllDMcnt
	sc.AllAS = fc.AllAS
	sc.AllAScnt = fc.AllAScnt
	return sc
}

//func (fc *FullCoverage) CNVCalling(sf *SimpleCoverage) {
//	k := "chr1"
//	v := fc.bedset.Starts[k]
//	f, err := os.Create("cnv.txt")
//	commonFunc.CheckErr(err)
//	defer f.Close()
//
//	fWriter := bufio.NewWriter(f)
//	defer fWriter.Flush()
//
//}

func (fc *FullCoverage) ShowrPM(sf *SimpleCoverage) {
	//messages := []string{"total", "hap1", "hap2"}
	f, err := os.Create(OutRPM)
	commonFunc.CheckErr(err)
	defer f.Close()
	fWriter := bufio.NewWriter(f)
	defer fWriter.Flush()

	fmt.Fprintf(fWriter, "[\n")
	chrIdx := 0
	chrTtl := len(fc.bedset.Starts)
	for k, v := range fc.bedset.Starts {
		for i := 0; i < len(v); i++ {
			fcA := fc.AllCnt[k][i]
			rPMs := make([][]int, 3)
			fmt.Fprintf(fWriter, "    ")
			for j := 0; j < 3; j++ {
				rPMs[j] = make([]int, len(fcA[j]))
				if j == 0 {
					fmt.Fprintf(fWriter, "[")
				}
				fmt.Fprintf(fWriter, "%d,%d,", sf.AllCnt[k][i][j],
					sf.AllCnt[k][i][j]-len(fcA[j]))
				barCnt := 0
				sz := len(fcA[j])
				for _, vv := range fcA[j] {
					//fmt.Println(kk, vv, barCnt)
					rPMs[j][barCnt] = vv
					barCnt++
				}
				sort.Ints(rPMs[j])
				fmt.Fprint(fWriter, "[")
				cnt := 0
				for _, v := range rPMs[j] {
					fmt.Fprint(fWriter, v)
					cnt++
					if cnt < sz {
						fmt.Fprint(fWriter, ",")
					}
				}
				fmt.Fprint(fWriter, "],")
				//if j < 2 {
				//	fmt.Fprint(fWriter, ",")
				//}
			}
			//fmt.Fprint(fWriter, ",", k, ",", i)
			fmt.Fprint(fWriter, "\"", k, "\",", fc.bedset.Starts[k][i],
				",", fc.bedset.Ends[k][i], "]")
			//fmt.Fprint(fWriter, "]")
			if chrIdx < chrTtl-1 || i < len(v)-1 {
				fmt.Fprint(fWriter, ",")
			}
			fmt.Fprint(fWriter, "\n")
		}
		chrIdx++
	}
	fmt.Fprintf(fWriter, "]")
}

func (fc *FullCoverage) ShowrSimplePM(sf *SimpleCoverage) {
	//messages := []string{"total", "hap1", "hap2"}
	f, err := os.Create(OutRPM)
	commonFunc.CheckErr(err)
	defer f.Close()
	fWriter := bufio.NewWriter(f)
	defer fWriter.Flush()

	//fmt.Fprintf(fWriter, "[\n")
	//chrIdx := 0
	//chrTtl := len(fc.bedset.Starts)
	//fmt.Fprint(fWriter, "chrom,start,end,exonID,ID,totalMol,totalMolNoRead,HP1Mol,HP1MolNoRead,HP2Mol,HP2MolNoRead\n")
	for k, v := range fc.bedset.Starts {
		for i := 0; i < len(v); i++ {
			exonID := fc.bedset.Exon[k][i]
			id := fc.bedset.ID[k][i]
			start := fc.bedset.Starts[k][i]
			end := fc.bedset.Ends[k][i]
			fcA := fc.AllCnt[k][i]
			fmt.Fprintf(fWriter, "%s,%d,%d,%s,%s", k, start, end, exonID, id)
			for j := 0; j < 3; j++ {
				fmt.Fprintf(fWriter, ",%d,%d", sf.AllCnt[k][i][j], sf.AllCnt[k][i][j]-len(fcA[j]))
			}
			fmt.Fprintf(fWriter, ",%d", fc.NumNadeRead[k][i])
			// write out DM, AS info
			if fc.AllDMcnt[k][i] > 0 {
				fmt.Fprintf(fWriter, ",%.4f", fc.AllDM[k][i]/float64(fc.AllDMcnt[k][i]))
			} else {
				fmt.Fprintf(fWriter, ",%.4f", -999.0)
			}
			if fc.AllAScnt[k][i] > 0 {
				fmt.Fprintf(fWriter, ",%.4f", fc.AllAS[k][i]/float64(fc.AllAScnt[k][i]))
			} else {
				fmt.Fprintf(fWriter, ",%.4f", -999.0)
			}
			fmt.Fprint(fWriter, "\n")
		}
	}
}

func (fc *FullCoverage) ShowrPMWithMID(allMolInfo *FullCoverage) {
	//messages := []string{"total", "hap1", "hap2"}
	f, err := os.Create(OutRPM)
	commonFunc.CheckErr(err)
	defer f.Close()
	fWriter := bufio.NewWriter(f)
	defer fWriter.Flush()

	f2, err2 := os.Create(CovStat)
	commonFunc.CheckErr(err2)
	defer f2.Close()
	fWriter2 := bufio.NewWriter(f2)
	defer fWriter2.Flush()

	OverallTargetSZ := 0
	OverallMolSeen := 0.0
	OverallMol := 0

	bs := fc.bedset
	for chrom, startsByChrom := range bs.Starts {
		for binIdx := 0; binIdx < len(startsByChrom); binIdx++ {
			readInfo := fc.AllCnt[chrom][binIdx]
			molInfo := allMolInfo.AllCnt[chrom][binIdx]
			fmt.Fprintf(fWriter, "%s,%d,%d,%s,%s,", chrom,
				bs.Starts[chrom][binIdx], bs.Ends[chrom][binIdx],
				bs.Exon[chrom][binIdx], bs.ID[chrom][binIdx])

			for hp := 0; hp < 3; hp++ {
				ttlMol := len(molInfo[hp])
				seenMol := len(readInfo[hp])
				fmt.Fprintf(fWriter, "%d,%d,", ttlMol, ttlMol-seenMol)
			}
			fmt.Fprintf(fWriter, "%d,", fc.NumNadeRead[chrom][binIdx])
			for hp := 0; hp < 3; hp++ {
				midsMol := make([]int, 0)
				for mid := range molInfo[hp] {
					midsMol = append(midsMol, mid)
				}
				sort.Ints(midsMol)
				for _, mid := range midsMol {
					fmt.Fprintf(fWriter, "%d:", mid)
				}
				fmt.Fprint(fWriter, ",")
				midsRead := make([]int, 0)
				for mid := range readInfo[hp] {
					midsRead = append(midsRead, mid)
				}
				sort.Ints(midsRead)
				for _, mid := range midsRead {
					fmt.Fprintf(fWriter, "%d:", mid)
				}
				fmt.Fprint(fWriter, ",")
			}
			// write out DM, AS info
			if fc.AllDMcnt[chrom][binIdx] > 0 {
				fmt.Fprintf(fWriter, "%.4f,", fc.AllDM[chrom][binIdx]/
					float64(fc.AllDMcnt[chrom][binIdx]))
			} else {
				fmt.Fprintf(fWriter, "%.4f,", -999.0)
			}
			if fc.AllAScnt[chrom][binIdx] > 0 {
				fmt.Fprintf(fWriter, "%.4f,", fc.AllAS[chrom][binIdx]/
					float64(fc.AllAScnt[chrom][binIdx]))
			} else {
				fmt.Fprintf(fWriter, "%.4f", -999.0)
			}
			fmt.Fprint(fWriter, "\n")

			/// CovStat
			exonID, _ := strconv.Atoi(bs.Exon[chrom][binIdx])
			if exonID%10 == 5 {
				OverallTargetSZ += (bs.Ends[chrom][binIdx] - bs.Starts[chrom][binIdx])
				OverallMolSeen += (float64(len(readInfo[0])) + 0.5*float64(fc.NumNadeRead[chrom][binIdx]))
				OverallMol += len(molInfo[0])
			}
		}
	}

	fmt.Fprint(fWriter2, "OverallMolecularSeen\tOverallMolecular\tOverallTargetSize\n")
	fmt.Fprintf(fWriter2, "%.1f\t%d\t%d\n", OverallMolSeen, OverallMol, OverallTargetSZ)
}
