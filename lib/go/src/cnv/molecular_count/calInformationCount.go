// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	//"os/exec"

	"cnv/commonFunc"
	"github.com/codegangsta/cli"
	//"github.com/gonum/stat"
)

var (
	Silence            bool
	BedFile            string
	FragmentFile       string
	PhasedFragmentFile string
	BamFile            string
	BamIdxFile         string
	OutFull            string
	Outstat            string
	Outprobe_ttl_cnv   string
	Outprobe_hp1_cnv   string
	Outprobe_hp2_cnv   string
	OutRPM             string
	FracPhased         string
	CovStat            string
	MAXCV              float64
	MINREADCNT         int
	CPU                int
	WGSMODE            bool
	BEDTYPE            int
	mapq               int
	MAPQ               byte
	OverlapThr         int

	//pathGetTargetCNV string = "/mnt/home/wei/longranger/CNV/get_target_cnv.sh"
)

func main() {
	baitSet := ReadBaits(BedFile)
	allMolInfo, fragS := baitSet.GetAllMolInfo(FragmentFile, PhasedFragmentFile)
	fullCoverage := ReadBamViaFetch(BamFile, baitSet, fragS, CPU)

	fullCoverage.ShowrPMWithMID(allMolInfo)
	os.Exit(0)
}

func Mean(nums map[int]int) float64 {
	if len(nums) == 0 {
		return -1.0
	}
	sum := 0
	for _, v := range nums {
		sum += v
	}
	return float64(sum) / float64(len(nums))
}

func CVEstimate(n, N int) float64 {
	if n == 0 || N == 0 {
		return -1.0
	}

	if N < n {
		fmt.Println("Wrong number of innerFrag smaller than innerRead")
		return -1.0
	}

	return math.Sqrt(1/float64(n) - 1/float64(N))

}

func writeOut(fragCov, readCov *SimpleCoverage) {
	readThr := [3]int{MINREADCNT, 2, 2}
	cvThr := [3]float64{MAXCV, 0.5, 0.5}

	outFileNames := map[string]string{
		"full": OutFull,
		"stat": Outstat,
		"ttl":  Outprobe_ttl_cnv,
		"hp1":  Outprobe_hp1_cnv,
		"hp2":  Outprobe_hp2_cnv,
	}
	outWriters := make(map[string]*bufio.Writer)
	for k, v := range outFileNames {
		f, err := os.Create(v)
		commonFunc.CheckErr(err)
		defer f.Close()
		w := bufio.NewWriter(f)
		defer w.Flush()
		outWriters[k] = w
	}

	bedWriters := map[int]*bufio.Writer{
		0: outWriters["ttl"],
		1: outWriters["hp1"],
		2: outWriters["hp2"],
	}

	// CNV, Read, Frag
	valid := [][]float64{make([]float64, 0), make([]float64, 0), make([]float64, 0)}

	for chrom, v := range fragCov.AllCnt {
		for i, _ := range v {
			allRead := readCov.AllCnt[chrom][i]
			innerRead := readCov.InnerCnt[chrom][i]
			allFrag := fragCov.AllCnt[chrom][i]
			innerFrag := fragCov.InnerCnt[chrom][i]

			cv := [3]float64{-1.0, -1.0, -1.0}
			est := [3]float64{-1.0, -1.0, -1.0}
			s := fragCov.bedset.Starts[chrom][i]
			e := fragCov.bedset.Ends[chrom][i]
			fmt.Fprintf(outWriters["full"], "%4s,%8d,%8d\t\t", chrom, s, e)
			for j := 0; j < 3; j++ {
				cv[j] = CVEstimate(innerRead[j], innerFrag[j])
				if innerRead[j] >= readThr[j] && cv[j] > 0.0 && cv[j] < cvThr[j] {
					est[j] = float64(allRead[j]) * float64(innerFrag[j]) / float64(innerRead[j])
					fmt.Fprintf(bedWriters[j], "%s\t%d\t%d\t%s\n", chrom, s, e, toString(est[j]))
				}
				fmt.Fprintf(outWriters["full"], "%4d,%4d,%4d,%4d:%8.2f\t\t",
					allRead[j], innerRead[j], allFrag[j], innerFrag[j], est[j])
			}
			fmt.Fprintf(outWriters["full"], "\n")
			if est[0] >= 0.0 {
				valid[0] = append(valid[0], est[0])
				valid[1] = append(valid[1], float64(allRead[0]))
				valid[2] = append(valid[2], float64(allFrag[0]))
			}
		}
	}

	fmt.Fprintf(outWriters["stat"], "method\tmean\tstd\tcv\t#probe\n")
	sz := len(valid[0])
	methods := [3]string{"NEW METHOD", "UNIQUE INSERT", "MOLECULAR"}
	for i := 0; i < 3; i++ {
		//mn, std := stat.MeanStdDev(valid[i], nil)
		mn, std := 1.0, 1.0
		fmt.Fprintf(outWriters["stat"], "%16s\t%.2f\t%.2f\t%.6f\t%d\n", methods[i], mn, std, std/mn, sz)
	}
}

func toString(est float64) string {
	if est >= 0.0 {
		return fmt.Sprintf("%.2f", est)
	} else {
		return "."
	}
}

func init() {
	app := cli.NewApp()
	app.Name = "10X CNV caller"
	app.Usage = "using molecular information to correct CNV callings"

	app.Flags = []cli.Flag{
		cli.BoolTFlag{
			Name:        "silenct, s",
			Usage:       "do not print the detailed information",
			Destination: &Silence,
		},
		cli.StringFlag{
			Name:        "bedfile, b",
			Value:       "/mnt/opt/meowmix_git/genometracks/hg19/cnv_stratification/pipeline_related/exome/v6_core_final_merged.bed",
			Usage:       "the input bed file",
			Destination: &BedFile,
		},
		cli.StringFlag{
			Name:        "fragment, f",
			Value:       "/mnt/home/wei/algoDev/longranger/CNV/17341_fragments.csv",
			Usage:       "the fragments csv file with positions and molecular id",
			Destination: &FragmentFile,
		},
		cli.StringFlag{
			Name:        "phasedfragment, p",
			Value:       "/mnt/home/wei/algoDev/longranger/CNV/17341_fragment_phasing.tsv",
			Usage:       "the phased fragments tsv file with positions and molecular id",
			Destination: &PhasedFragmentFile,
		},
		cli.StringFlag{
			Name:        "bam, m",
			Value:       "/mnt/home/wei/algoDev/longranger/CNV/17341_phased_possorted_bam.bam",
			Usage:       "the input bam file",
			Destination: &BamFile,
		},
		cli.StringFlag{
			Name:        "full",
			Value:       "cnv.tsv",
			Usage:       "output full resutls",
			Destination: &OutFull,
		},
		cli.StringFlag{
			Name:        "fracPhased",
			Value:       "fracPhased.csv",
			Usage:       "fraction of reads phased per each exon",
			Destination: &FracPhased,
		},
		cli.StringFlag{
			Name:        "stat",
			Value:       "stats.tsv",
			Usage:       "output stat resutls",
			Destination: &Outstat,
		},
		cli.StringFlag{
			Name:        "probe_ttl_cnv",
			Value:       "probe_ttl_cnvs.tsv",
			Usage:       "output probe_ttl_cnv resutls",
			Destination: &Outprobe_ttl_cnv,
		},
		cli.StringFlag{
			Name:        "probe_hp1_cnv",
			Value:       "probe_hp1_cnvs.tsv",
			Usage:       "output probe_hp1_cnv resutls",
			Destination: &Outprobe_hp1_cnv,
		},
		cli.StringFlag{
			Name:        "probe_hp2_cnv",
			Value:       "probe_hp2_cnvs.tsv",
			Usage:       "output probe_hp2_cnv resutls",
			Destination: &Outprobe_hp2_cnv,
		},
		cli.StringFlag{
			Name:        "rPM",
			Value:       "rPM.json",
			Usage:       "output rPM",
			Destination: &OutRPM,
		},
		cli.StringFlag{
			Name:        "covstat",
			Value:       "covstat.txt",
			Usage:       "output total size of tageted regions and total number of molecular",
			Destination: &CovStat,
		},
		cli.Float64Flag{
			Name:        "maxcv, c",
			Value:       0.333,
			Usage:       "maximum cv, below which the total DNA will be estimated",
			Destination: &MAXCV,
		},
		cli.IntFlag{
			Name:        "minreadcnt, n",
			Value:       4,
			Usage:       "min readcnt, above or equal to which the total DNA will be estimated",
			Destination: &MINREADCNT,
		},
		cli.IntFlag{
			Name:        "overlap",
			Value:       20,
			Usage:       "threshold beyond which a bait is regarded as overlapping",
			Destination: &OverlapThr,
		},
		cli.IntFlag{
			Name:        "fragendbuffer",
			Value:       500,
			Usage:       "threshold beyond which a bait is regarded as in the side of a fragment",
			Destination: &FragEndBuffer,
		},
		cli.BoolFlag{
			Name:        "testing, t",
			Usage:       "when true perform a quick run",
			Destination: &TESTING,
		},
		cli.IntFlag{
			Name:        "cpu",
			Value:       1,
			Usage:       "number of cpus used for the potential parallelization mode",
			Destination: &CPU,
		},
		cli.IntFlag{
			Name:        "mapq",
			Value:       30,
			Usage:       "lowest mapq value for a read to be considered",
			Destination: &mapq,
		},
		cli.BoolFlag{
			Name:        "wgs",
			Usage:       "whether the application is for WGS. Default WES or targeted sequencing",
			Destination: &WGSMODE,
		},
	}
	app.Action = func(c *cli.Context) {}

	app.Run(os.Args)
	BEDTYPE = -1
	MAPQ = (byte)(mapq)
	fmt.Println(BedFile)
	fmt.Println("WGSMODE", WGSMODE)
	fmt.Println("testing mode", TESTING)
	fmt.Println("number of fileds in the bed file ", BEDTYPE)
	fmt.Println("Silent", Silence)
}
