// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This file defines the "master" preprocessor functions that
 * compute a indexed-JSON file for loupe from the raw BED / VCF/ FA
 * data.
 * TODO: This file needs a lot of help with error handling. Many errors
 * are dropped or just turn into panics.
 */
package preprocessor

import (
	"crypto/md5"
	"encoding/hex"
	"log"
	. "loupe/formats"
	"os"
)

/*
 * This file produces the master concatenated .loupe file for the given data.
 * the VCF and BED files are mandatory and the rest may be empty strings.
 * THhis function works by creating a separate file for each "section" of
 * output and then concatonating all of the files together.
 * The final file at output_path will be contain several sections, each
 * of which is valid JSON:
 * A index, describing where the other sections are.
 * A summary section,
 * A VCF index section
 * A phase block track
 * A refseq index
 * a gene track
 * a exn track
 * a bam index
 * A "marker" section for easily alignment of binary readers
 * A bunch of jsonified-VCF sections, index via the VCF index section
 * A compressed "bam" section containing informatione extracted from the BAM file
 */
func PreprocessRunData(temp_dir_name string,
	vcf_path string,
	bed_path string,
	refseq_path string,
	bam_path string,
	structvar_path string,
	shortstructvar_path string,
	breakpoints_path string,
	targets_path string,
	fragments_path string,
	contig_path string,
	output_path string,
	external_summary_paths map[string]string) error {

	log.Printf("Storing temporary files in %v", temp_dir_name)

	writer, err := NewCompressedWriter(output_path, 4096)

	if err != nil {
		return err
	}

	index := make(map[string]LoupeSection)

	reference_hash := HashFiles([]string{bed_path, contig_path})
	log.Printf("Reference set hash code: %v", reference_hash)

	/* This is a little bit disgusting and perhaps could be improved...
	 * We create a new json file that will have the phase block histogram
	 * in it and then add that file to external_summary_paths so that it
	 * gets merged in with the rest of the summary data.
	 */
	phase_block_histogram_path := temp_dir_name + "/pb-histogram.json"
	external_summary_paths[phase_block_histogram_path] = "phase_block_histogram"

	log.Printf("Processing contig DB")
	index["contigs_db"] = BuildContigDB(contig_path, writer)

	/* Step 1: Load the VCF file */
	log.Printf("Loading Customer VCF file")
	var vcf_data []*SimpleVCFRow
	if vcf_path != "" {
		vcf_data, err = ReadVCFToArray(vcf_path)
		if err != nil {
			log.Printf("Cannot read VCF file: %v", err)
			return err
		}
	} else {
		vcf_data = make([]*SimpleVCFRow, 0, 0)
	}

	/* Step 2: Compute a histogram of phase block lengths.  This will be
	 * assembled into the summary section, later
	 */
	log.Printf("Computeing Phase block histogram")
	WritePBHistogram(phase_block_histogram_path, vcf_data)

	/* Step 3: Compute a nice gene index */
	log.Printf("Loading genes file")
	
	index["ref_annotations"], err = writer.WriteFromFile(bed_path)

	/* Step 4. Associate DBSNP/refseq data */
	log.Printf("Associating SNPs with RefSeq data")
	refseq_annotation, err := ComputeVCFAnnotation(refseq_path, vcf_data)

	if err != nil {
		log.Printf("Cannot load refseq data: %v", err)
		log.Printf("Continuing without refseq data")
	}
	/* Step 5: Compute the jsonified VCF data and index. This will also decorate SNPs
	 * with their genes and genes with the number of SNPs they contain
	 */
	log.Printf("Computing indexed JSON representation of VCF data")
	index["vcf_index"] = VCFtoJSON(vcf_data, writer, 512, nil, refseq_annotation)

	/* STep 6, write out the gene data */

	index["summary"] = WriteSummaryFile(writer,
		"Some loupe data",
		vcf_path,
		bed_path,
		external_summary_paths)

	/* Step 7 Build the phase summary file */
	log.Printf("Creating Phase-block summary file ")
	pbs, err := GeneratePhaseBlockSummary(nil, vcf_data)
	if err != nil {
		return err
	}

	index["phase_block_summary2"] = WriteTrackDataToFile(writer, pbs)
	if err != nil {
		return err
	}

	/* Step 8 Write out the refseq index */
	log.Printf("Creating RSID index")
	index["refseq_index"] = WriteRSIndexToFile(writer, refseq_annotation)

	/* Step 9 Write out the BAM data and index. This appends the BAM
	 * data to the VCF data so that offsets that are encoded into the indeces
	 * will line up nicely.
	 */

	log.Printf("Storing and indexing read data")
	if bam_path != "" {
		index["bam_data"], _ = BuildBAMIndexAndData(bam_path, writer, nil)
	} else {
		/* No BAM data? Write an innocuous bam index. */
		index["bam_data"], err = writer.WriteChunk([]byte("{\"Database\":{}, \"Index\":{}}"))
		if err != nil {
			panic(err)
		}

	}

	log.Printf("Loading structural variant calls")
	/* Step 10. WRite structural variant data tin .loupe file */
	if structvar_path != "" {
		var sv_track *BundledTrackData
		var err error
		if breakpoints_path != "" {
			sv_track, err = LoadSVsWithDetails(structvar_path, breakpoints_path)
		} else {
			sv_track, err = LoadSVs(structvar_path)
		}
		if err != nil {
			return err
		}
		index["sv_track"] = WriteTrackDataToFile(writer, sv_track)

	} else {
		index["sv_track"], err = writer.WriteChunk([]byte("{}"))
		if err != nil {
			panic(err)
		}
	}

	/* Step 11. Write short structural variant data in .loupe file */
	// todo ShortSV 
	log.Printf("Loading short structural variant calls")
	if shortstructvar_path != "" {
		var shortSV_track *BundledTrackData
		var err error
		if breakpoints_path != "" {
			shortSV_track, err = LoadShortSVsWithDetails(shortstructvar_path, breakpoints_path)
		} else {
			shortSV_track, err = LoadShortSVs(shortstructvar_path)
		}
		if err != nil {
			return err
		}
		index["shortSV_track"] = WriteTrackDataToFile(writer, shortSV_track)

	} else {
		index["shortSV_track"], err = writer.WriteChunk([]byte("{}"))
		if err != nil {
			panic(err)
		}
	}

	/* Step 12 write target set information */
	log.Printf("Loading target set")
	if targets_path != "" {
		targets_track, err := ReadGenericBedFile(targets_path, []string{})
		if err != nil {
			return err
		}
		_ = targets_track
		index["targets_track"] = BuildBlockedJSONIndex(writer, targets_track)
	}

	/* step 13. Write the preamble */
	log.Printf("Creating summary file ")

	err = WriteIndexFile(writer, index, reference_hash)
	if err != nil {
		return err
	}

	return nil
}

/*
 * Generate the contents of the index section
 */
func GenerateSummaryIndex(data []string) map[string]int {
	m := make(map[string]int)

	for index, value := range data {
		m[value] = index
	}

	return m
}

/*
 * Compute a md5 hash of a bunch or files. We use this to compute a reference
 * hash so that Loupe will only have to load the reference once, even when
 * loading multiple files.
 */
func HashFiles(files []string) string {
	hasher := md5.New()

	for i := 0; i < len(files); i++ {
		fp, err := os.Open(files[i])
		if err != nil {
			continue
		}
		for {
			buffer := make([]byte, 8192, 8192)
			bytes_read, err := fp.Read(buffer)
			if bytes_read > 0 {
				hasher.Write(buffer)
			}

			if err != nil {
				break
			}
		}
	}

	s := make([]byte, 0)

	return hex.EncodeToString(hasher.Sum(s))

}
