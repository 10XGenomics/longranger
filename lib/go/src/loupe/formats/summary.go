// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This file implements the file format for the "summary" part of loupe's data
 * file.  The summary will contain overall information about a particular loupe
 * data set. This would include:
 *  - date and name of the sequence
 *  - links to "full" datasets (VCF/bam/fa) files.
 *  - version information
 *  - sequence quality data
 *
 * Right now, this is essentially a placeholder.
 */

package formats

import (
	"encoding/json"
	"io/ioutil"
	"log"
	"strings"
)

func CombineJsonFiles(paths map[string]string) (interface{}, error) {
	all_data := make(map[string]interface{})

	for path, key := range paths {
		fdata, err := ioutil.ReadFile(path)
		if err != nil {
			log.Printf("Cannot read data from: %v: %v", path, err)
			continue
		}

		var json_data interface{}
		err = json.Unmarshal(fdata, &json_data)

		if err != nil {
			log.Printf("Cannot parse json from %v: %v", path, err)
			continue
		}
		all_data[key] = json_data
	}

	return all_data, nil
}

func WriteIndexFile(writer *CompressedWriter, index map[string]LoupeSection, hash string) error {
	index_data := make(map[string]interface{})

	index_data["index"] = index
	index_data["version"] = "2.4"
	index_data["reference"] = hash

	js, _ := json.Marshal(index_data)

	err := writer.WriteHeader(js)
	return err
}

/*
 * Write summary data to the given path.
 */
func WriteSummaryFile(writer *CompressedWriter,
	run_name string,
	vcf_path string,
	bed_path string,
	external_summaries map[string]string) LoupeSection {

	summary_data_iface, err := CombineJsonFiles(external_summaries)
	if err != nil {
		panic(err)
	}
	summary_data := summary_data_iface.(map[string]interface{})
	/* Current summary information is only a name and version */
	summary_data["name"] = run_name
	summary_data["vcf_source_path"] = vcf_path
	summary_data["bed_source_path"] = bed_path

	section, err := writer.WriteJSON(summary_data)

	if err != nil {
		panic(err)
	}
	return section
}

func ParseSummarySpecification(spec string) map[string]string {
	m := make(map[string]string)
	parts := strings.Split(spec, ",")
	for _, p := range parts {
		key_val := strings.Split(p, ":")

		key := key_val[0]
		val := key_val[1]

		m[key] = val
	}

	return m
}
