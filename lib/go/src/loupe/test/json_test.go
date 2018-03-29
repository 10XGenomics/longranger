// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

package test

import (
	. "loupe/formats"
	"testing"

	"encoding/json"
	//	"log"
)

func TestCombineJSON(t *testing.T) {
	j, _ := CombineJsonFiles(ParseSummarySpecification("inputs/combine_json_1.json:a,inputs/combine_json_2.json:b,inputs/combine_json_3.json:c,inputs/broken_json.json:d"))

	json.Marshal(j)
	//log.Printf("*********************")
	//log.Printf("J: %v", (string)(d))

	/* TODO: test this output */
}
