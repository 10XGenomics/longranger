// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package formats

/*
 This converts a "contig index" from the "fai" index in the reference data package
 to a JSON encoding that loupe understands.
*/
func BuildContigDB(ref_contigs_path string, output *CompressedWriter) LoupeSection {

	fields := []string{"ContigName", "ContigLength", "Offset", "a", "b"}
	contig_list := []map[string]string{}

	ReadGenericTsvFile(ref_contigs_path, fields, func(i int, f map[string]string) error {
		contig_list = append(contig_list, f)
		return nil
	})

	LoupeSection, err := output.WriteJSON(contig_list)

	if err != nil {
		panic(err)
	}
	return LoupeSection

}
