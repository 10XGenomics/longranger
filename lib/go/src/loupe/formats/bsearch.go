// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This implements code to binary search the contents of a loaded
 * VCF array for a SNP at a particular position
 */

package formats

/*
 * Compare sort-order of an VCF entry with what we're looking for.
 * We return less-than-zero, zero, or greater-than-zero depending
 * on how they compare.
 */
func VCFEntryCompare(entry *SimpleVCFRow, chr string, position int) int {
	/* First, compare the chromosomes */
	c := ChromosomeCMP(chr, entry.Chromosome)
	if c == 0 {
		/* If those match, compare the positions */
		return position - entry.Position
	}
	return c
}

/*
 * Bsearch a VCF array for chr+position. If chr+position does not occur in the
 * array, we return the item immediatly before it, or potentially zero.
 */
func BsearchVCFFile(vcf_array []*SimpleVCFRow, chr string, position int) int {

	where := 0

	/*
	 * Bsearch by setting bits, starting at the high bit. If setting the bit
	 * keeps is less-than our target position, keep it set. Otherwise,
	 * move to the next bit without setting it.
	 */

	/* XXX: This assumes 32 bit its and < 1<<31 entries in the VCF
	 * array!!! */

	for i := 31; i >= 0; i-- {
		test := where | 1<<uint32(i)
		if test >= len(vcf_array) {
			continue
		}

		vcf_entry := vcf_array[test]
		c := VCFEntryCompare(vcf_entry, chr, position)
		if c < 0 {
		} else if c > 0 {
			where = test
		} else {
			return test
		}
	}
	return where
}
