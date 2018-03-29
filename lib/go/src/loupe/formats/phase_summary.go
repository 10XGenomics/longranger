// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This module implements code to calculate a zoomed out view of phase block
 * data.  We store the data is a generic track with each entries attached
 * to an ExtraPhaseBlockData object containing phase-block specific information.
 */

package formats

import (
	"fmt"
)

/*
 * Phase block specific information beyond whats in GenericTrackData
 */
type ExtraPhaseBlockData struct {
	/* SNP counts */
	PhaseASNP   int
	PhaseBSNP   int
	UnPhasedSNP int
}

/*
 * This function decides of a row from a VCF file is compatible with a phase
 * block summary or not.
 * SNP's that share: chromosome, and phase id belong in the same summary
 */
func PBMatch(row *SimpleVCFRow, pbs *GenericTrackData, PhaseBlock int) bool {
	if pbs == nil || row == nil {
		return false
	}

	extra := pbs.Info.(*ExtraPhaseBlockData)

	/* Gene or Chromosome change always starts a now summary block */
	if pbs.Chromosome != row.Chromosome {
		return false
	}

	if row.PhaseId == 0 {
		return true
	}

	/* A phase block transition only causes a new summary block if there
	 * was at least one phased SNP on the old block.
	 */
	if PhaseBlock != row.PhaseId &&
		(extra.PhaseASNP != 0 || extra.PhaseBSNP != 0 || row.WellPhased) {
		return false
	}

	return true
}

/*
 Extend a PhaseBlockSummary struct to include another row of VCF data */
func ExtendPB(row *SimpleVCFRow, pbs *GenericTrackData) {
	extra := pbs.Info.(*ExtraPhaseBlockData)

	pbs.Stop = row.Position + 1
	if row.WellPhased {
		if row.Sequence[0].Sequence != row.Reference {
			extra.PhaseASNP++
		}
		if row.Sequence[1].Sequence != row.Reference {
			extra.PhaseBSNP++
		}
	} else {
		if row.Sequence[0].Sequence != row.Reference {
			extra.UnPhasedSNP++
		}
		if row.Sequence[1].Sequence != row.Reference {
			extra.UnPhasedSNP++
		}
	}
}

func HasMultiplePhasedSNPs(pbs *GenericTrackData) bool {
	extra := pbs.Info.(*ExtraPhaseBlockData)
	return (extra.PhaseASNP+extra.PhaseBSNP > 1)
}

func GeneratePhaseBlockSummary(gene_index *BundledTrackData, vcf_rows []*SimpleVCFRow) (*BundledTrackData, error) {

	/* Get an array of all rows in the VCF file */
	all_phase_block_data := make([]GenericTrackData, 0, 0)

	var current_pb *GenericTrackData
	current_phase_block_id := -1

	/* Greedily assemble phase block summaries from the rows of VCF data.
	 * Every time we transition to a new gene or a new phase block, start
	 * a new summary object.
	 */
	for _, row := range vcf_rows {

		/* Homozygous SNPs do not effect phasing at all */
		if row.Sequence[0].Sequence == row.Sequence[1].Sequence {
			continue
		}

		/* Unphased SNPs do not effect phasing at all */

		if !row.WellPhased {
			continue
		}

		/* Does this SNP belong in the same summary as the previous SNP?
		 */
		if PBMatch(row, current_pb, current_phase_block_id) {
			/* SNP is compatable. Extend the summary. */
			ExtendPB(row, current_pb)
		} else {

			/* Save it */
			if current_pb != nil && HasMultiplePhasedSNPs(current_pb) {
				/* XXX: This is a nasty kludge. This for loop
				 * starts with an uninitialized current_bp and
				 * we rely on falling-through to this code the
				 * first time the code runs to initialize a
				 * current_pb object.  Hence the nil check.
				 */
				all_phase_block_data = append(all_phase_block_data, *current_pb)
			}

			/* Create a new summary block and copy data into it */
			current_pb = new(GenericTrackData)
			current_pb.Chromosome = row.Chromosome
			current_pb.Start = row.Position
			current_pb.Stop = row.Position + 1
			current_pb.Name = fmt.Sprintf("%v", row.PhaseId)
			current_pb.Info = &ExtraPhaseBlockData{}
			/* Count SNP data from this row */
			ExtendPB(row, current_pb)
			current_phase_block_id = row.PhaseId
		}
	}
	if current_pb != nil {
		/* And don't forget to add the last one! */
		all_phase_block_data = append(all_phase_block_data, *current_pb)
	}

	return GenerateTrackIndex("phase_data", all_phase_block_data), nil

}
