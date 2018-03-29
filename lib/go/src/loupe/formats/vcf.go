// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

/*
 * This is a very simple parser for VCF files. It unpacks the "metadata" for
 * each variation as well as the first variant specimen and places them into
 * the structs below.  The ReadVCFWithCallback function provides efficient
 * iteration over all such structures in a single VCF file.
 */

package formats

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

const MAX_BARCODE_BYTES = 5
const MAX_GEM_GRP_BYTES = 2

/*
 * This represents a single 10X barcode with associated quality score.
 * Compressed to 8 bytes for memory efficiency.
 */
type TenXBarcode struct {
	Barcode  [MAX_BARCODE_BYTES]byte
	GemGroup [MAX_GEM_GRP_BYTES]byte
	Length   byte
}

/* Represent a very simple VCF file */
type SimpleVCFRow struct {
	/* The chromosome we're talking about */
	Chromosome string
	/* The offset into this chromosome */
	Position int
	/* The reference id specified in the file */
	Id string
	/* The reference nucleotide sequence specified in the file */
	Reference string
	/* The phased read of the "experimental" data */
	Sequence [2]PhasedAlternative
	/* Was the data above correctly phased or not */
	WellPhased bool
	/* The phase block the sequence is in*/
	PhaseId int

	CallQuality float64

	PhaseQuality float64
}

type PhasedAlternative struct {
	/* The experimental nucleiotide sequence */
	Sequence  string
	Frequency int
	/* The barcodes supporting this variant call */
	BarcodeEvidence []TenXBarcode
}

func RenderVCFSliceToJSON(rows []SimpleVCFRow) string {
	bytearray, err := json.Marshal(rows)
	if err != nil {
		panic(err)
	}
	return string(bytearray)
}

type VCFArray []*SimpleVCFRow

func TenXBarcodeInit(seq string) TenXBarcode {
	self := TenXBarcode{}
	self.SetBarcode(seq)
	return self
}

func (self *TenXBarcode) GetBarcode() string {
	barcode := SeqBytesToString(self.Barcode, int(self.Length))
	if self.GemGroup[0] == 0 {
		// Barcode has no gem group
		return barcode
	}

	// Support gem group string before AND after barcode for backwards compatibility
	if self.GemGroup[0] == '-' {
		return fmt.Sprintf("%s-%d", barcode, int(self.GemGroup[1]))
	} else if self.GemGroup[1] == '-' {
		return fmt.Sprintf("%d-%s", int(self.GemGroup[0]), barcode)
	} else {
		panic(fmt.Sprintf("TenXBarcode gem group field has invalid format %v", self.GemGroup))
	}
}

func (self *TenXBarcode) SetBarcode(seq string) {
	dash_offset := strings.IndexRune(seq, '-')
	if dash_offset > 0 {
		if IsSeq(seq[:dash_offset]) {
			// Barcode has format <seq>-<gem group>
			self.GemGroup[0] = '-'
			self.GemGroup[1] = ParseGemGroup(seq[dash_offset+1:])
			seq = seq[:dash_offset]
		} else if IsSeq(seq[dash_offset+1:]) {
			// Barcode has format <gem group>-<seq>
			self.GemGroup[0] = ParseGemGroup(seq[:dash_offset])
			self.GemGroup[1] = '-'
			seq = seq[dash_offset+1:]
		} else {
			panic(fmt.Sprintf("Barcode sequence %s does not contain barcode base string", seq))
		}
	} else {
		// Barcode has no gem group
		self.GemGroup[0] = 0
	}
	SeqStringToBytes(seq, &self.Barcode, &self.Length)
}

func (self TenXBarcode) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Barcode": self.GetBarcode(),
	})
}

func ParseGemGroup(str string) byte {
	gem_group, err := strconv.Atoi(str)
	if err != nil {
		panic(fmt.Sprintf("Gem group %s is not an integer", str))
	}
	if gem_group > 255 {
		panic(fmt.Sprintf("Gem group %d is greater than maximum value 255", gem_group))
	}
	return byte(gem_group)
}

func IsSeq(seq string) bool {
	for _, c := range seq {
		if c != 'A' && c != 'G' && c != 'C' && c != 'T' {
			return false
		}
	}
	return true
}

func SeqBytesToString(seq [MAX_BARCODE_BYTES]byte, l int) string {
	var val byte

	if l > 4*len(seq) {
		panic(fmt.Sprintf("Barcode sequence has length %d, greater than maximum length %d", l, 4*len(seq)))
	}

	c := make([]byte, l)
	for i, b := range seq {
		for j := 0; j < 4; j++ {
			idx := 4*i + j
			if idx < l {
				num := (b >> uint(2*j)) & 3
				switch num {
				case 0:
					val = 'A'
				case 1:
					val = 'G'
				case 2:
					val = 'C'
				case 3:
					val = 'T'
				}
				c[idx] = val
			} else {
				break
			}
		}
	}
	return string(c)
}

func SeqStringToBytes(seq string, b *[MAX_BARCODE_BYTES]byte, l *byte) {
	var val byte

	if len(seq) > 4*len(b) {
		panic(fmt.Sprintf("Barcode sequence %s is longer than max length %d", seq, 4*len(b)))
	}
	*l = byte(len(seq))

	for i, c := range seq {
		switch c {
		case 'A':
			val = 0
		case 'G':
			val = 1
		case 'C':
			val = 2
		case 'T':
			val = 3
		default:
			panic(fmt.Sprintf("Invalid character %c in DNA sequence string", c))
		}
		(*b)[i/4] |= val << uint(2*(i%4))
	}
}

/*
 * Compare the (correct) sorting order of two chromosomes. The order that this
 * yields is:
 * chr1, chr2, ... chr10, chr11, .. chr20, chr21, .. chrX, chrY.
 * The behavior is totally undefined if a and b are not in the format: chr{"X","Y",N,NN}.
 *
 * TODO: This is an ugly disaster. There has got to be a more pleasant approach.
 */
func ChromosomeCMP(a string, b string) int {
	if a < b {
		return -1
	}
	if a > b {
		return 1
	} else {
		return 0
	}
}

func (a VCFArray) Less(i, j int) bool {
	if a[i].Chromosome == a[j].Chromosome {
		return a[i].Position < a[j].Position
	} else {
		return (ChromosomeCMP(a[i].Chromosome, a[j].Chromosome)) < 0
	}
}

func (a VCFArray) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a VCFArray) Len() int {
	return len(a)
}

func ReadVCFWithCallback(path string,
	callback func(logical_row int,
		physical_row int,
		data *SimpleVCFRow) bool) error {

	data, err := ReadVCFToArray(path)
	if err != nil {
		return err
	}

	for i := 0; i < len(data); i++ {
		/*
		 * Oops! We just forgot what row the data was out.
		 */
		callback(i, 0, data[i])
	}

	return nil

}

func ReadVCFToArray(path string) ([]*SimpleVCFRow, error) {

	vcf_output_array := make([]*SimpleVCFRow, 0, 0)

	err := ReadVCFWithCallbackUnordered(path, func(logical int, physical int, data *SimpleVCFRow) bool {
		vcf_output_array = append(vcf_output_array, data)
		return true
	})

	if err != nil {
		return nil, err
	}

	sort.Sort(VCFArray(vcf_output_array))
	return vcf_output_array, nil
}

/*
 * Parse (very simply) a VCF file into SimpleVCFRow structs.
 * This parses a file that looks like:
 *
 * [chromosome] [offset] [reference-id] [reference-sequence] \
 * [alternative1,alternative2,...] [ignored float] [ignored string] \
 * [ignored string] [ignored string] [ALTERNATIVE/ALTERNATIVE]:.:.:.:.:\
 * [barcode;barcode;barcode,barcode,barcode]
 *
 * The data is parsed into the SimpleVCFRow struct defined above the callback
 * function is executed once per each row of the file.  Lines starting  with #
 * are ignored.
 *
 * BUG: This implementation will produce a row even if the alternatives
 * are listed as 0/0 which means there was no real variation called on that
 * row. However, other code (gene-index.go) expects ReadVCFWithCallback to
 * only return genuinly differences.
 *
 * |path| The path of the file to open |callback| a function to be called for
 * every record pulled from the file.
 */

func ReadVCFWithCallbackUnordered(path string,
	callback func(logical_row int,
		physical_row int,
		data *SimpleVCFRow) bool) error {

	fp, err := os.Open(path)
	defer fp.Close()

	/* What line number of the file are we on? */
	var line_number int
	/* How many records have we processed so far? (many lines are comments
	 * that we ignore*/
	var row_number int

	if err != nil {
		log.Printf("Cannot open %v: %v", path, err)
		return err
	}

	input_buffer := bufio.NewReader(fp)
	for {
		row := new(SimpleVCFRow)
		/* Grab the next line from the file */
		line, err := input_buffer.ReadString('\n')

		/* (probably an EOF condition. TODO: check for worse errors */
		if err != nil {
			break
		}
		line_number++
		if line_number%1000000 == 0 {
			log.Printf("VCF line: %v", line_number)
		}
		/* Skip comment lines */
		if line[0] == '#' {
			continue
		}
		var chromosome string
		var offset int
		var refId string
		var alternatives string
		var call_quality_string string
		var filter_string string
		var extras_string string
		var format_specifier string
		var reference string
		var variant_string string

		/* Parse data from the file */
		_, err = fmt.Sscanf(line, "%s%d%s%s%s%s%s%s%s%s",
			&chromosome,
			&offset,
			&refId,
			&reference,
			&alternatives,
			&call_quality_string,
			&filter_string,
			&extras_string,
			&format_specifier,
			&variant_string)

		if err != nil {
			log.Printf("Trouble parsing at line %v: %v",
				line_number,
				err)
			continue
		}
		row_number++

		/* Ignore filtered data */
		if filter_string != "." && filter_string != "PASS" && filter_string != "pass" {
			continue
		}

		row.Chromosome = NormalizeChromosomeName(chromosome)

		row.Position = offset
		row.Id = refId
		row.Reference = reference

		/* Try to grab the call quality. Leave it at zero if it doesn't
		 * parse
		 */
		fmt.Sscanf(call_quality_string, "%f", &row.CallQuality)

		/* Extract the specific keys out of the variant specifier that we care about, using
		 * the format_specifier to know where they are.
		 */
		ok, barcodes, genotype, phase_id_str, phase_quality, ref_frequency_str, frequency_string := ExtractVariantData(format_specifier, variant_string)

		if !ok {
			log.Printf("Trouble parsing per-variant data at line: %v: %v %v",
				line_number,
				format_specifier,
				variant_string)
		}

		var frequencies *[]string
		if frequency_string != "" {
			a := (strings.Split(frequency_string, ","))
			f := make([]string, len(a)+1)
			frequencies = &f
			f[0] = ref_frequency_str
			for j := 0; j < len(a); j++ {
				f[j+1] = a[j]
			}
		}
		row.PhaseQuality = phase_quality

		/*
		 * Extract variant information from the genotype field and
		 * match that against the reference and alternative fields to
		 * get the actual variant sequences.
		 * TODO: eww.... this is nasty!
		 */

		ok,
			s1,
			s2,
			wp,
			variant0_alternative_idx,
			variant1_alternative_idx := ExtractSequenceData(reference,
			alternatives,
			genotype)

		if !ok {
			continue;
		}

		row.Sequence[0].Sequence, row.Sequence[1].Sequence = s1, s2

		/*
		 * Extract the barcodes and attach them to the right sequence object.
		 */
		if barcodes != "" {
			barcodes_by_allele, ok := ParseBarcode(barcodes)
			if !ok {
				log.Printf("Can't parse barcode data: %v at: %v",
					barcodes,
					line_number)
				/* TODO XXX Do something here to mark the barcode data as bad*/

			} else {
				/* Copy barcode information for this variant into row object*/
				row.Sequence[0].BarcodeEvidence = barcodes_by_allele[variant0_alternative_idx]
				row.Sequence[1].BarcodeEvidence = barcodes_by_allele[variant1_alternative_idx]

				/* If we have supporting-read-frequencies, copy that into the row object,
				 * too.
				 */
				if frequencies != nil && len(*frequencies) > 0 {
					fmt.Sscanf((*frequencies)[variant0_alternative_idx],
						"%d",
						&row.Sequence[0].Frequency)
					fmt.Sscanf((*frequencies)[variant1_alternative_idx],
						"%d",
						&row.Sequence[1].Frequency)
				}
			}
		}

		/*
		 * Extract the phase id.
		 */
		row.WellPhased = false
		if phase_id_str != "" {
			var phase_id int
			_, err = fmt.Sscanf(phase_id_str, "%d", &phase_id)
			if err != nil {
				log.Printf("can't parse phase id data at %v: %v",
					line_number,
					phase_id_str)
				/* TODO: do something here to mark the phase data as bad */
			}
			row.PhaseId = phase_id
			if phase_id > 0 {
				row.WellPhased = wp
			}

		}

		callback(line_number, row_number, row)
	}
	return nil
}

/* This function takes a "variant format specifier" string and a variant string
 * and parses out the useful fields.
 * a variant format specifier is a string like AQ:PS:QQ:TR.... which is used
 * to find the right portions of the variant string which is also a series of
 * datums separated by colons.
 */
func ExtractVariantData(format_specifier string, variant string) (bool, string, string, string, float64, string, string) {
	var barcodes, genotype, phase_id, ref_frequency, frequencies string
	var phase_quality float64

	specs := strings.Split(format_specifier, ":")
	datas := strings.Split(variant, ":")

	for i := 0; i < len(specs) && i < len(datas); i++ {
		switch specs[i] {
		case ("GT"):
			genotype = datas[i]
		case ("PS"):
			phase_id = datas[i]
		case ("BX"):
			barcodes = datas[i]
		case ("RO"):
			ref_frequency = datas[i]
		case ("AO"):
			frequencies = datas[i]
		case ("PQ"):
			fmt.Sscanf(datas[i], "%f", &phase_quality)
		}
	}

	/* phase_id may be legit'ly blank for homozygous variants */
	if phase_id == "" {
		phase_id = "0"
	}

	if barcodes != "" && genotype != "" && phase_id != "" {
		return true, barcodes, genotype, phase_id, phase_quality, ref_frequency, frequencies
	} else {
		return false, barcodes, genotype, phase_id, phase_quality, ref_frequency, frequencies
	}
}

func ParseOneBarcodePhase(barcodes string) ([]TenXBarcode, bool) {
	b := make([]TenXBarcode, 0, 0)

	/* Separate out individual barcodes in this string. */
	barcode_str_array := strings.Split(barcodes, ";")

	for _, one_barcode := range barcode_str_array {
		if one_barcode == "" {
			continue
		}
		/* Separate the barcode from the quality score.
		 * TODO XXX BUG: This code is seriously incorrect! A barcode may be
		 * associated with multiple quality scores!
		 * This will return seriously wrong data for a barcode like
		 * ATAAG_12_55_09
		 */
		underbar_offset := strings.IndexRune(one_barcode, '_')
		if underbar_offset > 0 {
			var score int
			fmt.Sscanf(one_barcode[underbar_offset:len(one_barcode)], "%d", &score)
			/*
			 * TODO: Need to check for errors from Sscanf above and return
			 * false if it failed.
			 */

			/* Append this new barcode to the array.
			 * TODO: We may need to make this more efficient. This will
			 * cause a lot of reallocs and a lot of stack spilling.
			 */
			b = append(b, TenXBarcodeInit(one_barcode[0:underbar_offset]))
		} else {
			/*
			 * Can't understand this data. We return what we have and stop
			 */
			return b, false
		}
	}
	return b, true
}

/*
 * Parse a 10X barcode field from a VCF file.  This parses a string like
 * ATTTAAA_12;ACCATA_14,CCCCCCCC_99;GGGGGGGGG_0;ATATATATAT_43
 * into a pair of TenXBarcode arrays.  The "comma" splits the phase
 * and the "semicolon" separates the barcodes within that phase.
 */
func ParseBarcode(bc_string string) ([][]TenXBarcode, bool) {

	/* Split the two phases of barcode data */
	barcodes_per_phase := strings.Split(bc_string, ",")

	/* Only proceed if data looks ok so far */
	if len(barcodes_per_phase) < 1 {
		return nil, false
	}
	barcodes_for_alleles := make([][]TenXBarcode, len(barcodes_per_phase))

	ok := true
	for i, bc := range barcodes_per_phase {
		parsed_bc, parsed_ok := ParseOneBarcodePhase(bc)
		ok = ok && parsed_ok
		barcodes_for_alleles[i] = parsed_bc
	}
	/* Return success only if each phase was successful */
	return barcodes_for_alleles, ok
}

/*
 * This extracts data from the genotype field and returns the sequence information.
 * The genotype string is a string that looks like <integer>|<integer> or <integer>\<integer>
 * the index reference and alternatives string (0 = reference, 1 = first alternative, ....) and
 * the separator tells us if the call is well phased or not.
 */
func ExtractSequenceData(reference_string string, alternatives string, genotype string) (bool, string, string, bool, int, int) {
	var variant0_alternative_idx int
	var variant1_alternative_idx int
	var well_phased_chr rune
	if len(genotype) == 3 {
		well_phased_chr = rune(genotype[1])
		/* A . in the GT field means we're explicitly not making a call at all.
		 * So we just ignore those lines.
		 */
		if (genotype[0] == '.' || genotype[2] == '.') {
			return false, "", "", false, 0, 0
		}
		variant0_alternative_idx = int(genotype[0]) - int('0')
		variant1_alternative_idx = int(genotype[2]) - int('0')
	} else if len(genotype) == 1 {
		/* Ditto */
		if (genotype[0] == '.') {
			return false, "", "", false, 0, 0
		}
		variant0_alternative_idx = int(genotype[0]) - int('0')
		variant1_alternative_idx = variant0_alternative_idx
		well_phased_chr = '|'
	} else {
		log.Printf("BADNESS AT: %v", genotype)
		return false, "", "", false, 0, 0
	}

	alternatives_list := strings.Split(alternatives, ",")

	var sequence0 string
	var sequence1 string
	var wellphased bool

	/* TODO: We should bounds check variant_alternatieve_idx before
	 * using it */
	if variant0_alternative_idx == 0 {
		sequence0 = reference_string
	} else {

		sequence0 = alternatives_list[variant0_alternative_idx-1]
	}

	if variant1_alternative_idx == 0 {
		sequence1 = reference_string
	} else {
		sequence1 = alternatives_list[variant1_alternative_idx-1]
	}

	/* A "/" in the variant specification means that the variant
	 * was well phased.  A | means that it was not well phased.
	 * TODO: we should figure out how to grab and process the phase
	 * block, here.*/

	wellphased = well_phased_chr == '|'

	/* Ignore records wherein both hapotypes have the reference */
	if variant0_alternative_idx == 0 && variant1_alternative_idx == 0 {
		return false, sequence0, sequence1, wellphased, variant0_alternative_idx, variant1_alternative_idx
	}

	return true, sequence0, sequence1, wellphased, variant0_alternative_idx, variant1_alternative_idx
}
