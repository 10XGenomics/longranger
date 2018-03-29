package test

import "reads/validate_barcodes"
import "testing"
import "fmt"

func TestBCCorrection1(t *testing.T) {

	barcode_validator := validate_barcodes.NewBarcodeValidator("inputs/737K-april-2014.txt", 2.0, "inputs/bc_counts.json", 0.975, "1")

	barcode, correct := barcode_validator.ValidateBarcode("ACGCCAGTAGATTG\n", []byte{66, 66, 36, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, byte('\n')})
	//ACACCAGTAGATTG
	fmt.Println(barcode, correct)
	check(t, string(barcode) == "ACACCAGTAGATTG\n", "barcode should error correct to 'ACACCAGTAGATTG\\n'")
	check(t, correct, "barcode should be valid")
	barcode_validator2 := validate_barcodes.NewBarcodeValidator("inputs/737K-april-2014.txt", 2.0, "inputs/bc_counts2.json", 0.975, "1")
	barcode2, correct2 := barcode_validator2.ValidateBarcode("CTCGTCTCCACCCC\n", []byte("<'<B7''0<<7<'<\n"))
	fmt.Println(string(barcode2), correct2)
	check(t, string(barcode2) == "CTCGTCTCCACCTC\n", "barcode should error correct to 'CTCGTCTCCACCTC\\n'")
	check(t, correct2, "barcode should be valid")
	barcode3, correct3 := barcode_validator2.ValidateBarcode("ACGCCAGTAGATTG\n", []byte{66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, byte('\n')})
	check(t, string(barcode3) == "ACGCCAGTAGATTG\n", "barcodes shouldnt be corrected")
	check(t, !correct3, "should be false, not valid barcode")

}

func check(t *testing.T, test bool, err string) {
	if !test {
		t.Error(err)
	}
}
