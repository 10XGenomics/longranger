package formats

import (
	"unsafe"
)

/* Convert a 4-byte array to a float32 */
func BytesToFloat(b []byte) float32 {
	if b == nil {
		return 0.0
	} else {
		return *(*float32)(unsafe.Pointer(&b[0]))
	}
}

/*
 * Convert a little-endiat byte array to an integer.
 */
func BytesToInt(b []byte) int {

	if b == nil {
		return 0.0
	} else {
		var x int
		l := uint(len(b))
		for i := uint(0); i < l; i++ {
			x |= int(b[i]) << (8 * i)
		}
		return x
	}
}
