// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package formats

import (
	"strings"
)

func NormalizeChromosomeName(in string) string {

	if strings.HasPrefix(in, "chr") ||
		strings.HasPrefix(in, "CHR") {
		return in
	} else {
		return "chr" + in
	}

}
