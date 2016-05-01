/*
Package fastaindexer reads a fasta file and creates a index
*/

package main


import (
	"errors"
//	"fmt"
//	"io"
	"log"
	"os"
)



import (
	"github.com/sauloalgolang/fastareader/lib/fastaindex"
//	"github.com/sauloalgolang/fastareader/lib/fastatools"
)


// Error codes returned by failures to parse
var (
	ErrInternal   = errors.New("fastaindexer: internal error"  )
	ErrInvalidSeq = errors.New("fastaindexer: invalid sequence")
)


/*
main: checks if index exists, creating it otherwise, read index and create a go routine to read each sequence
*/
func main() {
	argsWithoutProg := os.Args[1:]

	if len(argsWithoutProg) != 1 {
		log.Println("no argument or too many arguments given")
		os.Exit(1)
	}

	filename        := argsWithoutProg[0]

	fastaindex.CreateFastaIndex(filename)
}
