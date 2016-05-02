/*
Package fastareader reads a fasta file using go routines
*/

package main

import (
	"errors"
	"fmt"
//	"io"
	"log"
	"os"
	"runtime"
)


import (
	"github.com/sauloalgolang/fastareader/lib/fastaindex"
	"github.com/sauloalgolang/fastareader/lib/fastatools"
)

// http://www.golangbootcamp.com/book/tricks_and_tips
// compile passing -ldflags "-X main.Build <build sha1>"
var Build string

// Error codes returned by failures to parse
var (
	ErrInternal   = errors.New("fastareader: internal error"  )
	ErrInvalidSeq = errors.New("fastareader: invalid sequence")
)


/*
check: This helper will streamline our error checks below.
src  : https://gobyexample.com/reading-files
input: e error
*/
func check(e error) {
        if e != nil {
                log.Fatal( e )
                panic(e)
        }
}



/*
main: checks if index exists, creating it otherwise, read index and create a go routine to read each sequence
*/
func main() {
	if Build != "" {
		log.Println("fastareader build:", Build)
	}


	argsWithoutProg := os.Args[1:]

        if len(argsWithoutProg) != 1 {
                log.Println("no argument or too many arguments given")
                os.Exit(1)
        }


        filename        := argsWithoutProg[0]

	f, err := os.Open(filename)
	check(err)
	defer f.Close()

	_, err = f.Stat()
	check(err)
	//log.Print(d)


	var numCPU = runtime.GOMAXPROCS(0)

	log.Println("numCPU",numCPU)


	idxName := filename + ".idx"
	if _, err := os.Stat(idxName); os.IsNotExist(err) {
		log.Println("Index does not exists. creating")
		fastaindex.CreateFastaIndex(filename)
	} else {
		log.Println("Index alread exists")
	}

	log.Println("Reading Index")

	idxData := fastaindex.ReadFastaIndex(filename)
	seqData := make([]*fastatools.SeqData, len(*idxData))

	for _, idx := range *idxData {
		idx.Print()
	}


	log.Println("Reading File")
	var limit  = make(chan int, numCPU)
	var waiter = make(chan int)
	for _, idx := range *idxData {
		//idx.Print()

		f := func (idx2 *fastaindex.IdxData) {
			seqd := fastatools.ReadFastaSeq(filename, idx2.SeqPos)
			log.Printf("RES: IDX: NAME '%s' ID %d SIZE %d POSITION %d FASTA: NAME '%s' SIZE %d\n", idx2.SeqName, idx2.SeqId, idx2.SeqSize, idx2.SeqPos, seqd.SeqName, seqd.Size())
			if (( idx2.SeqName != seqd.SeqName ) || (idx2.SeqSize != seqd.Size())) {
				log.Fatal(fmt.Sprintf("Sequence mismatch. expexted '%s', found '%s'. Expected size %d, found %d", idx2.SeqName, seqd.SeqName, idx2.SeqSize, seqd.Size()))
				os.Exit(1)
			}
			seqData[idx2.SeqId - 1] = seqd
			waiter <- 1
		}


		go func(w func(idx2 *fastaindex.IdxData), idx2 *fastaindex.IdxData) {
			limit <- 1
			w(idx2)
			<-limit
		}(f, idx)
	}

	log.Println("Waiting")

	for i := 1; i <= len(*idxData); i++ {
		<-waiter
	}

	log.Println("Done")

	for _, seq := range seqData {
		seq.Print()
	}
}
