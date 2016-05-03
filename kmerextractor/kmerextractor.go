/*
Package kmerextracter reads a fasta file using go routines
*/

package main

import (
	"errors"
	"fmt"
	"flag"
	"log"
	"os"
	"runtime"
)


import (
	"github.com/sauloalgolang/fastareader/lib/fastaindex"
	"github.com/sauloalgolang/fastareader/lib/fastatools"
	"github.com/sauloalgolang/fastareader/lib/kmertools"
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
                //panic(e)
        }
}


var filename string
var format   string
var kmerSize int
var threads  int
func init() {
	if Build != "" {
		log.Println("kmerextracter build:", Build)
	}



	flag.StringVar(&filename, "filename",      "", "input fasta")
	flag.StringVar(&format  , "format"  , "fasta", "format: fasta, csv, list" )
	flag.IntVar(   &kmerSize, "kmersize",       0, "kmer size"  )
	flag.IntVar(   &threads , "threads" ,       0, "number of threads. 0 for max"  )
	flag.Parse()



	if filename == "" {
		flag.PrintDefaults()
                log.Fatal("No input file given\n")
	}

	if _, err := os.Stat(filename); os.IsNotExist(err) {
		flag.PrintDefaults()
                log.Fatal("Input file '" + filename + "' does not exist\n")
	}

	fi, err := os.Open(filename)
	check(err)
	//defer fi.Close()

	st, err := fi.Stat()
	check(err)
	//log.Print(d)

	// http://stackoverflow.com/questions/8824571/golang-determining-whether-file-points-to-file-or-directory
	switch mode := st.Mode(); {
		case mode.IsDir():
			// do directory stuff
			flag.PrintDefaults()
			log.Fatal("File '"+filename+"' is a directory\n")
	}
	fi.Close()




	fmatch := false
	for _, fmt := range kmertools.AvailableFormats {
		if format == fmt {
			fmatch = true
			break
		}
	}

	if ! fmatch {
		flag.PrintDefaults()
		log.Println("Invalid format: '" + format)
		log.Println("Possibilities are:")
		for _, fmt := range kmertools.AvailableFormats {
			log.Println("\t"+fmt)
		}
		os.Exit(1)
	}



	if kmerSize == 0 {
		flag.PrintDefaults()
		log.Fatal("No kmer size set\n")
	}



	var numCPU = runtime.GOMAXPROCS(0)
	if threads < 0 {
		flag.PrintDefaults()
                log.Fatal("Number of threads (",threads,") must be greater or equalt to 0\n")
	} else
	if threads == 0 {
		threads = numCPU
	}

	log.Println("threads",threads)

	if threads > numCPU {
		log.Println("Number of threads (", threads, ") greater than number of CPUs (", numCPU, "). expect slow downs")
	}
}



/*
main: checks if index exists, creating it otherwise, read index and create a go routine to read each sequence
*/
func main() {
	log.Println("Reading Index")


	idxData         := fastaindex.ReadFastaIndexCreatingIfNotExists(filename)
	//seqData         := make([]*fastatools.SeqData, len(*idxData))


	for _, idx := range *idxData {
		idx.Print()
	}



	log.Println("Reading File")
	limit           := make(chan int, threads)
	waiter          := make(chan int         )
	data            := make(map[string]int   )

	for _, idx := range *idxData {
		//idx.Print()

		f := func (idx2 *fastaindex.IdxData) {
			seqd := fastatools.ReadFastaSeq(filename, idx2.SeqPos)

			log.Printf("RES: IDX: NAME '%s' ID %d SIZE %d POSITION %d FASTA: NAME '%s' SIZE %d\n", idx2.SeqName, idx2.SeqId, idx2.SeqSize, idx2.SeqPos, seqd.SeqName, seqd.Size())

			if (( idx2.SeqName != seqd.SeqName ) || (idx2.SeqSize != seqd.Size())) {
				log.Fatal(fmt.Sprintf("Sequence mismatch. expexted '%s', found '%s'. Expected size %d, found %d", idx2.SeqName, seqd.SeqName, idx2.SeqSize, seqd.Size()))
				os.Exit(1)
			}

			//seqData[idx2.SeqId - 1] = seqd
			kmertools.ExtractKmers(seqd, kmerSize, data)

			seqd.Sequence = make([]rune,0)

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


	log.Println("Done Counting")


	log.Println("Kmers", len(data))//, data)


	outFileName     := fmt.Sprintf("%s_%d.kmers.%s", filename, kmerSize, format)
	log.Println("Saving to", outFileName)
	kmertools.SaveKmers(outFileName, format, data)


	log.Println("Done")

	/*
	for _, seq := range seqData {
		seq.Print()
	}
	*/
}
