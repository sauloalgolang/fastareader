
/*
Package kmertools contains tools for kmer processing
*/

package kmertools


import (
//	"bufio"
//	"bytes"
	"errors"
	"log"
//	"os"
	"strings"
	"unicode"
)


import (
        "github.com/sauloalgolang/fastareader/lib/fastatools"
)

// Error codes returned by failures to parse
var (
	ErrInternal   = errors.New("kmertools: internal error" )
	ErrInvalidSeq = errors.New("kmertools: invalid fasta"  )
)

/*
check: This helper will streamline our error checks below.
src  : https://gobyexample.com/reading-files
input: e error
*/
func check(e error) {
	if e != nil {
		log.Fatal(e)
		panic(e)
	}
}


/*
src: https://github.com/golang/exp/blob/master/shootout/reverse-complement.go
*/
/*
var complement = [256]uint8{
	'A': 'T', 'a': 'T',
	'C': 'G', 'c': 'G',
	'G': 'C', 'g': 'C',
	'T': 'A', 't': 'A',
	'U': 'A', 'u': 'A',
	'M': 'K', 'm': 'K',
	'R': 'Y', 'r': 'Y',
	'W': 'W', 'w': 'W',
	'S': 'S', 's': 'S',
	'Y': 'R', 'y': 'R',
	'K': 'M', 'k': 'M',
	'V': 'B', 'v': 'B',
	'H': 'D', 'h': 'D',
	'D': 'H', 'd': 'H',
	'B': 'V', 'b': 'V',
	'N': 'N', 'n': 'N',
}
*/

func ToUpper(r []rune) {
    for i, v := range r {
        r[i] = unicode.ToUpper(v)
    }
}

func Map(f func(rune) rune, vs []rune) ([]rune) {
    vsm := make([]rune, len(vs))
    for i, v := range vs {
        vsm[i] = f(v)
    }
    return vsm
}

//https://golang.org/pkg/strings/#Map
func rcer(r rune) rune {
	switch {
		case r == 'A': return 'T'
		case r == 'C': return 'G'
		case r == 'G': return 'C'
		case r == 'T': return 'A'
	}
	return r
}

//http://golangcookbook.com/chapters/arrays/reverse/
func Reverse(sequence []rune) {
	for i, j := 0, len(sequence)-1; i < j; i, j = i+1, j-1 {
		sequence[i], sequence[j] = sequence[j], sequence[i]
	}
}

func IsInSlice(a rune, list []rune) bool {
    for _, b := range list {
        if b == a {
            return true
        }
    }
    return false
}

/*
ExtractKmers: Extract kmers
input       : seqd     *fastatools.SeqData
              kmerSize int
              data     map[string]int
*/
func ExtractKmers(seqd *fastatools.SeqData, kmerSize int, data map[string]int) {
	seqLen   := len(seqd.Sequence)

	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: LENGTH", seqLen, "KMER SIZE", kmerSize)

	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: UPPERING")
	ToUpper(seqd.Sequence)

	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: COMPLEMENTING")
	rc := Map(rcer, seqd.Sequence)
	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: REVERSING")
	Reverse(rc)

	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: FWD:", string(seqd.Sequence[0:10]), "-", string(seqd.Sequence[len(seqd.Sequence)-10:]))
	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: RCP:", string(           rc[0:10]), "-", string(           rc[len(rc           )-10:]))

	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: COUNTING")

	fwd     := ""
	rev     := ""

	fStart  := 0
	rStart  := 0

	fEnd    := 0
	rEnd    := 0

	for fStart = 0; fStart < (seqLen-kmerSize+1); fStart++ {
		if seqLen > 1000 {
			if (fStart != 0) && (fStart % (seqLen/10) == 0) {
				log.Printf("Start %12d %3.0f%%\n", fStart, ((float64(fStart)/float64(seqLen))*float64(100)))
			}
		}

		fEnd   = fStart + kmerSize
		fwd    = string(seqd.Sequence[fStart:fEnd])

		if ! strings.ContainsRune(fwd, 'N') {

		rStart = seqLen - fStart   - kmerSize
		rEnd   = rStart + kmerSize
		rev    = string(rc[           rStart:rEnd])

		/*
		log.Printf( "fStart %12d fEnd %12d fwd %s", fStart, fEnd, fwd)
		log.Printf( "rStart %12d rEnd %12d rev %s", rStart, rEnd, rev)
		log.Println("smaller", (fwd<rev))
		log.Println()
		*/

		if fwd<rev {
			data[fwd] += 1
		} else {
			data[rev] += 1
		}

		}
	}

	//log.Print(data)
	log.Println("EXTRACTING KMER FROM:", seqd.SeqName, " :: DONE :: KMERS", len(data))
}



/*
https://godoc.org/bitbucket.org/santucco/btree
https://godoc.org/github.com/bsm/rumcask/btree
https://godoc.org/github.com/cznic/b
https://godoc.org/github.com/droxer/binarysearch
https://godoc.org/github.com/droxer/btree
https://godoc.org/github.com/droxer/bloomfilter
https://godoc.org/github.com/droxer/mergesort
https://godoc.org/github.com/droxer/quicksort
https://godoc.org/github.com/gopkg/container/btree
https://godoc.org/github.com/Workiva/go-datastructures/bitarray
https://godoc.org/github.com/Workiva/go-datastructures/hashmap/fastinteger
https://godoc.org/github.com/Workiva/go-datastructures/rtree
https://godoc.org/github.com/Workiva/go-datastructures/rangetree
https://godoc.org/github.com/Workiva/go-datastructures/set
https://godoc.org/github.com/Workiva/go-datastructures/sort
https://godoc.org/github.com/Workiva/go-datastructures/trie/xfast
https://godoc.org/github.com/Workiva/go-datastructures/trie/yfast
*/
