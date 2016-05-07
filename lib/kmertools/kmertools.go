
/*
Package kmertools contains tools for kmer processing
*/

package kmertools


import (
//	"bufio"
//	"bytes"
	"errors"
	"fmt"
	"log"
	"os"
	"strings"
	"sync"
	"math"
//	"unicode"
)


import (
        "github.com/sauloalgolang/fastareader/lib/fastatools"
)

// Error codes returned by failures to parse
var (
	ErrInternal   = errors.New("kmertools: internal error" )
	ErrInvalidSeq = errors.New("kmertools: invalid fasta"  )
)


// Available output formats
var AvailableFormats = [3]string{ "fasta", "list", "csv" }


/*
check : This helper will streamline our error checks below.
src   : https://gobyexample.com/reading-files
inputs: e error
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





/*
Map    : runs a function to all elements of a slice, returning a new slice
inputs : f func(byte) byte
         vs []byte
outputs: []byte
*/
func Map(f func(byte) byte, vs []byte) ([]byte) {
    vsm := make([]byte, len(vs))
    for i, v := range vs {
        vsm[i] = f(v)
    }
    return vsm
}



/*
rcer   : completement a byte nucleotide, returning a byte
inputs : r byte
outputs: byte
src    : https://golang.org/pkg/strings/#Map
*/
func rcer(r byte) (byte) {
	switch {
		case r == 'A': return 'T'
		case r == 'a': return 'T'

		case r == 'C': return 'G'
		case r == 'c': return 'G'

		case r == 'G': return 'C'
		case r == 'g': return 'C'

		case r == 'T': return 'A'
		case r == 't': return 'A'

		case r == 'N': return 'N'
		case r == 'n': return 'N'
	}
	log.Panic("unknown byte", string(r), " ", r)
	return r
}



/*
Reverse: reverses a byte slice inplace
inputs : sequence []byte
src    : http://golangcookbook.com/chapters/arrays/reverse/
*/
func Reverse(sequence []byte) {
	for i, j := 0, len(sequence)-1; i < j; i, j = i+1, j-1 {
		sequence[i], sequence[j] = sequence[j], sequence[i]
	}
}


func ReverseComplement(src []byte) ([]byte) {
	dst := make([]byte, len(src))

	ls  := len(src)

	for i := 0; i < ls; i++ {
		dst[i] = rcer(src[ls-i-1])
	}

	return dst
}



/*
IsInSlice: checks whether a byte is present in a slice
inputs   : a byte
           list []byte
outputs  : bool
src      : http://stackoverflow.com/questions/15323767/does-golang-have-if-x-in-construct-similar-to-python
*/
func IsInSlice(a byte, list []byte) bool {
    for _, b := range list {
        if b == a {
            return true
        }
    }
    return false
}





// https://tour.golang.org/concurrency/9
// SafeCounter is safe to use concurrently.
type Data struct {
	v       map[string]int
	Total   uint64
	MaxSize int
	mux     sync.RWMutex
}

// http://stackoverflow.com/questions/4498998/how-to-initialize-members-in-go-struct
func (c *Data) New(kmerSize int) {
    c.MaxSize = int(math.Pow(4,float64(kmerSize))/2)
    c.v       = make(map[string]int, 0)//c.MaxSize)
}

// Inc increments the counter for the given key.
func (c *Data) Inc(key string) {
	c.mux.Lock()
	// Lock so only one goroutine at a time can access the map c.v.
	c.v[key]++
	c.Total++

	if len(c.v) > c.MaxSize {
		log.Panic(key)
	}

	c.mux.Unlock()

	if (c.Total % 10000000) == 0 {
		log.Println("Total Kmers:", c.Total, "Unique", len(c.v))
	}
}

// Value returns the current value of the counter for the given key.
func (c *Data) Value(key string) int {
	c.mux.Lock()
	// Lock so only one goroutine at a time can access the map c.v.
	defer c.mux.Unlock()
	return c.v[key]
}

func (c *Data) Len() int {
	c.mux.Lock()
	defer c.mux.Unlock()
	return len(c.v)
}

/*
SaveKmers: save data into a fasta in a given format
inputs   : outFileName string
           as          string
           data        map[string]int
*/
func (c *Data) SaveAs(outFileName string, as string) {
	outFileNameTmp := outFileName + ".tmp"

	//log.Println(c)

	fo, err        := os.Create(outFileNameTmp)
	check(err)
	defer func() {
		os.Remove(outFileNameTmp)
	}()
	defer fo.Close()

	i := 0
	for k, v := range c.v {
		i++

		//log.Println(i,v,k)

		if as == "fasta" {
			fmt.Fprintf(fo, ">%d count: %d\n%s\n\n", i, v, k)
		} else
		if as == "list"  {
			fmt.Fprintf(fo, "%s\n", k)
		} else
		if as == "csv"   {
			fmt.Fprintf(fo, "%s\t%d\n", k, v)
		}
	}

	fo.Close()

	os.Rename(outFileNameTmp, outFileName)
}



/*
ExtractKmers: Extract kmers, storing them in data
input       : seqd     *fastatools.SeqData
              kmerSize int
              data     map[string]int
*/

func ExtractKmers(seqd *fastatools.SeqData, kmerSize int, data *Data) {
	ExtractKmersFromSeq(seqd.Sequence, seqd.SeqName, kmerSize, data)
}




/*
ExtractKmers: Extract kmers, storing them in data
input       : seqd     *fastatools.SeqData
              kmerSize int
              data     map[string]int
*/
func ExtractKmersFromSeq(sequence []byte, seqName string, kmerSize int, data *Data) {
	if len(sequence) < kmerSize {
		return
	}

	seqLen  := len(sequence)

	rc      := ReverseComplement(sequence)

	fwd     := ""
	rev     := ""

	fStart  := 0
	rStart  := 0

	fEnd    := 0
	rEnd    := 0

	hadN    := true

	for fStart = 0; fStart < (seqLen-kmerSize+1); fStart++ {
		fEnd   = fStart + kmerSize
		fwd    = string(sequence[fStart:fEnd])

		if ( hadN ) {
			hadN = false
			for _, n := range fwd {
				if ! ( n == 'A' || n == 'C' || n == 'G' || n == 'T' ) {
					hadN = true
					break
				}
			}
			if hadN {
				continue
			}
		} else {
			n := fwd[kmerSize-1]
			if ! ( n == 'A' || n == 'C' || n == 'G' || n == 'T' ) {
				hadN = true
				continue
			}
		}

		rStart = seqLen - fStart   - kmerSize
		rEnd   = rStart + kmerSize
		rev    = string(rc[      rStart:rEnd])

		/*
		if p {
		log.Printf( "%s", sequence)
		log.Printf( "%s", rc      )
		log.Printf( "fStart %12d fEnd %12d len %d fwd %s", fStart, fEnd, len(fwd), fwd)
		log.Println( []byte(fwd) )
		log.Printf( "rStart %12d rEnd %12d len %d rev %s", rStart, rEnd, len(rev), rev)
		log.Println( []byte(rev) )
		log.Println("smaller", (fwd<rev))
		log.Println()
		}
		*/

		if fwd<rev {
			data.Inc(fwd)
		} else {
			data.Inc(rev)
		}
	}
}





func GetExtractKmersIterClbk( kmerSize int, data *Data ) func(*string)bool {
        log.Println("GetExtractKmersIterClbk")

	sequence := make([]byte, 0)
	seqName  := ""

	lineTot := 0
	lineSeq := 0

	ExtractKmersIterClbk := func(line *string)bool {
		lineTot++
                if (*line)[0] == '>' {
                        if seqName == "" { // first
                                seqName  = strings.TrimSpace((*line)[1:])
                                log.Println("Seq", seqName, "STARTING")
                                return true

                        } else { //next
                                log.Println("Seq", seqName, "DONE"    )
                                return false

                        }
                } else {
			lineSeq++
                        if len(*line) != 0 {
				//log.Println(len(*line), *line)

				sequence = append( sequence, []byte(*line)... )

				if len(sequence) >= kmerSize {
					ExtractKmersFromSeq(sequence, seqName, kmerSize, data)
					sequence = sequence[len(sequence)-kmerSize+1:]
				}

				//log.Println(len(sequence), string(sequence))

				//log.Println()
                        }

                        return true
                }
	}

	return ExtractKmersIterClbk
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
