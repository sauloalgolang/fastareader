/*
Package fastatools contains tools for fasta processing
*/

package fastatools


import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"log"
	"os"
	"strings"
)

// Error codes returned by failures to parse
var (
	ErrInternal   = errors.New("fastatools: internal error" )
	ErrInvalidSeq = errors.New("fastatools: invalid fasta"  )
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

type SeqData struct {
	SeqName  string
	Sequence []byte
}

func (seqd *SeqData) Size() (size int64) {
	return int64(len(seqd.Sequence))
}

func (seqd *SeqData) Print () {
	log.Printf("SeqData: NAME '%s' SIZE %d\n", seqd.SeqName, seqd.Size())
}

func (seqd *SeqData) SaveToFasta(fo *os.File) {
        fmt.Fprintf(fo, ">%s\n", seqd.SeqName)

	leng := len(seqd.Sequence)

	start := 0
	end   := 80
	sum   := 0

	for ; start < leng; start += 80 {
		end = start + 80

		if end >= leng {
			end = leng
		}

		if end != start {
			frag := string(seqd.Sequence[start:end])
        		fmt.Fprintf(fo, frag + "\n")
			sum  += len(frag)
		}
	}

	log.Println( "LENG", leng )
	log.Println( "SUM ", sum  )
}



/*
ReadFileLineByLine: reads line by line using callback
input             : fi   *os.File
                    clbk func(*string)(res bool)
*/
func ReadFileLineByLine(fi *os.File, clbk func(*string)(res bool)) {
	scanner  := bufio.NewScanner(fi)
	scanner.Split(bufio.ScanLines)

	for scanner.Scan() {
		line     := scanner.Text()
		res      := clbk(&line)
		if ! res {
			break
		}
  	}
}


/*
readSeqFromFasta: read a fasta Sequence inside a file
input           : file     *os.File
return          : sd       *SeqData
*/
func readSeqFromFasta(file *os.File) (sd *SeqData) {
	sd       = new(SeqData)

	var buffer bytes.Buffer

	processFastaLine := func( line *string ) ( res bool ){
		if (*line)[0] == byte('>') {
			if sd.SeqName == "" { // first
				sd.SeqName = strings.TrimSpace((*line)[1:])
				log.Println("Seq", sd.SeqName, "STARTING")
				return true

			} else { //next
				log.Println("Seq", sd.SeqName, "DONE"    )
				return false

			}
		} else {
			if len(*line) != 0 {
				buffer.WriteString(*line)
			}

			return true
		}
	}

	log.Println("Seq", sd.SeqName, "READING")
	ReadFileLineByLine(file, processFastaLine)

	log.Println("Seq", sd.SeqName, "CONVERTING")
	sd.Sequence = []byte(buffer.String())

	return sd
}

func OpenAndSeek(filename string, position int64) (file *os.File) {
        file, err := os.Open(filename)
	check(err)

        _, err  = file.Stat()
	check(err)

	//log.Println(d)

	_, err = file.Seek(position, 0)
	check(err)

	return file
}


/*
ReadFastaSeq: read a fasta Sequence inside a file
input       : filename string
              position string
return      : sd       *SeqData
*/
func ReadFastaSeq(filename string, position int64) (sd *SeqData) {
	file := OpenAndSeek(filename, position)
	defer file.Close()

	//log.Println("new positions", pos)

	sd   = readSeqFromFasta(file)

	sd.Print()

	return sd
}



//func(string)bool
func ReadFastaSeqIter(filename string, position int64, clbk func(*string)bool) {
	file := OpenAndSeek(filename, position)
	defer file.Close()

	log.Println("ReadFastaSeqIter :: filename:", filename, "position:", position)

	//kmertools.GetExtractKmersIterClbk( kmerSize, data))

	ReadFileLineByLine(file, clbk)

	log.Println("ReadFastaSeqIter :: filename:", filename, "position:", position, "DONE")
}


func GetPipeFastaBackClbk( of *os.File ) func(*string)bool {
	seqName               := ""

	pipeFastaBackIterClbk := func(line *string)bool {
                if (*line)[0] == '>' {
                        if seqName == "" { // first
                                seqName = strings.TrimSpace((*line)[1:])
                                log.Println("Seq", seqName, "STARTING")
				of.WriteString( string(*line) + "\n" )
                                return true

                        } else { //next
                                log.Println("Seq", seqName, "DONE"    )
                                return false

                        }
                } else {
                        if len(*line) != 0 {
				of.WriteString( string(*line) + "\n" )
                        }

                        return true
                }
	}

	return pipeFastaBackIterClbk
}
