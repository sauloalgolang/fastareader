/*
Package fastatools contains tools for fasta processing
*/

package fastatools


import (
	"bufio"
	"bytes"
	"errors"
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
	Sequence string
}

func (seqd *SeqData) Size() (size int64) {
	return int64(len(seqd.Sequence))
}

func (seqd *SeqData) Print () {
	log.Printf("SeqData: NAME '%s' SIZE %d\n", seqd.SeqName, seqd.Size())
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

	ReadFileLineByLine(file, processFastaLine)

	sd.Sequence = buffer.String()

	return sd
}


/*
ReadFastaSeq: read a fasta Sequence inside a file
input       : filename string
              position string
return      : sd       *SeqData
*/
func ReadFastaSeq(filename string, position int64) (sd *SeqData) {
        file, err := os.Open(filename)
	check(err)
	defer file.Close()

        _, err  = file.Stat()
	check(err)

	//log.Println(d)

	_, err = file.Seek(position, 0)
	check(err)

	//log.Println("new positions", pos)

	sd = readSeqFromFasta(file)

	sd.Print()

	return sd
}
