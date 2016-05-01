/*
Package fastatools contains tools for fasta processing
*/

package fastatools

import (
	"bufio"
	"bytes"
	"errors"
//	"fmt"
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
ReadFastaSeq: read a fasta Sequence inside a file
input       : filename string
              position string
return      : *SeqData
*/
func ReadFastaSeq(filename string, position int64) *SeqData {
        fi, err := os.Open(filename)
	check(err)
	defer fi.Close()

        _, err  = fi.Stat()
	check(err)

	//log.Println(d)

	pos, err := fi.Seek(position, 0)
	check(err)
	log.Println("new positions", pos)

	scanner  := bufio.NewScanner(fi)
	scanner.Split(bufio.ScanLines)

	sd       := new(SeqData)

	var buffer bytes.Buffer

	for scanner.Scan() {
		line     := scanner.Text()

		if line[0] == byte('>') {
			if sd.SeqName == "" { // first
				sd.SeqName = strings.TrimSpace(line[1:])
				log.Println("Seq", sd.SeqName)
			} else { //next
				log.Println("Seq", sd.SeqName, "DONE")
				break
			}
		} else {
			if len(line) != 0 {
				//fmt.Print(".")
				buffer.WriteString(line)
				//sd.Sequence += line
			}
		}
  	}

	sd.Sequence = buffer.String()

	sd.Print()

	return sd
}

