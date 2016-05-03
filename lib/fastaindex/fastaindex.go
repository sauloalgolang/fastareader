/*
Package fastaindex contains routines to work with fasta file indexes
*/

package fastaindex

import (
	"bufio"
	"errors"
	"fmt"
//	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

// Error codes returned by failures to parse
var (
        ErrInternal   = errors.New("fastareader: internal error"    )
        ErrInvalidFa  = errors.New("fastareader: invalid fasta file")
        ErrInvalidIdx = errors.New("fastareader: invalid index file")
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

type IdxData struct {
	SeqName string
	SeqId   int
	SeqSize int64
	SeqPos  int64
}

func (idx *IdxData) Print () {
	log.Printf("IDX: NAME '%s' ID %d SIZE %d POSITION %d\n", idx.SeqName, idx.SeqId, idx.SeqSize, idx.SeqPos)
}

func (idx *IdxData) Write (file *os.File) {
	fmt.Fprintf(file, "%s\t%d\t%d\t%d\n", idx.SeqName, idx.SeqId, idx.SeqSize, idx.SeqPos)
}

func (idx *IdxData) Read (line string) {
	cols        := strings.Split(line, "\t")

	//log.Println(cols, len(cols))

	idx.SeqName   =              cols[0]
	idx.SeqId  ,_ = strconv.Atoi(    cols[1]        )
	idx.SeqSize,_ = strconv.ParseInt(cols[2], 10, 64)
	idx.SeqPos ,_ = strconv.ParseInt(cols[3], 10, 64)
}



/*
CreateFastaIndex: index a fasta file
input           : filename string
creates         : filename.idx
*/
func CreateFastaIndex(filename string) {
        fi, err := os.Open(filename)
	check(err)
	defer fi.Close()

        d, err := fi.Stat()
	check(err)

	log.Println(d)

	// open output file
	idxName    := filename+".idx"
	idxNameTmp := idxName + ".tmp"
	fo, err    := os.Create(idxNameTmp)
	check(err)
	defer func() {
		os.Remove(idxNameTmp)
	}()
	defer fo.Close()

	scanner  := bufio.NewScanner(fi)
	scanner.Split(bufio.ScanLines)

	idx      := new(IdxData)
	position := 0

	for scanner.Scan() {
		line     := scanner.Text()
		position += len(line) + 1

		if line[0] == byte('>') {
			if idx.SeqName != "" {
				idx.Print()
				idx.Write(fo)
			}

			idx.SeqPos   = int64(position - len(line) - 1)
			idx.SeqName  = strings.TrimSpace(line[1:])
			idx.SeqId   += 1
			idx.SeqSize  = 0

		} else {
			if len(line) != 0 {
				idx.SeqSize += int64(len(line))
			}
		}
  	}

	idx.Print()

	fo.Close()
	os.Rename(idxNameTmp, idxName)
}



/*
ReadFastaIndex: reads a fasta index
input         : filename string
returns       : IdxData
*/
func ReadFastaIndex(filename string) *[]*IdxData {
	idxName    := filename+".idx"

	if _, err := os.Stat(idxName); os.IsNotExist(err) {
		log.Fatal("Index file", idxName, "does not exists")
		os.Exit(1)
	}

        fi, err := os.Open(idxName)
	check(err)
	defer fi.Close()

        d, err := fi.Stat()
	check(err)

	log.Println(d)

	data     := []*IdxData{}

	scanner  := bufio.NewScanner(fi)
	scanner.Split(bufio.ScanLines)

	for scanner.Scan() {
		line     := scanner.Text()

		//log.Println(line)

		idx      := new(IdxData)

		idx.Read(line)

		data      = append(data, idx)
	}

	return &data
}

