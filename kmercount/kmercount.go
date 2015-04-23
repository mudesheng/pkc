package kmercount

import (
	"fmt"
	"ga/bnt"
	"io"
	"log"
	"math"
	"os"
	"strconv"

	"compress/gzip"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
)

type KVT uint16 // kmer vaule type

func transUint64(ks alphabet.Letters) (kaddr uint64) {
	// fmt.Printf("len(ks):%d\n", len(ks))
	for _, v := range ks {
		kaddr <<= bnt.NumBitsInBase
		kaddr |= uint64(bnt.Base2Bnt[v])
	}

	return kaddr
}

func RCVaule(v uint64, size int) (rv uint64) {
	for i := 0; i < size; i++ {
		b := v & 0x3
		v >>= bnt.NumBitsInBase
		rv <<= bnt.NumBitsInBase
		rv |= b
	}

	return rv
}

func KC(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())
	faName := c.Parent().Flag("f").String()
	kmerlen, err := strconv.Atoi(c.Parent().Flag("K").String())
	fafp, err := os.Open(faName)
	if err != nil {
		log.Fatal("[kc] open fa input file failed")
	}
	defer fafp.Close()
	fagzip, err := gzip.NewReader(fafp)
	hasize := int64(math.Pow(float64(4), float64(kmerlen)))
	ha := make([]KVT, hasize, hasize)
	r := fasta.NewReader(fagzip, linear.NewSeq("", nil, alphabet.DNA))
	for {
		if s, err := r.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("failed read %v:%s\n", fafp, err)
			}
		} else {
			t := s.(*linear.Seq)
			tl := t.Slice().(alphabet.Letters)
			for i := 0; i < tl.Len()-kmerlen+1; i++ {
				v := transUint64(tl[i : i+kmerlen])
				ha[v]++
			}
		}
	}

	// statics hash table
	max := KVT(10000)
	ca := make([]int64, max, max)
	nozerosum := 0
	reverPairNoZeroSum := 0
	for i, v := range ha {
		if v < max-1 {
			ca[v]++
		} else {
			ca[max-1]++
		}
		if v > 0 {
			ri := RCVaule(uint64(i), kmerlen)
			// fmt.Printf("v:%v, rv:%v\n", i, ri)
			if ri > uint64(i) && ha[ri] > 0 {
				reverPairNoZeroSum += 2
			}

			nozerosum++
		}
	}
	// print output
	fmt.Printf("frequency\tcount\n")
	for i, v := range ca {
		fmt.Printf("%d\t%d\n", i, v)
	}
	fmt.Printf("nozerosum: %d\n", nozerosum)
	fmt.Printf("reverPairNoZeroSum: %d\n", reverPairNoZeroSum)
}
