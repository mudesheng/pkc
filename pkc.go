package main

import (
	// "fmt"
	"github.com/jwaldrip/odin/cli"
	"pkc/discretecosinetransform"
	"pkc/kmercount"
)

const Kmerdef = 10

/*
type Args struct {
	cfg    string
	kmer   int
	prefix string
	numCPU int
	cfSize int64
}*/

var app = cli.New("1.0.0", "kmer count from pacbio or other long reads", func(c cli.Command) {})

//var gaargs GAArgs

func init() {
	app.DefineStringFlag("f", "input.fa", "input long reads file")
	app.DefineIntFlag("K", Kmerdef, "kmer length")
	app.DefineStringFlag("p", "pkc", "prefix of the output file")
	//app.DefineIntFlag("t", 4, "number of CPU used")
	// ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	// {
	// 	ccf.DefineInt64Flag("S", 0, "the number of item cuckoofilter set")

	// }
	// cdbg := app.DefineSubCommand("kc", "kmer count", kmercount.KC)
	app.DefineSubCommand("kc", "kmer count", kmercount.KC)
	app.DefineSubCommand("dct", "Discrete Cosine Transform", discretecosinetransform.DCT)
	// {
	// 	cdbg.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	// }

}

func main() {
	app.Start()
}
