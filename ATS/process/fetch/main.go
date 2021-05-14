package main

import (
	"fmt"
	pb "github.com/cheggaaa/pb/v3"
	"github.com/voxelbrain/goptions"
	"os"
	"path/filepath"
	"sync"
)

// config is command line parameters
type config struct {
	Input  []string      `goptions:"-i, --input, obligatory, description='The input bam file'"`
	Output string        `goptions:"-o, --output, description='The output directory'"`
	UTR    string        `goptions:"-u, --utr, description='The path to utr file'"`
	Thread int           `goptions:"-t, --thread, description='how many thread to read and write bam file'"`
	Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

// defaultConfig provide default parameters
func defaultConfig() config {
	return config{Thread: 1}
}

func worker(wg *sync.WaitGroup, in chan *Region, bams []string, output string) {
	defer wg.Done()

	indexes, err := readIndex(bams)

	if err != nil {
		sugar.Fatal(err)
	}

	for {
		utr, ok := <-in

		if !ok {
			break
		}

		var data *Data
		for _, idx := range indexes {
			if data == nil {
				data = fetchBam(utr, idx.I, idx.R)
			} else {
				temp := fetchBam(utr, idx.I, idx.R)
				data.Add(temp)
				temp = nil
			}
		}

		err := data.ToJson(filepath.Join(output, utr.Index))

		if err != nil {
			sugar.Fatal(err)
		}

		data = nil
	}

	for _, i := range indexes {
		i.Close()
	}
}

// process is used to read and process bam file
func main() {
	conf := defaultConfig()
	goptions.ParseAndFail(&conf)
	setLogger(false, "")

	if len(conf.Input) == 0 {
		sugar.Fatal("bam file is mandatory. Please, provide one (-i|--input)")
	}

	if conf.Output == "" {
		sugar.Fatal("Please, provide one (-o|--output)")
	}

	if _, err := os.Stat(conf.Output); os.IsNotExist(err) {
		err = os.MkdirAll(conf.Output, os.ModePerm)
		if err != nil {
			sugar.Fatal(err)
		}
	}

	utr := loadUTR(conf.UTR)
	bar := pb.StartNew(len(utr))

	var wg sync.WaitGroup
	in := make(chan *Region)

	for i := 0; i < conf.Thread; i++ {
		go worker(&wg, in, conf.Input, conf.Output)
		wg.Add(1)
	}

	for i, u := range utr {
		u.Index = fmt.Sprintf("%v.json.gz", i)
		in <- u

		bar.Increment()
	}

	close(in)

	go func() {
		err := RegionToJson(utr, filepath.Join(conf.Output, "index.json.gz"))
		if err != nil {
			sugar.Fatal(err)
		}
	}()

	wg.Wait()
	bar.Finish()
}
