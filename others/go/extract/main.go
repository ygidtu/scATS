package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strconv"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	pb "github.com/cheggaaa/pb/v3"
	"github.com/pkg/errors"
	"github.com/voxelbrain/goptions"
)

// config is command line parameters
type config struct {
	Input  string        `goptions:"-i, --input, description='The input bam file'"`
	Output string        `goptions:"-o, --output, description='The output bam file'"`
	Number int           `goptions:"-n, --number, description='extract one record per ?'"`
	Thread int           `goptions:"-t, --thread, description='how many thread to read and write bam file'"`
	Debug  bool          `goptions:"-d, --debug, description='Debug level log'"`
	Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

// defaultConfig provide default parameters
func defaultConfig() config {
	return config{Number: 10, Debug: false, Thread: 10}
}

// createIndex
func createIndex(conf *config) error {
	path, err := exec.LookPath("samtools")
	if err != nil {
		return errors.Wrap(err, "install samtools first")
	}

	cmd := exec.Command(path, "index", "-@", strconv.Itoa(conf.Thread), conf.Output)
	err = cmd.Run()
	return err
}

// process is used to read and process bam file
func process(conf *config) error {
	outDir, err := filepath.Abs(filepath.Dir(conf.Output))

	if err != nil {
		return err
	}

	if _, err := os.Stat(outDir); os.IsNotExist(err) {
		if err := os.Mkdir(outDir, 0644); err != nil {
			return errors.Wrap(err, "failed to create output directory")
		}
	}

	stat, err := os.Stat(conf.Input)
	if os.IsNotExist(err) {
		return errors.Wrap(err, "input bam not exists")
	}

	ifh, err := os.Open(conf.Input)
	//Panic if something went wrong:
	if err != nil {
		return errors.Wrap(err, "faild to open input bam")
	}
	defer ifh.Close()ƒreplac

	//Create a new BAM reader with maximum
	//concurrency:
	bar := pb.Full.Start64(stat.Size())

	bamReader, err := bam.NewReader(bar.NewProxyReader(ifh), conf.Thread)
	if err != nil {
		return errors.Wrap(err, "failed to create bam reader")
	}
	defer bamReader.Close()

	// open output bam
	ofh, err := os.OpenFile(conf.Output, os.O_TRUNC|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return errors.Wrap(err, "failed to open output bam")
	}
	defer ofh.Close()

	w, err := bam.NewWriter(ofh, bamReader.Header(), conf.Thread)
	if err != nil {
		return errors.Wrap(err, "failed to create writer")
	}
	defer w.Close()

	bc := make(chan *sam.Record)

	go func() {
		for {
			rec, ok := <-bc
			if !ok {
				break
			}

			w.Write(rec)
		}
	}()

	count := 0
	for {
		count++
		rec, err := bamReader.Read()
		if err != nil {
			break
		}
		if count == 10 {
			bc <- rec
			count = 0
		}
	}
	close(bc)

	return nil
}

// Python used 8:44:30
func main() {
	conf := defaultConfig()
	goptions.ParseAndFail(&conf)
	setLogger(conf.Debug, "")

	sugar.Debugf("options: %v", conf)

	if conf.Input == "" {
		sugar.Fatal("bam file is mandatory. Please, provide one (-i|--input)")
	}

	if conf.Output == "" {
		sugar.Fatal("Please, provide one (-o|--output)")
	}

	sugar.Info("filtering bam")
	if err := process(&conf); err != nil {
		sugar.Fatal(err)
	}
	sugar.Info("creating index")
	if err := createIndex(&conf); err != nil {
		sugar.Fatal(err)
	}
}
