package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/rs/zerolog"
	"github.com/rs/zerolog/log"
	"github.com/schollz/progressbar/v3"
	"github.com/voxelbrain/goptions"
)

// SetLogger as name says
func SetLogger() zerolog.Logger {
	zerolog.SetGlobalLevel(zerolog.InfoLevel)

	output := zerolog.ConsoleWriter{
		Out:        os.Stderr,
		TimeFormat: time.RFC3339Nano,
		NoColor:    false,
	}
	output.FormatLevel = func(i interface{}) string {
		return strings.ToUpper(fmt.Sprintf("| %-6s |", i))
	}

	log.Logger = log.Output(output)

	return zerolog.New(output).With().Timestamp().Logger()
}

// zeros generate slice full of "0" by specific size
func zeros(size int) []string {
	res := make([]string, size)

	for i := 0; i < size; i++ {
		res[i] = "0"
	}
	return res
}

// format is used to format string
func format(inChan chan [][]string, outChan chan string, barcodes map[string]int, wg *sync.WaitGroup) {
	defer wg.Done()

	for {
		rows, ok := <-inChan

		if !ok {
			break
		}

		res := zeros(len(barcodes))
		rowId := ""
		for _, cols := range rows {
			if rowId == "" {
				rowId = cols[0]
			}

			if idx, ok := barcodes[cols[1]]; ok {
				res[idx] = cols[2]
			}
		}

		outChan <- fmt.Sprintf("%s,%s\n", rowId, strings.Join(res, ","))
	}
}

// write is save format string to file
func write(output string, barcodes map[string]int, outChan chan string, wg *sync.WaitGroup) {
	f, err := os.OpenFile(output, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Err(err).Msgf("failed to open %s", output)
	}

	gw := gzip.NewWriter(f)

	w := bufio.NewWriter(gw)

	// write header
	colnames := make([]string, len(barcodes))
	for i, j := range barcodes {
		colnames[j] = i
	}

	w.WriteString(",")
	w.WriteString(strings.Join(colnames, ","))
	w.WriteString("\n")

	// save content
	for {
		line, ok := <-outChan

		if !ok {
			break
		}

		w.WriteString(line)
	}

	w.Flush()
	gw.Flush()
	gw.Close()

	f.Close()
	wg.Done()
}

// openFileToRead as name says
func openFileToRead(input string, barText string) (*progressbar.ProgressBar, *bufio.Reader, *os.File) {
	// Try to format file
	stats, err := os.Stat(input)
	if os.IsNotExist(err) {
		log.Err(err).Msgf("%s not exists", input)
	}

	f, err := os.Open(input)
	if err != nil {
		log.Err(err).Msgf("failed to open %s", input)
	}

	bar := progressbar.DefaultBytes(
		stats.Size(),
		barText,
	)

	// wrap a gzip reader, if input file in gzipped
	reader := bufio.NewReader(f)
	if strings.HasSuffix(input, ".gz") {
		gr, err := gzip.NewReader(f)

		if err != nil {
			log.Err(err).Msgf("failed to open %s", input)
		}

		reader = bufio.NewReader(gr)
	}
	return bar, reader, f
}

// updateBar as name says
func updateBar(bar *progressbar.ProgressBar, already int64, f *os.File) int64 {
	curOffset, _ := f.Seek(0, os.SEEK_CUR)
	bar.Add(int(curOffset - already))
	return curOffset
}

// collectBarcode used to collect all barcodes from file
func collectBarcode(path string) (map[string]int, error) {
	barcodes := make(map[string]int)

	bar, reader, f := openFileToRead(path, "Barcoding")

	haveRead := int64(0) // use this to log the reading process
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return barcodes, errors.Wrapf(err, "failed to read from %s", path)
			}
		}

		haveRead = updateBar(bar, haveRead, f)

		fields := strings.Split(strings.TrimSpace(line), "\t")

		if _, ok := barcodes[fields[1]]; !ok {
			barcodes[fields[1]] = len(barcodes)
		}
	}

	bar.Finish()

	return barcodes, nil
}

// loadReference as name says
func loadReference(path string) (map[string]int, error) {
	barcodes := make(map[string]int)

	bar, reader, f := openFileToRead(path, "Loading ref")
	defer f.Close()

	haveRead := int64(0) // use this to log the reading process
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return barcodes, errors.Wrapf(err, "failed to read from %s", path)
			}
		}

		haveRead = updateBar(bar, haveRead, f)

		barcodes[strings.TrimSpace(line)] = len(barcodes)
	}

	bar.Finish()

	return barcodes, nil
}

// processCounts as name says
func processCounts(
	input, output string,
	thread, min int,
	barcodes map[string]int,
) map[string]int {
	inChan := make(chan [][]string)
	outChan := make(chan string)

	var wg sync.WaitGroup
	wg.Add(1)
	go write(output, barcodes, outChan, &wg)

	var wg1 sync.WaitGroup
	for i := 0; i < thread; i++ {
		wg1.Add(1)
		go format(inChan, outChan, barcodes, &wg1)
	}

	bar, reader, f := openFileToRead(input, "Formatting count")
	defer f.Close()
	lastRowId, sum := "", 0
	vals, validRowId := make([][]string, 0), make(map[string]int)

	haveRead := int64(0) // use this to log the reading process
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Err(err).Msgf("failed to read from %s", input)
			}
		}

		haveRead = updateBar(bar, haveRead, f)

		fields := strings.Split(strings.TrimSpace(line), "\t")

		if lastRowId != "" && lastRowId != fields[0] {
			if sum > min {
				inChan <- vals

				validRowId[lastRowId] = 0
			}

			vals = [][]string{}
			sum = 0
		}
		val, err := strconv.Atoi(fields[2])
		if err != nil {
			log.Err(err).Msgf("failed to convert %v", fields)
		}
		sum += val
		lastRowId = fields[0]
		vals = append(vals, fields)
	}

	if len(vals) > 0 {
		inChan <- vals
	}

	bar.Finish()

	close(inChan)
	wg1.Wait()
	close(outChan)
	wg.Wait()

	return validRowId
}

func processPSI(
	input, output string,
	thread int,
	barcodes, valid map[string]int,
) {
	inChan := make(chan [][]string)
	outChan := make(chan string)

	var wg sync.WaitGroup
	wg.Add(1)
	go write(output, barcodes, outChan, &wg)

	var wg1 sync.WaitGroup
	for i := 0; i < thread; i++ {
		wg1.Add(1)
		go format(inChan, outChan, barcodes, &wg1)
	}

	bar, reader, f := openFileToRead(input, "Formatting psi")
	defer f.Close()

	lastRowId := ""
	vals := make([][]string, 0)

	haveRead := int64(0) // use this to log the reading process
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Err(err).Msgf("failed to read from %s", input)
			}
		}

		haveRead = updateBar(bar, haveRead, f)

		fields := strings.Split(strings.TrimSpace(line), "\t")

		if lastRowId != "" && lastRowId != fields[0] {
			if _, ok := valid[lastRowId]; ok {
				inChan <- vals
			}

			vals = [][]string{}
		}
		lastRowId = fields[0]
		vals = append(vals, fields)
	}

	if _, ok := valid[lastRowId]; ok && len(vals) > 0 {
		inChan <- vals
	}

	bar.Finish()

	close(inChan)
	wg1.Wait()
	close(outChan)
	wg.Wait()
}

func main() {
	options := struct {
		Input      string        `goptions:"-i, --input, description='The path to count or psi mtx'"`
		Psi        string        `goptions:"-p, --psi, description='The path to psi mtx'"`
		Output     string        `goptions:"-o, --output, description='The path to output file'"`
		Thread     int           `goptions:"-t, --thread, description='How many threads to use'"`
		Min        int           `goptions:"-m, --min, description='The minimum counts of ATS to output'"`
		Referernce string        `goptions:"-r, --reference, description='List of cell barcodes to kept'"`
		Help       goptions.Help `goptions:"-h, --help, description='Show this help'"`
	}{
		Thread: 4,
	}
	goptions.ParseAndFail(&options)

	if len(os.Args) <= 1 {
		goptions.PrintHelp()
		os.Exit(0)
	}

	SetLogger()

	barcodes, err := loadReference(options.Referernce)
	if err != nil {
		log.Error().Err(err).Msg("")
	}

	if len(barcodes) == 0 {
		barcodes, err = collectBarcode(options.Input)
		if err != nil {
			log.Error().Err(err).Msg("")
		}
	}

	valid := processCounts(
		options.Input,
		fmt.Sprintf("%s.count.gz", options.Output),
		options.Thread,
		options.Min,
		barcodes,
	)

	processPSI(
		options.Psi,
		fmt.Sprintf("%s.psi.gz", options.Output),
		options.Thread,
		barcodes,
		valid,
	)
}
