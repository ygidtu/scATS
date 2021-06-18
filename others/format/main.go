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

// collectBarcode used to collect all barcodes from file
func collectBarcode(path string) (map[string]int, error) {
	barcodes := make(map[string]int)

	stats, err := os.Stat(path)
	if os.IsNotExist(err) {
		return barcodes, errors.New(path + " not exists")
	}

	f, err := os.Open(path)
	if err != nil {
		return barcodes, errors.Wrapf(err, "failed to open %s", path)
	}
	defer f.Close()

	bar := progressbar.DefaultBytes(
		stats.Size(),
		"barcoding",
	)

	reader := bufio.NewReader(f)

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return barcodes, errors.Wrapf(err, "failed to read from %s", path)
			}
		}

		bar.Add(len([]byte(line)))

		fields := strings.Split(strings.TrimSpace(line), "\t")

		if _, ok := barcodes[fields[1]]; !ok {
			barcodes[fields[1]] = len(barcodes)
		}
	}

	bar.Finish()

	return barcodes, nil
}

// zeros generate slice full of "0" by specific size
func zeros(size int) []string {
	res := make([]string, size)

	for i := 0; i < size; i++ {
		res[i] = "0"
	}
	return res
}

// write is save format string to file
func write(output string, barcodes map[string]int, outChan chan string) {
	f, err := os.OpenFile(output, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Err(err).Msgf("failed to open %s", output)
	}
	defer f.Close()

	gw := gzip.NewWriter(f)
	defer gw.Flush()
	defer gw.Close()

	w := bufio.NewWriter(gw)
	defer w.Flush()

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

			res[barcodes[cols[1]]] = cols[2]
		}

		outChan <- fmt.Sprintf("%s,%s\n", rowId, strings.Join(res, ","))
	}
}

func processCounts(input, output string, thread, min int, barcodes map[string]int) map[string]int {
	inChan := make(chan [][]string)
	outChan := make(chan string)

	go write(output, barcodes, outChan)

	var wg sync.WaitGroup
	for i := 0; i < thread; i++ {
		wg.Add(1)
		go format(inChan, outChan, barcodes, &wg)
	}

	// Try to format file
	stats, err := os.Stat(input)
	if os.IsNotExist(err) {
		log.Err(err).Msgf("%s not exists", input)
	}

	f, err := os.Open(input)
	if err != nil {
		log.Err(err).Msgf("failed to open %s", input)
	}
	defer f.Close()

	bar := progressbar.DefaultBytes(
		stats.Size(),
		"formatting count",
	)

	reader := bufio.NewReader(f)

	lastRowId := ""
	vals := make([][]string, 0)
	sum := 0
	validRowId := make(map[string]int)
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Err(err).Msgf("failed to read from %s", input)
			}
		}

		bar.Add(len([]byte(line)))

		fields := strings.Split(strings.TrimSpace(line), "\t")

		val, err := strconv.Atoi(fields[2])

		if err != nil {
			log.Err(err).Msgf("failed to convert %v", fields)
		}
		sum += val

		if lastRowId != "" && lastRowId != fields[0] {
			if sum > min {
				inChan <- vals

				validRowId[lastRowId] = 0
			}

			vals = [][]string{}
			sum = 0
		}
		lastRowId = fields[0]
		vals = append(vals, fields)
	}

	if len(vals) > 0 {
		inChan <- vals
	}

	bar.Finish()

	close(inChan)
	wg.Wait()
	close(outChan)

	return validRowId
}

func processPSI(input, output string, thread int, barcodes, valid map[string]int) {
	inChan := make(chan [][]string)
	outChan := make(chan string)

	go write(output, barcodes, outChan)

	var wg sync.WaitGroup
	for i := 0; i < thread; i++ {
		wg.Add(1)
		go format(inChan, outChan, barcodes, &wg)
	}

	// Try to format file
	stats, err := os.Stat(input)
	if os.IsNotExist(err) {
		log.Err(err).Msgf("%s not exists", input)
	}

	f, err := os.Open(input)
	if err != nil {
		log.Err(err).Msgf("failed to open %s", input)
	}
	defer f.Close()

	bar := progressbar.DefaultBytes(
		stats.Size(),
		"formatting psi",
	)

	reader := bufio.NewReader(f)

	lastRowId := ""
	vals := make([][]string, 0)

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Err(err).Msgf("failed to read from %s", input)
			}
		}

		bar.Add(len([]byte(line)))

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
	wg.Wait()
	close(outChan)
}

func main() {
	options := struct {
		Input  string        `goptions:"-i, --input, description='The path to count or psi mtx'"`
		Psi    string        `goptions:"-p, --psi, description='The path to psi mtx'"`
		Output string        `goptions:"-o, --output, description='The path to output file'"`
		Thread int           `goptions:"-t, --thread, description='How many threads to use'"`
		Min    int           `goptions:"-m, --min, description='The minimum counts of ATS to output'"`
		Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
	}{
		Thread: 4,
	}
	goptions.ParseAndFail(&options)

	if len(os.Args) <= 1 {
		goptions.PrintHelp()
		os.Exit(0)
	}

	SetLogger()

	barcodes, err := collectBarcode(options.Input)
	if err != nil {
		log.Error().Err(err).Msg("")
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
