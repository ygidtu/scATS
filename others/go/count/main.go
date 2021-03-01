package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	pb "github.com/cheggaaa/pb/v3"
	"github.com/pkg/errors"
	"github.com/voxelbrain/goptions"
)

// config is command line parameters
type config struct {
	Input  []string      `goptions:"-i, --input, obligatory, description='The input bam file'"`
	Bed    string        `goptions:"-b, --bed, description='The path to bed'"`
	Output string        `goptions:"-o, --output, description='The output csv'"`
	Thread int           `goptions:"-t, --thread, description='how many thread to read and write bam file'"`
	Debug  bool          `goptions:"-d, --debug, description='Debug level log'"`
	Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

// defaultConfig provide default parameters
func defaultConfig() config {
	return config{Debug: false, Thread: 10}
}

// Data as name says
type Data struct {
	Region string
	Data   map[string]int
}

// Bed as name says
type Bed struct {
	Chrom  string
	Start  int
	End    int
	Name   string
	Score  string
	Strand string
}

// BamReader as name says
type BamReader struct {
	Bam *bam.Reader
	Idx *bam.Index
}

type readPair struct {
	Bed    *Bed
	Reader *BamReader
	Key    string
}

// createIndex
func createIndex(path string) error {
	path, err := exec.LookPath("samtools")
	if err != nil {
		return errors.Wrap(err, "install samtools first")
	}

	file := fmt.Sprintf("%s.bai", path)

	if _, err := os.Stat(file); os.IsNotExist(err) {
		cmd := exec.Command(path, "index", path)
		err = cmd.Run()
		return err
	}
	return nil
}

// IsRead1 is true if this is read1
func IsRead1(r *sam.Record) bool {
	return r.Flags&sam.Read1 != 0
}

// IsRead2 true if this is read2
func IsRead2(r *sam.Record) bool {
	return r.Flags&sam.Read2 != 0
}

// IsReverse true if this is reversed
func IsReverse(r *sam.Record) bool {
	return r.Flags&sam.Reverse != 0
}

// Strand is function that calculate the strand based on read1/read2 and reverse
func Strand(r *sam.Record) string {

	if IsRead1(r) {
		if IsReverse(r) {
			return "-"
		}
		return "+"
	} else if IsRead2(r) {
		if IsReverse(r) {
			return "-"
		}
		return "+"
	}

	return "+"
}

func filterReads(record *sam.Record) bool {
	if record.Flags&sam.Unmapped != 0 {
		// sugar.Debugf("Unmapped: %d", record.Flags)
		return false
	}

	if record.Flags&sam.QCFail != 0 {
		// sugar.Debugf("QCFail: %d", record.Flags)
		return false
	}

	if record.Flags&sam.Secondary != 0 || record.Flags&sam.Supplementary != 0 {
		return false
	}

	if record.Flags&sam.Duplicate != 0 {
		return false
	}

	if _, ok := record.Tag([]byte("SA")); ok {
		return false
	}

	if tag, ok := record.Tag([]byte("NH")); ok {
		return tag.Value().(uint8) <= 1
	}

	return true
}

func create(data string) (*Bed, error) {
	dataArr := strings.Split(strings.TrimSpace(data), "\t")

	start, err := strconv.Atoi(dataArr[1])
	if err != nil {
		return nil, err
	}
	end, err := strconv.Atoi(dataArr[2])
	if err != nil {
		return nil, err
	}

	return &Bed{
		dataArr[0],
		start,
		end,
		dataArr[3],
		dataArr[4],
		dataArr[5],
	}, nil
}

func (b *Bed) String() string {
	return fmt.Sprintf("%s|%s:%d-%d:%s", b.Score, b.Chrom, b.Start, b.End, b.Strand)
}

func write(output string, res chan *Data, reader map[string]string) {

	f, err := os.OpenFile(output, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		sugar.Fatalf("failed to open %s: %s", output, err.Error())
		os.Exit(1)
	}

	w := bufio.NewWriter(f)
	defer w.Flush()
	defer f.Close()

	keys := make([]string, 0, 0)
	for key := range reader {
		keys = append(keys, key)
	}

	w.WriteString(fmt.Sprintf(",%s\n", strings.Join(keys, ",")))

	for {
		row, ok := <-res

		if !ok {
			break
		}

		res := make([]string, 0, 0)
		for _, key := range keys {
			if val, ok := row.Data[key]; ok {
				res = append(res, fmt.Sprintf("%v", val))
			} else {
				res = append(res, "0")
			}
		}

		w.WriteString(fmt.Sprintf("%s,%v\n", row.Region, strings.Join(res, ",")))
		w.Flush()
	}
}

func openBam(path string) (*BamReader, error) {
	createIndex(path)
	ifh, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer ifh.Close()

	reader, err := bam.NewReader(ifh, 1)
	if err != nil {
		return nil, err
	}
	defer reader.Close()

	idxF, err := os.Open(path + ".bai")
	if err != nil {
		return nil, err
	}
	defer idxF.Close()

	idx, err := bam.ReadIndex(idxF)
	if err != nil {
		return nil, err
	}
	idxF.Close()

	return &BamReader{reader, idx}, nil
}

func fetchBams(path string, region *Bed) int {

	count := 0
	createIndex(path)
	ifh, err := os.Open(path)
	if err != nil {
		return count
	}
	defer ifh.Close()

	reader, err := bam.NewReader(ifh, 1)
	if err != nil {
		return count
	}
	defer reader.Close()

	idxF, err := os.Open(path + ".bai")
	if err != nil {
		return count
	}
	defer idxF.Close()

	idx, err := bam.ReadIndex(idxF)
	if err != nil {
		return count
	}
	idxF.Close()

	for _, ref := range reader.Header().Refs() {
		if region.Chrom != "" && ref.Name() == region.Chrom {
			// sugar.Debug(region.String())
			if chunks, err := idx.Chunks(ref, region.Start, region.End); err == nil {
				// sugar.Debug(len(chunks))
				if iter, err := bam.NewIterator(reader, chunks); err == nil {
					for iter.Next() {
						rec := iter.Record()
						if !filterReads(rec) || Strand(rec) != region.Strand {
							continue
						}

						count++
					}
					// iter.Close()
				}
			} else {
				sugar.Debug(err)
			}
			break
		}
	}
	return count
}

func fetchBamsReader(r *BamReader, region *Bed) int {
	count := 0
	for _, ref := range r.Bam.Header().Refs() {
		if region.Chrom != "" && ref.Name() == region.Chrom {
			// sugar.Debug(region.String())
			if chunks, err := r.Idx.Chunks(ref, region.Start, region.End); err == nil {
				// sugar.Debug(len(chunks))
				if iter, err := bam.NewIterator(r.Bam, chunks); err == nil {
					for iter.Next() {
						rec := iter.Record()
						if !filterReads(rec) || Strand(rec) != region.Strand {
							continue
						}

						count++
					}
					// iter.Close()
				}
			} else {
				sugar.Debug(err)
			}
			break
		}
	}
	return count
}

// process is used to read and process bam file
func process(inputBam []string, inputFile string, thread int, output string) error {

	sugar.Info("open bam files")
	bamReader := make(map[string]string)
	// bamReader := make(map[string]*BamReader)
	for _, b := range inputBam {
		key := strings.ReplaceAll(filepath.Base(b), ".Aligned.sortedByCoord.out.bam", "")

		// r, err := openBam(b)
		// if err != nil {
		// 	return err
		// }
		// r.Key = key
		// bamReader[key] = r
		bamReader[key] = b
	}

	var wg sync.WaitGroup
	bc := make(chan *Bed)
	out := make(chan *Data)

	for i := 0; i < thread; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()

			go func() {
				for {
					region, ok := <-bc
					if !ok {
						break
					}

					res := make(map[string]int)
					for key, path := range bamReader {
						res[key] = fetchBams(path, region)
					}

					out <- &Data{region.String(), res}
				}

			}()
		}()
	}

	go write(output, out, bamReader)

	stats, err := os.Stat(inputFile)
	if os.IsNotExist(err) {
		return err
	}

	// start new bar
	bar := pb.Full.Start64(stats.Size())

	f, err := os.Open(inputFile)
	if err != nil {
		return err
	}

	defer f.Close()

	r := bufio.NewReader(bar.NewProxyReader(f))

	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
		}

		bed, err := create(line)
		if err != nil {
			sugar.Error(err)
		}

		bc <- bed
	}
	wg.Wait()

	return nil
}

// Python used 8:44:30
func main() {
	conf := defaultConfig()
	goptions.ParseAndFail(&conf)
	setLogger(conf.Debug, "")

	sugar.Debugf("options: %v", conf)

	if len(conf.Input) == 0 {
		sugar.Fatal("bam file is mandatory. Please, provide one (-i|--input)")
	}

	if conf.Output == "" {
		sugar.Fatal("Please, provide one (-o|--output)")
	}

	outDir, err := filepath.Abs(filepath.Dir(conf.Output))
	if err != nil {
		sugar.Fatal(err)
	}

	if _, err := os.Stat(outDir); os.IsNotExist(err) {
		if err := os.Mkdir(outDir, 0644); err != nil {
			sugar.Fatal(errors.Wrap(err, "failed to create output directory"))
		}
	}

	sugar.Info("filtering bam")
	if err := process(conf.Input, conf.Bed, conf.Thread, conf.Output); err != nil {
		sugar.Fatal(err)
	}
}
