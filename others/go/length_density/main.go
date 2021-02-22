package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
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
	Output string        `goptions:"-o, --output, description='The output bam file'"`
	Thread int           `goptions:"-t, --thread, description='how many thread to read and write bam file'"`
	Debug  bool          `goptions:"-d, --debug, description='Debug level log'"`
	Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

// defaultConfig provide default parameters
func defaultConfig() config {
	return config{Debug: false, Thread: 10}
}

// Record is used to log reads status
type Record struct {
	Count int
	Rec   *sam.Record
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

// IsRead1 is true if this is read1
func IsRead1(r *sam.Record) bool {
	return r.Flags&sam.Read1 != 0
}

// IsRead2 true if this is read2
func IsRead2(r *sam.Record) bool {
	return r.Flags&sam.Read2 != 0
}

// Unmap true if reads or mate is unmapped
func Unmap(r *sam.Record) bool {
	return r.Flags&sam.Unmapped != 0 || r.Flags&sam.MateUnmapped != 0
}

func write(output string, res chan string) {

	f, err := os.OpenFile(output, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		sugar.Fatalf("failed to open %s: %s", output, err.Error())
		os.Exit(1)
	}

	w := bufio.NewWriter(f)
	defer w.Flush()
	defer f.Close()

	for {
		row, ok := <-res

		if !ok {
			break
		}

		w.WriteString(fmt.Sprintf("%v\n", row))
	}
}

// process is used to read and process bam file
func process(input []string, thread int, output string) error {

	var wg sync.WaitGroup
	bc := make(chan string)
	out := make(chan string)

	for i := 0; i < thread; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()

			for {
				i, ok := <-bc

				if !ok {
					break
				}

				reads := make(map[string]*Record)

				ifh, err := os.Open(i)
				//Panic if something went wrong:
				if err != nil {
					sugar.Fatal(errors.Wrap(err, "faild to open input bam"))
				}
				defer ifh.Close()

				stat, _ := os.Stat(i)

				//Create a new BAM reader with maximum
				//concurrency:
				bar := pb.Full.Start64(stat.Size())

				bamReader, err := bam.NewReader(bar.NewProxyReader(ifh), thread)
				if err != nil {
					sugar.Fatal(errors.Wrap(err, "failed to create bam reader"))
				}
				defer bamReader.Close()

				/*
					改栈来处理这个会比较好

					先从index里获取生成每个染色体上每个位点的array
					每个单元格内为字典，以reads name存储数据
					1. 没有mate的reads则按照MatePos存储
					2. 能够对应上mate的reads输出，清理旧reads
					3. 必须向前回滚，清理旧有超出范围的reads

					这个可能兼顾了性能，也不会漏

					只要新的rec超过了旧rec的MatePos，
					就意味着没有paired上，删掉即可
				*/
				for {
					rec, err := bamReader.Read()

					if err != nil {
						break
					}

					if !sam.IsValidRecord(rec) {
						continue
					}

					if Unmap(rec) {
						continue
					}

					nh := int(rec.AuxFields.Get(sam.NewTag("NH")).Value().(uint8))
					if err != nil {
						continue
					}

					if nh > 1 {
						continue
					}

					if preRec, ok := reads[rec.Name]; ok {
						r1 := preRec.Rec
						r2 := rec
						if IsRead1(rec) {
							r1 = rec
							r2 = preRec.Rec
						}

						sites := []int{
							r1.Start(), r1.End(),
							r2.Start(), r2.End(),
						}
						sort.Ints(sites)

						out <- fmt.Sprintf(
							"%v,%v,%v,%v,%v,%v,%v,%v,%s",
							r1.Start(), r1.End(),
							r2.Start(), r2.End(),
							r1.Len(), r1.Len(),
							r1.MatePos-r1.End(),
							sites[3]-sites[0],
							strings.Split(filepath.Base(i), ".")[0],
						)

						delete(reads, rec.Name)
						// fmt.Printf("%v\n", len(reads))
					} else {
						reads[rec.Name] = &Record{0, rec}
					}

					del := make([]string, 0)
					for k, v := range reads {
						v.Count++

						if v.Count > 100 {
							if rec.Start() > v.Rec.MatePos {
								del = append(del, k)
							}
						}
					}

					for _, k := range del {
						delete(reads, k)
					}
				}
			}
		}()
	}

	go write(output, out)
	for _, i := range input {
		bc <- i
	}
	wg.Wait()
	close(bc)
	close(out)
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
	if err := process(conf.Input, conf.Thread, conf.Output); err != nil {
		sugar.Fatal(err)
	}
	sugar.Info("creating index")
	if err := createIndex(&conf); err != nil {
		sugar.Fatal(err)
	}
}
