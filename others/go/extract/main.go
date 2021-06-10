package main

import (
	"bufio"
	"compress/gzip"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	pb "github.com/cheggaaa/pb/v3"
	"github.com/pkg/errors"
	"github.com/voxelbrain/goptions"
	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
	"gopkg.in/natefinch/lumberjack.v2"
)

// config is command line parameters
type config struct {
	Input  string        `goptions:"-i, --input, description='The input bam file'"`
	Output string        `goptions:"-o, --output, description='The output bam file'"`
	Target string        `goptions:"-s, --select, description='Selected reads'"`
	Thread int           `goptions:"-t, --thread, description='how many thread to read and write bam file'"`
	Debug  bool          `goptions:"-d, --debug, description='Debug level log'"`
	Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

// defaultConfig provide default parameters
func defaultConfig() config {
	return config{Debug: false, Thread: 10}
}

var logger *zap.Logger
var sugar *zap.SugaredLogger

// NewEncoderConfig as name says
func NewEncoderConfig() zapcore.EncoderConfig {
	return zapcore.EncoderConfig{
		// Keys can be anything except the empty string.
		TimeKey:        "T",
		LevelKey:       "L",
		NameKey:        "N",
		CallerKey:      "C",
		MessageKey:     "M",
		StacktraceKey:  "S",
		LineEnding:     zapcore.DefaultLineEnding,
		EncodeLevel:    zapcore.CapitalColorLevelEncoder,
		EncodeTime:     TimeEncoder,
		EncodeDuration: zapcore.StringDurationEncoder,
		EncodeCaller:   zapcore.ShortCallerEncoder,
	}
}

// TimeEncoder as name says
func TimeEncoder(t time.Time, enc zapcore.PrimitiveArrayEncoder) {
	enc.AppendString(t.Format("2006-01-02 15:04:05.000"))
}

func setLogger(debug bool, log string) {

	level := zap.InfoLevel
	encoder := NewEncoderConfig()
	if debug {
		level = zap.DebugLevel
		encoder = zap.NewDevelopmentEncoderConfig()
	}

	core := zapcore.NewCore(
		zapcore.NewConsoleEncoder(encoder),
		zapcore.AddSync(os.Stderr),
		level,
	)

	if log != "" {
		w := zapcore.AddSync(&lumberjack.Logger{
			Filename:   log,
			MaxSize:    500, // megabytes
			MaxBackups: 3,
			MaxAge:     28, // days
		})

		core = zapcore.NewCore(
			zapcore.NewJSONEncoder(encoder),
			zapcore.AddSync(w),
			level,
		)
	}

	logger = zap.New(core, zap.AddCaller())
	defer logger.Sync()
	sugar = logger.Sugar()
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

// load tagets
func load(path string) (map[string]int, error) {
	required := map[string]int{}

	stat, err := os.Stat(path)
	if err != nil {
		return required, err
	}

	f, err := os.Open(path)
	if err != nil {
		return required, err
	}
	defer f.Close()

	bar := pb.Full.Start64(stat.Size())

	gr, err := gzip.NewReader(bar.NewProxyReader(f))
	if err != nil {
		return required, err
	}
	defer gr.Close()

	r := bufio.NewReader(gr)

	for {
		line, err := r.ReadString('\n')
		if err != nil {
			break
		}
		line = strings.TrimSpace(line)

		if line == "" {
			continue
		}

		required[line] = 0
	}

	bar.Finish()

	return required, nil
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

	required, err := load(conf.Target)
	if err != nil {
		return err
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
	defer ifh.Close()

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

	for {
		rec, err := bamReader.Read()
		if err != nil {
			break
		}
		if _, ok := required[rec.Name]; ok {
			bc <- rec
		}
	}
	close(bc)
	bar.Finish()
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
