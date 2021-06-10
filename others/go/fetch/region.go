package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	pb "github.com/cheggaaa/pb/v3"
)

type Region struct {
	Chromosome string `json:"chromosome"`
	Start      int    `json:"start"`
	End        int    `json:"end"`
	Strand     string `json:"strand"`
	Name       string `json:"name"`
	Id         string `json:"id"`
	Index      string `json:"index"`
}

// NewRegion create new Region pointer by bed format string
func NewRegion(data string) *Region {
	lineArr := strings.Split(strings.TrimSpace(data), "\t")

	start, err := strconv.Atoi(lineArr[1])
	if err != nil {
		sugar.Fatalf("failed to convert %s to int, %v", lineArr[1], err)
	}

	end, err := strconv.Atoi(lineArr[2])
	if err != nil {
		sugar.Fatalf("failed to convert %s to int, %v", lineArr[2], err)
	}

	return &Region{
		Chromosome: lineArr[0],
		Start:      start,
		End:        end,
		Strand:     lineArr[5],
		Name:       lineArr[4],
		Id:         lineArr[3],
	}
}

// Str convert Region to bed string
func (r *Region) Str() string {
	res := []string{
		r.Chromosome,
		fmt.Sprintf("%d", r.Start),
		fmt.Sprintf("%d", r.End),
		r.Id,
		r.Name,
		r.Strand,
	}
	return strings.Join(res, "\t")
}

// loadUTR as name says
func loadUTR(path string) []*Region {
	sugar.Infof("Reading %s", path)
	res := make([]*Region, 0, 0)

	stat, err := os.Stat(path)
	if err != nil {
		sugar.Fatalf("%s not exists", path)
	}

	f, err := os.Open(path)
	if err != nil {
		sugar.Fatalf("failed to open %s: %v", path, err)
	}

	bar := pb.Full.Start64(stat.Size())

	r := bufio.NewReader(bar.NewProxyReader(f))

	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			sugar.Fatalf("error while reading %s: %v", path, err)
		}

		res = append(res, NewRegion(line))
	}
	bar.Finish()
	return res
}

// RegionToJson
func RegionToJson(region []*Region, path string) error {
	data, err := json.MarshalIndent(region, "", "    ")
	if err != nil {
		return err
	}

	f, err := os.OpenFile(path, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		return err
	}
	defer f.Close()

	gf := gzip.NewWriter(f)
	defer gf.Close()

	fw := bufio.NewWriter(gf)
	defer fw.Flush()

	fw.Write(data)
	return nil
}
