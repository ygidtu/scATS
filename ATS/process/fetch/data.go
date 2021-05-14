package main

import (
	"bufio"
	"compress/gzip"
	jsoniter "github.com/json-iterator/go"
	"os"
)

var json = jsoniter.ConfigCompatibleWithStandardLibrary

// Data reads of single UTR
type Data struct {
	R1 []*Record
	R2 []*Record
}

// Merge UTR
func (d *Data) Add(data *Data) {
	d.R1 = append(d.R1, data.R1...)
	d.R2 = append(d.R2, data.R2...)
}

// ToJson convert struct to json
func (d *Data) ToJson(path string) error {
	data, err := json.MarshalIndent(d, "", "    ")
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
