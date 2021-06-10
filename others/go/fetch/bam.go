package main

import (
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
)

type BamEntry struct {
	F *os.File
	R *bam.Reader
	I *bam.Index
}

func (b *BamEntry) Close() {
	b.R.Close()
	b.F.Close()
}

func readIndex(path []string) ([]*BamEntry, error) {
	res := make([]*BamEntry, 0, 0)
	for _, p := range path {
		ifh, err := os.Open(p)
		//Panic if something went wrong:
		if err != nil {
			sugar.Fatalf("faild to open input bam: %v", err)
		}

		bamReader, err := bam.NewReader(ifh, 10)
		if err != nil {
			sugar.Fatalf("failed to create bam reader: %v", err)
		}

		idxF, err := os.Open(p + ".bai")
		if err != nil {
			return nil, err
		}

		idx, err := bam.ReadIndex(idxF)
		if err != nil {
			return nil, err
		}
		idxF.Close()

		res = append(res, &BamEntry{F: ifh, R: bamReader, I: idx})
	}
	return res, nil
}

// fetchIndex fetch bam index
func fetchIndex(r *bam.Reader, idx *bam.Index, chromosome string, start, end int) ([]bgzf.Chunk, error) {

	chunks := make([]bgzf.Chunk, 0, 0)

	for _, ref := range r.Header().Refs() {
		if ref.Name() == chromosome {
			tempChunks, err := idx.Chunks(ref, start, end)
			return tempChunks, err
		}
	}

	return chunks, nil
}

// fetchBam as name says
func fetchBam(utr *Region, idx *bam.Index, b *bam.Reader) *Data {

	res := &Data{R1: make([]*Record, 0, 0), R2: make([]*Record, 0, 0)}

	chunks, err := fetchIndex(b, idx, utr.Chromosome, utr.Start, utr.End)
	if err != nil {
		sugar.Errorf("failed to get %s: %v", utr.Str(), err)
		return res
	}

	iter, err := bam.NewIterator(b, chunks)
	if err != nil {
		sugar.Fatal(err)
	}
	defer iter.Close()

	tempRes := make(map[string][]*Record)
	for iter.Next() {
		rec := NewRecord(iter.Record())

		// Add more check to avoid out range
		if rec.Start > utr.End {
			break
		}

		if rec.End < utr.Start {
			continue
		}

		if rec.Valid() {
			if temp, ok := tempRes[rec.Name]; ok {
				temp = append(temp, rec)
				tempRes[rec.Name] = temp
			} else {
				tempRes[rec.Name] = []*Record{rec}
			}
		}
	}

	// sugar.Infof("%s: %d", utr.Str(), count)

	for _, values := range tempRes {
		if len(values) == 2 {
			for _, val := range values {
				if val.IsRead1 {
					res.R1 = append(res.R1, val)
				} else {
					res.R2 = append(res.R2, val)
				}
			}
		}
	}

	return res
}
