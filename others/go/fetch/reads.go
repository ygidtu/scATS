package main

import (
	"fmt"

	"github.com/biogo/hts/sam"
)

// Record is a wrap of
type Record struct {
	Name    string      `json:"name"`
	Chrom   string      `json:"chromosome"`
	Start   int         `json:"start"`
	End     int         `json:"end"`
	Strand  string      `json:"strand"`
	IsRead1 bool        `json:"is_read1"`
	Cigar   string      `json:"cigar"`
	Flags   sam.Flags   `json:"-"`
	Record  *sam.Record `json:"-"`
}

// NewRecord is function that create a pointer to record
func NewRecord(r *sam.Record) *Record {

	queryLength := 0
	for _, c := range r.Cigar {
		if c.Type() == sam.CigarMatch {
			queryLength += c.Len()
		}
	}

	rec := &Record{
		Name:    r.Name,
		Cigar:   r.Cigar.String(),
		Chrom:   r.Ref.Name(),
		Start:   r.Start(),
		End:     r.End(),
		IsRead1: IsRead1(r),
		Strand:  Strand(r),
		Flags:   r.Flags,
		Record:  r,
	}

	return rec
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
			return "+"
		}
		return "-"
	}
	if IsReverse(r) {
		return "-"
	}
	return "+"
}

func (r *Record) Valid() bool {
	if r.Flags&sam.Unmapped != 0 {
		// sugar.Debugf("Unmapped: %d", r.Flags)
		return false
	}

	if r.Flags&sam.QCFail != 0 {
		// sugar.Debugf("QCFail: %d", r.Flags)
		return false
	}

	if r.Flags&sam.Secondary != 0 || r.Flags&sam.Supplementary != 0 {
		return false
	}

	if r.Flags&sam.Duplicate != 0 {
		return false
	}

	//sugar.Infof("seq length: %d, min read length: %d", r.Seq.Length, conf.MinReadLength)

	// if r.Flags&sam.Paired != 0 &&
	// 	!(r.Flags == sam.Paired|sam.ProperPair|sam.MateReverse|sam.Read1 || // 163
	// 		r.Flags == sam.Paired|sam.ProperPair|sam.MateReverse|sam.Read2 || // 83
	// 		r.Flags == sam.Paired|sam.ProperPair|sam.Reverse|sam.Read2 || // 147
	// 		r.Flags == sam.Paired|sam.ProperPair|sam.Reverse|sam.Read1) { // 99
	// 	return false
	// }

	// if _, ok := r.Record.Tag([]byte("SA")); ok {
	// 	return false
	// }

	if tag, ok := r.Record.Tag([]byte("NH")); ok {
		return tag.Value().(uint8) <= 1
	}

	return true
}

func (r *Record) Str() string {
	return fmt.Sprintf("%s\t%d\t%d\t%s\t%s", r.Chrom, r.Start, r.End, r.Strand, r.Cigar)
}
