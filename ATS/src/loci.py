#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Dedifned genomic loci object
"""

from typing import List, Optional

import pysam


class Region(object):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the genomic regions
    """

    __slots__ = ["chromosome", "start", "end", "strand"]

    def __init__(self, chromosome: str, start: int, end: int, strand: str):
        start, end = int(start), int(end)
        assert strand in ("+", "-"), "strand must be + or -"
        assert end >= start, "end must bigger than start, current -> start: {}, end: {}".format(start, end)

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __str__(self) -> str:
        return "{}:{}-{}:{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand
        )

    def __hash__(self):
        return hash(self.__str__())

    def __gt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome
        
        if self.start != self.start:
            return self.start > other.start

        if self.end != self.end:
            return self.end > other.end
        
        return self.strand > other.strand

    def __lt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome
        
        if self.start != self.start:
            return self.start < other.start

        if self.end != self.end:
            return self.end < other.end
        
        return self.strand < other.strand

    def __and__(self, other) -> bool:
        u"""
        override & for overlap checking
        """
        if self.chromosome == other.chromosome and self.strand == other.strand:
            return self.start < other.end and self.end > other.start
        return False

    def __add__(self, other) -> Optional:
        u"""
        override + for merging regions

        :return Region or None: None means there is no any overlap, Region is the merged region
        """
        if self & other:
            return Region(
                chromosome = self.chromosome,
                strand = self.strand,
                start = min(self.start, other.start),
                end = max(self.end, other.end)
            )

    def is_cover(self, other, tolerance: int = 3) -> bool:
        u"""
        check whether this region is covered another region
        :params tolerance: how many mismatch is allowed
        """
        if self.chromosome != other.chromosome or self.strand != other.strand:
            return False

        return self.start - tolerance <= other.start <= other.end <= self.end + tolerance
        

class GenomicLoci(Region):

    def __init__(
        self, chromosome: str, 
        start: int, end: int, 
        strand: str, gtf_line: str
    ):
        super(GenomicLoci, self).__init__(chromosome, start, end, strand)
        self.gtf_line = gtf_line

class GTF(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the records in BED format
    """
    __slots__ = ["attrs", "source", "ids"]
    __id_label__ = ["gene_id", "gene_name", "transcript_id", "transcript_name", "exon_id"]

    def __init__(
        self, 
        chromosome: str, 
        start: int, end: int, 
        strand: str, 
        source: str,
        attrs: dict = None
    ):
        super(GTF, self).__init__(chromosome, start, end, strand)
        self.attrs = attrs
        self.source = source

        self.ids = {}
        for i in self.__id_label__:
            self.ids[i] = [self.attrs[i]] if i in self.attrs.keys() else []

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end, self.strand, self.source))

    def __add__(self, other):
        if self & other:
            gtf = GTF(
                chromosome=self.chromosome,
                start=min(self.start, other.start),
                end=max(self.end, other.end),
                strand=self.strand,
                attrs = self.attrs,
                source = self.source
            )

            for i in self.__id_label__:
                self.ids[i] += other.ids[i]
                self.ids[i] = list(set(self.ids[i]))
            return gtf

    @classmethod
    def decode_attrs(cls, attrs: str) -> dict:
        u"""
        Decode attrs from gtf files

        :param attrs: gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1";
        :return: {gene_id: "ENSG00000223972"; gene_version: "5", gene_name: "DDX11L1"}
        """
        res = {}

        for i in attrs.split(";"):
            i = i.strip().replace('"', "").replace(";", "")
            values = i.split(" ")

            if len(values) > 1:
                res[values[0].strip()] = values[1].strip()

        return res

    @classmethod
    def create(cls, record: str):
        u"""
        Create GTF object from gtf record
        :param record: 1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; 
        """

        record = record.strip().split("\t")

        source = record[2]
        if "RNA" in record[2]:
            source = "transcript"

        return cls(
            chromosome=record[0], 
            start=int(record[3]), 
            end=int(record[4]),
            strand=record[6],
            attrs = cls.decode_attrs(record[8]) if len(record) >= 8 else {},
            source=source
        )

    @property
    def transcript_id(self) -> str:
        return ",".join(sorted(self.ids["transcript_id"]))

    @property
    def gene_id(self) -> str:
        return ",".join(sorted(self.ids["gene_id"]))

    @property
    def exon_id(self) -> str:
        return ",".join(sorted(self.ids["exon_id"]))

    @property
    def transcript_name(self) -> str:
        return ",".join(sorted(self.ids["transcript_name"]))

    @property
    def gene_name(self) -> str:
        return ",".join(sorted(self.ids["gene_name"]))

    @property
    def bed(self) -> str:
        r_id = ",".join([x for x in [self.gene_id, self.transcript_id, self.exon_id] if x])
        r_name = ",".join([x for x in [self.gene_name, self.transcript_name] if x])

        return "{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chromosome,
            self.start,
            self.end,
            r_id,
            r_name,
            self.strand
        )


class BED(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the records in BED format
    """

    __slots__ = ("name", "id")

    def __init__(
        self, 
        chromosome: str, 
        start: int, 
        end: int, 
        strand: str, 
        name: str, 
        record_id: str
    ):
        super(BED, self).__init__(chromosome, start, end, strand)
        self.name = name
        self.id = record_id

    def __str__(self) -> str:
        return "{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.id,
            self.name,
            self.strand
        )

    def __len__(self) -> int:
        return self.end - self.start

    @classmethod
    def create(cls, bed: str):
        u"""
        Create Region obj from bed record

        :param bed: the str of bed record
        """
        bed = bed.strip().split("\t")

        if len(bed) not in [3, 4, 6]:
            raise TypeError("Invalid columns, 3, 4 or 6 is required, current: {}".format(bed))

        strand = "+"
        name = "."
        record_id = "."
        
        if len(bed) == 4:
            strand = bed[3]

        if len(bed) == 6:
            record_id = bed[3]
            name = bed[4]
            strand = bed[5]

        return cls(
            chromosome=bed[0],
            start=int(bed[1]),
            end=int(bed[2]),
            record_id=record_id,
            name=name,
            strand=strand
        )

    def to_str(self) -> str:
        return "{}:{}-{}:{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand
        )


class Reads(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the reads and it's junctions
    """

    __slots__ = ["name", "is_read1", "cigar", "kwargs"]

    def __init__(
        self, 
        ref: str, 
        start: int, 
        end: int, 
        name: str,
        strand: str,
        is_read1: bool, 
        cigar = None,
        **kwargs
    ):
        u"""
        init the reads object

        :param ref: reference id, chr1, chr2 etc.
        :param start: left-most start
        :param end: right-most site
        :param name: reads name
        :param strand: strand
        :param is_read1: as name says
        :param cigar: cigar tuples, detailes see pysam.AlignedSegment
        :param paired: SE -> None, PE -> pointer to paired reads 
        """

        super(Reads, self).__init__(ref, start, end, strand)

        self.name = name
        self.is_read1 = is_read1
        self.cigar = cigar
        self.kwargs = kwargs

    @classmethod
    def create(cls, record: pysam.AlignedSegment, skip_qc: bool = False):
        u"""
        Create Reads obj from pysam.AlignedSegment

        :param record: as type
        :param skip_qc: skip QC filtering
        :return if qc is enabled and the records failed qc, then return None
        """
        if not skip_qc:
            if record.is_unmapped or record.is_qcfail or not record.is_proper_pair:
                return None

            try:
                if record.get_tag("NH") > 1:
                    return None
            except ValueError:
                pass

        return cls(
            ref=record.reference_name,
            start=record.reference_start + 1,
            end = record.reference_end + 1,
            strand = cls.__determine_strand__(record),
            name = record.query_name,
            is_read1 = record.is_read1,
            cigar = record.cigarstring
        )

    @staticmethod
    def __determine_strand__(record: pysam.AlignedSegment) -> str:
        u"""
        determine the strand from pysam.AlignedSegment
        :param record: pysam.AlignedSegment
        """

        if record.is_read1:
            return "-" if record.is_reverse else "+"
        
        if record.is_read2:
            return "+" if record.is_reverse else "-"

        if not record.is_paired:
            return "-" if record.is_reverse else "+"
        return "*"

    @property
    def exon_parts(self) -> List[int]:
        u"""
        generate exon parts from cigar

        M	BAM_CMATCH	0
        I	BAM_CINS	1
        D	BAM_CDEL	2
        N	BAM_CREF_SKIP	3
        S	BAM_CSOFT_CLIP	4
        H	BAM_CHARD_CLIP	5
        P	BAM_CPAD	6
        =	BAM_CEQUAL	7
        X	BAM_CDIFF	8
        B	BAM_CBACK	9

        :return list of int
        """
        pos = self.start
        pos_list = [pos]

        skipped_code = ("I", "D", "S", "H")
        
        j = ""
        for i in self.cigar:
            try:
                int(i)
                j += i
            except ValueError:
                j = int(j)
                if i not in skipped_code:
                    pos += j

                if i == "N":
                    pos_list.append(pos - j)
                    pos_list.append(pos - 1)

                j = ""
        
        pos_list.append(self.end)

        return pos_list


if __name__ == '__main__':
    print(Region.create_from_bed("chr1\t1\t22\t.\t.\t+\tasfd"))
    pass
