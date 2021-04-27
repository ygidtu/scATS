#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Dedifned genomic loci object
"""

from typing import Optional


class Region(object):

    def __init__(self, chromosome: str, start: int, end: int, strand: str):
        assert strand in ("+", "-"), "strand must be + or -"
        assert end >= start, f"end must bigger than start, current -> start: {start}, end: {end}"

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}:{self.strand}"

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

    @classmethod
    def create_from_bed(cls, bed: str):
        u"""
        Create Region obj from bed record

        :params bed: the str of bed record
        """
        bed = bed.strip().split("\t")

        if len(bed) == 3:
            return Region(chromosome=bed[0], start=int(bed[1]), end=int(bed[2]), strand = "+")
        
        if len(bed) == 4:
            return Region(chromosome=bed[0], start=int(bed[1]), end=int(bed[2]), strand = bed[3])
        if len(bed) == 6:
            return Region(chromosome=bed[0], start=int(bed[1]), end=int(bed[2]), strand = bed[5])

        raise TypeError(f"Invalid columns, 3, 4 or 6 is required, current: {bed}")


class GTF(Region):

    __id_label__ = ["gene_id", "gene_name", "transcript_id", "transcript_name", "exon_id"]

    def __init__(
        self, 
        chromosome: str, 
        start: int, end: int, 
        strand: int, 
        source: str,
        attrs: dict = None
    ):
        self.region = Region(chromosome, start, end, strand)
        self.attrs = attrs
        self.source = source

        self.ids = {}
        for i in self.__id_label__:
            self.ids[i] = [self.attrs[i]] if i in self.attrs.keys() else []

    def __hash__(self):
        return hash((self.region, self.source))

    def __gt__(self, other):
        return self.region > other.region

    def __lt__(self, other):
        return self.region < other.region

    def __and__(self, other):
        return self.region & other.region

    def __add__(self, other):
        region = self.region + other.region
        if not region:
            return region

        gtf = GTF(
            chromosome=region.chromosome,
            start=region.start,
            end=region.end,
            strand=region.strand,
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

        :params attrs: gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1";
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
        :params record: 1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; 
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
    def chromosome(self) -> str:
        return self.region.chromosome
    
    @property
    def start(self) -> int:
        return self.region.start

    @property
    def end(self) -> int:
        return self.region.end

    @property
    def strand(self) -> int:
        return self.region.strand

    @property
    def transcript_id(self) -> str:
        return ",".join(sorted(self.ids["transcript_id"]))

    @property
    def gene_id(self) -> str:
        return ",".join(sorted(self.ids["gene_id"]))

    @property
    def exon_id(self) -> str:
        return ",".join(sorted(self.ids["exon_id")])

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

        return f"{self.region.chromosome}\t{self.region.start}\t{self.region.end}\t{r_id}\t{r_name}\t{self.region.strand}"


if __name__ == '__main__':
    print(Region.create_from_bed("chr1\t1\t22\t.\t.\t+\tasfd"))
    pass