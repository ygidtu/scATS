#!/usr/bin/env python3
# -*- coding:utf-8 -*-

u"""
Created & Modified at 2020.12.15

This scripts is used to filter the  ATS results of julia version ATSmix

standard
    1. if there is any first exons in reads region, at least 1
    2. if there is any junctions around this first exons
    3. the counts of exons needs to big enough
"""
import logging
import os
import re
import typing

from glob import glob
from multiprocessing import Pool
from subprocess import check_call

import click
import pysam

from rich.logging import RichHandler
from tqdm import tqdm


logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)

log = logging.getLogger("rich")


class GenomicRegions(object):
    u"""
    utils related to genomic regions
    """

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __str__(self) -> str:
        return "{}:{}-{}".format(self.chrom, self.start, self.end)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))

    def __gt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom > other.chrom

        if self.start != other.start:
            return self.start > other.start

        return self.end > other.end

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end

    def covered(self, other):
        u"""
        whether self is totally coveraged other
        """
        return self.chrom == other.chrom and self.start <= other.start and self.end >= other.end


class GTF(object):
    u"""
    class to record first exons information
    """

    def __init__(self, path: str):
        u"""
        init
        """
        self.first_exons = self.load_gtf(path)

    @staticmethod
    def load_gtf(path: str) -> typing.Dict:
        u"""
        load first exons from gtf
        """

        def decode_attrs(string: str):
            string = string.replace('"', '')
            attrs = string.split("; ")

            res = {}
            for i in attrs:
                i = [x.strip() for x in i.split()]
                if len(i) > 1:
                    res[i[0]] = i[1]

            return res
        
        log.info("reading {}".format(path))
        exons = {}   # transcript id -> list of exons
        with open(path) as r:
            for line in r:
                if line.startswith("#"):
                    continue
                
                lines = line.strip().split("\t")

                if lines[2] == "exon":
                    site = GenomicRegions(lines[0], int(lines[3]), int(lines[4]))
                    attrs = decode_attrs(lines[8])

                    t_id = "{}:{}".format(attrs["transcript_id"], lines[6])
                    temp = exons.get(t_id, [])
                    temp.append(site)
                    exons[t_id] = temp

        log.info("extract first exons from {} transcripts".format(len(exons)))
        first_exons = {}
        for t_id, exons in tqdm(exons.items(), desc="finding first exons"):
            exons = sorted(exons)
            e = exons[0] if t_id.endswith("+") else exons[-1]
            temp = first_exons.get(e.chrom, set())
            temp.add(e)
            first_exons[e.chrom] = temp

        log.info("first exons found".format(len(first_exons)))
        return {x: sorted(y) for x, y in first_exons.items()}

    def get(self, chrom: str) -> typing.List[GenomicRegions]:
        return self.first_exons.get(chrom, [])


def extract_from_cigar_string(record: pysam.AlignedSegment, mode:str=""):
    u"""
    extract junctions from cigar string

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

    IGV only skip S and I

    :param record: SAMSegment
    :return:
    """
    pos = record.reference_start + 1
    pos_list = [pos]

    skipped_code = (1, 2, 4, 5)
    
    modes = {
        "star": (2, 5, 6),
        "igv": (1, 4),
        "none": ()
    }

    skipped_code = modes.get(mode, skipped_code)

    for i, j in record.cigartuples:

        if i not in skipped_code:
            pos += j

        if i == 3:
            pos_list.append(pos - j)
            pos_list.append(pos - 1)
    
    pos_list.append(record.reference_end + 1)

    return pos_list


class Site(object):
    u"""
    Class to restore the ATS site information
    """

    def __init__(self, chrom:str, utr_start: int, utr_end: int, reads_start: int, reads_end: int, bic: float):
        u"""
        init this class
        """
        self.chrom = chrom
        self.utr_start = utr_start
        self.utr_end = utr_end
        self.reads = GenomicRegions(chrom, reads_start, reads_end)
        self.bic = bic

        self.source = ""
        self.covered_exons = []

    @property
    def reads_start(self):
        return self.reads.start

    @property
    def reads_end(self):
        return self.reads.end

    def __str__(self):
        return self.source

    @classmethod
    def create_from_string(cls, string: str):

        line = string.strip().split()

        reads_starts = [int(x) for x in line[1].split(",")]
        reads_ends = [int(x) for x in line[2].split(",")]

        try:
            bic = float(line[-1])
        except ValueError:
            return None

        site = re.split(r"[:-]", line[0])
        
        s = Site(
            chrom = site[0], 
            utr_start = int(site[1]),
            utr_end = int(site[2]),
            reads_start = reads_starts[0],
            reads_end = reads_ends[-1],
            bic = bic
        )
        s.source = string.strip()
        return s

    def check_first_exons(self, gtf:GTF) -> bool:
        u"""
        check whether this region covered any first exons
        """
        for i in gtf.get(self.chrom):
            if self.reads.covered(i):
                self.covered_exons.append(i)
            elif i > self.reads:
                break

        return len(self.covered_exons) > 0

    def check_junction(self, bam: str, min_counts: int = 10) -> bool:
        u"""
        check the junction of covered exons
        """
        with pysam.AlignmentFile(bam) as r:
            for i in set(self.covered_exons):
                counts = {}
                for rec in r.fetch(i.chrom, i.start, i.end):
                    if rec.is_qcfail or rec.is_unmapped:
                        continue

                    if not rec.cigarstring:
                        continue

                    if rec.has_tag("NH"):
                        if rec.get_tag("NH") != 1:
                            continue

                    junc = extract_from_cigar_string(rec)

                    for idx in range(0, len(junc), 2):
                        if junc[idx] <= i.end and junc[idx + 1] >= i.start:
                            return True
                            # temp =  GenomicRegions(rec.reference_name, junc[idx], junc[idx + 1])
                            # temp_counts = counts.get(temp, 0) + 1

                            # if temp_counts > min_counts:
                            #     return True

                            # counts[temp] = temp_counts
        return False                        


def __check__(args):
    site, gtf, bam = args

    if site.check_first_exons(gtf) and site.check_junction(bam):
        return site
    return None


@click.command()
@click.option('-p', '--path', type=click.Path(exists=True), help="Path to ATS output")
@click.option('-b', '--bam', type=click.Path(exists=True), help="Path to bam file")
@click.option('-g', '--gtf', type=click.Path(exists=True), help="Path to reference gtf file")
@click.option('-o', '--output', type=click.Path(), help="Path to output file")
def main(path: str, bam: str, gtf: str, output: str):
    u"""
    main
    """

    gtf = GTF(gtf)

    sites = []

    with open(path) as r:
        for line in tqdm(r):
            site = Site.create_from_string(line)

            if site:
                sites.append([site, gtf, bam])
    
    with Pool(10) as pool:
        res = list(tqdm(pool.imap(__check__, sites), total=len(sites)))
    
    with open(output, "w+") as w:
        for site in res:
            if site:
                w.write(str(site) +  "\n")


# def main(input_dir:str, output:str, expression:str):
#     gtf="/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.gtf"
#     cellranger="/mnt/raid64/Covid19_Gravida/cellranger"

#     for i in glob(os.path.join(input_dir, "*_R1.bed")):
#         key = os.path.basename(i)
#         o = os.path.join(output, key)
#         key = key.replace("_R1.bed", "")

#         print(i, key)
#         try:
#             check_call(f"julia filter.jl -a {i} -o {o} -e {expression}/{key}_gene.counts -g {gtf} -m 20", shell=True)
#         except Exception as err:
#             print(err)


if __name__ == '__main__':
    # main(
    #     "/mnt/raid64/Covid19_Gravida/apamix/R1",
    #     "/mnt/raid64/Covid19_Gravida/apamix/filtered/",
    #     "/mnt/raid64/Covid19_Gravida/apamix/expression/"
    # )
    main()
