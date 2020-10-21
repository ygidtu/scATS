#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.08.21 By Zhang Yiming

This script is used to extract 
"""

import gzip
import os
import typing
from multiprocessing import Pool, cpu_count

import click
from tqdm import tqdm

import pysam
from loguru import logger
from utils.bed6 import Bed6

"""
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
"""

CIGAR = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]


class Info(object):

    def __init__(self, relative_pos: int, start:  int, end: int, length: int, paired: str, unknown: int = 0):
        u"""
        init
        """
        self.relative_pos = relative_pos
        self.length = length
        self.unknown = unknown
        self.paired = paired
        self.start = start
        self.end = end

    def __str__(self) -> str:
        u"""
        convert this object to string
        """
        return "{}\t{}\t{}\t{}\t{}\t{}".format(
            self.relative_pos, self.length,
            self.unknown, self.paired,
            "NA", "NA"
        )


class RelativeInfo(object):
    u"""
    Class handle the information of relative to utr from read
    """

    def __init__(self, utr: Bed6):
        u"""
        init
        """
        self.utr = utr
        self.__reads__ = {}

    @classmethod
    def header(cls)  -> str:
        return "\t".join([
            "BC", "Chrom", "UTRStart", "UTREnd",
            "ReadStart", "ReadEnd", "ReadLabel",
            "StartLocInUTR", "LenInUTR",
            "LenInPA", "LenPA", "PASite",
        ])

    def __str__(self) -> str:
        u"""
        convert this object to string
        """
        res = []

        for bc, reads in self.__reads__.items():
            for i in reads:
                res.append("\t".join([
                    bc, self.utr.chrom, 
                    str(self.utr.start), str(self.utr.end),
                    str(i.start), str(i.end), str(i.paired),
                    str(i.relative_pos), str(min(self.utr.end, i.end) - max(self.utr.start, i.start)),
                    "NA", "NA", "NA"
                ]))
        return "\n".join(res)
 
    def add(self, reads: pysam.AlignedSegment):
        u"""
        add reads to this utr
        """

        # if this reads is paired-end, then determine the strand by read1
        # if reads.is_paired:
        #     if reads.is_read1:
        #         strand = reads.is_reverse
        #     else:
        #         strand = reads.mate_is_reverse
        # else:
        #     strand = reads.is_reverse

        # strand = "-" if strand else "+"

        utr_site = self.utr.end # if strand == "+" else self.utr.start
        start_site = reads.reference_start

        if reads.reference_start <= utr_site <= reads.reference_end:
            bc = reads.get_tag("CB")

            if bc not in self.__reads__.keys():
                self.__reads__[bc] = []

            self.__reads__[bc].append(Info(
                relative_pos=start_site - utr_site,
                start=reads.reference_start,
                end=reads.reference_end,
                length=reads.reference_length,
                unknown=-1,
                paired="{}_{}".format(
                    "pair" if reads.is_paired else "unpair",
                    "unmap" if reads.is_unmapped else "map"
                )
            ))


def filter(read: pysam.AlignedSegment) -> bool:
    u"""
    check whether the read is qualified
    """
    if read.has_tag("NH") and read.get_tag("NH") > 1:
        return False

    if read.is_duplicate or read.is_qcfail: # reads.reference_start
        return False
    
    if not read.has_tag("CB"):
        return False

    return True


def process_bam(args) -> typing.List[RelativeInfo]:
    u"""
    Extract information of relative path of utr from bam
    """
    bam_path, utr = args
    logger.info("{} processing {}".format(os.getpid(), utr[0].chrom))
    res = []
    with pysam.AlignmentFile(bam_path) as r:
        for u in utr:
            
            reads = r.fetch(u.chrom, u.start, u.end)
            u = RelativeInfo(u)

            for read in reads:
                if not filter(read):
                    continue
                u.add(read)

            res.append(u)

    return res
        

@click.command()
@click.option(
    '-u',
    '--utr',
    type=click.Path(exists=True),
    help='The path to bed file of utr.',
    required=True
)
@click.option(
    '-b',
    '--bam',
    type=click.Path(exists=True),
    help='The gtf (gz or not) file for preparing utr and intron region.',
    required=True
)
@click.option(
    '-o',
    '--output',
    type=str,
    help ='The path to output file'
)
@click.option(
    '-p',
    '--process',
    type=click.IntRange(1, cpu_count(), clamp=True),
    help='The number of cpu to use',
    default=1
)
def extract(
    utr, bam, output, process
    ):

    logger.info("load from utr")

    utr_data = {}
    with gzip.open(utr, "rt") if utr.endswith("gz") else open(utr, "r") as r:
        for line in tqdm(r):
            lines = line.strip().split("\t")

            try:
                temp = Bed6(
                    chrom=lines[0],
                    start=int(lines[1]),
                    end=int(lines[2]),
                    strand=lines[3]
                )
                temp_lst = utr_data.get(temp.chrom, [])
                temp_lst.append(temp)
                utr_data[temp.chrom] = temp_lst
            except ValueError as err:
                logger.error(err)
                logger.error(line)
                exit(1)

    with pysam.AlignmentFile(bam) as r:
        cmds = [[bam, utr_data[x]] for x in set(utr_data.keys()) & set(r.references)]
        # cmds = [[bam, utr_data[str(x)]] for x in [6, 8, 9, 14, "Y"]]

    with Pool(process) as p:
        res = list(tqdm(p.imap_unordered(process_bam, cmds), total=len(cmds)))
            
    p.close()
    p.join()

    outdir = os.path.dirname(os.path.abspath(output))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with open(output, "w+") as w:
        w.write(RelativeInfo.header() + "\n")
        for r in tqdm(res, desc="Writing"):
            for u in r:
                w.write(str(u) + "\n")
