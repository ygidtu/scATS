#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Preprocess of UTR bed and BAM, to get the data for model
"""
import gzip
import os
import re

from typing import Dict, List

import pysam

from rich.progress import Progress, track

try:
    from src.loci import BED, Reads
except ImportError:
    from loci import BED, Reads


def load_ats(path: str) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output ifle
    :return dict: keys -> utr region; values -> list of splice region
    """
    beds = {}
    header = True
    with open(path) as r:
        for line in r:
            if header:
                header = False
                continue

            line = line.strip().split("\t")
            
            site = re.split(r"[:-]", line[0])
            # print(site)
            strand = site[-1]

            if strand != "+":
                strand = "-"
            chrom, start_pos, end_pos = site[0], int(site[1]), int(site[2])

            utr = BED(chrom, start_pos, end_pos, strand, name = line[2], record_id = line[2])
            alpha = line[3].split(",")

            if len(alpha) > 1:
                if utr not in beds.keys():
                    beds[utr] = []
                try:
                    for x in alpha:
                        if x != "":
                            s = int(float(x))
                            
                            if strand == "+":
                                beds[utr].append(BED(chrom, s, s + 1, strand, line[0], str(len(beds) + 1))) 
                            else:
                                beds[utr].append(BED(chrom, s - 1, s, strand, line[0], str(len(beds) + 1)))
                except Exception as e:
                    print(e)
                    pass

    return beds


def load_utr(path: str) -> List[BED]:
    u"""
    Load extracted UTR from bed file
    :param path: path to bed file
    :return list of BED objects
    """

    res = []

    with Progress() as progress:
        task_id = progress.add_task(f"Reading {os.path.basename(path)}", total=os.path.getsize(path))

        with gzip.open(path, "rt") if path.endswith("gz") else open(path) as r:
            for line in r:
                progress.update(task_id, advance=len(str.encode(line)))

                res.append(BED.create(line))

    return res


def load_reads(bams: List[str], region: BED) -> List[Reads]:
    u"""
    Load reads from bams files
    :param bams: list of bam path
    :param region: Region in BED object
    :return list of Reads
    """
    res = []
    for bam in bams:
        with pysam.AlignmentFile(bam) as r:
            for record in r.fetch(region.chromosome, region.start, region.end):
                record = Reads.create(record)
                if record and record.strand == region.strand:
                    res.append(record)

    return res


def load_paired_reads(bam: List[str], region: BED) -> Dict:
    u"""
    Load paired reads, keys -> R1; values -> R2
    """
    res = {}
    for b in bam:
        with pysam.AlignmentFile(b) as r:
            for rec in r.fetch(region.chromosome, region.start, region.end):

                if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                    continue

                if rec.is_read1:
                    mate = r.mate(rec)
                    r1 = Reads.create(rec)
                    r2 = Reads.create(mate)

                    if r1 and r2:
                        res[r1] = r2

    return res


if __name__ == '__main__':
    pass
