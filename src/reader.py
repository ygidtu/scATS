#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Preprocess of UTR bed and BAM, to get the data for model
"""
import gzip
import pickle
import os
import random
import re
from typing import Dict, List

import pysam

from src.loci import BED, Reads
from src.progress import custom_progress


def load_ats(path: str, min_ats: int = 1) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output file
    :param min_ats: the minimum number of ATS locate in single UTR
    :return dict: keys -> utr region; values -> list of splice region
    """

    beds = {}
    header = True
    progress = custom_progress(io = True)
    with progress:
        task_id = progress.add_task(
            f"Reading... ", total=os.path.getsize(path))

        with open(path) as r:
            for curr_idx, line in enumerate(r):
                if header:
                    header = False
                    continue
                progress.update(task_id, advance=len(str.encode(line)))
                line = line.strip().split("\t")

                site = re.split(r"[:-]", line[0])

                strand = site[-1]

                if strand != "+":
                    strand = "-"
                chrom, start_pos, end_pos = site[0], int(site[1]), int(site[2])

                utr = BED(
                    chrom, start_pos, end_pos, strand,
                    name="",
                    record_id=""
                )

                alpha = line[3].split(",")
                if len(alpha) > min_ats:
                    if utr not in beds.keys():
                        beds[utr] = set()
                    try:
                        for x in alpha:
                            if x != "":
                                x = int(float(x))

                                beds[utr].add(
                                    BED(chrom, x - 1, x, strand, line[0], str(len(beds) + 1)))
                    except Exception as e:
                        print(e)
                        pass
    return beds


def load_utr(path: str, utr_length: int = 1000, debug: bool = False) -> List[BED]:
    u"""
    Load extracted UTR from bed file
    :param path: path to bed file
    :return list of BED objects
    """

    res = []

    progress = custom_progress(io = True)

    with progress:
        task_id = progress.add_task(
            f"Reading...  ", total=os.path.getsize(path))

        with gzip.open(path, "rt") if path.endswith("gz") else open(path) as r:
            for line in r:

                if debug and len(res) > 50:
                    return res
                progress.update(task_id, advance=len(str.encode(line)))

                b = BED.create(line)
                if len(b) > utr_length:
                    continue
                res.append(b)

    if not debug:
        random.seed(42)
        random.shuffle(res)
    return res


def __get_strand__(read: pysam.AlignedSegment) -> str:
    u"""
    determine the reads strand

    :params read: 
    """

    if read.is_paired:
        if read.is_read1 and read.is_reverse:
            return "-"
        elif read.is_read2 and not read.is_reverse:
            return "-"
        return "+"
 
    return "-" if read.is_reverse else "+"
   

def load_barcodes(path: str) -> dict:
    u"""
    load required barcodes

    :params path
    """
    res = {}

    if not os.path.exists(path):
        return res

    with open(path) as r:
        for line in r:
            line = line.strip()

            key = line[:min(3, len(line))]

            temp = res.get(key, set())
            temp.add(line)
            res[key] = temp

    return res


def __is_barcode_exists__(barcodes: dict, rec: pysam.AlignedSegment) -> bool:
    u"""
    check whether this read contains required barcode
    :params barcodes: a collection of required barcodes
    :params rec: 
    """
    if not rec.has_tag("CB"):
        return False
    
    cb = rec.get_tag("CB")

    return cb[:min(3, len(cb))] in barcodes and cb in barcodes[cb[:min(3, len(cb))]]


def load_reads(bam: List[str], region: BED, barcode):
    u"""
    Load reads, keys -> R1; values -> R2
    Only both R1

    @2021.06.10 remove dict from this function, use list and sort to increase function speed, with disadvantage like more memory usage

    :params bam: list of bam files
    :params region:
    :return generator: generate r1 and r2
    """
    for b in bam: 
        r1s, r2s = [], []

        r = pysam.AlignmentFile(b) if isinstance(b, str) else b

        # fetch with until_eof is faster on large bam file according to
        # https://pysam.readthedocs.io/en/latest/faq.html

        for rec in r.fetch(region.chromosome, region.start, region.end, until_eof=True):
            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                continue

            if __get_strand__(rec) != region.strand:
                continue
            
            # only use the required barcodes for analysis
            if barcode[b]:
                if not __is_barcode_exists__(barcode[b], rec):
                    continue

            if rec.is_read1:
                r1s.append(rec)
            else:
                r2s.append(rec)

        if isinstance(b, str):
            r.close()

        r1s = sorted(r1s, key=lambda x: x.query_name)
        r2s = sorted(r2s, key=lambda x: x.query_name)

        i, j = 0, 0
        while i < len(r1s) and j < len(r2s):
            r1 = r1s[i]
            r2 = r2s[j]

            if r1.query_name < r2.query_name:
                i += 1
            elif r1.query_name > r2.query_name:
                j += 1
            else:
                r1 = Reads.create(r1)
                r2 = Reads.create(r2)
                if r1 and r2 and region.is_cover(r1, 0) and region.is_cover(r2, 0):
                    yield r1, r2

                i += 1


def check_bam(path: str) -> bool:
    try:
        with pysam.AlignmentFile(path) as r:
            pass
    except Exception as err:
        return False
    return True


if __name__ == '__main__':
    pass
