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


def load_ats(path: str, julia: bool = False) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output ifle
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
                    name=line[2] if not julia else line[0],
                    record_id=line[2] if not julia else str(curr_idx)
                )

                alpha = line[4].split(",") if not julia else line[2].split(",")
                if len(alpha) > 1:
                    if utr not in beds.keys():
                        beds[utr] = []
                    try:
                        for x in alpha:
                            if x != "":
                                x = int(float(x))
                                s = utr.start + x if strand == "+" else utr.end - x

                                beds[utr].append(
                                    BED(chrom, s - 1, s, strand, line[0], str(len(beds) + 1)))
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


def load_reads(bam: List[str], region: BED):
    u"""
    Load reads, keys -> R1; values -> R2
    Only both R1

    @2021.06.10 remove dict from this function, use list and sort to increase function speed, with disadvantage like more memory usage

    :params bam: list of bam files
    :params region:
    :return generator: generate r1 and r2
    """
    res = {}
    for b in bam:
        r1s, r2s = [], []

        r = pysam.AlignmentFile(b) if isinstance(b, str) else b

        # fetch with until_eof is faster on large bam file according to
        # https://pysam.readthedocs.io/en/latest/faq.html

        for rec in r.fetch(region.chromosome, region.start, region.end, until_eof=True):
            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
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
