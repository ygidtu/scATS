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
from rich.progress import Progress

from src.loci import BED, Reads


def load_ats(path: str, julia: bool = False) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output ifle
    :return dict: keys -> utr region; values -> list of splice region
    """
    
    beds = {}
    header = True
    with Progress() as progress:
        task_id = progress.add_task(f"Reading {os.path.basename(path)}", total=os.path.getsize(path))

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
                    name = line[2] if not julia else line[0], 
                    record_id = line[2] if not julia else str(curr_idx)
                )

                alpha = line[3].split(",") if not julia else line[2].split(",")

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

def load_utr(path: str, utr_length: int = 1000, debug: bool = False) -> List[BED]:
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


# @profile
def load_reads(bam: List[str], region: BED) -> Dict:
    u"""
    Load reads, keys -> R1; values -> R2
    Only both R1

    :params bam: list of bam files
    :params region:
    """
    res = {}
    for b in bam:
        paired = {}
        r = pysam.AlignmentFile(b) if isinstance(b, str) else b

        # fetch with until_eof is faster on large bam file according to 
        # https://pysam.readthedocs.io/en/latest/faq.html
        for rec in r.fetch(region.chromosome, region.start, region.end, until_eof=True):
            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                continue

            if rec.is_read1:
                r1 = Reads.create(rec)
                if rec.query_name not in paired.keys():
                    paired[rec.query_name] = r1
                else:
                    r2 = paired.pop(rec.query_name)
                    if r1 and r2:
                        res[r1] = r2
            else:
                r2 = Reads.create(rec)

                if rec.query_name not in paired.keys():
                    paired[rec.query_name] = r2
                else:
                    r1 = paired.pop(rec.query_name)
                    if r1 and r2:
                        res[r1] = r2

        if isinstance(b, str):
            r.close()

    return {i: j for i, j in res.items() if region.is_cover(i, 0) and region.is_cover(j, 0)}


def check_bam(path: str) -> bool:
    try:
        with pysam.AlignmentFile(path) as r:
            pass
    except Exception as err:
        return False
    return True


class Index(object):

    def __init__(self, path: str):
        self.path = path
        self.utr = self.__load_utr__()

    def __len__(self):
        return len(self.utr)

    def __load_utr__(self):
        with open(os.path.join(self.path, "index.pkl"), "rb") as r:
            return pickle.load(r)

    def get(self, idx: int):
        if idx < len(self.utr):
            utr = self.utr[idx]

            with open(os.path.join(self.path, f"{idx}.pkl"), "rb") as r:
                values = pickle.load(r)

            return utr, values
        return None, None


if __name__ == '__main__':
    pass
