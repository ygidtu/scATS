#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Preprocess of UTR bed and BAM, to get the data for model
"""
import gzip
import os

from typing import List

import pysam

from rich.progress import Progress, track


from src.loci import BED, Reads



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


def load_reads(bams: List[str], utr: BED) -> List[Reads]:
    u"""
    Load reads from bams files
    :param bams: list of bam path
    :param utr: UTR in BED object
    :return list of Reads
    """
    res = []
    for bam in bams:
        with pysam.AlignmentFile(bam) as r:
            for record in r.fetch(utr.chromosome, utr.start, utr.end):
                record = Reads.create(record)
                if record and record.strand == utr.strand:
                    res.append(record)

    return res


if __name__ == '__main__':
    pass
