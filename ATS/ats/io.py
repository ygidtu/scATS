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


from src.loci import BED



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


def load_reads(bam: List[str], utr: BED) -> :

