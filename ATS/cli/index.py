#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.05.14 by Zhang
Create reads index from bam to accelerate AtsModel using GO Core
"""

import os
from multiprocessing import cpu_count, Pool
from rich.progress import track
import pysam
import pickle
import click
from logger import log
from src.reader import load_utr, load_reads

__dir__ = os.path.dirname(os.path.abspath(__file__))
__exe__ = os.path.join(__dir__, "../process/fetch/afe")


def check_pickle(path: str) -> bool:
    try:
        with open(path, "rb") as r:
            pickle.load(r)
        return True
    except Exception:
        return False


def run(args):
    o, u, bams = args

    if not check_pickle(o):
        return 

    try:
        res = load_reads(bams, u)

        with open(o, "wb+") as w:
            pickle.dump(res, w)
    except ValueError as err:
        log.error(err)

        with open(o, "wb+") as w:
            pickle.dump({}, w)
    

@click.command()
@click.option(
    "-u", "--utr",
    type=click.Path(exists = True),
    required=True,
    help=""" The path to utr file, bed format. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
)
@click.option(
    "-p", "--processes",
    type=click.IntRange(1,cpu_count()),
    default = 1,
    help=""" How many cpu to use. """
)
@click.argument("bams", nargs = -1, type=click.Path(exists=True), required=True)
def index(bams, utr, processes: int, output: str):
    u"""
    Index is used to reduce the IO time, It‘s helpful for the parameter testing
    """
    
    os.makedirs(output, exist_ok=True)

    utr = load_utr(utr)

    with open(os.path.join(output, "index.pkl"), "wb+") as w:
        pickle.dump(utr, w)

    cmds = [[os.path.join(output, f"{idx}.pkl"), u, bams] for idx, u in enumerate(utr)]

    with Pool(processes) as p:
        list(track(p.imap_unordered(run, cmds), total = len(cmds)))
    

if __name__ == '__main__':
    pass
