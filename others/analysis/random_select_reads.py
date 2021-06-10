#!/usr/bin/env python3
# -*- coding:utf -8 -*-
u"""
Create smaller bam files
"""
import logging
import os
import gzip
import random
from multiprocessing import Pool
from subprocess import check_call
from shutil import rmtree

import click

from rich.logging import RichHandler

log = logging.getLogger(("rich"))

def init_logger(level = "NOTSET"):
    FORMAT = "%(message)s"
    logging.basicConfig(
        level=level, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
    )
    global log
    log = logging.getLogger(("rich"))


def create_name(bam, output):
    temp = f"{output}.temp"

    os.makedirs(temp, exist_ok=True)

    check_call("samtools view %s | pv | awk '{print $1}' | sort -u -T %s | gzip > %s" % (bam, temp, output), shell=True)

    if os.path.exists(temp):
        rmtree(temp)


def process(args):
    bam, output = args
    log.info(f"start processing {bam}")
    name = f"{bam}.name.gz"

    if not os.path.exists(name):
        create_name(bam, name)
    
    if os.path.getsize(name) > 0 and not os.path.exists(f"{bam}.name_sel.gz"):
        names = set()
        with gzip.open(name, "rt") as r:
            with gzip.open(f"{bam}.name_sel.gz", "wt+") as w:
                for line in r:
                    names.add(line)

                    if len(names) > 10000:
                        random.seed(42)
                        names = random.sample(names, len(names) // 10)
                        for i in names:
                            w.write(i + "\n")
                        names = set()

                if names:
                    random.seed(42)
                    names = random.sample(names, len(names) // 10)
                    for i in names:
                        w.write(i + "\n")
                    names = set()

    if os.path.exists(f"{bam}.name_sel.gz") and not os.path.exists(output):
        # with gzip.open(f"{bam}.name_sel.gz", "rt") as r:
        #     names = r.readlines()

        # log.info(f"random select {len(names)} from {bam}")
        # with pysam.AlignmentFile(bam) as r:
        #     with pysam.AlignmentFile(output, "wb", template = r) as w:
        #         for rec in r:
        #             if rec.query_name in names:
        #                 w.write(rec)

        check_call(f"{os.path.dirname(os.path.abspath(__file__))}/go/extract/afe -i {bam} -s {bam}.name_sel.gz -o {output} -t 20", shell=True)
        log.info(f"finished with {bam}")


@click.command()
@click.option("-p", "--processes", type = int, default = 10)
@click.option("-o", "--output", type = click.Path())
@click.argument("bams", nargs = -1, type=click.Path(exists=True), required=True)
def main(processes: int, bams, output: str):
    os.makedirs(output, exist_ok = True)
    init_logger()
    cmds = [[b, os.path.join(output, os.path.basename(b))] for b in bams]
    
    with Pool(processes) as p:
        p.map(process, cmds)


if __name__ == '__main__':
    main()
