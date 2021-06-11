#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.05.10 by Zhang

Functions used to count ats and calculate psi
"""
import click
from multiprocessing import cpu_count

from process.postprocess import count, psi


@click.command()
@click.option(
    "-i", "--ats",
    type=click.Path(exists=True),
    help=""" The path to inferred ats sites. """
)
@click.option(
    "-b", "--bam",
    type=click.Path(exists=True),
    help=""" The file contains path to bams. """
)
@click.option(
    "--delimiter",
    type=str,
    default="\t",
    help="The delimiter of input bam list"
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    help=""" The path to output utr file, bed format. """
)
@click.option(
    "-p", "--processes",
    type=click.IntRange(1, cpu_count()),
    default=1,
    help=""" How many cpu to use. """
)
def postprocess(ats: str, output: str, bam: str, delimiter: str, processes: int):
    u"""
    Postprocess: count ats and calculate psi
    \f
    """
    bams = {}
    with open(bam) as r:
        for line in r:
            line = line.strip().split(delimiter)

            if line:
                bams[line[0]] = None if len(line) < 2 else line[1]

    count(bams, ats = ats, output = f"{output}.count", processes = processes)
    psi(f"{output}.count", f"{output}.psi")


if __name__ == '__main__':
    pass
