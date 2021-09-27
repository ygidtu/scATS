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
@click.option(
    "-c", "--compress",
    is_flag=True,
    type=click.BOOL,
    help=""" Wheter to save in gzip format. """
)
@click.option(
    "--bulk",
    is_flag=True,
    type=click.BOOL,
    help=""" Wheter the input bam is Nanopore or PacBio. """
)
def postprocess(
    ats: str, output: str, 
    bam: str, delimiter: str, 
    processes: int,compress: bool,
    bulk: bool
):
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

    c = f"{output}.count.gz" if compress else f"{output}.count"
    p = f"{output}.psi.gz" if compress else f"{output}.psi"

    count(bams, ats = ats, output = c, processes = processes, bulk = bulk)
    psi(c, p)


if __name__ == '__main__':
    pass
