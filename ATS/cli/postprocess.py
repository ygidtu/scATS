#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.05.10 by Zhang

Functions used to count ats and calculate psi
"""
import click

from process.postprocess import count, psi

@click.command()
@click.option(
    "-i", "--ats",
    type=click.Path(exists = True),
    help = """ The path to inferred ats sites. """
)
@click.option(
    "-b", "--bam",
    type=int,
    default = 500,
    help=""" The file contains path to bams. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    help=""" The path to output utr file, bed format. """
)
def postprocess(ats: str, output: str, bam: int):
    u"""
    Postprocess: count ats and calculate psi
    \f
    """

    count(bams, ats, processes, f"{output}.count")
    psi(f"{output}.count", f"{output}.psi")

                 

if __name__ == '__main__':
    pass