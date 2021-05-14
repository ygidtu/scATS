#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.05.14 by Zhang
Create reads index from bam to accelerate AtsModel using GO Core
"""

import os
from subprocess import check_call
from multiprocessing import cpu_count

import click

__dir__ = os.path.dirname(os.path.abspath(__file__))
__exe__ = os.path.join(__dir__, "../process/fetch/afe")



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
def index(bams, utr: str, processes: int, output: str):
    cmd = f"{__exe__} -u {os.path.abspath(utr)} -o {os.path.abspath(output)} -t {processes}"
    for i in bams:
        cmd = f"{cmd} -i {os.path.abspath(i)}"

    check_call(cmd, shell=True)


if __name__ == '__main__':
    pass
