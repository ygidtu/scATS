#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang at 2021.06.16

This is used to calculate the co-expression ATS
"""

import click

from process.correlation import corr


@click.command()
@click.option(
    "-i", "--input",
    type=click.Path(exists=True),
    required=True,
    help=""" The path to counts. Sorted by first column. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
)
@click.option(
    "-d", "--distance",
    type=int,
    default=1000,
    help=""" The maxmimum distance between two UTRs. """
)
@click.option(
    "--spearman",
    is_flag = True,
    help=""" Use spearman correlation instead of pearson. """
)

def coexp(input: str, output: str, distance: int, spearman: bool):
    u"""
    Co-expression

    Note: The input count file must be sorted by first column.
    """
    corr(input, output, distance=distance, pearson = not spearman)

if __name__ == '__main__':
    pass
