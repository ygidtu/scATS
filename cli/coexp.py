#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang at 2021.06.16

This is used to calculate the co-expression ATS
"""

import click

from core.correlation import corr


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
    "-b", "--barcode",
    type=click.Path(),
    help=""" The list of barcodes to use. """
)
@click.option(
    "--spearman",
    is_flag = True,
    help=""" Use spearman correlation instead of pearson. """
)
@click.option(
    "-g", "--group-info",
    type=click.Path(),
    help=""" Path to file contains group information, two columns required, 1st barcodes, 2nd group name. """
)
def coexp(input: str, output: str, distance: int, spearman: bool, barcode: str, group_info: str):
    u"""
    Co-expression

    Note: The input count file must be sorted by first column.
    """
    corr(input, output, distance=distance, pearson=not spearman, barcode=barcode, group_info=group_info)


if __name__ == '__main__':
    pass
