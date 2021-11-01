#!/usr/bin/envpython3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.27 by Zhang

Contians all the parameters and command line params handler
"""
import click

from process.preprocess import process


@click.command()
@click.option(
    "-g", "--gtf",
    type=click.Path(exists=True),
    help=""" The path to input gtf file. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    help=""" The path to output utr file, bed format. """
)
@click.option(
    "-l", "--utr-length",
    type=int,
    default=1000,
    help=""" The radius of UTR. """
)
def preprocess(gtf: str, output: str, utr_length: int):
    u"""
    Preprocess: extract UTR from gtf file
    \f
    """
    process(gtf=gtf, output=output, span=utr_length // 2)


if __name__ == '__main__':
    preprocess()
