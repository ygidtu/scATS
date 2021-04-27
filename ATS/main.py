#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Main function 
"""
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.use('Agg')

from cli.cli import cli

if __name__ == '__main__':
    cli()
    pass


