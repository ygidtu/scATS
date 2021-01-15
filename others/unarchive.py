#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.01.07
"""
import os

from glob import glob
from shutil import rmtree

import patoolib

from rich import print


def main(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    for f in glob(os.path.join(input_dir, "*.tar")):
        o = os.path.join(output_dir, os.path.basename(f).replace(".tar", ""))
        if not os.path.exists(o):
            try:
                patoolib.extract_archive(f, outdir=output_dir)
            except Exception as err:
                print(err)
                if os.path.exists(o):
                    rmtree(o)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)