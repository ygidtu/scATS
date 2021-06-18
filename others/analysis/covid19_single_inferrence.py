#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.06.17 by Zhang
"""
import os
from glob import glob
from subprocess import check_call


def main(input_bam: str, utr: str, output: str):
    u"""
    :params input_bam: path to bam directory
    :params utr: path to utr
    :params output: path to output directory
    """

    os.makedirs(output, exist_ok=True)

    for b in glob(os.path.join(input_bam, "*.bam")):
        key = os.path.basename(b).split(".")[0]
        print(key)
        check_call(f"python /mnt/raid61/Personal_data/zhangyiming/code/afe/main.py ats --utr {utr} --utr-length 1000 -p 20 --output {output}/{key} {b}", shell=True)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
