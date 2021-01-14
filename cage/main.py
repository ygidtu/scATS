#!/usr/bin/env python3
# -*- utr-8 -*-
u"""
Created at 2021.01.06
"""
import os
import re

from glob import glob
from multiprocessing import Pool
from subprocess import check_call

__dir__ = os.path.abspath(os.path.dirname(__file__))


def call(cmd):
    check_call(cmd, shell=True)


def main(input_dir: str, output: str, n_jobs: int = 12):
    os.makedirs(output, exist_ok = True)

    bams = glob(f"{input_dir}/*/outs/possorted_genome_bam.bam")

    cmds, cmds1 = [], []
    for b in bams:
        c = os.path.join(os.path.dirname(b), "filtered_feature_bc_matrix/barcodes.tsv.gz")
        key = re.sub(r"[^A-Z]", "", os.path.basename(os.path.dirname(os.path.dirname(b))))

        # b = f"/mnt/raid64/Covid19_Gravida/apamix/bam/{key}.bam"

        cmds.append(f"python {__dir__}/generate_ctss.py {b} {output}/{key}_R1.bed")
        cmds1.append(f"python {__dir__}/bc_umi.4.py {b} {c} {output}/peak_region.txt {output}/{key}.csv.gz {n_jobs}")
        
    # with Pool(len(cmds)) as p:
    #     p.map(call, cmds)
    #     # check_call(, shell=True)
    
    # call(f"Rscript {__dir__}/cage.R {output}")

    for c in cmds1:
        print(c)

        # if os.path.exists(c.split()[-2]):
        #     continue

        call(c)
    

if __name__ == '__main__':


    main(
        "/mnt/raid64/Covid19_Gravida/cellranger/",
        "/mnt/raid64/ATS/Personal/zhangyiming/cage"
    )


