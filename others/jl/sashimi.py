#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os

from glob import glob
from subprocess import check_call

__dir__ = os.path.abspath(os.path.dirname(__file__))


def main(input_dir:str, output:str):
    gtf="/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.gtf"
    cellranger="/mnt/raid64/Covid19_Gravida/cellranger"

    for i in glob(os.path.join(input_dir, "*_R2.bed")):
        key = os.path.basename(i)
        key = key.replace("_R2.bed", "")
        
        print(i, key)
        bam = glob(os.path.join(cellranger, key + "*", "outs/possorted_genome_bam.bam"))
        bam=bam[0]
        try:
            check_call(f"julia {__dir__}/sashimi.jl -i {i} -o {output}/{key} -b {bam} -r {gtf} -p 10", shell=True, timeout=200)
        except Exception as err:
            print(err)


if __name__ == '__main__':
    main(
        "/mnt/raid64/Covid19_Gravida/apamix/R2",
        "/mnt/raid64/Covid19_Gravida/apamix/R2/sashimi",
    )

