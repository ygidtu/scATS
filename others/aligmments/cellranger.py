#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.01.06
"""
import os
import re
from glob import glob
from subprocess import check_call

from rich import print


REFERENCE = {
    "hs": "/mnt/raid64/ATS/cellranger/ref/Homo_sapiens/",
    "mm": "/mnt/raid64/ATS/cellranger/ref/Mus_musculus/"
}


def collect_files(input_dir:str, softlinks: str):
    os.makedirs(softlinks, exist_ok = True)
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith("fq.gz") or f.endswith("fastq.gz"):
                o = os.path.join(softlinks, f)
                if not os.path.exists(o):
                    os.symlink(os.path.join(parent, f), o)

    fs = set()
    for f in glob(os.path.join(softlinks, "*.fastq.gz")):
        key = re.sub(r"_S1_L\d{3}_[RI]\d_\d{3}.f(ast)?q.gz", "", os.path.basename(f)) 
        fs.add(key)

    return fs


def main(input_dir: str, softlinks: str, output: str):
    # cellranger count --id LF-2-1 --fastqs=/mnt/raid64/Covid19_Gravida/raw_data/filter/LF-2-1 --transcriptome=/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens --localcores=20 --localmem=80 --sample=LF-2-1

    fs = collect_files(input_dir, softlinks)

    print(fs)
    os.makedirs(output, exist_ok=True)
    for f in fs:
        
        try:
            ref = REFERENCE[f.split("_")[2]]
            cmd = f"cellranger count --id {f} --fastqs={softlinks} --transcriptome={ref} --localcores=20 --localmem=80 --sample={f}"
            if f.split("_")[2] == "hs":
                continue
            print(cmd)
            check_call(cmd, shell=True, cwd=os.path.join(output, f.split("_")[2]))
        except Exception as err:
            print(f)
            print(err)
        pass


if __name__ == '__main__':
    main(
        "/mnt/raid64/ATS/Personal/zhangyiming/data",
        "/mnt/raid64/ATS/alignments/cellranger/softlinks",
        "/mnt/raid64/ATS/alignments/cellranger/"
    )

    # set bg_rgb=[1,1,1]
    # space cmyk
    # set orthoscopic, on   # turns it on
    # set ray_trace_fog=0
    # set depth_cue=0
    # util.performance(0)
    # set antialias, 4
    # png /Users/zhangyiming/Downloads/nCOV/software/500_532_N.png, width=10cm, dpi=3000, ray=1
    
    