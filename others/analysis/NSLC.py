#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Create at 2021.03.27
make sashimi plot of NSCLC ATS
"""
import os

from multiprocessing import Pool
from subprocess import check_call

import pandas as pd

from rich import print
from tqdm import tqdm


data = pd.read_csv("/mnt/raid61/Personal_data/zhangyiming/code/afe/others/lusc.csv")

data.head()


gtf = "/mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz"
bam = "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam.tsv"
# "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/nslc_bam.tsv"
output = "/mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/sashimi/tissue"
file_use = {
    "Class": [
        "/mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/celltype1_bc.txt",
        "/mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/cellclass1_order.txt"
    ],
    "Cell": [
        "/mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/celltype_bc.txt",
        "/mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/cellclass_order.txt"
    ]
}


def decode_site(site, span = 100):
    site = site.split(":")

    sites = site[1].split("-")

    return f"{site[0].replace('chr', '')}:{int(sites[0]) - span}:{int(sites[1]) + span}"


cmds = []
for _, row in data.iterrows():
    for key, bo in file_use.items():
        o = os.path.join(output, key)
        os.makedirs(o, exist_ok=True)

        junc = decode_site(row['hg38'])

        if row['geneSymbol'] != "ARHGAP24":
            continue

        cmd = f"sashimiplot junc --gtf {gtf} --bam {bam} --sj 100 --junc {junc} --fileout {o}/{row['geneSymbol']}_PBMC.pdf --trackline {row['site']} --ie 1,1  --ps RF --ssm R1 " # --bc {bo[0]} --co {bo[1]}"

        cmds.append(cmd)


def call(cmd):
    with open(os.devnull, "w+") as w:
        check_call(cmd, shell=True, stdout = w, stderr = w)


with Pool(10) as p:
    list(tqdm(p.imap(call, cmds), total=len(cmds)))

