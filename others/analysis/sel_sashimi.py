import os
import pandas as pd
from subprocess import check_call
from multiprocessing import Pool
import numpy as np


data = pd.read_excel("/mnt/raid64/ATS/Personal/zhangyiming/quant/sel_site.xlsx")
output = "/mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel"

os.makedirs(output, exist_ok=True)

cmds = []
for _, row in data.iterrows():

    cmd = f"sashimiplot junc --gtf {row['gtf']} --bam {row['bam']} --sj 100 --junc {row['chrom']}:{row['start']}:{row['end']} --fileout {output}/{row['id']}.pdf --trackline {row['line']} --ie 1,1  --ps RF --ssm R1"
    # print(row)
    if row['bc'] != "aa":
        cmd = cmd + f" --bc {row['bc']} --co {row['order']}"
    cmds.append(cmd)


def call(cmd):
    # print(cmd)
    with open(os.devnull, "w") as w:
        try:
            check_call(cmd, shell=True, stdout = w, stderr = w)
        except Exception as err:
                print(err)
                print(cmd)


with Pool(10) as p:
    # p.map(call, cmds)
    for i in cmds:
        print(i)

# sashimiplot junc --gtf /mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz --bam /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam2.tsv --sj 100 --junc X:65667610:65667720 --fileout /mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel/MSN_mono.pdf --trackline 65667652,65667711 --ie 1,1  --ps RF --ssm R1 --bc /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_bc.txt --co /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_order.txt
# sashimiplot junc --gtf /mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz --bam /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/nslc_bam.tsv --sj 100 --junc 3:81761600:81761800 --fileout /mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel/nslc_2_1.pdf --trackline 81761739,81761642,81761542 --ie 1,1  --ps RF --ssm R1 --bc /mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/celltype_bc.txt --co /mnt/raid64/ATS/Personal/zhangyiming/quant/NSLC/cellclass_order.txt

# sashimiplot junc --gtf /mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz --bam /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam4.tsv --sj 100 --junc 11:110430070:110430300 --fileout /mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel/FDX1_mono.pdf --trackline 110430104,110430166 --ie 1,1  --ps RF --ssm R1 --bc /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_bc.txt --co /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_order.txt


# sashimiplot junc --gtf /mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz --bam /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam2.tsv --sj 100 --junc 1:152036970:152037025 --fileout /mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel/S100A11_mono.pdf --trackline 152036972,152037021 --ie 1,1  --ps RF --ssm R1 --bc /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_bc.txt --co /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/mono_order.txt

# sashimiplot junc --gtf /mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz --bam /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam2.tsv --sj 100 --junc 13:46168460:46168600 --fileout /mnt/raid64/ATS/Personal/zhangyiming/quant/sashimi_sel/LCP1_covid.pdf --trackline 46168577,46168494 --ie 1,1  --ps RF --ssm R1 --bc /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/cellclass_bc.txt --co /mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/cellclass_order.txt