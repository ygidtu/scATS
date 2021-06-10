#!/usr/bin/env python3
# -*- coding:utf-8
import os
from subprocess import check_call

from glob import glob


jl = "/mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/quant.jl"
i = "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_mm1.txt"
fs = glob("/mnt/raid64/ATS/alignments/cellranger/mm/*pbmc*/outs/possorted_genome_bam.bam")
output = "/mnt/raid64/ATS/Personal/zhangyiming/quant/mm"

for f in fs:
    key = os.path.basename(os.path.dirname(os.path.dirname(f)))

    check_call(f"julia -t 10 {jl} -i {i} -c {os.path.dirname(f)} -o {output}/{key}.quant", shell=True)

jl = "/mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/psi.jl"
g = "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_mm.postprocess"
for f in fs:
    key = os.path.basename(os.path.dirname(os.path.dirname(f)))

    check_call(f"julia -t 10 {jl} -i {output}/{key}.quant -g {g} -o {output}/{key}.psi", shell=True)
