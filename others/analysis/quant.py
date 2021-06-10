import os
from glob import glob

from subprocess import check_call


bams = [os.path.realpath(x) for x in glob("/mnt/raid64/ATS/Personal/zhangyiming/bams/*.bam") if "posso" in os.path.realpath(x)]

output = "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/"
os.makedirs(output, exist_ok = True)


for b in bams:
    b = os.path.dirname(b)
    key = os.path.basename(os.path.dirname(b))
    key = key.split("-")[0]
    print(key)
    o = f"{output}/{key}.quant"
    if os.path.exists(o):
        continue
    check_call(f"julia -t 10 /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/quant.jl -i /mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_cage_R.txt -c {b} -o {o}", shell=True)

