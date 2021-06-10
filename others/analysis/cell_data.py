import os
from glob import glob
from subprocess import check_call
import random

bams = glob("/mnt/raid64/Chen_Cell_2016/alignment/STAR/EGAD0000100267*/*.sorted*.bam")

data = {}

for b in bams:
    key = os.path.basename(os.path.dirname(b))
    temp = data.get(key, [])
    temp.append(b)
    data[key] = temp

sel_bams = []

for val in data.values():
    val = sorted(val)
    random.seed(42)
    sel_bams += random.sample(val, 10)

ref = "/mnt/raid61/Personal_data/zhangyiming/ePSI/ref/genecode.v30lift37.gene_utr1.bed"
# ref/genecode.v30lift37.gene_utr.bed
# 
check_call(f"julia /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/run.jl  -i {ref} -o /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_new.txt --using-R 12 --cage-mode {' '.join(sel_bams)}", shell=True)


sel_bams = []

for val in data.values():
    val = sorted(val)
    random.seed(42)
    sel_bams += random.sample(val, 3)

try:
    check_call(f"julia -t 12 /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/quant1.jl  -i /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200.postprocess.bed -o /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200_sample.quant {' '.join(sel_bams)}", shell=True)
except Exception as err:
    pass


try:
    check_call(f"julia -t 12 /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/quant1.jl  -i /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200.postprocess_cage.bed -o /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200.quant {' '.join(bams)}", shell=True)
except Exception as err:
    pass

try:
    check_call(f"/mnt/raid61/Personal_data/zhangyiming/code/afe/others/go/count/count -t 12 -b /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200.postprocess.bed -o /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall_200_go.quant {' '.join(['-i ' + x for x in bams])}", shell=True)
except Exception as err:
    pass


# /mnt/data5/zhangyiming/

#     Mono = get(psi[["EGAD00001002674"]])
#     Neu = get(psi[["EGAD00001002675"]])
#     TCell = get(psi[["EGAD00001002671"]])

# julia /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/sashimi.jl -i /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall.txt -b /mnt/raid61/Personal_data/zhangyiming/code/afe/tests/cell_bam.tsv -r /mnt/raid61/Personal_data/zhangyiming/ePSI/ref/gencode.v30lift37.annotation.sorted.gtf.gz -o /mnt/raid61/Personal_data/zhangyiming/ePSI/ats/sashimi

ref = "/mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall.txt"
target = "/mnt/raid61/Personal_data/zhangyiming/ePSI/ats/overall1.txt"

ref = "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_mm.txt"
target = "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_mm1.txt"


ref = "/mnt/raid64/ATS/Personal/zhangyiming/nslc_ATS.bed"
target = "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_nslc.txt"

with open(target, "w+") as w:
    with open(ref) as r:
        header = None
        for line in r:
            if line.startswith("OrderedCollections"):
                line = line.strip().split("(")[1]
                line = line.replace(")", "")

                temp = {}

                for t in line.split('","'):
                    t = t.replace('"', '')
                    t = t.split(" => ")
                    temp[t[0]] = t[1]
                
                w.write("\t".join([temp[x] for x in header]) + "\n")
            else:
                header = line.strip().split("\t")
                w.write("\t".join(header) + "\n")


