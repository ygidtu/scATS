import os
from glob import glob
from subprocess import check_call

indir = "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi//celltype_bc"
order = "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/n_order.txt"
jl = "/mnt/raid61/Personal_data/zhangyiming/code/afe/others/sashimi.jl"
markers = "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/markers_between_stages_psi.rds"
output = "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/sashimi/stage"
gtf = "/mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.sorted.gtf.gz"
bam = "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam3.tsv"

fs = glob(os.path.join(indir, "*_bc.txt"))


for i in fs:
    key = os.path.basename(i).replace("_bc.txt", "")

    o = os.path.join(output, key)

    check_call(f"julia -t 10 {jl} -i {markers} -g {gtf} --bc {i} -b {bam} --order {order} --cell {key} -e 50 -n 30 -o {o}", shell=True)
