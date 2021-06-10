import os
from subprocess import check_call
from multiprocessing import Pool
from tqdm import tqdm
from glob import glob
import pysam

path = "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/test"
output = "/mnt/raid64/ATS/Personal/zhangyiming/STRT_QC"
qualimap = "/mnt/raid64/ATS/Personal/zhangyiming/software/qualimap_v2.2.1/qualimap"
gtf = "/mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.gtf"
fq = "/mnt/raid64/ATS/rawdata/STRT/PRJNA394919/"


# files = []
# with open(path) as r:
#     for line in r:
#         files.append(line.strip().split())


# os.makedirs(output, exist_ok = True)


# def qc(args):
#     sample, path = args

#     outdir = os.path.join(output, sample)
#     check_call(f"{qualimap} rnaseq -bam {path} -gtf {gtf} -outdir {outdir}", shell=True)



# with Pool(1) as p:
#     list(tqdm(p.imap(qc, files), total = len(files)))


fastq = {os.path.basename(x).split("_")[0]: x for x in glob(os.path.join(fq, "*.fastq.gz"))}

# print(fastq)

with open("stats.csv", "w+") as w:
    for b in glob(os.path.join(output, "*/rnaseq_qc_results.txt")):
        print(b)
        key = os.path.basename(os.path.dirname(b))

        total = 0
        non_uniq = 0
        with open(b) as r:
            for line in r:
                if "reads aligned" in line:
                    line = line.strip().split("=")[-1].strip().replace(',', '')
                    total = int(line)

                if "non-unique alignments" in line:
                    line = line.strip().split("=")[-1].strip().replace(',', '')
                    non_uniq = int(line)
        
        fqc = 0
        with pysam.FastxFile(fastq[key]) as r:
            for line in r:
                fqc += 1
        w.write(f"{key},{fqc},{total},{non_uniq}\n")
            

