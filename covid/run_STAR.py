import os
import re
from glob import glob
from subprocess import check_call

import pandas as pd

from tqdm import tqdm


fq = "/mnt/raid64/BetaCoV_sc/Coivd19_amit/fq"
SOFTLINKS = "/mnt/raid64/Covid19_Gravida/cellranger/softlinks"
INDEX = "/mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/star_index"
HISAT = "/mnt/data8/zhangyiming/software/hisat2-2.2.1"
HISAT_INDEX = "/mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/hisat2_index/hisat2_index"
META = "/mnt/raid64/BetaCoV_sc/Coivd19_amit/meta.tsv"
OUTPUT = "/mnt/raid64/Covid19_Gravida/STAR/"
BAMS = "/mnt/raid64/ATS/Personal/zhangyiming/bams/"
GENOME = "/mnt/raid62/BetaCoV/Person/tangchao/data/reference/MN908947.genome.fasta"


def main():
    meta = pd.read_csv(META, header = None, sep = "\t")

    sample = []
    for _, row in meta.iterrows():
        key = row[0].split(":")[1].split(" (")[0]
        key = key.split(", ")[-1].strip()

        fs = glob(os.path.join(fq, key + "*.fastq.gz"))

        row[1] = row[1].strip().strip("%")
        
        if len(fs) < 1:
            fs = glob(os.path.join(fq, row[1] + "*.fastq.gz"))

            if len(fs) < 1:
                fs = glob(os.path.join(fq, row[1], "*.fastq.gz"))

        sample.append([key, row[1], row[2], fs])

    for row in tqdm(sample):
        check_call(f"STAR --runThreadN 12 --outSAMtype BAM SortedByCoordinate --outBAMcompression 9 --genomeDir {INDEX} --readFilesCommand zcat --readFilesIn {' '.join(sorted(row[3]))} --outFileNamePrefix {OUTPUT}/{row[0]}_{row[2]}_", shell=True)


def covid19():
    bams = [x for x in glob(os.path.join(BAMS, "*.bam")) if os.path.islink(x)]

    fqs = {}
    for x in glob(os.path.join(SOFTLINKS, "*.fastq.gz")):
        key = os.path.basename(x).split("_")[0].split("-")[0]
        temp = fqs.get(key, [])
        temp.append(x)
        fqs[key] = temp

    for x in tqdm(bams):
        ckey = os.path.basename(x).replace(".bam", "")
        okey = re.sub(r"[\d-]", "", os.path.basename(os.path.dirname(os.path.dirname(os.path.realpath(x)))))

        l = os.path.join(OUTPUT, ckey + "_Log.progress.out")
        if os.path.exists(l):
            with open(l) as r:
                for line in r:
                    if "ALL DONE!" in line:
                        continue
        print(x)
        check_call(f"STAR --runThreadN 12 --outSAMtype BAM SortedByCoordinate --outBAMcompression 9 --genomeDir {INDEX} --readFilesCommand zcat --readFilesIn {' '.join(sorted(fqs[okey]))} --outFileNamePrefix {OUTPUT}/{ckey}_", shell=True)


def hisat2():
    if not os.path.exists(os.path.dirname(HISAT_INDEX)):
        check_call(f"{HISAT}/hisat2-build -p 10 {GENOME} {HISAT_INDEX}", shell=True)

    bams = [x for x in glob(os.path.join(BAMS, "*.bam")) if os.path.islink(x)]
    fqs = {}
    for x in glob(os.path.join(SOFTLINKS, "*.fastq.gz")):
        key = os.path.basename(x).split("_")[0].split("-")[0]
        temp = fqs.get(key, [])
        temp.append(x)
        fqs[key] = temp

    for x in tqdm(bams):
        ckey = os.path.basename(x).replace(".bam", "")
        okey = re.sub(r"[\d-]", "", os.path.basename(os.path.dirname(os.path.dirname(os.path.realpath(x)))))

        l = os.path.join(OUTPUT, ckey + "_Log.progress.out")
        if os.path.exists(l):
            with open(l) as r:
                for line in r:
                    if "ALL DONE!" in line:
                        continue
        # print(x)
        f = sorted(fqs[okey])
        print(f)
        cmd = f"{HISAT}/hisat2 -p 10 -x {HISAT_INDEX} -1 {f[0]} -2 {f[1]} | samtools view -b | samtools sort > {OUTPUT}/{ckey}.bam && samtools index {OUTPUT}/{ckey}.bam"
        
        check_call(cmd, shell=True)
 


if __name__ == '__main__':
    from fire import Fire
    Fire()


