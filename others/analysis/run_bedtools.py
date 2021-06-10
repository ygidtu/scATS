import os
from glob import glob

from subprocess import check_call

from multiprocessing import Pool
from tqdm import tqdm


def call(cmd):
    check_call(cmd, shell=True)


def main(input_dir, output, n_jobs=10):
    os.makedirs(output, exist_ok=True)

    bams = []
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f == "possorted_genome_bam.bam":
                bams.append(f"bedtools bamtobed -i {os.path.join(parent, f)} -bedpe -split 1> {os.path.join(output, os.path.basename(os.path.dirname(parent)) + '.bed')} 2> /dev/null")

    with Pool(n_jobs) as p:
        list(tqdm(p.imap(call, bams), total = len(bams)))


if __name__ == '__main__':
    from fire import Fire

    Fire(main)
