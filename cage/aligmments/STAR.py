#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.01.07
"""

import os
import re

from glob import glob
from shutil import rmtree
from subprocess import check_call

from rich import print


REFERENCE = {
    "human": "/mnt/raid61/Personal_data/zhangyiming/genome/Homo_sapiens/STAR_index",
    "mouse": "/mnt/raid61/Personal_data/zhangyiming/genome/Mus_musculus/STAR_index"
}

class SampleInfo(object):

    def __init__(self, ref: str, fid: str):
        self.ref = ref
        self.id = fid
        self.fs = []
        self.n_jobs = 8

    def add(self, path: str):
        if os.path.exists(path):
            self.fs.append(os.path.abspath(path))

    def is_finish(self, cwd: None) -> bool:
        name = f"{self.id}_Log.progress.out"
        if cwd:
            name = os.path.join(cwd, name)
        
        if not os.path.exists(name):
            return False
        
        with open(name) as r:
            for line in r:
                line = line.strip()

                if "ALL DONE!" in line:
                    return True
        return False       

    def __hash__(self):
        return hash(self.id)

    def __str__(self) -> str:

        if len(self.fs) > 2:
            fs = [x for x in self.fs if re.search(r"[12].f(ast)?q.gz", x)]
        else:
            fs = self.fs

        return f"STAR --runThreadN {self.n_jobs} --outSAMtype BAM SortedByCoordinate --outBAMcompression 9 --genomeDir {self.ref} --readFilesCommand zcat --readFilesIn {' '.join(sorted(fs))} --outFileNamePrefix {self.id}_"


def call(sample: SampleInfo, cwd = None):
    check_call(str(sample), cwd=cwd, shell=True)


def collect_files(input_dir: str):

    fs = {}
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith("fastq.gz") and  "RNA" in f:
                key = f.split("_")[0]

                if "human" in f or "Homo_sapiens" in f:
                    ref = REFERENCE["human"]
                else:
                    ref = REFERENCE["mouse"]

                temp = SampleInfo(ref = ref, fid = key)

                if key not in fs.keys():
                    fs[key] = temp
                
                fs[key].add(os.path.join(parent, f))

    return sorted(fs.values(), key=lambda x: x.id)
    
def main(input_dir: str, output_dir: str, n_jobs: int = 10):
    os.makedirs(output_dir, exist_ok=True)

    fs = collect_files(input_dir)

    for f in fs:
        print(f)
        f.n_jobs = n_jobs

        if not f.is_finish(output_dir):
            try:
                call(f, output_dir)
            except Exception as err:
                print(err)

                for i in glob(os.path.join(output_dir, f.id + "*")):
                    if os.path.isfile(i):
                        os.remove(i)
                    else:
                        rmtree(i)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
