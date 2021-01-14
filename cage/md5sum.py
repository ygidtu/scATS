#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.01.07
"""
import hashlib
import json
import os

from multiprocessing import Pool

from tqdm import tqdm


def md5sum(filename: str):
    with open(filename, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(8192)

    return {"filename": filename, "md5": file_hash.hexdigest()}


def main(input_dir: str, n_jobs: int = 10):
    fs = []
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith("fq.gz") or f.endswith("fastq.gz"):
                fs.append(os.path.join(parent, f))

    with Pool(n_jobs) as p:
        res = list(tqdm(p.imap(md5sum, fs), total=len(fs)))

    with open(os.path.join(input_dir, "md5sum.json"), "w+") as w:
        json.dump(res, w, indent = 4)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
