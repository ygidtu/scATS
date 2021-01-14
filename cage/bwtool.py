#!/usr/bin/env python3
# -*- utf-8 -*-
u"""
Created at 2021.01.08
"""
import gzip
import os

from glob import glob
from subprocess import check_call, check_output
from multiprocessing import Pool

import numpy as np
import pyBigWig

from tqdm import tqdm


bwtool = "/mnt/raid61/Personal_data/zhangyiming/software/bwtool/bwtool"


def format_peaks(path: str):
    sites = []

    r = gzip.open(path, "rt") if path.endswith("gz") else open(path)
    for line in r:
        line = line.replace(":", "\t").replace(";", "\t")
        line = line.split()

        if "-" in line[1]:
            site = [int(x) for x in line[1].split('-')]
            
            if site[1] - site[0] > 200:
                continue
            
            # site = site[1] if line[2] == "-" else site[0]
            site = int((site[0] + site[1]) /  2)
        else:
            site = int((int(line[1]) + int(line[2])) / 2)
        sites.append([line[0], site])
        # w.write("{}\t{}\t{}\t.\t0\t{}\n".format(
        #     line[0], site, site, line[2]
        # ))
    r.close()
    return sites


def load(args):
    bw, sites, expand = args
    bw = pyBigWig.open(bw)

    temp = np.zeros(expand * 2)
    total = 0
    for site in sites:
        try:
            temp += np.nan_to_num(np.array(bw.values(site[0], site[1] - expand, site[1] + expand)))
            total += 1
        except Exception as err:
            # print(site, err)
            try:
                temp += np.nan_to_num(np.array(bw.values("chr" + site[0], site[1] - expand, site[1] + expand)))
                total += 1
            except Exception as err:
                # print(err)
                continue
    # print(temp)
    return temp, total

    
def main(input_dir: str, peaks: str, output: str, expand: int  = 2000, n_jobs: int = 10):
    
    bws = {}
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.lower().endswith("bigwig") or f.lower().endswith("bw"):
                bws[os.path.basename(parent)] = os.path.join(parent, f)

    sites = format_peaks(peaks)
    res  = {}
    bk =  len(sites) // n_jobs
    for key, bw in tqdm(bws.items()):
        if key == "CTCF":
            continue
        # bw = pyBigWig.open(bw)
        
        total =  0
        temp = np.zeros(expand * 2)
        cmds = [[bw, sites[x: x+bk], expand] for x in range(0, len(sites), bk)]
        with Pool(n_jobs) as p:
            for val, n in list(tqdm(p.imap(load, cmds), total=len(cmds), desc=key)):
                if val is not None:
                    # print(t)
                    temp += val
                    total += n
        # print(temp)
        # print(total)

        if total > 0:
            temp = temp / total

        # for r in tqdm(sites):
        #     try:
        #         temp += np.array(bw.values(r[0], max(0, r[1] - expand), r[1] + expand))
        #     except Exception as err:
        #         try:
        #             temp += np.array(bw.values("chr" + r[0], max(0, r[1] - expand), r[1] + expand))
        #         except Exception as err:
        #             continue

        # temp = temp / len(sites)
        res[key] = temp

    with open(output, "w+") as w:
        keys = sorted(res.keys())

        w.write("\t".join(keys) + "\n")

        for i in range(expand * 2):
            temp = [str(res[key][i]) for key in keys]
            w.write("\t".join(temp) + "\n")


if __name__ == '__main__':
    from fire import Fire
    Fire(main)

