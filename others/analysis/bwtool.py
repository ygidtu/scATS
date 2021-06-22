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
    r.close()
    return sites


def load(args):
    bw, sites, expand = args
    try:
        bw = pyBigWig.open(bw)
    except Exception as err:
        print(err)
        print(args[0])
        exit(1)

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
    return temp, total

    
def peaks(input_dir: str, peaks: str, output: str, expand: int  = 2000, n_jobs: int = 10):
    u"""
    Count the histone markers signal around peaks

    :param input_dir path to dirctory contains bigwig files, path/histone names/sample.bigwig
    :param peaks path to peaks file that output by cage.R
    :param output path to output file
    :param expand how many bp to expand based on the center of peaks
    :param n_jobs how many processes to use
    """
    output = os.path.abspath(output)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    bws_ss = {"+": {}, "-": {}, "ns": {}}
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.lower().endswith("bigwig") or f.lower().endswith("bw"):
                key = os.path.basename(parent)
                if key in bws_ss.keys():
                    bws_ss[key][os.path.basename(os.path.dirname(parent))] = os.path.join(parent, f)
                else:
                    bws_ss["ns"][key] = os.path.join(parent, f)

    sites = format_peaks(peaks)
    res  = {}
    bk =  len(sites) // n_jobs

    for ss, bws in bws_ss.items():
        if not bws:
            continue

        for key, bw in tqdm(bws.items(), desc = ss):
            total =  0
            temp = np.zeros(expand * 2)
            cmds = [[bw, sites[x: x+bk], expand] for x in range(0, len(sites), bk)]
            with Pool(n_jobs) as p:
                for val, n in list(tqdm(p.imap(load, cmds), total=len(cmds), desc=key)):
                    if val is not None:
                        temp += val
                        total += n

            if total > 0:
                temp = temp / total
            res[key] = temp

        o = output + "_" + (ss if ss != "ns" else "")
        with open(o, "w+") as w:
            keys = sorted(res.keys())

            w.write("\t".join(keys) + "\n")

            for i in range(expand * 2):
                temp = [str(res[key][i]) for key in keys]
                w.write("\t".join(temp) + "\n")
        

def format_model(path: str, shift: int = 0, ws_t: float = 0.0):
    sites = []

    with open(path) as r:
        for line in r:
            if line.startswith("utr"):
                continue
            line = line.strip().split("\t")

            if not line[5] and line[5] != "NA":
                continue

            site = line[0].split(":")

            try:
                ss = [int(float(x)) for x in line[4].split(",")]
                ws = [float(x) for x in line[1].split(",")]

                ss = [x for x, y in zip(ss, ws) if y > ws_t]

                # if len(ss) < 2:
                #     continue

                if len(ss) % 2 == 0:
                    mid = len(ss) // 2
                    ss[:mid] = [x + shift for x in ss[:mid]]
                    ss[mid:] = [x - shift for x in ss[mid:]]
                else:
                    mid = len(ss) // 2
                    ss[:(mid-1)] = [x + shift for x in ss[:(mid-1)]]
                    ss[(mid+1): ] = [x - shift for x in ss[(mid+1):]]

                for i in ss:
                    sites.append([site[0], i])
            except ValueError:
                continue

    return sites


def model(input_dir: str, model: str, output: str, expand: int  = 2000, n_jobs: int = 10, shift: int = 0):
    u"""
    Count the histone markers signal around peaks

    :param input_dir path to dirctory contains bigwig files, path/histone names/sample.bigwig
    :param model path to output file of ATSmix
    :param output path to output file
    :param expand how many bp to expand based on the center of peaks
    :param n_jobs how many processes to use
    """
    output = os.path.abspath(output)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    bws = {}
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f.lower().endswith("bigwig") or f.lower().endswith("bw"):
                bws[os.path.basename(parent)] = os.path.join(parent, f)

    sites = format_model(model, shift = shift)
    res  = {}
    bk =  len(sites) // n_jobs
    for key, bw in tqdm(bws.items()):

        total =  0
        temp = np.zeros(expand * 2)
        cmds = [[bw, sites[x: x+bk], expand] for x in range(0, len(sites), bk)]
        with Pool(n_jobs) as p:
            for val, n in list(tqdm(p.imap(load, cmds), total=len(cmds), desc=key)):
                if val is not None:
                    temp += val
                    total += n

        if total > 0:
            temp = temp / total
        res[key] = temp

    with open(output, "w+") as w:
        keys = sorted(res.keys())

        w.write("\t".join(keys) + "\n")

        for i in range(expand * 2):
            temp = [str(res[key][i]) for key in keys]
            w.write("\t".join(temp) + "\n")



if __name__ == '__main__':
    from fire import Fire
    Fire({"peaks": peaks, "model": model})
