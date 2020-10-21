#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2020.09.10 by Zhang Yiming

This script is used to test the connect of python scripts and R
"""

import os
import json
import pickle
import re

from multiprocessing import Pool
from subprocess import check_call

import numpy as np

from tqdm import tqdm

from apamix.apamix import APA

from loguru import logger

logger.add('apamix.log',
            rotation='10 MB',
            colorize=True,
            level="INFO")


__dir__ = os.path.abspath(os.path.dirname(__file__))
ERROR = os.path.join(__dir__, "errors_cupy")
os.makedirs(ERROR, exist_ok=True)


class AFELoc:

    def __init__(self, chrom: str, start: int, end: int, utr_start: int, utr_end: int, start_in_utr: int, len_in_utr: int, barcode: str = ""):
        u"""
        init this object
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.utr_start = utr_start
        self.utr_end = utr_end
        self.start_in_utr = start_in_utr
        self.len_in_utr = len_in_utr
        self.barcode = barcode

    @classmethod
    def create(cls, data: str):
        u"""
        Class method
        create AFELoc object based on \t seperated string
        """
        info = data.strip().split()
        if not info:
            return None
        return AFELoc(
            chrom=info[1], start=int(info[4]), end=int(info[5]),
            utr_start=int(info[2]), utr_end=int(info[3]),
            start_in_utr=int(info[7]), len_in_utr=int(info[8]),
            barcode=info[0]
        )

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))

    def __str__(self) -> str:
        return f"{self.chrom}\t{self.utr_start}\t{self.utr_end}\t{self.start}\t{self.end}\t{self.start_in_utr}\t{self.len_in_utr}"

    def to_key(self) -> str:
        return f"{self.chrom}:{self.utr_start}-{self.utr_end}"

# @profile
def __single_process__(args):
    # same_start, same_end, key, start_in_utr = args
    res = {}
    error = []

    for key, common_records, start, end in args:

        L = int(np.mean([x.utr_end - x.utr_start for x in common_records]))

        r1_utr_st_arr, polya_len_arr = np.full((1, len(common_records)), L - 500)[0], np.full((1, len(common_records)), np.NaN)[0]
        r1_len_arr = []

        output_key = ""
        for x in common_records:
            r1_len_arr.append(x.len_in_utr)
            output_key += f"{x.start},{x.end},"

        try:
            test = APA(n_max_apa=5,
                n_min_apa=1,
                r1_utr_st_arr=r1_utr_st_arr,
                r1_len_arr=r1_len_arr,
                r2_len_arr=polya_len_arr,
                polya_len_arr=polya_len_arr,
                pa_site_arr=polya_len_arr,
                utr_len=np.max(r1_utr_st_arr) + np.max(r1_len_arr) + 300,
                verbose=False
            )
            # logger.debug('init infer')
        except Exception as err:
            error.append({
                "r1_utr_st_arr": r1_utr_st_arr, "polya_len_arr": polya_len_arr,
                "r1_len_arr": r1_len_arr, "chrom": key, "start": start, "end": end,
                "error": str(err)
            })
      
            continue

        a = test.inference()

        if a and all([x is not None for x in a.values()]):
            res[f"{key}:{output_key.strip(',')}"] = a

    return res


def __assign_tasks__(args):
    data = []
    same_start, same_end, key, start_in_utr = args
    for start, records in same_start[key].items():
        if len(records) <= 1:
            continue
        
        ends = set()
        for record in records:
            if record.end in ends:
                continue

            ends.add(record.end)

            same_end_records = same_end[key].get(record.end, [])

            common_records = [x for x in list(set(records) & set(same_end_records)) if x.start_in_utr > start_in_utr]
            if len(common_records) <= 1:
                continue

            data.append([key, common_records, start, record.end])
    return data


#    main('../../tests/extract/C141.txt', '../../tests/extract/C141_pypy.json')
    
def main(input_file:str='../../tests/extract/C141.txt', output_file: str='../../tests/extract/C141_pypy.json', start_in_utr: int=-500, n_jobs=5):
    data = {}
    with open(input_file) as r:
        for idx, line in enumerate(r):
            if idx:
                temp = AFELoc.create(line)

                if not temp:
                    continue
           
                # if idx > 10000:
                #     break

                if temp.barcode not in data.keys():
                    data[temp.barcode] = {
                        "start": {},
                        "end": {}
                    }

                chrom = temp.chrom

                temp_same_start = data[temp.barcode]["start"].get(chrom, {})
                temp_same_end = data[temp.barcode]["end"].get(chrom, {})

                temp_start = temp_same_start.get(temp.utr_start, set())
                temp_start.add(temp)
                temp_same_start[temp.utr_start] = temp_start

                temp_end = temp_same_end.get(temp.utr_end, set())
                temp_end.add(temp)
                temp_same_end[temp.utr_end] = temp_end

                data[temp.barcode]["start"][chrom] = temp_same_start
                data[temp.barcode]["end"][chrom] = temp_same_end

    temp_output = output_file + "_temp"
    os.makedirs(temp_output, exist_ok = True)

    for bc, temp_data in data.items():
        temp_json = os.path.join(temp_output, bc + ".json")

        if os.path.exists(temp_json):            
            try:
                with open(temp_json) as r:
                    for line in r:
                        continue
                continue
            except EOFError:
                pass

        res, data_mp, cmds = {}, [], []
  
        same_start, same_end = temp_data["start"], temp_data["end"]
        for key in tqdm(set(same_start.keys()) & set(same_end.keys()), "Manage data"):
            cmds.append([same_start, same_end, key, start_in_utr])

        with Pool(n_jobs) as p:
            temp = list(tqdm(p.imap_unordered(__assign_tasks__, cmds), total=len(cmds)))
        p.close()
        p.join()
    
        for t in temp:
            data_mp += t

        cmds = []
        bk = len(data_mp) // n_jobs

        # print(data)
        if bk <= 1:
            cmds.append(data_mp)
        else:
            for i in range(0, len(data_mp), bk):
                cmds.append(data_mp[i: i+bk])

        with Pool(n_jobs) as p:
            temp = list(tqdm(p.imap_unordered(__single_process__, cmds), total=len(cmds)))
        p.close()
        p.join()

        for t in temp:
            res.update(t)

        if not res:
            continue

        with open(temp_json, "w+") as w:
            json.dump(res, w, indent=4)

    # with open(output_file, "w+") as w:
    #     header = ["sites"]
    #     v = list(res.values())
    #     if v:
    #         v = v[0]
    #     else:
    #         exit(0)
    #     header += list(v.keys())

    #     w.write("\t".join(header) + "\n")

    #     for key, vals in res.items():
    #         row = [key]
    #         for  i in header[1:]:
    #             val = vals[i]
    #             if isinstance(val, float):
    #                 row.append(str(val))
    #             else:
    #                 row.append(",".join([str(x) for x in val]))
    #         w.write("\t".join(row) + "\n")


def sashimi(in_json, output, bam, gtf="/mnt/raid61/Personal_data/zhangyiming/genome/Homo_sapiens/Homo_sapiens.GRCh38.93.gtf"):
    SCRIPT = "/mnt/raid61/Personal_data/zhangyiming/code/pysashimi/main.py"

    with open(in_json) as r:
        data = json.load(r)
    
    os.makedirs(output, exist_ok=True)
    for key, value in data.items():
        key = re.split(r"[:,]", key)
        line = ",".join(key[1:])
        keys = sorted([int(x) for x in key[1:]])
        keys[0] = max(0, int(keys[0]) - 100)
        keys[-1] = int(keys[-1]) + 100
        key = f"{key[0]}:{keys[0]}-{keys[-1]}"
        check_call(f"python3 {SCRIPT} normal -b {bam} -g {gtf} -o {output}/{key}.pdf -e {key}:+ --title {value['bic']} -p 10 --indicator-lines {line}", shell=True)


def batch(input_dir, gtf, output_dir="./", n_jobs=10):
    utr = os.path.join(output_dir, "ref_utr.bed")
    if not os.path.exists(utr):
        print("prepare utr bed")
        check_call(f"python {__dir__}/main.py prepare -g {gtf} -p {output_dir}/ref_", shell=True)

    fs = []
    for parent, _, files in os.walk(input_dir):
        for f in files:
            if f == "possorted_genome_bam.bam":
                key = ""
                for i, j in enumerate(parent.split("/")):
                    if j == "outs":
                        key = parent.split("/")[i-1]
                        break
                print(key)
                fs.append(os.path.join(output_dir, "extract", key + ".txt"))
                check_call(f"python {__dir__}/main.py extract -u {utr} -b {os.path.join(parent, f)} -p {n_jobs} -o {fs[-1]}", shell=True)
                main(fs[-1], os.path.join(output_dir, key + ".txt"), -500, n_jobs)


if __name__ == '__main__':
    # try :
    #     import profile
    # except:
    #     import cProfile as profile

    # profile.run("main('../../tests/extract/C141.txt', '../../tests/extract/C141_cupy.json')")
 
    from fire import Fire
    Fire()