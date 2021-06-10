#!/usr/bin/env python3
# -*- utf-8 -*-
u"""
Stats R1 length, R2 length, fragment size and reads region
"""
import os
import json

from glob import glob
from multiprocessing import Pool

import pysam

from tqdm import tqdm


def call(args):
    bam, output = args
    count = 0
    with open(output, "w+") as w:
        with pysam.AlignmentFile(bam) as r:
            for rec in r:
                if not rec.is_paired:
                    continue
                
                try:
                    if rec.is_read1:
                        rec2 = r.mate(rec)

                        sites = [rec.reference_start, rec.reference_end, rec2.reference_start, rec2.reference_end]

                        w.write(f"{rec.reference_length},{rec.query_alignment_length},{rec2.reference_length},{rec2.query_alignment_length},{rec2.reference_start-rec.reference_end},{max(sites) - min(sites)}\n")

                        count += 1
                except Exception as err:
                    # print(err)
                    continue
                
                
                if count > 100000:
                    break
            

def main(input_dir: str, output: str):

    bams = [x for x in glob(os.path.join(input_dir, "*.bam")) if "pos" in os.path.realpath(x)]

    # with open(output, "w+") as w:
    #     json.dump({os.path.realpath(x): x for x in bams}, w , indent=4)

    keys = {x: os.path.basename(os.path.dirname(os.path.dirname(os.path.realpath(x)))).split("-")[0] for x in bams}

    os.makedirs(output, exist_ok = True)

    cmds = []
    for x, y in keys.items():
        cmds.append([
            # f"/mnt/raid64/Covid19_Gravida/apamix/bam/{y}.bam", 
            x,
            os.path.join(output, os.path.basename(x).replace(".bam", ".csv"))
        ])

    # os.path.basename(x).replace(".bam", "")

    with Pool(len(cmds)) as p:
        list(tqdm(p.imap(call, cmds), total = len(cmds)))


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
