#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Functions used to count the single cell level expression and calculate PSI
"""
import pysam

from multiprocessing import Process, Queue
from typing import Dict, List

from rich.progress import track

from src.reader import check_bam, load_ats
from src.loci import BED


class Bam(object):
    u"""
    Handle the bam related processing
    """

    def __init__(self, path: str, label: str = None):
        u"""
        init this class with bam path and file label
        """
        if not check_bam(path):
            log.error(f"{bam} is not a valid bam file: {err}")
            exit(1)
        self.path = path

        self.label = os.path.basename(path) if label is None else label
    
    def reads(self, region: BED, cell_tag: str = "CB", umi_tag: str = "UB"):
        u"""
        generator to get all cell barcode and umi barcode from specific region
        :params region: target region
        :params cell_tag: the cell barcode tag in 10x bam
        :params umi_tag: the umi tag in 10x bam
        :return cell barcode and umi
        """
        with pysam.AlignmentFile(self.path) as r:
            for record in r.fetch(region.chromosome, region.start, region.end):
                if record.is_qcfail or record.is_unmapped:
                    continue

                if record.has_tag("NH") and record.get_tag("NH") > 1:
                    continue
                    
                if record.has_tag(cell_tag) and record.has_tag(umi_tag):
                    yield record.get_tag(cell_tag), record.get_tag(umi_tag)


def count_consumer(bam_files: List[Bam], input_queue: Queue, output_queue: Queue):
    u"""
    Count ATS
    """
    while True:
        data = input_queue.get()
        utr, regions = data

        res = {}
        for r in regions:
            row_id = f"{utr.to_str()}_{}"

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                for cb, ub in b:
                    col_id = f"{cb}_{b.label}"

                    if col_id not in res[row_id].keys():
                        res[row_id][col_id] = set()
                    
                    res[row_id][col_id].add(ub)
        
        output_queue.put(res)


def count(bams: Dict, ats: str, processes: int = 1, output: str) :
    u"""
    count 
    """
    bam_files = []
    for i, j in bams.items():
        bam_files.append(Bam(i, j))

    ats = load_ats(ats)

    # row_ids = []
    # col_ids = set()

    res = {}
    for utr, regions in track(ats.items()):
        for r in regions:
            row_id = f"{utr.to_str()}_{}"
            # row_ids.append(row_id)

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                for cb, ub in b:
                    col_id = f"{cb}_{b.label}"

                    # col_ids.add(col_id)

                    if col_id not in res[row_id].keys():
                        res[row_id][col_id] = set()
                    
                    res[row_id][col_id].add(ub)

    with open(output, "w+") as w:
        for row, data in res.items():
            for col, val in data.items():
                if len(val) > 0:
                    w.write(f"{row}\t{col}\t{len(val)}\n")

    return res


def psi(mtx: str, output: str):
    u"""
    Calculate PSI base on count matrix
    :params mtx: count matrix
    """

    # cell -> utr -> sum
    summurize = {}
    res = {}
    with open(mtx) as r:
        for line in r:
            line = line.split()

            col_id = line[1]
            utr = line[0].split("_")[0]
            row_id = line[0]

            # sum
            if col_id not in summurize.keys():
                summurize[col_id] = {}
            
            if utr not in summurize[col_id].keys():
                summurize[col_id][utr] = 0

            summurize[col_id][utr] += int(line[2])

            # collect data
            if row_id not in res.keys():
                res[row_id] = {}

            res[row_id][col_id] = int(line[2])

    with open(output, "w+") as w:
        for row, data in res.items():
            for col, val in col.items():
                total = summurize[col][row.split("_")[0]]
                w.write(f"{row}\t{col}\t{val / total}\n")


if __name__ == '__main__':
    pass
