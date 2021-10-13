#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Functions used to count the single cell level expression and calculate PSI
"""
import gzip

from multiprocessing import Process, Queue
from typing import Dict, List

import pysam

from src.loci import BED
from src.logger import log
from src.progress import custom_progress
from src.reader import check_bam, load_ats, load_gtf
from src.expression import Expr


class Bam(object):
    u"""
    Handle the bam related processing
    """

    def __init__(self, path: str, label: str = None):
        u"""
        init this class with bam path and file label
        """
        if not check_bam(path):
            log.error(f"{path} is not a valid bam file")
            exit(1)
        self.path = path

        self.label = "" if label is None else label
    
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

    def reads_bulk(self, region: BED):
        u"""
        generator to get all cell barcode and umi barcode from specific region
        :params region: target region
        :return the nubmer of reads
        """
        count = 0
        with pysam.AlignmentFile(self.path) as r:
            for record in r.fetch(region.chromosome, region.start, region.end):
                if record.is_qcfail or record.is_unmapped:
                    continue

                if record.has_tag("NH") and record.get_tag("NH") > 1:
                    continue
                    
                count += 1
        return count


def count_consumer(bam_files: List[Bam], input_queue: Queue, output_queue: Queue, gtf: bool):
    u"""
    Count ATS
    """
    while True:
        data = input_queue.get()
        utr, regions = data

        res = {}
        for r in regions:
            row_id = f"{utr.to_str()}_{r.to_str()}" if not gtf else f"{utr}-{r.name}"

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                for cb, ub in b.reads(region = r):
                    col_id = f"{cb}_{b.label}" if b.label else cb

                    if col_id not in res[row_id].keys():
                        res[row_id][col_id] = set()
                    
                    res[row_id][col_id].add(ub)
        
        output_queue.put(res)


def count_bulk_consumer(bam_files: List[Bam], input_queue: Queue, output_queue: Queue, gtf: bool):
    u"""
    Count ATS from bulk bam
    """
    while True:
        data = input_queue.get()
        utr, regions = data

        res = {}
        for r in regions:
            row_id = f"{utr.to_str()}_{r.to_str()}" if not gtf else f"{utr}-{r.name}"

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                res[row_id][b.label] = b.reads_bulk(r)
        
        output_queue.put(res)


def count(bams: Dict, output: str, ats: str, processes: int = 1, bulk: bool = False) :
    u"""
    count 
    """
    bam_files = []
    for i, j in bams.items():
        bam_files.append(Bam(i, j))

    gtf = ats.endswith("gtf")

    ats = load_gtf(ats) if gtf else load_ats(ats)

    input_queue = Queue()
    output_queue = Queue()

    # generate consumers
    consumers = []
    for _ in range(processes):
        p = Process(
            target=count_consumer if not bulk else count_bulk_consumer,
            args=(
                bam_files,
                input_queue,
                output_queue,
                gtf,
            )
        )
        p.daemon = True
        p.start()
        consumers.append(p)

    for utr, regions in ats.items():
        input_queue.put([utr, regions])

    progress = custom_progress()

    with progress:
        w = gzip.open(output, "wt+") if output.endswith(".gz") else open(output, "w+")
        
        task = progress.add_task("Counting...", total=len(ats))

        while not progress.finished:
            res = output_queue.get()

            for row, data in res.items():
                for col, val in data.items():
                    if not bulk:
                        val = len(val)
                        
                    if val:
                        w.write(f"{row}\t{col}\t{val}\n")
                    w.flush()
            progress.update(task, advance=1)

        w.close()


def psi(mtx: str, output: str):
    u"""
    Calculate PSI base on count matrix
    :params mtx: count matrix
    """

    # cell -> utr -> sum
    expr = Expr.get_psi(mtx)
        
    w = gzip.open(output, "wt+") if output.endswith(".gz") else open(output, "w+")
    for i in expr:
        w.write("\t".join([str(x) for x in i]) + "\n")
    
    w.close()



if __name__ == '__main__':
    pass
