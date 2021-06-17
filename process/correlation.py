#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang at 2021.06.16

Core function to calculate co-expression
"""
import os

from typing import List
from scipy.stats import pearsonr, spearmanr

from src.loci import Region
from src.progress import custom_progress


class UTRIndex(object):
    u"""
    as name says
    """

    def __init__(self, row_id: str):
        self.utr = Region.create(row_id.split("_")[0])
        self.start_idx = 0
        self.end_idx = 0

    def __str__(self):
        return str(self.utr)

    def __len__(self):
        return self.end_idx - self.start_idx

    @property
    def chromosome(self) -> str:
        return self.utr.chromosome

    @property
    def start(self) -> int:
        return self.utr.start

    @property
    def end(self) -> int:
        return self.utr.end


class RowIds(object):
    u"""
    class to handle the expression row ids
    """

    def __init__(self):
        self.row_ids = []
        self.utrs = []
        self.expr = []
        self.barcodes = set()

        self.__row_id_idx__ = 0
        self.__last_row_id__ = None
        self.__last_utr__ = None

    def __len__(self):
        return len(self.utrs)

    def add_row_ids(self, row_id: str):
        u"""
        add new row_id to this object
        """
        if self.__last_row_id__ != row_id:
            self.row_ids.append(row_id)
            self.__row_id_idx__ += 1
            self.__last_row_id__ = row_id

        utr = row_id.split("_")[0]

        if self.__last_utr__ != utr:
            if self.__last_utr__:
                self.utrs[-1].end_idx = self.__row_id_idx__ - 1
            
            utrIdx = UTRIndex(row_id)
            utrIdx.start_idx = self.__row_id_idx__ - 1
            self.utrs.append(utrIdx)

            self.__last_utr__ = utr

    def add_expr(self, col_id: str, value: int):
        if self.__row_id_idx__ >= len(self.expr):
            self.expr.append({})
        self.expr[self.__row_id_idx__ - 1][col_id] = value
        self.barcodes.add(col_id)

    def close(self):
        self.utrs[-1].end_idx = self.__row_id_idx__ - 1

    def sort(self):
        u"""
        sort list of utrs
        """
        self.utrs = sorted(self.utrs, key=lambda x: x.utr)

    @staticmethod
    def __create_region_from_row_id__(row_id: str):
        u"""
        as name says

        :params row_id: 1:100-200:+_1:150-151:+
        """
        reg = Region.create(row_id.split("_")[-1])
        reg.original = row_id
        return reg

    def get_row_ids(self, utr: UTRIndex) -> dict:
        u"""
        return row_id with it's index
        """
        res = {}
        for i in range(utr.start_idx, utr.end_idx):
            res[self.__create_region_from_row_id__(self.row_ids[i])] = i

        return res

    def get_expr(self, index: int):
        return [self.expr[index].get(x, 0) for x in self.barcodes]


def __correction_between_ats__(
    first_ats: dict, 
    second_ats: dict, 
    distance: int,
    expr: RowIds,
    corr_func
):
    u"""
    as namy says

    :params first_ats: list of row_id
    :params second_ats: list of row_id
    :params distance: the maxmimum distance between two UTR
    :params expr: all the expression data
    :params barcodes: all the barcodes
    :params corr_func: 
    """
    atslist1 = sorted(first_ats.keys())
    atslist2 = sorted(second_ats.keys())

    i, j = 0, 0 # backup index, current index
    res = []
    for ats1 in atslist1:

        for j in range(i, len(atslist2)):
            ats2 = atslist2[j]

            if ats1.start - ats2.end > distance:
                i = -1
                break
            elif ats2.start - ats1.end > distance:
                continue
            else:
                if i < 0:
                    i = j
                ats1_expr = expr.get_expr(first_ats[ats1])
                ats2_expr = expr.get_expr(second_ats[ats2])

                r, p = corr_func(ats1_expr, ats2_expr)

                res.append(f"{ats1}\t{ats2}\t{r}\t{p}")

    return res


def corr(mtx: str, output: str, distance: int = 1000, pearson: bool = True):
    u"""
    calculate correction of ATS count

    :params mtx: path to count file
    :params output: path to output file
    :params distance: the maxmimum distance between two UTR
    :params pearson: whether to calculate pearson or spearman correction
    """

    corr_func = pearsonr if pearson else spearmanr

    expr = RowIds()

    progress = custom_progress(io = True)
    with progress:
        task_id = progress.add_task("Loading... ", total = os.path.getsize(mtx))
        
        with open(mtx) as r:
            for idx, line in enumerate(r):
                progress.update(task_id, advance=len(str.encode(line)))
                line = line.split()

                col_id = line[1]
                row_id = line[0]

                expr.add_row_ids(row_id)
                expr.add_expr(col_id, value = int(line[2]))
    expr.close()
    # sort and get all UTRs
    expr.sort()

    progress = custom_progress()

    with progress:
        progress.add_task("Computing... ", total=len(expr))

        with open(output, "w+") as w:
            for i in range(len(expr.utrs)):
                last_utr = expr.utrs[i]

                for j in range(i+1, len(expr.utrs)):
                    curr_utr = expr.utrs[j]

                    # if utr is too far away, skip
                    if curr_utr.chromosome != last_utr.chromosome or curr_utr.start - last_utr.end > distance:
                        break

                    for r in __correction_between_ats__(
                        first_ats=expr.get_row_ids(last_utr),
                        second_ats=expr.get_row_ids(curr_utr),
                        distance=distance,
                        expr=expr,
                        corr_func=corr_func
                    ):
                        w.write(r + "\n")

                progress.update(task_id, advance=1)


if __name__ == '__main__':
    pass