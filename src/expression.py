#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Created by Zhang at 2021.06.16

Contains class related to expression matrix
"""
import os

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


class Expr(object):
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

    @classmethod
    def create(cls, mtx: str):
        u"""
        """
        expr = cls()
        progress = custom_progress(io = True)
        with progress:
            task_id = progress.add_task("Loading... ", total = os.path.getsize(mtx))
            
            with open(mtx) as r:
                for line in r:
                    progress.update(task_id, advance=len(str.encode(line)))
                    line = line.split()

                    col_id = line[1]
                    row_id = line[0]

                    expr.add_row_ids(row_id)
                    expr.add_expr(col_id, value = int(line[2]))

        expr.close()
        # sort and get all UTRs
        expr.sort()
        return expr

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
        u"""
        get list of expression values
        """
        return [self.expr[index].get(x, 0) for x in self.barcodes]

    @classmethod
    def get_psi(cls, mtx: str):
        u"""
        """

        last_utr = None

        summaries = {}
        expr = []

        progress = custom_progress(io = True)
        with progress:
            task_id = progress.add_task("PSI... ", total = os.path.getsize(mtx))
            
            with open(mtx) as r:
                for line in r:
                    progress.update(task_id, advance=len(str.encode(line)))
                    line = line.split()

                    row_id, col_id, value = line[0],line[1], int(line[2])
                    utr = row_id.split("_")[0]

                    if utr != last_utr:
                        if last_utr:
                            for row in expr:
                                row[-1] = row[-1] / summaries[row[1]]
                                yield row
                            summaries = {}
                            expr = []

                        last_utr = utr

                    summaries[col_id] = value + summaries.get(col_id, 0)
                    expr.append([row_id, col_id, value])

        for row in expr:
            row[-1] = row[-1] / summaries[row[1]]
            yield row

  

if __name__ == '__main__':
    pass
