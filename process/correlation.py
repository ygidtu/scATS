#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang at 2021.06.16

Core function to calculate co-expression
"""
import os

from scipy.stats import pearsonr, spearmanr

from src.expression import Expr
from src.progress import custom_progress


def __correction_between_ats__(
    first_ats: dict, 
    second_ats: dict, 
    distance: int,
    expr: Expr,
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

    expr = Expr()

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

    progress = custom_progress()

    with progress:
        progress.add_task("Computing... ", total=len(expr))

        with open(output, "w+") as w:
            for i in range(len(expr)):
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