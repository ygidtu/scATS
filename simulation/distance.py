#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Stats the distance between inferred sites and benckmarks
"""
import matplotlib

matplotlib.use("Agg")

import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from rich.progress import Progress


def load_tss(tss):
    data = {}

    with open(tss) as r:
        for line in r:
            _, st, en,  name, _, strand = line.strip().split("\t")

            data[name] =  int(st) if strand == "+"  else int(en)
    return data


def load_results(path):
    header = None
    data = {}
    with open(path) as r:
        for line in r:
            line = line.split()
            if not header:
                header = line
                continue
            
            temp = {i: j for i, j in zip(header, line)}

            key = "{}_{}".format(temp.get("gene_name"), temp.get("transcript_name"))
            sites = [int(x) for x in temp.get("infered_sites", "").split(",")]
            data[key] = sites + data.get(key, [])

    return data


def main(tss, output: str, *files):
    u"""stats the distance between inferred sites and simulated sites.
    
    Keyword arguments:
    tss -- path to tss.bed generated by simulation.py
    output -- path to output directory
    files -- path to result files by ats model.
    Return: return_description
    """
    os.makedirs(output, exist_ok=True)
    print("Loading tss")
    tss = load_tss(tss)

    with Progress() as progress:
        task1 = progress.add_task("[red]Processing...", total = len(files))
        res = {
            "file": [],
            "gene": [],
            "dist": []
        }
        for f in files:
            print(f)
            name = os.path.basename(f)
            data = load_results(f)
            # print(list(tss.keys())[:5])
            # print(list(data.keys())[:5])
            keys = sorted(set(tss.keys()) & set(data.keys()))

            task2 = progress.add_task("[cyan]Comparing... ", total=len(keys))
            for key in keys:
                for i in data[key]:
                    res["file"].append(name)
                    res["gene"].append(key)
                    res["dist"].append(tss[key] - i)


                progress.update(task2, advance=1)

            progress.update(task1, advance=1)

    # pandas and draw
    df = pd.DataFrame(res)

    _, ax = plt.subplots(figsize=(12, 6))
    sns.kdeplot(data = df, x="dist", hue="file", ax=ax)

    for x in [-100, 100]:
        ax.axvline(x=x, color='grey', linestyle='--')
    plt.savefig(os.path.join(output, "dist.png"), dpi = 300)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
