#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Train a RandomForestClassifier
"""

import gzip
import pickle
import random
from multiprocessing import Pool

import click
import numpy as np
import pyBigWig
from pybedtools import BedTool
from rich.progress import track
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

from src.logger import log


class Bed(object):
    def __init__(self, chrom: str, start: int, end: int, tag: int = 0):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.data = {}
        self.tag = tag

    @classmethod
    def create(cls, record: str):

        elements = record.split()
        return cls(chrom=elements[0], start=int(elements[1]), end=int(elements[2]), tag=1)

    @classmethod
    def convert(cls, record):

        return cls(chrom=record.chrom, start=record.start, end=record.end, tag=0)

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    def __hash__(self) -> int:
        return hash((self.chrom, self.start, self.end))

    def __eq__(self, other) -> bool:
        return self.chrom == other.chrom and \
               self.start == other.start and \
               self.end == other.end

    def __lt__(self, other) -> bool:
        if self.chrom != other.chrom:
            return self.chrom < other.chrom

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end

    def __gt__(self, other) -> bool:
        if self.chrom != other.chrom:
            return self.chrom > other.chrom

        if self.start != other.start:
            return self.start > other.start

        return self.end > other.end

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    def load_signal(self, chroms: str, span: int = 500, key: str = None):

        try:
            self.data[key] = chroms[self.chrom][(self.center - span):(self.center + span)]
        except Exception as err:
            log.debug(err)
            self.data[key] = np.zeros(span * 2)

        if not self.data[key]:
            self.data[key] = np.zeros(span * 2)

    def check(self, path) -> bool:
        with pyBigWig.open(path) as r:
            return r.chroms(self.chrom) is not None


def load(args):
    path, chrom = args

    with pyBigWig.open(path) as r:
        for x, y in r.chroms().items():
            if x == chrom:
                return {x: r.values(x, 0, y)}

    return {}


def smaller(a, b, span=500):
    if a.chrom != b.chrom:
        return a.chrom < b.chrom
    return b.start - a.end >= span


def bigger(a, b, span=500):
    if a.chrom != b.chrom:
        return a.chrom > b.chrom

    return a.start - b.end > span


CLASSIFIERS = {
    "knn": {
        "classifier": KNeighborsClassifier,
        "params": {
            "n_neighbors": [3, 4, 5, 6, 7, 8, 9, 10],
            "weights": ["uniform", "distance"],
            "algorithm": ["ball_tree", "kd_tree", "brute"],
            "leaf_size": [1, 2],
            "metric": ["minkowski", "euclidean", "manhattan", "chebyshev", "wminkowski", "seuclidean", "mahalanobis"]
        }
    },
    "rf": {
        "classifier": RandomForestClassifier,
        "params": {
            'n_estimators': [10, 20, 30, 40, 50, 75, 100, 200, 300, 400, 500],
            'max_features': ['sqrt', 'log2'],
            'criterion': ['gini', 'entropy', "log_loss"]
        }
    },
    "svc": {
        "classifier": SVC,
        "params": {
            "kernel": ["linear", "poly", "rbf", "sigmoid", "precomputed"],
            "degree": [3, 4, 5, 6],
            "gamma": ["scale", "auto"],
            "decision_function_shape": ["ovo", "ovr"]
        }
    },
    "dt": {
        "classifier": DecisionTreeClassifier,
        "params": {
            'criterion': ['gini', 'entropy', 'log_loss'],
            "splitter": ["best", "random"],
            "max_features": ["sqrt", "log2"]
        }
    }
}


@click.command()
@click.option('-i', '--atac', type=click.Path(exists=True), required=True, help='The path to ATAC bigwig file.')
@click.option('-c', '--cage', type=click.Path(exists=True), required=True, help='The path to CAGE peaks.')
@click.option('-o', '--output', type=click.Path(), required=True, help='The path to output file.')
@click.option('-g', '--genome', type=click.Path(), required=True,
              help='the genome version of CAGE peaks, eg: hg19, hg38, mm10.')
@click.option('-n', '--n-jobs', type=int, default=1, help='The number of processes to use.')
@click.option('-s', '--sample-size', type=float, default=1, help='0~1 to control the data size used for training.')
@click.option('-t', '--test-size', type=float, default=.3, help='the size of test dataset.')
@click.option('--classifier', type=click.Choice(["knn", "rf", "svc", "dt"]), default="rf", help='the default classifier')
def train(atac: str, cage: str, output: str,
          genome: str = "mm10", n_jobs: int = 1,
          sample_size: float = 1, test_size: float = .3, classifier: str="rf"):
    u"""
    Train machine learning classifier.

    :param atac: path to ATAC bigwig file
    :param cage: path to CAGE peaks
    :param output: path to output path of model
    :param genome: the genome version of CAGE peaks, eg: hg19, hg38, mm10
    :param n_jobs: the number of processes used to train model
    :param sample_size: 0~1 to control the data size used for training
    :param test_size: the size of test dataset
    """
    seed = 42
    log.info(f"Loading ATAC signals from {atac}")
    chroms = {}
    with pyBigWig.open(atac) as r:
        with Pool(n_jobs) as p:
            for x in p.map(load, [[atac, x] for x in r.chroms()]):
                chroms.update(x)

    log.info(f"Loading CAGE peaks from {cage}")
    cages = set()
    if cage.endswith(".gz"):
        r = gzip.open(cage, "rt")
    else:
        r = open(cage, "r")

    for line in r:
        b = Bed.create(line)
        b.tag = 1

        if b.chrom in chroms.keys():
            cages.add(b)
    r.close()

    log.info(f"Generating random peaks for {genome}")
    x = BedTool()
    y = x.random(l=1, n=len(cages), genome=genome, seed=seed)
    rand = [Bed.convert(x) for x in y if x.chrom in chroms.keys()]

    log.info(f"Merging peaks")
    res = []
    cages = sorted(cages)
    rand = sorted(rand)

    i, j = 0, 0
    while i < len(cages) and j < len(rand):
        bi = cages[i]
        bj = rand[j]

        if smaller(bi, bj):
            res.append(bi)
            i += 1
        elif bigger(bi, bj):
            res.append(bj)
            j += 1
        else:
            j += 1

    for b in cages[i:]:
        res.append(b)

    for b in rand[j:]:
        res.append(b)

    if 0 < sample_size < 1:
        random.seed(seed)
        res = random.sample(res, int(len(res) * sample_size))

    log.info("Prepare sample data")
    for r in track(res):
        r.load_signal(chroms, key="1")

    X, y = [], []
    for r in track(res):
        if len(r.data["1"]) == 1000:
            data = np.array(r.data["1"], dtype="float")
            data[np.isnan(data)] = 0
            X.append(np.array(data))
            y.append(r.tag)

    log.info("Grid search for best parameters")
    X_train, X_test, y_train, y_test = train_test_split(np.array(X), np.array(y), test_size=test_size)

    accuracy_score = -1
    classifier, param_grid = CLASSIFIERS[classifier]["classifier"], CLASSIFIERS[classifier]["params"]

    rfc = classifier(n_jobs=n_jobs)
    cv_rfc = GridSearchCV(estimator=rfc, param_grid=param_grid, cv=5)
    cv_rfc.fit(X_train, y_train)
    log.info(f"Best params: {cv_rfc.best_params_}")

    log.info("Generating classifier using best parameters")
    rfc = classifier(**cv_rfc.best_params_)
    rfc.fit(X_train, y_train)
    y_pred = rfc.predict(X_test)

    if metrics.accuracy_score(y_test, y_pred) > accuracy_score:
        log.info("Accuracy:", metrics.accuracy_score(y_test, y_pred))

        with open(output, "wb+") as w:
            pickle.dump(rfc, w)


if __name__ == '__main__':
    train()
