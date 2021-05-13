#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.08 by Zhang
This is used for debugging
"""

import random
from math import e, log
from multiprocessing import Process

import click
import matplotlib
import numpy as np
import scipy.stats as stats
from KDEpy import FFTKDE
from pyitlib import discrete_random_variable as drv
from rich import print
from src.reader import load_reads
from collections import Counter

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read(args):
    bams, i = args
    reads = load_reads(bams, i)


class Consumer(Process):
    
    def __init__(self, bams, queue, out, **kwargs):
        super().__init__()
        self.bams = bams
        self.queue = queue
        self.out = out
        self.kwargs = kwargs

    def run(self):

        while True:
            i = self.queue.get()

            res = load_reads(self.bams, i)
            self.out.put(res)


    def finish(self):
        [b.close() for b in self.bams if not isinstance(b, str)]


def __entropy__(labels, base=None):
    n_labels = len(labels)

    if n_labels <= 1:
        return 0

    probs = labels / n_labels
    ent = 0.

    # Compute entropy
    base = e if base is None else int(base)

    for i in probs:
        if i > 0:
            ent -= i * log(i, base)

    return ent


def entropy1(mtx, base=None):
    """ Computes entropy of label distribution. """

    return np.array([__entropy__(x, base=base) for x in mtx])


@click.command()
def test():
    
    def normpdf(x, mu, sigma, log = False):
        u = (np.array(x)-mu)/np.abs(sigma)
        y = (1/(np.sqrt(2*np.pi)*np.abs(sigma)))*np.exp(-u*u/2)
        return np.log(y) if log else y

    print(normpdf([20, 10, 15], 20, 10, log = True))
    print(stats.norm(20, 10).logpdf([20, 10, 15]))

    # gaussian kde
    random.seed(42)
    st_arr = random.sample(range(0, 100), 50)
    
    x_arr = np.arange(-100, 200)  # extend to include peaks in 0 or L-1
    kernel = stats.gaussian_kde(st_arr)
    y1 = kernel.pdf(x_arr)

    kernel = FFTKDE(kernel = "gaussian", bw='silverman').fit(st_arr)
    x2, y2 = kernel.evaluate()


    # plt.plot(x_arr, st_arr, label = "raw")
    plt.plot(x_arr, y1, label = "scipy")
    plt.plot(x2, y2, label = "KDE")
    plt.legend()
    plt.savefig("test.png", dpi = 600)

    ## entropy
    Z = np.array([[0, 1, 2], [3, 4, 5]])

    print(stats.entropy(Z[0], base = 10))
    print(__entropy__(Z[0], base = 10))


if __name__ == '__main__':
    pass
