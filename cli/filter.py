#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle

from multiprocessing import Pool
from typing import List

import click
import numpy as np
import pyBigWig

from rich.progress import track
from scipy.stats import gaussian_kde, pearsonr

from src.loci import Region
from src.logger import log
from src.reader import load_reads


def func(x, a, b, c):  # Gaussian peak
    return a * np.exp(-0.5 * np.power((x - b) / c, 2.0))


def calculate_gaussian_density_similarity(x: np.array, y: np.array):
    if len(x) < 3:
        return float('nan'), [x, y, "lt3"]
    x = x - min(x)
    # estimate initial parameters from the data
    # curve fit the test data
    """
    try:
        fitted_parameters, p_cov = curve_fit(func, x, y, initial_parameters, maxfev=500000)
    except RuntimeError as err:
        log.debug(err)
        return float("nan"), [x, y, "0", "0"]
    model_predictions = func(x, *fitted_parameters)
    """
    try:
        kde = gaussian_kde(np.repeat(x, y.astype(np.int)))
        model_predictions1 = kde.pdf(x)
    except Exception as err:
        log.info(err)
        return float("nan"), [x, y, "0"]

    """
    abs_error = model_predictions - y
    squared_errors = np.square(abs_error)
    mean_squared_errors = np.mean(squared_errors)
    root_mean_squared_error = np.sqrt(mean_squared_errors)
    r_squared = 1.0 - (np.var(abs_error) / np.var(y))
    """
    return pearsonr(model_predictions1, y)[0], [x, y, model_predictions1]


class ATS:
    u"""summary_line

    utr	reference_id	reference_name	number_of_reads	inferred_sites	alpha_arr	beta_arr	ws	L
    5:122704569-122705569:-	ENSG00000251538,ENST00000514657,ENSE00002071527	LINC02201,LINC02201-203	13	122704940,
    122704839	629,730	5,5	0.07566976284989525,0.9218149110366607,0.002515326113444209	1000
    """

    def __init__(self,
                 utr: str,
                 reference_id: str,
                 reference_name: str,
                 number_of_reads: int,
                 inferred_sites: List[int],
                 alpha_arr: List[int],
                 beta_arr: List[int],
                 ws_arr: List[float],
                 L: int
                 ):
        self.utr = Region.create(utr)
        self.reference_id = reference_id
        self.reference_name = reference_name
        self.number_of_reads = number_of_reads
        self.inferred_sites = inferred_sites
        self.alpha_arr = alpha_arr
        self.beta_arr = beta_arr
        self.ws = ws_arr
        self.L = L
        self.predict = []
        self.similarity = []
        self.data = []

    @property
    def chrom(self) -> str:
        return self.utr.chromosome

    @property
    def strand(self) -> str:
        return self.utr.strand

    def __str__(self) -> str:
        res = []
        for i in [
            self.utr, self.reference_id, self.reference_name,
            self.number_of_reads, self.inferred_sites, self.alpha_arr,
            self.beta_arr, self.ws, self.L, self.predict, self.similarity
        ]:
            if isinstance(i, list):
                if not i:
                    i = "na"
                else:
                    i = ",".join([str(x) for x in i])
            else:
                i = str(i)

            res.append(i)

        return "\t".join(res)

    @classmethod
    def create(cls, record: str):
        u"""

        """
        records = record.strip().split("\t")
        return cls(
            utr=records[0],
            reference_id=records[1],
            reference_name=records[2],
            number_of_reads=int(records[3], ),
            inferred_sites=[int(x) for x in records[4].split(",")],
            alpha_arr=[int(x) for x in records[5].split(",")],
            beta_arr=[int(x) for x in records[6].split(",")],
            ws_arr=[float(x) for x in records[7].split(",")],
            L=int(records[8])
        )

    def load(self, bams: List[str], site: int):
        vals = {}
        start = site if self.utr.strand == "+" else site - 1
        end = site + 1 if self.utr.strand == "+" else site

        for r1, r2 in load_reads(
                bams,
                Region(chromosome=self.utr.chromosome, start=start, end=end, strand=self.strand),
                barcode={}, remove_duplicate_umi=True, inside_region=False
        ):
            if r2 is None:
                r2 = r1

            if r2 is None:
                continue

            ss = r1.start if self.utr.strand == "+" else r2.end
            vals[ss] = vals.get(ss, 0) + 1

        return np.array(list(vals.keys())), np.array(list(vals.values()))

    def check(self, bigwig: str, classifier, bams: List[str]):
        u"""
        check the whether the ATS is valid

        :param bigwig: path to atac bigwig file
        :param classifier: the machine learning classifier
        :param bams: path to bam files to calculate the gaussian density similarity
        """
        with pyBigWig.open(bigwig) as r:
            chrom = self.chrom
            if r.chroms(chrom) is None:
                chrom = "chr" + self.chrom if not self.chrom.startswith("chr") else self.chrom.lstrip("chr")

                if r.chroms(chrom) is None:
                    chrom = ""

            if r.chroms(chrom) is not None:
                for i in range(len(self.inferred_sites)):
                    if classifier:
                        start, end = self.inferred_sites[i] - 500, self.inferred_sites[i] + 500
                        y_pred = classifier.predict([r.values(chrom, start, end)])
                        self.predict.append("true" if y_pred[0] == 1 else "false")

                    if bams:
                        x, y = self.load(bams, self.inferred_sites[i])
                        s, val = calculate_gaussian_density_similarity(x, y)
                        self.similarity.append(s)
                        self.data.append(val)


def call(args):
    ats, bigwig, classifier, bams = args
    ats.check(bigwig, classifier, bams)

    if ats.inferred_sites:
        return ats
    return None


@click.command()
@click.option('-i', '--infile', type=click.Path(exists=True), required=True, help='The path to ats file.')
@click.option('-b', '--bigwig', type=click.Path(exists=True), required=True, help='The path to ATAC bigwig file.')
@click.option('-o', '--output', type=click.Path(), required=True, help='The path to output file.')
@click.option('-c', '--classifier', type=click.Path(), help='The path to classifier file.')
@click.option('-n', '--n-jobs', type=int, default=1, help='The number of processes to use.')
@click.option('-d', '--debug', type=click.Path(), help='The number of processes to use.')
@click.argument("bams", nargs=-1, type=click.Path(exists=True))
def filter(infile: str, bigwig: str, output: str, classifier: str, n_jobs: int, debug: str,  bams: List[str]):
    u"""
    Filter ATSs by machine learning classifier or gaussian density similarity.
    """

    if not classifier and not bams:
        log.error("Please provide classifier or bam files")

    rfc = None
    if classifier:
        log.info("Load classifier from %s" % classifier)
        with open(classifier, "rb") as r:
            rfc = pickle.load(r)

    log.info("Load ATSs")
    header = None
    values = []
    with open(infile) as r:
        for line in r:
            if not header:
                header = line.strip()
            else:
                values.append(ATS.create(line))

            # if debug and len(values) > 100:
            #     break

    log.info("Predict")

    w1 = None
    if debug:
        w1 = open(debug, "w+")

    with open(output, "w+") as w:
        w.write(f"{header}\tpredict\tsimilarity\n")

        with Pool(n_jobs) as p:
            for row in list(track(
                    p.imap_unordered(call, [[ats, bigwig, rfc, bams] for ats in values]),
                    total=len(values), description="Filtering...")):
                if row:
                    w.write(f"{row}\n")

                    if w1:
                        val = []
                        for data in row.data:
                            if isinstance(data, str):
                                cols = data
                            else:
                                cols = [",".join([str(i) for i in x]) for x in data]
                            cols.append(str(row.utr))
                            val.append("\t".join(cols))
                        w1.write("\n".join(val) + "\n")
    if w1:
        w1.close()


if __name__ == '__main__':
    filter()

