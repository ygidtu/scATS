#!/usr/bin/envpython3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Contians all the parameters and command line params handler
"""
import math

from typing import Optional, List

import click
import pysam

from ats.io import load_utr
from logger import log, init_logger


class ATSParams(object):
    u"""
    This a a param handler for ATS inferrence
    """

    def __init__(
        self,
        utr: str,
        bam: List[str], 
        n_max_ats: int = 5, 
        n_min_ats: int = 1,
        utr_length: int = 2000,
        mu_f: int = 300,
        sigma_f: int = 50,
        min_ws: float = 0.01,
        max_unif_ws: float = 0.1,
        max_beta: int = 50,
        fixed_inference_flag: bool = False
    ):
        u"""
        init this function
        
        :params bam: path to bam files or list of bams
        """
        log.info("Load UTR")
        self.utr = load_utr(utr)
        self.bam = self.check_path(bam)
        self.n_max_ats = n_max_ats
        self.n_min_ats = n_min_ats
        self.utr_length = utr_length
        self.muf_f = mu_f
        self.sigma_f = sigma_f
        self.min_ws = min_ws
        self.min_ws = min_ws
        self.max_unif_ws = max_unif_ws
        self.max_beta = max_beta
        self.fixed_inference_flag = fixed_inference_flag
        self.st_arr = []

    @staticmethod
    def check_path(bams: List[str]) -> List[str]:
        u"""
        check input bam files
        """
        for bam in bams:
            try:
                with pysam.AlignmentFile(bam) as r:
                    pass
            except Exception as err:
                log.error(f"{bam} is not a valid bam file: {err}")
                exit(err)
        return bams

    def __str__(self) -> str:
        res = []
        for i in dir(self):
            if i.startswith("__"):
                continue

            if "function" in str(type(getattr(self, i))):
                continue
            
            res.append(f"- {i}: {getattr(self, i)}")

        return "\n".join(res)


@click.command()
@click.option(
    "--utr",
    type=click.Path(exists = True),
    help=""" The path to utr file, bed format """
)
@click.option(
    "--n-max-ats",
    type=click.IntRange(1, math.inf),
    default = 5,
    help=""" The maximum number of ATSs in same UTR """
)
@click.option(
    "--n-min-ats",
    type=click.IntRange(0, math.inf),
    default = 0,
    help=""" The minimum number of ATSs in same UTR. """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The length of UTR """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The estimate length of gene """
)
@click.option(
    "--mu-f",
    type=int,
    default = 300,
    help=""" The mean of fragment length. """
)   
@click.option(
    "--sigma-f",
    type=int,
    default = 50,
    help=""" The standard deviation of fragment length. """
)
@click.option(
    "--min-ws",
    type=float,
    default = 0.01,
    help=""" The minimum weight of ATSs. """
)
@click.option(
    "--max-unif-ws",
    type=float,
    default = 0.1,
    help=""" The maximum weight of uniform component. """
)
@click.option(
    "--max-beta",
    type=int,
    default = 50,
    help=""" The maximum std for ATSs. """
)
@click.option(
    "--fixed-inference",
    is_flag=True,
    type=click.BOOL,
    help=""" Inference with fixed parameters. """
)
@click.option(
    "-d", "--debug",
    is_flag=True,
    type=click.BOOL,
    help=""" Enable debug mode. """
)
@click.argument("bams", nargs = -1, type=click.Path(exists=True))
def inference(
    utr: str,
    n_max_ats: int, 
    n_min_ats: int,
    utr_length: int,
    mu_f: int,
    sigma_f: int,
    min_ws: float,
    max_unif_ws: float,
    max_beta: int,
    fixed_inference: bool,
    debug: bool, bams: List[str],
):
    u"""
    Inference
    \f

    :param debug: enable debug mode
    """

    init_logger("DEBUG" if debug else "INFO")

    if not bams:
        log.error("")

    params = ATSParams(
        utr = utr,
        bam = bams,
        n_max_ats = n_max_ats, 
        n_min_ats = n_min_ats,
        utr_length = utr_length,
        mu_f = mu_f,
        sigma_f = sigma_f,
        min_ws = min_ws,
        max_unif_ws = max_unif_ws,
        max_beta = max_beta,
        fixed_inference_flag = fixed_inference
    )

    print(len(params.utr))
    pass



if __name__ == '__main__':
    inference()

