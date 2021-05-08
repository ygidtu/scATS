#!/usr/bin/envpython3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Contians all the parameters and command line params handler
"""
import gzip
import json
import math
import os
from multiprocessing import Process, Queue, cpu_count
from typing import List, Optional

import click
import pysam
from ats.core import AtsModel, Parameters
from ats.reader import load_reads, load_utr
from logger import init_logger, log
# from rich import print
from rich.progress import Progress


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
        self.mu_f = mu_f
        self.sigma_f = sigma_f
        self.min_ws = min_ws
        self.min_ws = min_ws
        self.max_unif_ws = max_unif_ws
        self.max_beta = max_beta
        self.fixed_inference_flag = fixed_inference_flag

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

    def __iter__(self):
        for i in range(len(self)):
            yield i

    def __len__(self):
        return len(self.utr)

    @staticmethod
    def __format_reads_to_relative__(reads: List, utr) -> List[int]:
        u"""
        format list of reads into start site array
        """
        st_arr = []

        utr_site = utr.start if utr.strand == "+" else utr.end

        for r in reads:
            if utr.start <= r.start <= r.end <= utr.end:
                site = r.start if utr.strand == "+" else r.end
                st_arr.append(site - utr_site if utr.strand == "+" else utr_site - site)

        return st_arr

    def keys(self) -> List[str]:
        res = ["utr", "reference_id", "reference_name", "infered_sites"]
        res += Parameters.keys()
        return res

    def get_model(self, idx: int):
        u"""
        get model by index
        :param idx: the idx of utr
        """
        if idx < len(self.utr):
            # only use R1
            reads = load_reads(self.bam, self.utr[idx])
            reads = list(reads.keys())
            st_arr = self.__format_reads_to_relative__(reads, self.utr[idx])
            if len(st_arr) <= 1:
                return None
  
            m = AtsModel(
                n_max_ats=self.n_max_ats,
                n_min_ats=self.n_min_ats,
                st_arr=st_arr,
                utr_length=self.utr_length,
                mu_f=self.mu_f,
                sigma_f=self.sigma_f,
                min_ws=self.min_ws,
                max_unif_ws=self.max_unif_ws,
                max_beta=self.max_beta,
                fixed_inference_flag=self.fixed_inference_flag
            )

            return m
        return None

    def format_res(self, idx: int, res: Parameters) -> str:
        u"""
        as name says format ATS model results to meaningful str
        """
        utr = self.utr[idx]
        site = utr.start if utr.strand == "+" else utr.end

        sites = [str(site + x if utr.strand == "+" else site - x) for x in res.alpha_arr]

        data = [
            f"{utr.chromosome}:{utr.start}-{utr.end}:{utr.strand}",
            utr.id,
            utr.name,
            ",".join(sites),
            res.to_res_str()
        ]
        return "\t".join(data)

    def run(self, idx: int, m: AtsModel) -> Optional[str]:
        u"""
        Factory function to execute the ATS model and format results
        """
        if m:
            res = m.run()
            if res:
                return self.format_res(idx, res)
        return None


def consumer(input_queue: Queue, output_queue: Queue, error_queue: Queue, params: ATSParams):
    u"""
    Multiprocessing consumer to perform the ATS core function
    :param input_queue: multiprocessing.Queue to get the index
    :param output_queue: multiprocessing.Queue to return the results
    :param params: the parameters for ATS model
    """

    while True:
        idx = input_queue.get()

        if idx is None:
            log.debug(f"{os.getpid()} existing")
            break
        log.debug(f"{os.getpid()} processing {idx}")
        m = params.get_model(idx)
        res = None
        try:
            res = params.run(idx, m)
        except Exception as err:
            output_queue.put(None)
            error_queue.put(m.dumps())
            log.error(err)
        finally:
            # input_queue.task_done()
            output_queue.put(res)


@click.command()
@click.option(
    "--utr",
    type=click.Path(exists = True),
    required=True,
    help=""" The path to utr file, bed format. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
)
@click.option(
    "--n-max-ats",
    type=click.IntRange(1, math.inf),
    default = 5,
    help=""" The maximum number of ATSs in same UTR. """
)
@click.option(
    "--n-min-ats",
    type=click.IntRange(1, math.inf),
    default = 1,
    help=""" The minimum number of ATSs in same UTR. """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The length of UTR. """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The estimate length of gene. """
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
    default = True,
    type=click.BOOL,
    help=""" Inference with fixed parameters. """
)
@click.option(
    "-d", "--debug",
    is_flag=True,
    type=click.BOOL,
    help=""" Enable debug mode to get more debugging information. """
)
@click.option(
    "-p", "--processes",
    type=click.IntRange(1,cpu_count()),
    default = 1,
    help=""" How many cpu to use. """
)
@click.argument("bams", nargs = -1, type=click.Path(exists=True), required=True)
def inference(
    utr: str,
    output: str,
    n_max_ats: int, 
    n_min_ats: int,
    utr_length: int,
    mu_f: int,
    sigma_f: int,
    min_ws: float,
    max_unif_ws: float,
    max_beta: int,
    fixed_inference: bool,
    processes: int,
    debug: bool, 
    bams: List[str],
):
    u"""
    Inference
    \f

    :param debug: enable debug mode
    """

    init_logger("DEBUG" if debug else "INFO")

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

    # init queues
    input_queue = Queue()
    output_queue = Queue()
    error_queue = Queue()

    # generate consumers
    consumers = []
    for _ in range(processes):
        p = Process(
            target=consumer, 
            args=(
                input_queue, 
                output_queue,
                error_queue,
                params,
            )
        )
        p.daemon = True
        p.start()
        consumers.append(p)

    # producer to assign task
    for i in params:
        input_queue.put(i)
        # join to wait consumers finished
        # input_queue.join()
    
    with Progress() as progress:
        task = progress.add_task("Computing...", total=len(params))
        with open(output, "w+") as w:
            header = '\t'.join(params.keys())
            w.write(f"{header}\n")

            while not progress.finished:
                res = output_queue.get(block=True, timeout=None)
                if res:
                    w.write(f"{res}\n")

                progress.update(task, advance = 1)

    # kept the error data for further debugging
    errors = []
    while not error_queue.empty():
        res = error_queue.get()
        errors.append(res)
    if errors:
        with gzip.open(f"{output}.error_data.json.gz", "wt+") as w:
            json.dump(errors, w, indent = 4)

    log.info("DONE")


if __name__ == '__main__':
    inference()

