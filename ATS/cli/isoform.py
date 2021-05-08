#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.06 by Zhang
"""
from multiprocessing import cpu_count, Process, Queue
from typing import List
import gzip
import click
from ats.isoform import GTFUtils, assign_isoform
from ats.reader import load_ats
from logger import init_logger, log
from rich import print
from rich.progress import Progress


def consumer(
    input_queue: Queue, 
    output_queue: Queue, 
    gtf:GTFUtils, 
    bams: List[str], 
    mu_f: int, 
    sigma_f: int, 
    min_frag_length: int
):
    u"""
    consumer
    """
    while True:
        utr, region = input_queue.get()
        res = []
        
        try:
            iso_tbl = gtf.read_transcripts(utr)

            if not iso_tbl:
                output_queue.put(res)
                continue
            
            iso_wins_list = iso_tbl.winList

            if not iso_wins_list:
                output_queue.put(res)
                continue

            iso_tbl.set_bams(bams)
            reads = iso_tbl.reads(utr)
            
            r1_wins_list = []
            r2_wins_list = []
            frag_inds = []
            for r1, r2 in reads.items():
                assign1 = iso_tbl.assign(r1)
                assign2 = iso_tbl.assign(r2)

                for a in set(assign1) & set(assign2):
                    r1_wins_list.append(iso_tbl.reads_to_relative(r1, return_winlist = True))
                    r2_wins_list.append(iso_tbl.reads_to_relative(r2, return_winlist = True))
                    frag_inds.append(a)

            for r in region:
                ats_pos = (r.start if r.strand == "+" else r.end) - iso_tbl.gene.start
                iso_ws = assign_isoform(
                    ats_pos, 
                    iso_wins_list, 
                    r1_wins_list, 
                    r2_wins_list, 
                    frag_inds,
                    mu_f, 
                    sigma_f, 
                    min_frag_len=min_frag_length
                )

                ws = []
                ids = []
                gene_id = iso_tbl.get(0)
                for i, j in enumerate(iso_ws):
                    if j > 0:
                        ws.append(str(j))
                        ids.append(iso_tbl.get(i))

                if ws:
                    res.append(f"{r.to_str()}\t{gene_id}\t{','.join(ids)}\t{','.join(ws)}\n")
                else:
                    res.append(f"{r.to_str()}\t{gene_id}\t.\t.\n")
        except Exception as err:
            log.exception(err)
            print(utr)
        output_queue.put(res)


@click.command()
@click.option(
    "-i", "--ats",
    type=click.Path(exists = True),
    required=True,
    help=""" The path to utr file, bed format. """
)
@click.option(
    "-g", "--gtf",
    type=click.Path(exists = True),
    required=True,
    help=""" The path to reference gtf file. """
)
@click.option(
    "-o", "--output", 
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
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
    "--min-frag-length",
    type=int,
    default = 200,
    help=""" The minimum fragment length. """
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
def isoform(
    ats: str,
    gtf: str,
    output: str,
    mu_f: int,
    sigma_f: int,
    min_frag_length: int,
    processes: int,
    debug: bool, 
    bams: List[str],
):
    u"""
    Infer isoforms
    \f

    :param debug: enable debug mode
    """

    init_logger("DEBUG" if debug else "INFO")

    gtf = GTFUtils(gtf)
    ats = load_ats(ats)

    input_queue = Queue()
    output_queue = Queue()

    # generate consumers
    consumers = []
    for _ in range(processes):
        p = Process(
            target=consumer, 
            args=(
                input_queue, 
                output_queue,
                gtf,
                bams,
                mu_f,
                sigma_f,
                min_frag_length,
            )
        )
        p.daemon = True
        p.start()
        consumers.append(p)

    # producer to assign task
    for utr, region in ats.items():
        input_queue.put([utr, region])

    with Progress() as progress:
        task = progress.add_task("Computing...", total=len(ats))

        with gzip.open(output, "wt+") if output.endswith(".gz") else open(output, "w+") as w:
            w.write("ats\tgene_id\ttranscript_ids\tws\n")
            
            while not progress.finished:
                res = output_queue.get(block=True, timeout=None)
                if res:
                    w.write("\n".join(res))

                progress.update(task, advance = 1)
    log.info("DONE")


if __name__ == '__main__':
    isoform()
