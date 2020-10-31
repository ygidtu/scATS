#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.10.26 by Zhang Yiming

Due to BioAlignments.jl currently seems cannot read bam by chromosomes, 
therefore using pysam for this 
"""
import os

from multiprocessing import Pool, cpu_count

import click
import pysam

from tqdm import tqdm


def filter(rec) -> bool:
    u"""
    filter low quality reads
    """
    if rec.has_tag("NH") and rec.get_tag("NH") != 1:
        return False
    
    if rec.is_duplicate and rec.is_qcfail:
        return False

    return rec.has_tag("CB") and rec.has_tag("UB")


def determine_strand(rec) -> str:
    u"""
    as function name says
    """
    strand = "*"
    if rec.is_read1:
        strand = "-" if rec.is_reverse else "+"
    else:
        strand = "+" if rec.is_reverse else  "-"

    return strand


def process_bam(args):
    u"""
    read and count the different start sites by cell and gene
    """
    path, ref = args
    start_sites = {}
    with pysam.AlignmentFile(path, "rb") as r:
        for rec in r.fetch(ref):
            if not filter(rec): continue

            strand = determine_strand(rec)
            if  strand == "*": continue

            start_pos = rec.reference_start if  strand == "+" else rec.reference_end

            cb = f"{rec.get_tag('CB')}\t{rec.get_tag('UB')}"
            temp_sites = start_sites.get(cb, {})
            temp_sites[start_pos] = temp_sites.get(start_pos, 0) + 1
            start_sites[cb] = temp_sites
    
    return start_sites


@click.command()
@click.option(
    '-i',
    '--input-bam',
    type=click.Path(exists=True),
    help='The path to bam file.',
    required=True
)
@click.option(
    '-o',
    '--output',
    type=str,
    help ='The path to output file'
)
@click.option(
    '-p',
    '--process',
    type=click.IntRange(1, cpu_count(), clamp=True),
    help='The number of cpu to use',
    default=1
)
def main(
    input_bam:str, process:int, output:str
):
    u"""
    """
    out_dir = os.path.dirname(os.path.abspath(output))
    os.makedirs(out_dir, exist_ok=True)

    with pysam.AlignmentFile(input_bam) as r:
        cmds = [(input_bam, r) for r in r.references]
    
    with Pool(process) as p:
        res = list(tqdm(p.imap(process_bam, cmds), total=len(cmds)))

    with open(output, "w+") as w:
        for r in tqdm(res):
            for cb, sites in r.items():
                temp = [f"{i}:{j}" for i, j in sites.items()]
                w.write(f"{cb}\t{';'.join(temp)}\n")



if __name__ == '__main__':
    main()
