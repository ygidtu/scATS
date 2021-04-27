#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Functions used to extract utr from gtf file
"""
import os
import gzip

from rich.progress import Progress, track

from src.region import GTF


def process(gtf: str, output: str, span: int = 500):
    u"""
    extract utr from gtf file
    :params gtf: path to gtf file
    :params output: path to output file
    :params span: the radius of UTR
    :return None
    """
    output = os.path.abspath(output)
    os.makedirs(os.path.dirname(output), exist_ok = True)

    transcripts = {}

    # read gtf
    with Progress() as progress:
        task_id = progress.add_task(f"Processing {os.path.basename(gtf)}", total=os.path.getsize(gtf))

        with gzip.open(gtf, "rt") if gtf.endswith("gz") else open(gtf, "r") as r:
            for line in r:
                progress.update(task_id, advance=len(str.encode(line)))

                if line.startswith("#"):
                    continue

                rec = GTF.create(line)

                if rec.source == "exon":
                    if rec.transcript_id not in transcripts.keys():
                        transcripts[rec.transcript_id] = []

                    transcripts[rec.transcript_id].append(rec)

    # get first exons
    first_exons = []
    for exons in track(transcripts.values(), description=f"Collecting first exons"):
        exons = sorted(exons)
        exon = exons[0]
        site = exon.start
        if exons[0].strand == "-":
            exon = exons[-1]
            site = exon.end
        
        exon = GTF(
            chromosome=exon.chromosome,
            start=max(site - span, 1),
            end=site + span,
            strand=exon.strand,
            attrs=exon.attrs,
            source = exon.source
        )
        first_exons.append(exon)
    first_exons = sorted(first_exons)

    # merging
    with gzip.open(output, "wt+") if output.endswith("gz") else open(output, "w+") as w:
        curr_exon = first_exons[0]
        for i in track(range(1, len(first_exons)), description=f"Writing {os.path.basename(output)}"):
            if curr_exon & first_exons[i]:
                curr_exon = curr_exon + first_exons[i]
            else:
                w.write(f"{curr_exon.bed}\n") 
                curr_exon = first_exons[i]
        w.write(f"{curr_exon.bed}\n")
                        

if __name__ == '__main__':
    process(
        "/mnt/raid61/Personal_data/zhangyiming/genome/Homo_sapiens/Homo_sapiens.GRCh38.93.gtf", 
        "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/blood_ctss.bed1"
    )
