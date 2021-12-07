#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import gzip
import os
import re
from multiprocessing import Process, Queue

import random
import numpy as np
import pyfaidx
import pysam
from rich.progress import (BarColumn, Progress, TextColumn, TimeElapsedColumn,
                           TimeRemainingColumn, TransferSpeedColumn, track)


def custom_progress(io: bool = False):
    if io:
        return Progress(
            "[progress.description]{task.description}",
            BarColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TextColumn("| Remaining:"),
            TimeRemainingColumn(),
            TextColumn("| Speed: "),
            TransferSpeedColumn()
        )
    return Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}% ({task.completed}/{task.total})",
        TextColumn("| Elapsed:"),
        TimeElapsedColumn(),
        TextColumn("| Remaining:"),
        TimeRemainingColumn(),
    )


fa = None
NOISE_WS = 0
MAX_READS = 1000
MIN_READS = 50
READS_LEN = 150


def __decode_attr__(attrs):
    data = {}
    for line in attrs.split(";"):
        if line:
            key, val = line.split(' "')
            key = re.sub(r'[";\s]', '', key)
            val = re.sub(r'[";\s]', '', val)
            data[key] = val
    return data


def generate_tss(gtf, bed, total=None):
    count = 0
    data = []
    generated = set()
    with open(bed, "w+") as w, \
            open(bed.replace("bed", "gtf"), "w+") as wg:
        with gzip.open(gtf, "rt") if gtf.endswith("gz") else open(gtf, "r") as r:
            for line in r:
                if line.startswith("#"):
                    continue

                wg.write(line)
                line = line.split()
                if line[2] == "transcript":
                    site = int(line[3]) if line[6] == "+" else int(line[4])
                    site = [site, site + 1 if line[6] == "+" else site - 1]
                    site = sorted(site)

                    attrs = __decode_attr__(" ".join(line[8:]))

                    gn = attrs.get("gene_name", attrs.get("Parent"))
                    tn = attrs.get("transcript_name", attrs.get("Name"))

                    line = f"{line[0]}\t{site[0]}\t{site[1]}\t{gn}_{tn}\t{count}\t{line[6]}\n".replace(
                        '"', '').replace(";", '')

                    site = f"{line[0]}\t{site[0]}\t{site[1]}\t{line[6]}"

                    if site in generated:
                        continue
                    generated.add(site)

                    data.append(line)
                    w.write(line)

                    count += 1

                    if total and count >= total:
                        break
    return data


def generate_reads_length(
    K=50, seed=42,
    len_mu_base=250, len_mu_rand=100,
    len_sd_base=20, len_sd_rand=20
):
    np.random.seed(seed)
    frag_len_mu = len_mu_base + int(len_mu_rand * np.random.rand())
    np.random.seed(seed)
    frag_len_sd = len_sd_base + int(len_sd_rand * np.random.rand())
    np.random.seed(seed)
    return np.random.normal(frag_len_mu, frag_len_sd, K)


def generate_normal_reads(
    K=50, mu=1000, seed=42,
    len_mu_base=250, len_mu_rand=100,
    len_sd_base=20, len_sd_rand=20,
    strand = "+"
):
    random.seed(seed)
    sigma = random.choices(range(30, 50, 5), k=1)

    np.random.seed(seed)
    r1_arr = np.random.normal(mu, sigma, K)

    if strand != "+":
        r1_arr -= READS_LEN

    length = generate_reads_length(
        K, seed=seed,
        len_mu_base=len_mu_base,
        len_mu_rand=len_mu_rand,
        len_sd_base=len_sd_base,
        len_sd_rand=len_sd_rand
    )

    r2_arr = r1_arr + [x // 2 if x < READS_LEN * 2 else x - READS_LEN for x in length]

    r1_arr = r1_arr.astype(int)
    r2_arr = r2_arr.astype(int)

    return [[x, x + READS_LEN - 1] for x in r1_arr], [[x - READS_LEN + 1, x] for x in r2_arr]


def generate_uniform_reads(
    K=50, seed=42, mu=1000,
    len_mu_base=250, len_mu_rand=100,
    len_sd_base=20, len_sd_rand=20, 
    strand = "+"
):
    length = generate_reads_length(
        K, seed=seed,
        len_mu_base=len_mu_base,
        len_mu_rand=len_mu_rand,
        len_sd_base=len_sd_base,
        len_sd_rand=len_sd_rand
    )

    np.random.seed(seed)
    r1_arr = np.random.uniform(mu - length, mu + length)

    if strand != "+":
        r1_arr -= READS_LEN

    r2_arr = r1_arr + length

    r1_arr = r1_arr.astype(int)
    r2_arr = r2_arr.astype(int)

    return [[x, x + READS_LEN - 1] for x in r1_arr], [[x - READS_LEN + 1, x] for x in r2_arr]


def generate_reads(
    K=50, seed=42, mu=1000,
    len_mu_base=250, len_mu_rand=100,
    len_sd_base=20, len_sd_rand=20,
    strand = "+"
):

    num_of_noise = int(K * NOISE_WS)

    r1_arr, r2_arr = generate_normal_reads(
        K - num_of_noise, seed=seed, mu=mu,
        len_mu_base=len_mu_base,
        len_mu_rand=len_mu_rand,
        len_sd_base=len_sd_base,
        len_sd_rand=len_sd_rand,
        strand = strand
    )
    noise1_arr, noise2_arr = generate_uniform_reads(
        num_of_noise, seed=seed, mu=mu,
        len_mu_base=len_mu_base,
        len_mu_rand=len_mu_rand,
        len_sd_base=len_sd_base,
        len_sd_rand=len_sd_rand,
        strand = strand
    )

    return sorted(r1_arr + noise1_arr), sorted(r2_arr + noise2_arr)


def create_reads(name, seq, r1_start, r2_start, ref_id=0, flag=67, tags=None):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag

    a.reference_id = ref_id
    a.reference_start = r1_start if flag in (67, 83) else r2_start
    a.next_reference_id = ref_id
    a.next_reference_start = r2_start if flag not in (67, 83) else r1_start

    a.mapping_quality = 20
    a.cigar = ((0, len(seq)),)
    a.template_length = len(seq)
    a.query_qualities = pysam.qualitystring_to_array('I'*len(seq))
    a.tags = [("NM", 1), ("RG", "L1"), ("NH", 1)] + tags

    return a


def generate_header(fa):
    res = {
        'HD': {'VN': '1.0'},
        'SQ': [{"SN": i, "LN": j.rlen} for i, j in fa.index.items()]
    }
    ref_ids = {i["SN"]: idx for idx, i in enumerate(res["SQ"])}
    return res, ref_ids


def consumer(inQueue, outQueue, ref_ids):

    while True:
        line = inQueue.get()

        if not line:
            break
        chrom, st, en, meta, seed, strand = line.strip().split("\t")

        gn, tn = meta.split("_")

        np.random.seed(int(seed))
        r1_arr, r2_arr = generate_reads(
            K=np.random.randint(MIN_READS, MAX_READS),
            seed=int(seed),
            mu=int(st) if strand == "+" else int(en),
            strand = strand
        )

        reads = []
        count = 1
        for i, j in zip(r1_arr, r2_arr):
            
            st1, en1 = i
            st2, en2 = j
            seq1 = fa.fetch(str(chrom), st1, en1)
            seq1 = seq1.reverse.complement.seq if strand == "+" else seq1.seq

            seq2 = fa.fetch(str(chrom), st2, en2)
            seq2 = seq2.reverse.complement.seq if strand == "+" else seq2.seq

            reads.append([
                f"{meta}_{count}", seq1, 
                st1, # if strand == "+" else en1, 
                st2, # if strand == "+" else en2,
                67 if strand == "+" else 83, 
                ref_ids[str(chrom)],
                [
                    ("GN", gn), ("TN", tn),
                    ("TS", int(st) if strand == "+" else int(en))
                ]
            ])

            reads.append([
                f"{meta}_{count}", seq2,
                st1, # if strand == "+" else en1, 
                st2, # if strand == "+" else en2,
                147 if strand == "+" else 131, 
                ref_ids[str(chrom)],
                [
                    ("GN", gn), ("TN", tn),
                    ("TS", int(st) if strand == "+" else int(en))
                ]
            ])
            count += 1

        outQueue.put(reads)


def generate_bam(gtf, fasta, bed="tss.bed", bam="simu.bam", total=None, n_jobs = 10):
    u"""generate testing files based on gtf and fasta files.
    
    Keyword arguments:
    gtf -- path to reference gtf file
    fasta -- path to genome fasta file
    bed -- path to output tss.bed. default: ./tss.bed
    bam -- path to output simulated bam file. default: ./simu.bam
    total -- how many transcripts used to generate simulated data. default: 10000
    n_jobs -- how many processes to use. default: 10
    Return: None
    """
    
    global fa
    fa = pyfaidx.Faidx(fasta)

    data = generate_tss(gtf, bed, total=total)

    inQueue, outQueue = Queue(), Queue()

    header, ref_ids = generate_header(fa)

    # generate consumers
    consumers = []
    for _ in range(n_jobs):
        p = Process(
            target=consumer,
            args=(inQueue, outQueue, ref_ids,)
        )
        p.daemon = True
        p.start()
        consumers.append(p)

    for line in track(data):
        inQueue.put(line)

    progress = custom_progress()
    temp_bam = bam + ".temp"
    with progress:
        task = progress.add_task("Computing...", total=len(data))

        with pysam.AlignmentFile(temp_bam, "wb", header=header) as outf:
            while not progress.finished:
                res = outQueue.get(block=True, timeout=None)
                if res:
                    for r in res:
                        outf.write(
                            create_reads(
                                name=r[0],
                                seq=r[1],
                                r1_start=r[2],
                                r2_start=r[3],
                                ref_id=r[5],
                                flag=r[4],
                                tags=r[6])
                        )

                progress.update(task, advance=1)

    print("sorting")
    pysam.sort("-o", bam, temp_bam)

    print("indexing")
    pysam.index(bam)

    if os.path.exists(temp_bam):
        os.remove(temp_bam)


if __name__ == '__main__':
    from fire import Fire
    Fire(generate_bam)
