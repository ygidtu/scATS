#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys

from argparse import ArgumentParser

import editdistance
import pysam

from tqdm import tqdm



def main(input_bam: str, output_bam: str, sequence: str, threshold: int = 3):
    u"""
    main 

    :param input_bam: the path to input bam file
    :param output_bam: the path to output bam
    :param sequence: target sequence
    :return None
    """
    output_bam = os.path.abspath(output_bam)

    outdir = os.path.dirname(output_bam)
    os.makedirs(outdir, exist_ok = True)

    with pysam.AlignmentFile(input_bam) as r:
        with pysam.AlignmentFile(output_bam, "wb", template = r) as w:
            for record in tqdm(r, desc = "Reading"):
                """
                name = record.query_name
                name = name.split("_")
                record.set_tag("UB", name[1])
                record.set_tag("CB", name[0] + "-1")
                w.write(record)
                """

                # using Levenshtein distance to calculate the number of mismatch between reference sequence and target sequence
                ref = record.query_sequence

                for i in range(0, len(ref), len(sequence)):
                    if editdistance.eval(ref[i:i+len(sequence)], sequence) < threshold:
                        w.write(record)
                        break
    pysam.index(output_bam)


if __name__ == '__main__':
    parser = ArgumentParser(description="Filter bam by ref sequence")

    parser.add_argument("-i", "--input", help = "The path to input bam file", type = str, required = True)
    parser.add_argument("-o", "--output", help = "The path to output bam file", type = str, required = True)
    parser.add_argument("-s", "--seq", help = "The reference sequence", type = str, required = True)
    parser.add_argument("-t", "--threshold", help = "The filter threshold", type = int, default = 3)

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(0)

    try:
        args = parser.parse_args(sys.argv[1:])
    except Exception as err:
        print(err)
        parser.print_help()
        exit(1)

    main(
        input_bam = args.input,
        output_bam = args.output,
        sequence = args.seq,
        threshold = args.threshold
    )
    