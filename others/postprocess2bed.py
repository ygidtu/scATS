#!/usr/bin/env python3
#  -*- coding:utf-8 -*-
u"""
Created at 2021.02.25
"""

import os


def main(input_file, output_file):

    data = []

    with open(input_file) as r:
        for line in r:
            line = line.split()

            site = line[1].split(":")
            site[1] = [int(x) for x in site[1].split("-")]

            data.append([site[0], int(site[1][0]), int(site[1][1]), line[0], line[2], site[2]])

    data = sorted(data, key = lambda x: (x[0], x[1], x[2]))

    with open(output_file, "w+") as w:
        for line in data:
            w.write('\t'.join([str(x) for x in line]) + "\n")


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
