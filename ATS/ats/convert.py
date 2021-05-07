#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.07 by Zhang

This script contains the objects and functions related to the coordinate convertion
"""

from typing import Dict, List, Optional

import pysam

try:
    from ats.reader import load_paired_reads
    from src.loci import BED, Reads
except ImportError:
    from loci import BED, Reads
    from reader import load_paired_reads


class Coordinate(object):
    u"""
    object to handle the coordanite convertion of isoforms from reference file
    """

    def __init__(self, gene: BED, isoforms:  Dict):
        u"""
        init
        :params gene: gene region
        :params isoforms: transcript id -> list of exons
        """

        self.gene = gene
        self.isoforms = isoforms
        self.ids = self.__generate_isoform_idx__(isoforms)
        self.bams = []
    
    def set_bams(self, bams: List[str]):
        u"""
        check input bam files
        """
        for bam in bams:
            try:
                with pysam.AlignmentFile(bam) as r:
                    pass
            except Exception as err:
                raise FileNotFoundError(f"{bam} is not a valid bam file: {err}")
                
        self.bams = bams

    @staticmethod
    def __generate_isoform_idx__(isoforms: Dict):
        u"""
        assign each isoform a unique index id
        """
        data = []
        for iso, exons in isoforms.items():
            data.append([iso, exons])

        data = sorted(data, key=lambda x: [x[1][0].start, x[1][-1].end])

        return [x[0] for x in data]

    @property
    def relative(self) -> List:
        u"""
        convert absolute coordinate to relative
        :return [[id, start, end], [id, start, end]]
        """

        res = [[0, 0, len(self.gene)]]

        for idx, iso in enumerate(self.ids):
            idx += 1
            exons = self.isoforms[iso]

            for e in exons:
                res.append([idx, e.start - self.gene.start, e.end - self.gene.start])

        return res

    def get(self, idx: int) -> Optional[BED]:
        u"""
        get transcript id by index
        """
        if idx == 0:
            return self.gene.id
        
        if idx > len(self.ids):
            return None

        return self.ids[idx - 1]

    def reads(self, region: BED) -> Dict:
        u"""
        load reads and convert and assign isoforms
        """
        return load_paired_reads(self.bams, region)


if __name__ == '__main__':
    from rich import print
    gene = BED("1", 1000, 4000, "+", "", "")

    isos = {
        "iso1": [
            BED("1", 1000, 1150, "+", "", ""),
            BED("1", 1116, 2100, "+", "", ""),
            BED("1", 1823, 3100, "+", "", "")
        ],
        "iso2": [
            BED("1", 1102, 3020, "+", "", ""),
            BED("1", 1203, 3300, "+", "", ""),
            BED("1", 2604, 3560, "+", "", "")
        ],
        "iso3": [
            BED("1", 3405, 4000, "+", "", ""),
            BED("1", 3706, 4000, "+", "", ""),
            BED("1", 3843, 4000, "+", "", "")
        ]
    }

    rel = ReferenceIsoform(gene, isos)

    print(rel.ids)
    print(rel.relative)
    print(rel.get(4))

