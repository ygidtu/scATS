#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import gzip
import os
import pickle
import sys

import click
from tqdm import tqdm

import cli
from loguru import logger
from utils.bed6 import Bed6, BedUtils
from utils.isoform import Transcripts
from utils.utils import window


@click.command()
@click.option(
    "-g",
    '--gtf',
    type=click.Path(exists=True),
    help='The gtf (gz or not) file for preparing utr and intron region.',
    required=True
)
@click.option(
    "-p",
    '--prefix',
    type=str,
    default='genes',
    help ='The prefix of output bed file.'
)
@click.option(
    "-d",
    '--distance',
    type=int,
    default=0,
    help ='The distance between utr to merge.'
)
def prepare(
    gtf,
    prefix,
    distance
    ):
    if not all([gtf]):
        cli(['prepare', '--help'])
        sys.exit(1)
    
    if distance < 0:
        distance = 0

    utr_outfile = prefix + '_utr.bed'
    intron_outfile = prefix + '_intron.bed'
    gtf_prefix = os.path.splitext(os.path.abspath(gtf))[0]

    gene_db = '{}.pkl'.format(gtf_prefix)
    
    if os.path.exists(gene_db):
        logger.info("load gene_db")
        # if exsit, just read it
        struc = pickle.load(open(gene_db, 'rb'))
    else:
        logger.info("create gene_db")
        # if not, prepare the transcript information
        struc = Transcripts(gtf).structure()
        pickle.dump(struc, open(gene_db, 'wb'))

    # prepare BedTool obj for exon line from gtf file
    logger.info("prepare BedTool obj for exon line from gtf file")
    exon_from_gtf = []
    with gzip.open(gtf, "rt") if gtf.endswith("gz") else open(gtf, "r") as r:
        for line in tqdm(r):
            if line.startswith("#"):
                continue
            line = str(line).split('\t')
            bio_type = line[2]
            if bio_type != 'exon':
                continue
            exon_from_gtf.append(Bed6(
                chrom=line[0],
                start=int(line[3]),
                end=int(line[4]),
                name=".",
                score=".",
                strand=line[6],
            ))

    intron_lst = set()
    # only keep these biotype, keep same with 10X, but not same to MCA.
    mt_lst = set()
    utr_lst = set()
    genes_use = {'antisense', 'lincRNA', 'protein_coding'}
    logger.info("prepare utr from gene_db")
    # 按照基因 -> 转录本 -> 外显子的层次顺序遍历所有外显子信息
    for gid, tids in struc.items():
        for tid, info in tids.items():
            if info['gene_biotype'] not in genes_use:
                continue
            
            ## 此处获取3' UTR区域
            ## +：则获取最后的exon
            ## -：则获取最开始的exon
            if info['strand'] == '-':
                st, en, = info['exon'][-1]
            else:
                st, en, = info['exon'][0]

            line_info = Bed6(
                chrom=str(info['chrom']), 
                start=int(st),
                end=int(en),
                name=';'.join([str(info['gene_biotype']), gid, tid, str(info['gene_name'])]),
                score=".",
                strand=str(info['strand']),
            )
            if str(info['chrom']).upper() == 'MT':
                mt_lst.add(line_info)
            else:
                # utr_out.write(str(line_info) + "\n")
                utr_lst.add(line_info)

            # intron information
            try:
                exon_lst = info['exon']
                exon_lst = sorted(exon_lst, key=lambda x:[x[0], x[1]])
                if len(exon_lst) == 1:
                    continue
                else:
                    for two_exons in window(exon_lst, 2):
                        try:
                            intron_lst.add(Bed6(
                                chrom=str(info['chrom']), 
                                start=int(two_exons[0][1] + 1), 
                                end=int(two_exons[1][0] - 1), 
                                name=';'.join([
                                        str(info['chrom']),
                                        gid, tid,
                                        str(info['gene_name'])
                                ]), 
                                score=".",
                                strand=str(info['strand']),
                            ))
                        except ValueError as err:
                            logger.error(err)
                            logger.error(info)
                            logger.error(two_exons)
                            exit(err)
            except KeyError:
                continue

    """
    sort and merge utr region information
    """
    logger.info("sort and merge utr region information")
    utr_lst = BedUtils.self_merge(list(utr_lst), distance=distance, strandness=True)
    utr_lst += BedUtils.self_merge(list(mt_lst), distance=0, strandness=True)

    with open(utr_outfile, "w+") as w:
        for i in tqdm(utr_lst, desc="Write utr"):
            w.write(i.get_bed4() + "\n")

    """
    prepare intronic region:
     - remove any overlap with utr region which output from last step
     - subtract any overlap with exon region
    """
    logger.info("prepare intronic region")
    intron_lst  = BedUtils.self_merge(list(intron_lst), distance=0, strandness=True)

    intron_lst = BedUtils.subtract(intron_lst, list(exon_from_gtf))

    with open(intron_outfile, "w+") as w:
        for i in tqdm(intron_lst, desc="Write intron"):
            w.write(i.get_bed4() + "\n")
