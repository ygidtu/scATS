#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import pysam
import pickle
import numpy as np

import pandas as pd
from loguru import logger
from collections import defaultdict

from apamix.apamix import APA
from utils.utils import dotdict, dict_list
from utils.bam import cigar_support, collapse_intron, check_strand



def run(arg):
    line, bamfile, cb_df, outdir, verbose = arg
    region_name = line
    chrom, left_site, right_site, strand = line.split('\t')
    # logger.info(f'processing {region_name}')
    # to check whether success
    pickle_tmp = f'{outdir}/tmp/{chrom}_{left_site}_{right_site}_{strand}.pickle'

    if os.path.exists(pickle_tmp) and os.stat(pickle_tmp).st_size != 0:
        logger.info(f'cache found, return cache results. {pickle_tmp}')
        pickle_in = open(pickle_tmp, 'rb')
        res = pickle.load(pickle_in)
        res = res.drop('cb',axis=1).fillna(0).astype('int64')
        pickle_in.close()
        return res
    
    pickle_out = open(pickle_tmp, 'wb')
    bamfh = pysam.AlignmentFile(bamfile,'rb', check_sq=False)

    apa_reads = dotdict({
            'r2_len_arr': [],
            'r2_utr_st_arr': [],
            'r1_len_arr': [],
            'polya_len_arr': [],
            'pa_site_arr': [],
            'cb': [],
            'umi': []
        })

    left_site = int(left_site)
    right_site = int(right_site)


    """
    Fix this line, because lambda function can't picklize and
    we can't parallel process this.
    """

    # res = defaultdict(lambda :defaultdict(list))
    res = defaultdict(dict_list)

    for line in bamfh.fetch(chrom,int(left_site) - 1, int(right_site) + 1):
        if not check_strand(line, strand):
            continue

        start_site = line.reference_start + 1
        end_site = line.reference_end + 1
        query_name = line.query_name

        if start_site < left_site and end_site > right_site:
            continue

        if line.is_read1:
            res[query_name]['r1'].append(line)
        else:
            res[query_name]['r2'].append(line)

    for k,v in res.items():
        label = '_'.join(sorted(v.keys()))
        if label == 'r1':
            continue
        elif label == 'r2':
            read2 = v[label]
            if len(read2) == 1:
                read2 = read2[0]
            else:
                # if there were multiple reads, continue
                continue

            r2_relative_start = collapse_intron(read2)
            if r2_relative_start > right_site and strand == '+':
                continue

            if r2_relative_start < left_site and strand == '-':
                continue

            pa_support = cigar_support(read2, strand)

            try:
                dt = len(read2.get_tag('DT'))
            except KeyError:
                dt = 0

            if strand == '+':
                r2_st = r2_relative_start - left_site +  1
                if pa_support == 'yes':
                    pa_site = read2.reference_end + 1 - left_site
                else:
                    pa_site = np.NaN
            else:
                r2_st = right_site - r2_relative_start +  1
                if pa_support == 'yes':
                    pa_site = right_site - read2.reference_start + 1
                else:
                    pa_site = np.NaN
            r1_dt = np.NaN
        else:
            r1 = v['r1']
            r2 = v['r2']
            if len(r1) != 1 or len(r2) != 1:
                continue
            read1 = r1[0]
            read2 = r2[0]

            """
            TODO
            Check the max gap, must below the max insert length?
            """

            if max(read1.reference_end, read2.reference_end) - \
                min(read1.reference_start, read2.reference_start) > 1000:
                continue

            cigar_info = read1.cigar


            r2_relative_start = collapse_intron(read2)

            if r2_relative_start > right_site and strand == '+':
                continue

            if r2_relative_start < left_site and strand == '-':
                continue

            pa_support = cigar_support(read2, strand)

            try:
                dt = len(read2.get_tag('DT'))
                if read2.is_reverse:
                    if cigar_info[0][0] == 4:
                        r1_dt = int(cigar_info[0][-1]) + dt
                    else:
                        r1_dt = dt
                else:
                    if cigar_info[-1][0] == 4:
                        r1_dt = int(cigar_info[-1][-1]) + dt
                    else:
                        r1_dt = dt

            except KeyError:
                continue
            if strand == '+':
                pa_site = read1.reference_end + 1 - left_site
                r2_st = r2_relative_start - left_site +  1
            else:
                pa_site = right_site - read1.reference_start + 1
                r2_st = right_site - r2_relative_start +  1
        try:
            # for 10X
            cb = read2.get_tag('CB')
            umi = read2.get_tag('UB')
            apa_reads.r2_utr_st_arr.append(r2_st)
            apa_reads.r2_len_arr.append(read2.query_alignment_length)
            apa_reads.r1_len_arr.append(dt)
            apa_reads.polya_len_arr.append(np.NaN)
            apa_reads.pa_site_arr.append(np.Nan)
            apa_reads.cb.append(cb)
            apa_reads.umi.append(umi)
        except KeyError:
            continue


    if len(apa_reads.r2_utr_st_arr) < 50:
        return None
    elif len(apa_reads.r2_utr_st_arr) > 1000000:
        # skip the huge region
        huge_tmp = f'{outdir}/huge/{chrom}_{left_site}_{right_site}_{strand}'
        huge_in = open(huge_tmp, 'w')
        huge_in.write('\n')
        huge_in.close()
        return None
    else:
        pass

    apa_reads = pd.DataFrame(apa_reads)

    relative_offset = 150
    apa_reads.r2_utr_st_arr = apa_reads.r2_utr_st_arr + relative_offset
    apa_reads.pa_site_arr = apa_reads.pa_site_arr + relative_offset
    utr_l = max([max(apa_reads.r1_utr_st_arr) + max(apa_reads.r1_len_arr) + 300, 1000])
        
    test = APA(
        n_max_apa=5,
        n_min_apa=1,
        r1_utr_st_arr=apa_reads.r1_utr_st_arr,
        r1_len_arr=apa_reads.r1_len_arr,
        r2_len_arr=apa_reads.r2_len_arr,
        polya_len_arr=apa_reads.polya_len_arr,
        pa_site_arr=apa_reads.pa_site_arr,
        utr_len = utr_l,
        cb=apa_reads.cb,
        umi=apa_reads.umi,
        region_name=region_name,
        verbose=False
    )

    cb = apa_reads.cb
    umi = apa_reads.umi
    del apa_reads

    logger.debug(f'Processing {chrom}, {left_site}, {right_site}, {strand}')

    res = test.inference()
    del test

    df = pd.concat([pd.Series(res.label, name='label'), cb, umi], axis=1)
    df.columns = ['label', 'cb', 'umi']


    ignore_label = len(res.alpha_arr) # remove last noise component

    # Correct pA sites from relative position to genome position
    if strand == '+':
        alpha_arr = left_site + res.alpha_arr - relative_offset - 1
    else:
        alpha_arr = right_site - res.alpha_arr + relative_offset + 1
    alpha_arr = [f'{chrom}:{int(alpha)}:{int(beta)}:{strand}' for alpha, beta in zip(alpha_arr, res.beta_arr.values.tolist())]

    # final_res = defaultdict(int_dic)
    df = df[df.label != ignore_label]
    df = df.groupby(['cb','umi'])['label'].max().reset_index()
    # TODO, add collapse UMI with one distane toleration
    df = df.groupby(['cb','label'])['umi'].size().reset_index()# calculate umi counts
    df = df.pivot_table(index=['cb'], columns='label', values='umi').reset_index()# make pivot table
    rename_dict = dict(zip(range(len(alpha_arr)), alpha_arr))
    df.rename(columns = rename_dict, inplace = True)
    res = pd.merge(cb_df, df, on = ['cb'], how='left')

    pickle.dump(res, pickle_out)
    pickle_out.close()

    res = res.drop('cb',axis=1).fillna(0).astype('int64')
    return res


def wraper_process(arg):
    return run(arg)
