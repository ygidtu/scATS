import numpy as np
from ATS_core import AtsModel, est_density, Parameters, plot_est_vs_real
from ATS_core import Window, WinList, assign_isoform, map_iso_to_gene, build_tree, query
from typing import List
from csv import reader
import re


class PreprocData:
    pass


def readcsv(filename):
    # read csv file as a list of lists
    with open(filename, 'r') as read_obj:
        # pass the file sobject to reader() to get the reader object
        csv_reader = reader(read_obj, delimiter=' ')
        # Pass reader object to list() to get a list of lists
        list_of_rows = list(csv_reader)
        list_of_rows = [[int(i) for i in l] for l in list_of_rows]
        return list_of_rows


def read_location(filename) -> List[WinList]:
    reslist = []
    with open(filename, 'r') as fh:
        for i, line in enumerate(fh):
            res = re.split(r'\t', line)
            res = [int(e) for e in res]
            tmpwinlist = WinList()
            for st, en in zip(res[::2], res[1::2]):
                tmpwinlist.append(Window(st, en))
            reslist.append(tmpwinlist)
    return reslist


# preprocessing data
def preproc_data() -> PreprocData:
    # prepare data
    filename = 'data/iso.txt'  # isoform information, columns (iso_id, start, end), iso_id=0 is the whole gene
    iso_tbl = readcsv(filename)

    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([Window(row[1], row[2]) for row in iso_tbl])
    all_win_on_gene = all_wins.split()
    # print(all_win_on_gene)

    # map all isoforms to the gene (atomic non-overlapping windows)
    iso_inds = list(set([row[0] for row in iso_tbl]))
    n_iso = len(iso_inds)
    iso_wins_on_gene_list = []
    iso_wins_list = []
    for i in range(n_iso):
        tmp_wins = WinList([Window(row[1], row[2]) for row in iso_tbl if row[0] == i])
        tmp_wins.sort()
        iso_wins_list.append(tmp_wins)
        iso_wins_on_gene_list.append(map_iso_to_gene(tmp_wins, all_win_on_gene))

    # build the tree
    tree = build_tree(iso_wins_on_gene_list)

    # get the positions of r1 and r2 on gene
    filename = 'data/r1_on_gene.txt'  # read 1 information, (start1, end1), (start2, end2), ...
    r1_list = read_location(filename)  # list of WinList

    filename = 'data/r2_on_gene.txt'  # read 2 information, (start1, end1), (start2, end2), ...
    r2_list = read_location(filename)  # list of WinList

    n_frag = len(r1_list)
    assert n_frag == len(r2_list)

    frag_label = [0 for _ in range(n_frag)]  # isoform index of fragments, 0 is the whole gene


    # # query the position of r1 and r2 on different isoforms
    # r1st = np.zeros((n_frag, n_iso), dtype=int)
    # r1en = np.zeros((n_frag, n_iso), dtype=int)
    # r2st = np.zeros((n_frag, n_iso), dtype=int)
    # r2en = np.zeros((n_frag, n_iso), dtype=int)
    # st_arr = np.zeros(n_frag, dtype=int)

    # for i in range(n_frag):
    #     # print(i)
    #     st_arr[i] = r1_list[i][0].start
    #     r1_wins = query(tree, r1_list[i], 0)
    #     r2_wins = query(tree, r2_list[i], 0)
    #     for j, (w1, w2) in enumerate(zip(r1_wins, r2_wins)):
    #         tmpflag = w1.is_empty() or w2.is_empty() or w1.start < 0 or w2.start < 0
    #         if tmpflag:
    #             continue
    #         r1st[i, j], r1en[i, j], r2st[i, j], r2en[i, j] = w1.start, w1.end, w2.start, w2.end

    # prepare results
    data = PreprocData()
    # data.st_arr = st_arr
    # data.r1st = r1st
    # data.r1en = r1en
    # data.r2st = r2st
    # data.r2en = r2en

    data.st_arr = np.array([winlist[0].start for winlist in r1_list])

    data.iso_wins_on_gene_list = iso_wins_on_gene_list
    data.iso_wins_list = iso_wins_list
    data.tree = tree
    data.n_frag = n_frag
    data.n_iso = n_iso
    data.gene_length = all_win_on_gene[-1].end

    data.r1_list = r1_list
    data.r2_list = r2_list
    data.frag_label = frag_label

    data.n_max_ats = 5
    data.n_min_ats = 1
    data.mu_f = 500
    data.sigma_f = 30
    data.min_ws = 0.05
    data.max_unif_ws = 0.1

    return data


# infer ATS positions
def infer_ats(data: PreprocData) -> Parameters:
    atsmix = AtsModel(n_max_ats=data.n_max_ats, n_min_ats=data.n_min_ats, st_arr=data.st_arr,
                      gene_length=data.gene_length,
                      mu_f=data.mu_f,  # fragment length mean
                      sigma_f=data.sigma_f,  # fragment length standard deviation

                      # pa site information
                      min_ws=data.min_ws,  # minimum weight of ATS site
                      max_unif_ws=data.max_unif_ws,

                      # inference with fixed parameters
                      fixed_inference_flag=False,

                      debug=False)
    para = atsmix.run()
    return para


# estimate the proportion of different isoforms given an ATS
def infer_isoform(data: PreprocData, ats_pos, valid_frag_inds):
    mu_frag = data.mu_f
    sd_frag = data.sigma_f
    min_frag_len = 200

    r1_wins_list = [data.r1_list[i] for i in valid_frag_inds]
    r2_wins_list = [data.r2_list[i] for i in valid_frag_inds]
    read_labels = [data.frag_label[i] for i in valid_frag_inds]

    iso_wins_list = data.iso_wins_list

    iso_ws = assign_isoform(ats_pos, iso_wins_list, r1_wins_list, r2_wins_list, read_labels, mu_frag, sd_frag, min_frag_len=min_frag_len)
    iso_ws = np.around(iso_ws[1:], decimals=3)  # remove weights for the first, which is the whole gene
    novel_flag = False
    if all(iso_ws < 1e-8):
        novel_flag = True
        print(f'ats_pos={ats_pos} potentially support a novel isoform')
    print(f'ats_pos={ats_pos} iso_weights={iso_ws}')
    return novel_flag, iso_ws


if __name__ == "__main__":
    data = preproc_data()
    para = infer_ats(data)
    label_arr = para.label_arr

    res = [('AtsIndex', 'AtsPos', 'NovelFlag', 'IsoformWeights')]

    for k, ats_pos in enumerate(para.alpha_arr):
        frag_inds = [i for i, label in enumerate(para.label_arr) if label == k]
        novel_flag, iso_ws = infer_isoform(data, ats_pos, frag_inds)
        res.append((k, ats_pos, novel_flag, iso_ws))
    print("\n", res)
