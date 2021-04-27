#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Modified at 2021.04.25 by Zhang

Core ATS model
"""
import math

from __future__ import annotations
from functools import lru_cache, reduce
from typing import List, Tuple

# setup matplotlib backend
import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt


from scipy import stats
from scipy.signal import find_peaks

from logger import log

#################### CODE STRUCTURE #######################
# Part 1: Isoform mapping, tree construction
# Part 2: isoform processing and assignment
# Part 3: ATS mixture

################### Part 1: Isoform mapping #####################################

class Window:
    DEL_VAL = 1 << 31

    def __init__(self, start: int = DEL_VAL, end: int = DEL_VAL) -> Window:
        """
        :type start: int
        :type end: int
        """
        self.start = start
        self.end = end

    def __key(self) -> Tuple[int, int]:
        return self.start, self.end

    def __hash__(self):
        return hash(self.__key())

    def is_empty(self):
        return self.start == self.end

    def __bool__(self):
        return self.start != Window.DEL_VAL

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f'{self.start}  {self.end}  {len(self)}'

    #################### window vs window -> bool #####################
    def __lshift__(self, other):
        return self.end <= other.start

    def __lt__(self, other):
        return self.start < other.start < self.end < other.end

    def __le__(self, other):
        return other.start <= self.start and self.end <= other.end

    def __eq__(self, other):
        if isinstance(other, Window):
            return self.start == other.start and self.end == other.end
        else:
            return False

    def __ge__(self, other):
        return self.start <= other.start and other.end <= self.end

    def __gt__(self, other):
        return other.start < self.start < other.end < self.end

    def __rshift__(self, other):
        return other.end <= self.start

    # overlap
    def __ne__(self, other):
        return not (other.end <= self.start or self.end <= other.start)

    # adjacent
    def adj(self, other):
        return self.end == other.start or self.start == other.end

    #################### window vs point operation #####################
    # if a point falls within the window
    def __contains__(self, point):
        if self.is_empty():
            return False
        return self.start <= point < self.end

    def shift_start(self, start):
        if self.start == Window.DEL_VAL:
            return self
        else:
            return Window(start, start + self.end - self.start)


class WinList(list):
    def __init__(self, *args, is_sorted=False):
        super().__init__(*args)
        self.is_sorted = is_sorted

    def __str__(self):
        return "\n".join(['start  end  length'] + [str(w) for w in self.__iter__()])

    def append(self, __object: Window) -> None:
        super().append(__object)
        self.is_sorted = False

    def rm_empty_win(self) -> WinList:
        return WinList([w for w in self if w], is_sorted=self.is_sorted)

    def sort(self) -> None:
        if not self.is_sorted:
            super().sort(key=lambda win: (win.start, win.end))
            self.is_sorted = True

    def rm_duplicate(self) -> WinList:
        res = WinList(set(self))
        res.sort()
        return res

    # merge overlapping windows
    def merge(self) -> WinList:
        winlist = self.rm_empty_win()
        assert len(winlist) > 0
        if not winlist.is_sorted:
            winlist.sort()
        res_list = WinList([winlist[0]])
        for win in winlist:
            curr_win = res_list[-1]  # note that windows are sorted first by start, then by end
            if curr_win.start <= win.start <= curr_win.end:
                res_list[-1] = Window(curr_win.start, win.end)
            else:
                res_list.append(win)
        res_list.is_sorted = True
        return res_list

    # split the winList into tiny non-overlapping intervals, including regions that are not covered
    # e.g. [(0,3), (2,4), (6,8)] => [(0,2), (2,3), (3,4), (4,6), (6,8)]
    def split(self):
        winlist = self.rm_empty_win()
        if len(winlist) == 0:
            return winlist
        if not winlist.is_sorted:
            winlist.sort()
        boarder = set()
        for win in winlist:
            if not win:
                continue
            boarder.add(win.start)
            boarder.add(win.end)
        boarder_arr = sorted(list(boarder))
        winlist = [Window(i, j) for i, j in zip(boarder_arr, boarder_arr[1:])]
        return WinList(winlist, is_sorted=True)

    # return the left most position of the windows list
    # mainly for building trees, note that there may be Empty windows in this list, non-empty windows are sorted
    def get_left(self):
        for win in self:
            if win.start != Window.DEL_VAL:
                return win.start
        return Window.DEL_VAL

    # return the right most position of the windows list
    def get_right(self):
        for win in reversed(self):
            if win.end != Window.DEL_VAL:
                return win.end
        return Window.DEL_VAL

    def get_range(self):
        return Window(self.get_left(), self.get_right())


class TreeNode:
    def __init__(self, left=None, right=None, winlist: WinList = WinList()):
        """
        :type left: TreeNode
        :type right: TreeNode
        :type winlist: a list of windows, each
        """
        self.left = left
        self.right = right
        self.winlist = winlist  # ith window corresponds to range of this tree node


def map_iso_to_gene(iso_list: WinList, gene_list: WinList) -> WinList:
    """
    :param iso_list: exon windows of the isoform on gene
    :param gene_list: smallest windows across isoforms by intersecting all exons of all isoforms on the gene
    :return: a WinList (same size as gene_list), with indexes relative to the isoform (local index)
    """
    assert iso_list.is_sorted
    assert gene_list.is_sorted
    res_win_list = WinList()
    curr_exon_j = 0  # jth exon on isoform
    local_iso_index = 0
    for win in gene_list:
        if curr_exon_j >= len(iso_list):
            res_win_list.append(Window())
            continue
        if win <= iso_list[curr_exon_j]:
            res_win_list.append(win.shift_start(local_iso_index))
            local_iso_index += len(win)
            if iso_list[curr_exon_j].end == win.end:
                curr_exon_j += 1
        else:
            res_win_list.append(Window())
    return res_win_list


def build_tree(isowins_list: List[WinList]):
    """
    :param isowins_list: each element is an aligned WinList of the same window across different isoforms
    :return: a tree
    """

    def dfs(i_left, i_right):
        root = TreeNode()
        if i_left == i_right - 1:
            root.winlist = alignwins_list[i_left]
            root.left = None
            root.right = None
        else:
            mid = (i_left + i_right) // 2
            tmplist = [WinList(isowins[i_left:i_right]).get_range() for isowins in isowins_list]
            root.winlist = WinList(tmplist)
            root.left = dfs(i_left, mid)
            root.right = dfs(mid, i_right)
        return root

    # each element is an aligned WinList of the same window across different isoforms
    alignwins_list = [WinList(wins) for wins in zip(*isowins_list)]
    return dfs(0, len(alignwins_list))


def concatenate_winlist(winlist1: WinList, winlist2: WinList) -> WinList:
    res_list = [concatenate_windows(w1, w2) for w1, w2 in zip(winlist1, winlist2)]
    return WinList(res_list)


def concatenate_windows(win_left: Window, win_right: Window) -> Window:
    if win_left and win_right and win_left.end == win_right.start:
        return Window(win_left.start, win_right.end)
    else:
        return Window()


def query(root: TreeNode, query_winlist: WinList, iso_index: int) -> WinList:
    """
    map the query_winlist to isoform, get the relative index of the query window on all isoforms
    :param root: root of the tree
    :param query_winlist: winlist of a read mapped to the "iso_index"th isoform
    :param iso_index: index of given isoform, used for locating the window
    :return: mapped windows on all isoforms
    """

    @lru_cache(maxsize=1024)
    def query_window(root: TreeNode, qwin: Window, i_iso: int) -> WinList:
        if root is None or root.winlist[i_iso].is_empty():
            return none_winlist
        if not qwin <= root.winlist[i_iso]:
            return none_winlist
        if root.left is None and root.right is None:  # leaf case
            local_st = qwin.start - root.winlist[i_iso].start
            local_win = Window(local_st, local_st + len(qwin))
            res_winlist = WinList()
            for win in root.winlist:
                if not win:
                    res_winlist.append(Window())
                else:
                    res_winlist.append(local_win.shift_start(win.start + local_st))
            return res_winlist
        # root is inner node, one child might be intron
        if root.left.winlist[i_iso].end != Window.DEL_VAL:
            mid = root.left.winlist[i_iso].end
        else:
            mid = root.right.winlist[i_iso].start

        if qwin.end <= mid:
            return query_window(root.left, qwin, i_iso)
        elif qwin.start >= mid:
            return query_window(root.right, qwin, i_iso)
        else:
            qwin_left = Window(qwin.start, mid)
            qwin_right = Window(mid, qwin.end)
            left_winlist = query_window(root.left, qwin_left, i_iso)
            right_winlist = query_window(root.right, qwin_right, i_iso)

            return WinList([concatenate_windows(lw, rw) for lw, rw in zip(left_winlist, right_winlist)])

    n_isoform = len(root.winlist)
    none_winlist = WinList([Window() for _ in range(n_isoform)])

    res_list = [query_window(root, win, iso_index) for win in query_winlist]
    res_winlist = reduce(concatenate_winlist, res_list)

    return res_winlist

# map a window to gene
def map_to_gene(root: TreeNode, query_winlist: WinList, iso_index: int) -> WinList:
    """
    map the query_winlist to gene, get the relative index of the query window on the gene
    :param root: root of the tree, the first exon must be the full gene
    :param query_winlist: a WinList on the "iso_index"th isoform
    :param iso_index: index of given isoform, used for locating the window
    :return: mapped windows on the gene
    """

    @lru_cache(maxsize=1024)
    def map_window(root: TreeNode, qwin: Window, i_iso: int) -> WinList:
        if root is None or root.winlist[i_iso].is_empty():
            return none_list
        if not qwin <= root.winlist[i_iso]:
            return none_list
        if root.left is None and root.right is None:  # leaf case
            local_st = qwin.start - root.winlist[i_iso].start
            local_win = Window(local_st, local_st + len(qwin))
            gene_win = root.winlist[0]
            res_winlist = WinList( [local_win.shift_start(gene_win.start + local_st), ] )
            return res_winlist
        # root is inner node, one child might be intron
        if root.left.winlist[i_iso].end != Window.DEL_VAL:
            mid = root.left.winlist[i_iso].end
        else:
            mid = root.right.winlist[i_iso].start

        if qwin.end <= mid:
            return map_window(root.left, qwin, i_iso)
        elif qwin.start >= mid:
            return map_window(root.right, qwin, i_iso)
        else:
            qwin_left = Window(qwin.start, mid)
            qwin_right = Window(mid, qwin.end)
            left_list = map_window(root.left, qwin_left, i_iso)
            right_list = map_window(root.right, qwin_right, i_iso)
            return left_list+right_list

    none_list = WinList()

    res_list = [map_window(root, win, iso_index) for win in query_winlist]
    res_winlist = reduce(lambda a, b: a + b, res_list)

    return res_winlist


################### Part 2: isoform processing and assignment #####################################
# modify exons of isoform such that ATS falls into an exon
def proc_isoform(ats_pos: int, iso_wins: WinList) -> WinList:
    """
    :param iso_wins_list: list of isoforms, 0th element is gene range
    :return: a list of isoforms with exons modified
    """
    # step 1: check if [max(0, ats-50), ats+50] falls in exon
    # step 1.1: if yes => no special treatment
    # step 1.2: if no  => merge [max(0, ats-50), min(ats+200, exon_start)] with the nearest first exon
    ats_win = Window(ats_pos, ats_pos + 1)
    i = 0
    while iso_wins[i] << ats_win:
        i += 1
    if i > 0:  # ATS is after the first exon
        return None
    first_exon_ind = i  # first downstream exon

    ats_win = Window(max(0, ats_pos - 50), ats_pos + 50)
    if ats_win <= iso_wins[first_exon_ind]:
        return iso_wins

    expanded_ats_win = Window(max(0, ats_pos - 50), min(ats_pos + 200, iso_wins[first_exon_ind].start + 1))

    if expanded_ats_win << iso_wins[first_exon_ind]:  # expanded ats window must overlap with isoform
        return None

    iso_wins.append(expanded_ats_win)
    iso_wins.is_sorted = False
    return iso_wins.merge()  # sorted wins


# calculate relative positions of reads on each isoform
def get_read_relative_pos(iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList],
                          read_labels):
    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([win for iso_wins in iso_wins_list for win in iso_wins])
    all_win_on_gene = all_wins.split()

    # map all isoforms to the gene (atomic non-overlapping windows)
    atom_isowins_list = []
    for iso_wins in iso_wins_list:
        atom_isowins_list.append(map_iso_to_gene(iso_wins, all_win_on_gene))

    # build the tree
    tree = build_tree(atom_isowins_list)

    # get relative position of R1 and R2 on the gene
    n_frag = len(read_labels)
    n_iso = len(iso_wins_list)
    read_st_mat = np.zeros((n_frag, n_iso), dtype=int)
    read_en_mat = np.zeros((n_frag, n_iso), dtype=int)

    for i in range(n_frag):
        r1_wins = r1_wins_list[i]
        r2_wins = r2_wins_list[i]
        res1 = query(tree, r1_wins, read_labels[i])
        res2 = query(tree, r2_wins, read_labels[i])
        for j in range(n_iso):
            if res1[j] and res2[j]:
                read_st_mat[i, j] = res1[j].start
                read_en_mat[i, j] = res2[j].end

    return read_st_mat, read_en_mat


# given an ATS position, calculate the posterior distribution of different isoforms
# if return all zeros, it hints this may be a novel ATS
def assign_isoform(ats_pos: int, iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList],
                   read_labels, mu_frag=400, sd_frag=30, min_frag_len=150):
    """
    :param ats_pos: position of given ATS
    :param iso_wins_list: a list of isoform, each is a WinList of Windows (exons).
            the first must be the full gene Window(0, utr_length)
    :param r1_wins_list: r1 winlist, a consecutive read can map to multiple windows on the gene
    :param r2_wins_list: r2 winlist
    :param read_labels: isoform index of each pair-read, 0 is the full gene
    :param mu_frag: mean of fragment length, e.g. 400
    :param sd_frag: sd of fragment length, e.g. 30
    :param min_frag_len: minimum fragment length, fragments shorter than this are equally assigned to all isoforms
    :return: posterior probability of different isoforms, e.g. np.array([0, 0.3, 0.7, 0])
            the first exon will always have weight 0
    """
    # extend ATS and merge with exons

    # 0th element of iso_wins_list is the whole gene window i.e. [Window(0,utr_length),]
    # return None if ats_pos is not near the first exon
    iso_wins_list = [proc_isoform(ats_pos, iso_wins) for iso_wins in iso_wins_list]
    n_orig_iso = len(iso_wins_list)
    valid_iso_inds = np.array([i for i, w in enumerate(iso_wins_list) if w])

    # in case ats_pos is not close to the start position of any isoform, i.e. ats_pos only compatible with the whole gene
    if len(valid_iso_inds) == 1:
        return np.zeros(n_orig_iso, dtype='float')

    # remove empty isoforms
    iso_wins_list = [w for i, w in enumerate(iso_wins_list) if w]

    # change read labels from all isoforms to selected isoforms
    label_map = {k: i for i, k in enumerate(valid_iso_inds)}
    read_labels = [label_map.get(k, 0) for k in read_labels]  # unselected isoforms will be mapped to the first exon

    # get relative positions of the reads on each isoform
    read_st_mat, read_en_mat = get_read_relative_pos(iso_wins_list, r1_wins_list, r2_wins_list, read_labels)

    # assign each pair-end read to different isoforms
    n_frag, n_iso = read_st_mat.shape
    read_len_mat = read_en_mat - read_st_mat
    post_prob_mat = np.zeros((n_frag, n_iso), dtype='float')
    post_prob_mat[:, 1:] = cal_post_prob(mu_frag, sd_frag, read_len_mat[:, 1:])  # all-zero rows in read_len_mat are properly handled

    # special treatment to short fragments, overwrite previous calculation
    # if the fragment is less then the minimum fragment size (150bp),
    # they should be assigned to all valid isoforms equally
    # min_frag_len = 150
    flag_arr1 = np.sum(read_len_mat[:, 1:], axis=1) > 1         # only need to be treated if mapped to multiple isoforms
    flag_arr2 = np.max(read_len_mat[:, 1:], axis=1) < min_frag_len   # fragment length less than min_frag_len
    short_frag_inds = np.logical_and(flag_arr1, flag_arr2)
    tmpmat = read_len_mat[short_frag_inds, 1:] > 0
    post_prob_mat[short_frag_inds, 1:] = tmpmat / np.sum(tmpmat, axis=1)[:, np.newaxis]

    iso_post_prob = np.sum(post_prob_mat, axis=0)
    iso_post_prob = iso_post_prob / np.sum(iso_post_prob)  # 对应于该ATS的各isoform比例

    res_post_prob = np.zeros(n_orig_iso, dtype='float')
    res_post_prob[valid_iso_inds] = iso_post_prob

    return res_post_prob


# calculate the posterior probability of each DNA fragment on different isoforms
def cal_post_prob(mu_frag_size: int, sd_frag_size: int, len_iso_mat: np.ndarray):
    """
    :param mu_frag_size: mean of fragment size
    :param sd_frag_size: std of fragment size
    :param len_iso_mat: n_frag x n_iso numpy matrix, each element is the length of DNA fragment on each isoform
    :return: n_frag x n_iso numpy matrix, each element is the posterior of a given fragment on an isoform
    """

    def norm_pdf(x, mean, sd):
        prob_density = (1 / (np.sqrt(2 * np.pi) * sd)) * np.exp(-0.5 * ((x - mean) / sd) ** 2)
        return prob_density

    max_arr = np.max(len_iso_mat, axis=1)
    val_inds = max_arr > 0

    prob_mat = np.zeros(len_iso_mat.shape, dtype="float")
    prob_mat[val_inds, :] = norm_pdf(len_iso_mat[val_inds, :], mu_frag_size, sd_frag_size)
    prob_mat[len_iso_mat == 0] = 0  # len=0 means reads cannot be mapped to isoform

    prob_mat[val_inds, :] = prob_mat[val_inds, :] / np.sum(prob_mat[val_inds, :], 1)[:, np.newaxis]

    return prob_mat


################### Part 3: ATS mixture #####################################

# calculate the posterior probability of each DNA fragment on different isoforms
def cal_post_prob(mu_frag_size: int, sd_frag_size: int, len_iso_mat: np.ndarray):
    """
    :param mu_frag_size: mean of fragment size
    :param sd_frag_size: std of fragment size
    :param len_iso_mat: n_frag x n_iso numpy matrix, each element is the length of DNA fragment on each isoform
    :return: n_frag x n_iso numpy matrix, each element is the posterior of a given fragment on an isoform
    """

    def norm_pdf(x, mean, sd):
        prob_density = (1 / (np.sqrt(2 * np.pi) * sd)) * np.exp(-0.5 * ((x - mean) / sd) ** 2)
        return prob_density

    max_arr = np.max(len_iso_mat, axis=1)
    val_inds = max_arr > 0

    prob_mat = np.zeros(len_iso_mat.shape, dtype="float")
    prob_mat[val_inds, :] = norm_pdf(len_iso_mat[val_inds, :], mu_frag_size, sd_frag_size)
    prob_mat[len_iso_mat == 0] = 0

    prob_mat[val_inds, :] = prob_mat[val_inds, :] / np.sum(prob_mat[val_inds, :], 1)[:, np.newaxis]

    return prob_mat

# calculate estimated density, last component is uniform component
def est_density(para, x_arr):
    K = para.K
    y_arr = np.zeros(len(x_arr))
    for k in range(K):
        y_arr += para.ws[k] * stats.norm(loc=para.alpha_arr[k], scale=para.beta_arr[k]).pdf(x_arr)
    y_arr += para.ws[K] * 1 / para.L
    return y_arr


# plot density given by the parameters
def plot_para(para, x_arr=None, line_style='-', color=None, label=None):
    if x_arr is None:
        x_arr = np.arange(para.L + 200)
    y_arr = est_density(para, x_arr)
    alpha_inds = np.searchsorted(x_arr, para.alpha_arr)

    plt.plot(x_arr, y_arr, linestyle=line_style, label=label, color=color)
    plt.vlines(para.alpha_arr, ymin=0, ymax=y_arr[alpha_inds], linestyle=line_style, color=color)


# plot estimated result
def plot_est_vs_real(est_para, real_para):
    """
    plot the estimated result versus
    :param est_para: estimated parameters
    :param real_para: ground truth parameters
    :return:
    """
    x_arr = np.arange(est_para.L + 200)
    pred_y_arr = est_density(est_para, x_arr)
    real_y_arr = est_density(real_para, x_arr)

    plt.style.use('ggplot')
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    plot_para(est_para, x_arr=x_arr, line_style='--', color=colors[0], label='pred')
    plot_para(real_para, x_arr=x_arr, line_style=':', color=colors[1], label='real')
    plt.legend(loc='best')

    plt.show()


# class for all parameters
class Parameters:
    def __init__(self, title='', alpha_arr=None, beta_arr=None, ws=None, L=None):
        self.title = title
        self.alpha_arr = alpha_arr
        self.beta_arr = beta_arr
        self.ws = ws
        self.K = len(self.alpha_arr)
        self.L = L

    def __str__(self):
        outstr = '-' * 10 + f'{self.title} K={self.K}' + '-' * 10 + '\n'
        outstr += f'K={self.K} L={self.L} Last component is uniform component.\n'
        outstr += f'alpha_arr={self.alpha_arr}\n'
        outstr += f'beta_arr={self.beta_arr}\n'
        outstr += f'ws={np.around(self.ws, decimals=2)}\n'
        if hasattr(self, 'bic'):
            outstr += f'bic={np.around(self.bic, decimals=2)}\n'
        outstr += '-' * 30 + '\n'
        return outstr


# ats mixture model inference
# using a class
class AtsModel(object):
    def __init__(
        self,
        n_max_ats,
        n_min_ats=1,
        st_arr=None,
        utr_length=2000,
        mu_f=300,
        sigma_f=50,
        min_ws=0.01,
        max_unif_ws=0.1,
        max_beta=50,
        fixed_inference_flag=False
    ):

        self.n_max_ats = n_max_ats  # maximum number of ATS sites
        self.n_min_ats = n_min_ats  # minimum number of ATS sites

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        self.st_arr = st_arr  # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR

        self.L = utr_length  # length of UTR region

        # fragment size information
        self.mu_f = mu_f  # fragment length mean
        self.sigma_f = sigma_f  # fragment length standard deviation

        # pa site information
        self.min_ws = min_ws  # minimum weight of ATS site
        self.max_unif_ws = max_unif_ws  # maximum weight of uniform component
        self.max_beta = max_beta  # maximum std for ATS site

        # inference with fixed parameters
        self.fixed_inference_flag = fixed_inference_flag

        # fragment number
        self.n_frag = len(self.st_arr)
        assert self.n_frag == len(self.st_arr) f"self.n_frag ({self.n_frag}) != len(self.st_arr) (l{len(self.st_arr)});"

        # initialize and other info
        self.step_size = 5
        self.nround = 50
        self.unif_log_lik = None
        self.predef_beta_arr = None

        self.pos_infinite = float("inf")
        self.neg_infinite = float('-inf')

    def cal_z_k(self, para, k, log_zmat):
        # K = len(ws) - 1  # last component is uniform component
        ws = para.ws
        alpha_arr = para.alpha_arr
        beta_arr = para.beta_arr
        if k < para.K:
            log_zmat[:, k] = np.log(ws[k]) + self.lik_l_ab(self.st_arr, alpha_arr[k], beta_arr[k], log=True)
        else:
            log_zmat[:, k] = np.log(ws[k]) + self.unif_log_lik
        return log_zmat

    # Z is log likelihood
    @staticmethod
    def norm_z(Z):
        Z = np.exp(Z - np.max(Z, axis=1, keepdims=True))
        Z = Z / np.sum(Z, axis=1, keepdims=True)
        return Z

    # maximize ws given Z
    def maximize_ws(self, Z):
        ws = np.sum(Z, axis=0) / Z.shape[0]
        if ws[-1] > self.max_unif_ws:
            ws[:-1] = (1 - self.max_unif_ws) * ws[:-1] / np.sum(ws[:-1])
            ws[-1] = self.max_unif_ws
        return ws

    def mstep(self, para, Z, k):
        alpha_arr = para.alpha_arr
        beta_arr = para.beta_arr

        tmp_sumk = np.sum(Z[:, k])
        # avoid division by zero
        if tmp_sumk < 1e-8:
            Z[:, k] += 1e-8
            Z = AtsModel.norm_z(Z)
            tmp_sumk = np.sum(Z[:, k])

        para.ws = self.maximize_ws(Z)

        alpha_arr[k] = np.sum(Z[:, k] * self.st_arr) / tmp_sumk

        tmp_beta = math.sqrt(np.sum(Z[:, k] * ((self.st_arr - alpha_arr[k]) ** 2)) / tmp_sumk)

        idx = np.searchsorted(self.predef_beta_arr, tmp_beta, side='left')
        if idx == len(self.predef_beta_arr):
            beta_arr[k] = self.predef_beta_arr[idx - 1]
        elif idx > 0 and self.predef_beta_arr[idx] - tmp_beta >= tmp_beta - self.predef_beta_arr[idx - 1]:
            beta_arr[k] = self.predef_beta_arr[idx - 1]
        else:
            beta_arr[k] = self.predef_beta_arr[idx]

        return para

    # mstep when alpha_arr and beta_arr are fixed
    def mstep_fixed(self, para, Z, k):
        # avoid division by zero
        if np.sum(Z[:, k]) < 1e-8:
            Z[:, k] += 1e-8
            Z = AtsModel.norm_z(Z)

        para.ws = self.maximize_ws(Z)

        return para

    @staticmethod
    def elbo(log_zmat, Z):
        lb = AtsModel.exp_log_lik(log_zmat, Z) + np.sum(stats.entropy(Z, axis=1))
        return lb

    # calculate the expected log joint likelihood
    @staticmethod
    def exp_log_lik(log_zmat, Z):
        return np.sum(Z[Z != 0] * log_zmat[Z != 0])

    # uniform component likelihood
    def lik_f0(self, log=False):
        if log:
            return -math.log(self.L)
        else:
            return 1 / self.L

    # p(l|alpha, beta)
    @staticmethod
    def lik_l_ab(l_arr, alpha, beta, log=False):
        if log:
            return stats.norm(loc=alpha, scale=beta).logpdf(l_arr)
        else:
            return stats.norm(loc=alpha, scale=beta).pdf(l_arr)

    # generate random k such that each K is a group and no consecutive elements are the same
    @staticmethod
    def gen_k_arr(K, n):
        def _gen(K):
            ii = 0
            last_ind = -1
            arr = np.random.permutation(K)
            while True:
                if ii % K == 0:
                    np.random.shuffle(arr)
                    if arr[0] == last_ind:
                        tmpi = np.random.choice(K - 1) + 1
                        arr[0], arr[tmpi] = arr[tmpi], arr[0]
                    ii = 0
                    last_ind == arr[-1]
                yield arr[ii]
                ii += 1

        if K == 0 or K == 1:
            return np.zeros(n, dtype='int')
        ite = _gen(K)
        res = []
        for _ in range(n):
            res.append(next(ite))
        return np.array(res, dtype='int')

    @staticmethod
    def cal_bic(log_zmat, Z):
        N, K = Z.shape
        res = -2 * AtsModel.exp_log_lik(log_zmat, Z) + (3 * K + 1) * np.log(N)  # the smaller bic, the better model
        return res

    def fixed_inference(self, para):
        para.ws = self.init_ws(len(para.alpha_arr))
        res = self.em_algo(para, fixed_inference_flag=True)
        return res

    # perform inference for K components
    def em_algo(self, para, fixed_inference_flag=False):
        lb = self.neg_infinite
        lb_arr = []
        N = self.n_frag
        K = para.K

        k_arr = self.gen_k_arr(K, self.nround)

        log_zmat = np.zeros((N, K + 1))
        for k in range(K + 1):
            log_zmat = self.cal_z_k(para, k, log_zmat)

        for i in range(self.nround):
            log.debug('iteration=', i + 1, '  lb=', lb)

            # E-Step
            log_zmat = self.cal_z_k(para, k_arr[i], log_zmat)

            Z = self.norm_z(log_zmat)

            if fixed_inference_flag:
                para = self.mstep_fixed(para, Z, k_arr[i])
            else:
                para = self.mstep(para, Z, k_arr[i])

            lb_new = self.elbo(log_zmat, Z)
            lb_arr.append(lb_new)

            if np.abs(lb_new - lb) < np.abs(1e-6 * lb):
                break
            else:
                lb = lb_new

        if i == self.nround:
            log.debug(f'Run all {i + 1} iterations. lb={lb}')
        else:
            log.debug(f'Converge in  {i + 1} iterations. lb={lb}')

        bic = AtsModel.cal_bic(log_zmat, Z)
        log.debug("bic=", bic)
        log.debug('estimated ws:  ', np.around(para.ws, decimals=2))
        log.debug("estimated alpha: ", np.around(para.alpha_arr, decimals=2))
        log.debug("estimated beta: ", np.around(para.beta_arr, decimals=2))

        # nd = len(lb_arr)
        # if nd >= 3:
        #     plt.plot(list(range(nd - 3)), lb_arr[3:nd])
        #     plt.show()

        # sorted_inds = sorted(range(len(alpha_arr)), key=lambda k: alpha_arr[k])
        sorted_inds = np.argsort(para.alpha_arr)
        para.alpha_arr = para.alpha_arr[sorted_inds]
        para.alpha_arr = np.rint(para.alpha_arr).astype('int')  # round to nearest integer
        para.beta_arr = para.beta_arr[sorted_inds]
        para.ws[0:K] = para.ws[sorted_inds]

        if not fixed_inference_flag:
            para.title = 'Estimated parameters'
        para.bic = bic
        para.lb_arr = lb_arr

        return para

    def sample_alpha(self, n_ats):
        # bw = (5 * self.step_size) / np.std(self.st_arr)
        # kernel = stats.gaussian_kde(self.st_arr, bw_method=bw)
        kernel = stats.gaussian_kde(self.st_arr)
        x_arr = np.arange(-100, self.L + 100)  # extend to include peaks in 0 or L-1
        y_arr = kernel.pdf(x_arr)
        peak_inds, _ = find_peaks(y_arr)
        peaks = x_arr[peak_inds]
        peaks_ws = y_arr[peak_inds] / sum(y_arr[peak_inds])

        if n_ats <= len(peaks):
            return np.random.choice(peaks, size=n_ats, replace=False, p=peaks_ws)
        else:
            mu = np.random.choice(peaks, size=n_ats - len(peaks), replace=True, p=peaks_ws)
            mu = np.sort(np.concatenate((peaks, mu)))
            shift = np.rint(np.random.normal(loc=0, scale=5 * self.step_size, size=n_ats))
            return mu + shift

    def init_ws(self, n_ats):
        ws = np.random.uniform(size=(n_ats + 1))
        ws = ws / sum(ws)
        if ws[-1] > self.max_unif_ws:
            ws[:-1] = ws[:-1] * (1 - self.max_unif_ws)
            ws[-1] = self.max_unif_ws
        return ws

    def init_para(self, n_ats):
        # alpha_arr = np.random.choice(self.st_arr, size=n_ats, replace=True)
        alpha_arr = self.sample_alpha(n_ats)
        beta_arr = np.random.choice(self.predef_beta_arr, size=n_ats, replace=True)
        ws = self.init_ws(n_ats)

        para = Parameters(title='Initial parameters', alpha_arr=alpha_arr, beta_arr=beta_arr, ws=ws, L=self.L)
   
        log.debug(para)

        return para

    # remove components with weight less than min_ws
    def rm_component(self, para):
        rm_inds = [i for i in range(para.K) if para.ws[i] < self.min_ws]
        if len(rm_inds) == 0:
            return para

        log.warn(f'Remove components {rm_inds} with weight less than min_ws={self.min_ws}.')
        keep_inds = np.array([i for i in range(para.K) if not para.ws[i] < self.min_ws])
        para.alpha_arr = para.alpha_arr[keep_inds]
        para.beta_arr = para.beta_arr[keep_inds]
        para.K = len(keep_inds)
        para.ws = None
        para = self.fixed_inference(para)
        return para

    def em_optim0(self, n_ats):
        n_trial = 5
        lb_arr = np.full(n_trial, self.neg_infinite)
        bic_arr = np.full(n_trial, self.pos_infinite)
        res_list = list()

        for i in range(n_trial):

            log.debug(f'K={n_ats} | i_trial={i + 1} | n_trial={n_trial}')
            para = self.init_para(n_ats)

            log.debug(para)

            res_list.append(self.em_algo(para))

            lb_arr[i] = res_list[i].lb_arr[-1]
            bic_arr[i] = res_list[i].bic

        min_ind = np.argmin(bic_arr)
        res = res_list[min_ind]

        res.title = 'Estimated Parameters'
        log.debug(res)

        return res

    def get_label(self, para, st_arr=None):
        if st_arr is None:
            st_arr = self.st_arr
        N = len(st_arr)
        K = para.K
        log_zmat = np.zeros((N, K + 1), dtype='float')
        for k in range(K + 1):
            log_zmat = self.cal_z_k(para, k, log_zmat)
        Z = self.norm_z(log_zmat)
        label_arr = np.argmax(Z, axis=1)
        return label_arr


    def run(self):
        if self.max_beta < self.step_size:
            raise Exception("max_beta=" + str(self.max_beta) + " step_size=" + str(self.step_size) +
                            ", max_beta has to be greater than step_size!")

        self.predef_beta_arr = np.arange(self.step_size, self.max_beta, self.step_size)

        n_ats_trial = self.n_max_ats - self.n_min_ats + 1
        bic_arr = np.full(n_ats_trial, self.pos_infinite)
        res_list = list()

        self.unif_log_lik = self.lik_f0(log=True)

        for i, n_ats in enumerate(range(self.n_max_ats, self.n_min_ats - 1, -1)):
            # print()
            # print(20 * '*' + ' k = ' + str(n_ats) + ' ' + 20 * '*')
            res = self.em_optim0(n_ats)
            res_list.append(res)
            bic_arr[i] = res.bic

        min_ind = np.argmin(bic_arr)
        res = res_list[min_ind]

        res = self.rm_component(res)
        res.label_arr = self.get_label(res)

        res.title = f'Final Result'
        log.debug(res)

        return res

