#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Modified at 2021.04.25 by Zhang

Core ATS model
"""
from __future__ import annotations

import math
import json

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

    @classmethod
    def keys(cls) -> List[str]:
        res = []
        for i in cls.__init__.__code__.co_varnames:
            if i.startswith("_") or i == "title" or i == "self":
                continue

            res.append(str(i))
        return res

    def to_res_str(self):
        res = []
        for i in self.keys():
            val = getattr(self, i)

            if "ndarray" in str(type(val)) or "list" in str(type(val)):
                val = ",".join([str(x) for x in val])
            else:
                val = str(val)
            res.append(val)
        return "\t".join(res)


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

        self.L = utr_length  # length of UTR region

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        self.st_arr = [x for x in st_arr if 0 <= x < utr_length]  # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR

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
        assert self.n_frag == len(self.st_arr), f"self.n_frag ({self.n_frag}) != len(self.st_arr) (l{len(self.st_arr)});"

        # initialize and other info
        self.step_size = 5
        self.nround = 50
        self.unif_log_lik = None
        self.predef_beta_arr = None

        self.pos_infinite = float("inf")
        self.neg_infinite = float('-inf')

    def dumps(self) -> Dict:
        res = {}
        for i in self.__init__.__code__.co_varnames:
            if i == "self":
                continue
            res[str(i)] = getattr(self, "L" if i == "utr_length" else i)
        return res

    def dump(self, path: str):
        u"""
        dump this model to json
        :param path: path to json
        """

        with open(path, "w+") as w:
            json.dump(self.dumps(), w, indent = 4)

    @classmethod
    def load(cls, path: str):
        u"""
        restore model from json file
        """
        with open(path) as r:
            data = json.load(r)
        
        return cls(**data)


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
        u"""

        """
        alpha_arr = para.alpha_arr
        beta_arr = para.beta_arr

        tmp_sumk = np.sum(Z[:, k])
        # avoid division by zero
        if tmp_sumk < 1e-8:
            Z[:, k] += 1e-8
            Z = AtsModel.norm_z(Z)
            tmp_sumk = np.sum(Z[:, k])

        para.ws = self.maximize_ws(Z)

        # log.debug(f"tmp_sumk={tmp_sumk}; k={k}; Z.shape={Z.shape}")
        u"""
        2020.05.06 Here, 
        para.alpha_arr is an empry list, but try to set value by index
        """
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
        u"""
        Call mstep and msted_fixed
        """
        lb = self.neg_infinite
        lb_arr = []
        N = self.n_frag
        K = para.K

        k_arr = self.gen_k_arr(K, self.nround)

        log_zmat = np.zeros((N, K + 1))
        for k in range(K + 1):
            log_zmat = self.cal_z_k(para, k, log_zmat)

        for i in range(self.nround):
            log.debug(f'iteration={i + 1}  lb={lb}')

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
        log.debug(f"bic={bic}")
        log.debug(f'estimated ws: {np.around(para.ws, decimals=2)}')
        log.debug(f"estimated alpha: {np.around(para.alpha_arr, decimals=2)}")
        log.debug(f"estimated beta: {np.around(para.beta_arr, decimals=2)}")

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

        # log.debug(f"n_ats={n_ats}; len(peaks)={len(peaks)}")

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
        # log.debug(f"alpha_arr inited: {alpha_arr}")
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

    def run(self) -> Parameters:
        if not self.st_arr:
            return None
        if self.max_beta < self.step_size:
            raise Exception("max_beta=" + str(self.max_beta) + " step_size=" + str(self.step_size) +
                            ", max_beta has to be greater than step_size!")

        self.predef_beta_arr = np.arange(self.step_size, self.max_beta, self.step_size)

        n_ats_trial = self.n_max_ats - self.n_min_ats + 1
        bic_arr = np.full(n_ats_trial, self.pos_infinite)
        res_list = list()

        self.unif_log_lik = self.lik_f0(log=True)

        for i, n_ats in enumerate(range(self.n_max_ats, self.n_min_ats - 1, -1)):
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


if __name__ == '__main__':
    pass