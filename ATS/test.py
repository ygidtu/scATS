#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.08 by Zhang
This is used for debugging
"""
from __future__ import annotations
from functools import lru_cache, reduce
from typing import List, Tuple

import math
from scipy import stats
from scipy.signal import find_peaks
import numpy as np

from datetime import datetime

from cli.ats import ATSParams
from ats.core import AtsModel as newModel
from ats.core import Parameters

from src.reader import load_reads, load_utr

import warnings
warnings.filterwarnings("ignore")


# ats mixture model inference
# using a class
class AtsModel(object):
    def __init__(self,
        n_max_ats,
        n_min_ats=1,
        st_arr=None,
        utr_length=2000,
        mu_f=300,
        sigma_f=50,
        min_ws=0.01,
        max_unif_ws=0.1,
        max_beta=50,
        fixed_inference_flag=False,
        debug = False
):
        self.n_max_ats = n_max_ats  # maximum number of ATS sites
        self.n_min_ats = n_min_ats  # minimum number of ATS sites

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        self.st_arr = [x for x in st_arr if 0 <= x < utr_length]   # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
        self.L = utr_length  # length of UTR region

        # pa site information
        self.min_ws = min_ws  # minimum weight of ATS site
        self.max_unif_ws = max_unif_ws  # maximum weight of uniform component
        self.max_beta = max_beta  # maximum std for ATS site
        self.mu_f =mu_f
        self.sigma_f=sigma_f

        # inference with fixed parameters
        self.fixed_inference_flag = fixed_inference_flag

        # single end mode
        self.debug = debug

        # fragment number
        self.n_frag = len(self.st_arr)
        assert self.n_frag == len(self.st_arr)

        # initialize and other info
        self.step_size = 5
        self.nround = 50
        self.unif_log_lik = None
        self.predef_beta_arr = None

        self.pos_infinite = float("inf")
        self.neg_infinite = float('-inf')

    #@profile
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
    #@profile
    def norm_z(Z):
        Z = np.exp(Z - np.max(Z, axis=1, keepdims=True))
        Z = Z / np.sum(Z, axis=1, keepdims=True)
        return Z

    # maximize ws given Z
    #@profile
    def maximize_ws(self, Z):
        ws = np.sum(Z, axis=0) / Z.shape[0]
        if ws[-1] > self.max_unif_ws:
            ws[:-1] = (1 - self.max_unif_ws) * ws[:-1] / np.sum(ws[:-1])
            ws[-1] = self.max_unif_ws
        return ws

    #@profile
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
    #@profile
    def mstep_fixed(self, para, Z, k):
        # avoid division by zero
        if np.sum(Z[:, k]) < 1e-8:
            Z[:, k] += 1e-8
            Z = AtsModel.norm_z(Z)

        para.ws = self.maximize_ws(Z)

        return para

    @staticmethod
    #@profile
    def elbo(log_zmat, Z):
        lb = AtsModel.exp_log_lik(log_zmat, Z) + np.sum(stats.entropy(Z, axis=1))
        return lb

    # calculate the expected log joint likelihood
    @staticmethod
    #@profile
    def exp_log_lik(log_zmat, Z):
        return np.sum(Z[Z != 0] * log_zmat[Z != 0])

    # uniform component likelihood
    #@profile
    def lik_f0(self, log=False):
        if log:
            return -math.log(self.L)
        else:
            return 1 / self.L

    # p(l|alpha, beta)
    @staticmethod
    #@profile
    def lik_l_ab(l_arr, alpha, beta, log=False):
        if log:
            return stats.norm(loc=alpha, scale=beta).logpdf(l_arr)
        else:
            return stats.norm(loc=alpha, scale=beta).pdf(l_arr)

    # generate random k such that each K is a group and no consecutive elements are the same
    @staticmethod
    #@profile
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
    #@profile
    def cal_bic(log_zmat, Z):
        N, K = Z.shape
        res = -2 * AtsModel.exp_log_lik(log_zmat, Z) + (3 * K + 1) * np.log(N)  # the smaller bic, the better model
        return res
        
    #@profile
    def fixed_inference(self, para):
        para.ws = self.init_ws(len(para.alpha_arr))
        res = self.em_algo(para, fixed_inference_flag=True)
        return res

    # perform inference for K components
    #@profile
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
            if self.debug:
                print('iteration=', i + 1, '  lb=', lb)

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

        if self.debug:
            if i == self.nround:
                print(f'Run all {i + 1} iterations. lb={lb}')
            else:
                print(f'Converge in  {i + 1} iterations. lb={lb}')

        bic = AtsModel.cal_bic(log_zmat, Z)
        if self.debug:
            print("bic=", bic)
            print('estimated ws:  ', np.around(para.ws, decimals=2))
            print("estimated alpha: ", np.around(para.alpha_arr, decimals=2))
            print("estimated beta: ", np.around(para.beta_arr, decimals=2))

        if self.debug:
            nd = len(lb_arr)
            if nd >= 3:
                plt.plot(list(range(nd - 3)), lb_arr[3:nd])
                plt.show()

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

    #@profile
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

    #@profile
    def init_ws(self, n_ats):
        ws = np.random.uniform(size=(n_ats + 1))
        ws = ws / sum(ws)
        if ws[-1] > self.max_unif_ws:
            ws[:-1] = ws[:-1] * (1 - self.max_unif_ws)
            ws[-1] = self.max_unif_ws
        return ws

    #@profile
    def init_para(self, n_ats):
        # alpha_arr = np.random.choice(self.st_arr, size=n_ats, replace=True)
        alpha_arr = self.sample_alpha(n_ats)
        beta_arr = np.random.choice(self.predef_beta_arr, size=n_ats, replace=True)
        ws = self.init_ws(n_ats)

        para = Parameters(title='Initial parameters', alpha_arr=alpha_arr, beta_arr=beta_arr, ws=ws, L=self.L)
        if self.debug:
            print(para)

        return para

    # remove components with weight less than min_ws
    #@profile
    def rm_component(self, para):
        rm_inds = [i for i in range(para.K) if para.ws[i] < self.min_ws]
        if len(rm_inds) == 0:
            return para

        print(f'Remove components {rm_inds} with weight less than min_ws={self.min_ws}.')
        keep_inds = np.array([i for i in range(para.K) if not para.ws[i] < self.min_ws])
        para.alpha_arr = para.alpha_arr[keep_inds]
        para.beta_arr = para.beta_arr[keep_inds]
        para.K = len(keep_inds)
        para.ws = None
        para = self.fixed_inference(para)
        return para

    #@profile
    def em_optim0(self, n_ats):
        n_trial = 5
        lb_arr = np.full(n_trial, self.neg_infinite)
        bic_arr = np.full(n_trial, self.pos_infinite)
        res_list = list()

        for i in range(n_trial):
            if self.debug:
                print('-----------------K=', n_ats, ' | ', 'i_trial=', i + 1, ' | n_trial=', n_trial, ' -------------')
            para = self.init_para(n_ats)

            if self.debug:
                print(para)

            res_list.append(self.em_algo(para))

            lb_arr[i] = res_list[i].lb_arr[-1]
            bic_arr[i] = res_list[i].bic

        min_ind = np.argmin(bic_arr)
        res = res_list[min_ind]

        res.title = 'Estimated Parameters'
        print(res)

        return res

    #@profile
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

    #@profile
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
        print(res)

        return res


def sample(utr, bams):
    ats = load_utr(utr)

    with open("data/reads.txt", "w+") as w:
        w.write(f"R1_START\tR1_END\tR2_START\tR2_END\tUTR_IDX\n")
        for idx, i in enumerate(ats[:50]):
            reads = load_reads(bams, i)

            for r1, r2 in reads.items():
                site = i.start if i.strand == "+" else i.end

                data = [r1.start - site, r1.end - site, r2.start - site, r2.end - site]
                if all([0<= x < 1000 for x in data]):
                    data.append(idx)
                    w.write("\t".join([str(x) for x in data]) + "\n")
        

def test():
    utr = "data/utr.bed"
    bams = ["/mnt/raid64/ATS/Personal/zhangyiming/bams/select_bam/NCovM1.bam"]

    sample(utr, bams)

    # params = ATSParams(
    #     utr = utr,
    #     bam = bams,
    #     n_max_ats = 5, 
    #     n_min_ats = 1,
    #     utr_length = 1000,
    #     debug = False
    # )
    
    # begin = datetime.now()
    # res = []
    # for idx in range(100):
    #     m = params.get_model(idx, model = AtsModel)
    #     try:
    #         temp = params.run(idx, m)
    #         if temp:
    #             res.append(temp)
    #     except Exception as err:
    #         continue

    # print(len(res))
    # print("old model: ", datetime.now() - begin)

    # begin = datetime.now()
    # res = []
    # for idx in range(100):
    #     m = params.get_model(idx, model = newModel)
    #     try:
    #         temp = params.run(idx, m)
    #         if temp:
    #             res.append(temp)
    #     except Exception as err:
    #         continue
    # print(len(res))
    # print("new model: ", datetime.now() - begin)
 


if __name__ == '__main__':
    test()
    pass
