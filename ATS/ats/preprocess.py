#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Preprocess of UTR bed and BAM, to get the data for model
"""
import gzip
import os

import pysam


def load_utr(path: str) -> 