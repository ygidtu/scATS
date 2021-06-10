#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Setup logger
"""

import logging

from rich.logging import RichHandler

log = logging.getLogger(("rich"))


def init_logger(level="NOTSET"):
    FORMAT = "%(message)s"
    logging.basicConfig(
        level=level,
        format=FORMAT,
        datefmt="[%Y-%m-%d %H:%M:%S]",
        handlers=[RichHandler()]
    )
    global log
    log = logging.getLogger(("rich"))


if __name__ == '__main__':
    pass
