import os
from glob import glob
from datetime import datetime

import numpy as np

import pyximport; pyximport.install()
import ats.entropy as cent


def main():
    fs = glob("testfile/*")

    for i in fs:
        # begin = datetime.now()
        data = np.loadtxt(i)
        cent.entropy(data)
        # print(datetime.now() - begin)


if __name__ == '__main__':
    main()
