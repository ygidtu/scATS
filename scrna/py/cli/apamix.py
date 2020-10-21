import os
import signal
import sys
import time
from multiprocessing import Pool, set_start_method

import click
import pandas as pd
import scipy.io
import scipy.sparse
import tqdm

from apamix.inference import wraper_process
from loguru import logger
from utils.utils import dict_list

logger.add('apamix.log',
            rotation='10 MB',
            colorize=True,
            level="INFO")



@click.command()
@click.option(
    '--bed',
    type=str,
    help='The target regions (bed format) used for mix inference',
    required=True
    )
@click.option(
    '--bam',
    type=str,
    help='The bam file (sorted and indexed)',
    required=True
    )
@click.option(
    '--out',
    '-o',
    type=str,
    help='The output path',
    required=True
    )
@click.option(
    '--cores',
    type=int,
    help='Num (cores) of region are infering at once',
    default=1
    )
@click.option(
    '--cb',
    type=str,
    help='The cell barcode file, one cb for one line.',
    required=True
    )
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose mode'
)
def apamix(
    bed,
    bam,
    out,
    cb,
    cores,
    verbose
    ):
    if not all([bed, bam, out, cb]):
        cli(['apamix', '--help'])
        sys.exit(1)

    if not os.path.exists(os.path.join(out, 'tmp')):
        os.makedirs(os.path.join(out, 'tmp'))
        os.makedirs(os.path.join(out, 'huge'))

    target_region = open(bed, 'r')
    res_lst = []
    cb_df = pd.read_csv(cb, names=['cb'])
    # cb_df.cb = list(map(lambda x : x.split('-')[0], cb_df.cb.values))
    
    peak_lst = []

    for i in target_region:
        if not i:
            continue
        peak_lst.append(i.strip())

    target_region.close()

    # pbar = tqdm.tqdm(total=len(peak_lst),desc="Progress", ncols=100, bar_format='{l_bar}{bar}|')
    # def update(*a):
    #     pbar.update()
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # set_start_method('spawn', force=True)
    pool = Pool(processes=cores)
    signal.signal(signal.SIGINT, original_sigint_handler)

    res_lst = []
    try:
        for x in range(len(peak_lst)):
            arg = [peak_lst[x], bam, cb_df, out, verbose]
            # res_lst.append(wraper_process(arg))
            res_lst.append(pool.apply_async(wraper_process, (arg,)))
            # res_lst.append(pool.apply_async(wraper_process, (arg,), callback=update))
            # time.sleep(5)

    except KeyboardInterrupt:
        logger.info('Caught KeyboardInterrupt, terminating workers')
        pool.terminate()
        sys.exit('KeyboardInterrupt')

    pool.close()
    pool.join()


    logger.info('Concating your final sheet')
    final = pd.concat([res.get() for res in res_lst], sort=False, axis=1)
    pa_sites = '\n'.join(final.columns.values.tolist())
    with open(f'{out}/pasite.txt', 'w') as fh:
        fh.write(pa_sites + '\n')

    final = scipy.sparse.csr_matrix(final)
    scipy.io.mmwrite(f"{out}/mm", final)
    logger.info('All done')
