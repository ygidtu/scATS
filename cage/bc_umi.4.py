import os
import sys
import pysam
import pandas as pd
import numpy as np
import pickle
import datetime
import random
import string
import shutil

from multiprocessing import Pool
from multiprocessing.pool import MaybeEncodingError

from tqdm import tqdm
from collections import defaultdict
from loguru import logger

def run(args):
	bamfile, chrom, s_t, strand, cb_df, temp_dir = args

	bam_fh = pysam.AlignmentFile(bamfile,'rb', check_sq=False)

	bc_pool = set(cb_df.cb.values.tolist())
	collapse_dict = defaultdict(set)
	left, right = map(int, s_t.split('-'))

	tmp_lst = []
	index_count = 0
	for line in bam_fh.fetch(chrom, left, right):

		# check read1 and quality
		if not line.is_read1 or line.mapping_quality != 255:
			continue

		#  check strand
		strand_from_read = '-' if line.is_reverse else '+'
		if strand != strand_from_read:
			continue

		start_site = line.reference_end if line.is_reverse else line.reference_start

		if not left <= start_site <= right:
			continue
		try:
			current_cb = line.get_tag('CB')
			current_umi = line.get_tag('UB')
		except KeyError as e:
			continue
		# cb_list.append(current_cb)
		if current_cb not in bc_pool:
			continue
		index_count += 1
		if index_count % 100000 == 0:
			suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
			random_string = ''.join(random.choice(string.ascii_lowercase) for i in range(16))
			tmpfile_name = f'{chrom}_{s_t}_{suffix}_{random_string}'
			tmp_lst.append(tmpfile_name)

			tmpfile = open(tmpfile_name, 'wb')
			pickle.dump(collapse_dict, tmpfile)
			tmpfile.close()
			collapse_dict = defaultdict(set)

		collapse_dict[current_cb].add(current_umi)

	if not collapse_dict:
		return None
	if tmp_lst:
		for tmpfile_name in tmp_lst:
			tmpfile = open(tmpfile_name, 'rb')
			tmpdict = pickle.load(tmpfile)
			tmpfile.close()
			collapse_dict.update(tmpdict)
			os.remove(tmpfile_name)

	df = pd.concat(
		[pd.Series(list(collapse_dict.keys()), name='cb', dtype=str),
		pd.Series(map(len,collapse_dict.values()), name=f'{chrom}:{left}-{right}:{strand}', dtype=np.int64)],
		axis=1
		)

	res = pd.merge(cb_df, df, on=['cb'], how='left')
	res = res.drop('cb', axis=1).fillna(0).astype('int64')
	res = res.transpose()

	peak_cache = f'{temp_dir}/{chrom}_{s_t}_{strand}'
	with open(peak_cache, 'wb') as fh:
		pickle.dump(res, fh)

	return peak_cache


def main(bamfile: str, cb_file: str, peakfile: str, outfile: str, n_jobs: int = 12, region_size: int = 200):
	cb_df = pd.read_csv(cb_file, names=['cb'])

	temp_dir = outfile + ".temp"
	if os.path.exists(temp_dir):
		return
		# shutil.rmtree(temp_dir, ignore_errors=True)

	os.makedirs(temp_dir)

	pool = Pool(processes=n_jobs)
	error_record = open(f'{outfile}.toolong.log', 'w')
	peak_region = defaultdict(list)
	with open(peakfile) as ph:
		for peak in tqdm(ph):
			chrom, s_t, strand = peak.strip().split('\t')
			chrom = 'MT' if chrom == 'M' else chrom
			if np.diff(list(map(int, s_t.split('-'))))[0] > region_size:
				error_record.write(f'{chrom}\t{s_t}\t{strand}\n')
				continue
			peak_region[chrom].append((s_t, strand))
	error_record.close()

	md, hd='w', True
	count_inf = 0
	res_lst = []
	try:
		for chrom, info in peak_region.items():
			for s_t, strand in info:
				arg = [bamfile, chrom, s_t, strand, cb_df, temp_dir]
				res_lst.append(pool.apply_async(run, (arg,)))
	except KeyboardInterrupt:
		logger.info('Caught KeyboardInterrupt, terminating workers')
		pool.terminate()
		sys.exit('KeyboardInterrupt')

	md, hd='w', True
	for cache_file in tqdm(res_lst):
		try:
			cache_file = cache_file.get()
		except Exception as e:
			logger.debug(f'error found, {e}')
			continue

		if not cache_file:
			continue
		with open(cache_file, 'rb') as cache_h:
			df = pickle.load(cache_h)
			if hd:
				df.columns = cb_df.cb.values.tolist()
			df.to_csv(outfile, mode=md, header=hd, compression='gzip')
			md, hd='a', False

	logger.info('All done')



if __name__ == '__main__':
	from fire import Fire
	Fire(main)
