import sys
import pysam
import pandas as pd

from multiprocessing import Pool
from multiprocessing.pool import MaybeEncodingError

from collections import defaultdict
from loguru import logger

def run(args):
	bamfile, chrom, s_t, strand, cb_df = args

	bam_fh = pysam.AlignmentFile(bamfile,'rb', check_sq=False)
	cb_list = []
	umi_list = []

	left, right = map(int, s_t.split('-'))
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
		cb_list.append(current_cb)
		umi_list.append(current_umi)

	if umi_list is None or umi_list is None:
		return None

	df = pd.concat([pd.Series(cb_list, name='cb', dtype=str), 
		pd.Series(umi_list, name= 'umi', dtype=str), 
		pd.Series([f'{chrom}:{left}-{right}:{strand}'] * len(umi_list), name='peak', dtype=str)],
		axis=1)

	df = df.groupby(['cb', 'peak'])['umi'].size().reset_index()
	df = df.pivot_table(index=['cb'], columns='peak', values='umi').reset_index()# make pivot table
	res = pd.merge(cb_df, df, on=['cb'], how='left')
	res = res.drop('cb', axis=1).fillna(0).astype('int64')
	res = res.transpose()

	return res


def main(args):
	bamfile, cb_file, peakfile, outfile = args
	cb_df = pd.read_csv(cb_file, names=['cb'])
	
	pool = Pool(processes=12)

	peak_region = defaultdict(list)
	with open(peakfile) as ph:
		for peak in ph:
			chrom, s_t, strand = peak.strip().split('\t')
			chrom = 'MT' if chrom == 'M' else chrom
			peak_region[chrom].append((s_t, strand))


	md, hd='w', True
	count_inf = 0
	res_lst = []
	try:
		for chrom, info in peak_region.items():
			for s_t, strand in info:
				arg = [bamfile, chrom, s_t, strand, cb_df]
				res_lst.append(pool.apply_async(run, (arg,)))
	except KeyboardInterrupt:
		logger.info('Caught KeyboardInterrupt, terminating workers')
		pool.terminate()
		sys.exit('KeyboardInterrupt')

	md, hd='w', True
	for df in res_lst:
		try:
			df = df.get()
		except Exception as e:
			logger.debug(f'error found, {e}')
			continue

		if not isinstance(df, pd.DataFrame):
			continue

		if hd:
			df.columns = cb_df.cb.values.tolist()
		df.to_csv(outfile, mode=md, header=hd, compression='gzip')
		md, hd='a', False
	logger.info('All done')



if __name__ == '__main__':
	main(sys.argv[1:])
