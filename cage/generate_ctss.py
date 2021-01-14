import sys
import pysam
from collections import defaultdict


def main(args):
	bamfile, outfile = args
	bam_fh = pysam.AlignmentFile(bamfile)
	freq_dict = defaultdict(int)
	strand_dict = defaultdict(str)

	count_inf = 0
	for line in bam_fh.fetch():
		count_inf += 1

		if count_inf % 100000 == 0:
			print(f'processed {count_inf} line')
		if not line.is_read1:
			continue

		chrom_id = bam_fh.get_reference_name(line.reference_id)

		if chrom_id == 'MT':
			chrom_id = 'M'

		strand = '-' if line.is_reverse else '+'
		leftsite = line.reference_end if line.is_reverse else line.reference_start
		rightsite = leftsite + 1

		tss_id = f'{chrom_id}\t{leftsite}\t{rightsite}\t.'

		freq_dict[tss_id] += 1

		if not strand_dict[tss_id]:
			strand_dict[tss_id] = strand

	with open(outfile, 'w') as fh:
		for tss_id, freq in freq_dict.items():
			fh.write(f'chr{tss_id}\t{freq}\t{strand_dict[tss_id]}\n')


if __name__ == '__main__':
	main(sys.argv[1:])
