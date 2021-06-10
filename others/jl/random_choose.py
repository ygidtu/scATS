import random
import pysam
from rich import print
from tqdm import tqdm


def main(input_bam, output_bam):
    print("get total records")

    with pysam.AlignmentFile(input_bam) as r:
        with pysam.AlignmentFile(output_bam, "wb", template = r) as w:
            for idx, rec in tqdm(enumerate(r)):
                if idx % 10 == 0:
                    w.write(rec)


if __name__ ==  '__main__':
    from fire import Fire
    Fire(main)