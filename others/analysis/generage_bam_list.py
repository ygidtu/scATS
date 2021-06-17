import os
from glob import glob

def main(indir="/mnt/raid64/ATS/Personal/zhangyiming/bams/"):
    for i in glob(os.path.join(indir, "*.bam")):
        rp = os.path.realpath(i)
        if rp.endswith("possorted_genome_bam.bam"):
            key = os.path.basename(os.path.dirname(os.path.dirname(rp)))
            key = key.split("-")[0]

            print(f"{i}\t{key}")


if __name__ == "__main__":
    main()
