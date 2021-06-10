
import os

from glob import glob
from multiprocessing import Pool
from tqdm import tqdm


def call(args):
    b, output = args

    key = os.path.basename(b).replace(".bed", "")

    with open(output, "w+") as w:
        with open(b) as r:
            for line in r:
                line = line.strip().split()
                try:
                    sites = [int(line[x]) for x in [1, 2, 4, 5]]

                    if any([x == -1 for x in sites]):
                        continue
                    
                    w.write(f"{sites[1] - sites[0]},{sites[3] - sites[2]},{sites[2] - sites[1]},{max(sites) - min(sites)},{key}\n")
                except Exception:
                    continue


def main(bed, output):
    os.makedirs(output, exist_ok=True)
    bed = glob(os.path.join(bed, "*.bed"))

    cmds = [[b, os.path.join(output, os.path.basename(b).replace(".bed", ".csv"))] for b in bed]

    with Pool(len(cmds)) as p:
        list(tqdm(p.imap(call, cmds), total=len(cmds)))


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
