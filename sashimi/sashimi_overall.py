import os
import re
import signal
import sys
from multiprocessing import Pool
from subprocess import CalledProcessError, check_call

from rich import print
from rich.progress import (BarColumn, Progress, TextColumn, TimeElapsedColumn,
                           TimeRemainingColumn, TransferSpeedColumn, track)


def custom_progress(io: bool = False):
    if io:
        return Progress(
            "[progress.description]{task.description}",
            BarColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TextColumn("| Remaining:"),
            TimeRemainingColumn(),
            TextColumn("| Speed: "),
            TransferSpeedColumn()
        )
    return Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}% ({task.completed}/{task.total})",
        TextColumn("| Elapsed:"),
        TimeElapsedColumn(),
        TextColumn("| Remaining:"),
        TimeRemainingColumn(),
    )


GTF = {
    "human": "/mnt/raid64/ATS/alignments/cellranger/ref/Homo_sapiens/genes/genes.gtf", "mouse": "/mnt/raid64/ATS/alignments/cellranger/ref/Mus_musculus/genes/genes.gtf",
    "gfp": "/mnt/raid64/ATS/alignments/cellranger/ref/refdata-gex-mm10-2020-A/genes/genes.gtf"
}

SCRIPT = "python /mnt/raid61/Personal_data/zhangyiming/code/pysashimi/main.py"


class Event(object):

    def __init__(self, region: str, sites: str, name: str):
        self.sites = sorted([int(x) for x in sites.split(",")])
        chrom, sites, strand = region.split("_")[-1].split(":")
        self.utr = sorted([int(x) for x in sites.split("-")])
        self.chrom = chrom
        self.strand = strand
        self.name = name
        self.ws = []

    def __len__(self):
        return len(self.sites)

    def __hash__(self):
        return hash(self.name)

    def __add__(self, other):
        self.utr += other.utr
        self.utr.sort()
        self.sites += other.sites
        self.sites.sort()

    def region(self, span: int = 0):
        return f"{self.chrom}:{self.utr[0] - span}-{self.utr[-1] + span}:{self.strand}"

    def __str__(self):
        return ",".join([str(x) for x, y in zip(self.sites, self.ws) if x > 0])

    @classmethod
    def create(cls, data: dict):
        e = cls(data["utr"], data["infered_sites"], data["gene_name"])

        for i in data["ws"].split(","):
            e.ws.append(float(i))

        return e


def format_event(idx: str, span: int = 100):
    chrom, sites, strand = idx.split("_")[-1].split(":")
    sites = sites.split("-")

    chrom, sites, strand = idx.split("_")[0].split(":")
    sites = [int(x) for x in sites.split("-")]

    return "{}:{}-{}:{}".format(chrom, sites[0] - span, sites[-1] + span, strand)


def __decode_attr__(attrs):
    data = {}
    for line in attrs.split(";"):
        if line:
            key, val = line.split(' "')
            key = re.sub(r'[";\s]', '', key)
            val = re.sub(r'[";\s]', '', val)
            data[key] = val
    return data


def load_gtf(path: str):
    res = {}
    progress = custom_progress(io=True)

    with progress:
        task = progress.add_task("Loading...", total=os.path.getsize(path))
        with open(path) as r:
            for line in r:
                progress.update(task, advance=len(str.encode(line)))
                if line.startswith("#"):
                    continue

                line = line.split()

                if line[2] == "transcript":
                    attrs = __decode_attr__(" ".join(line[8:]))

                    key = attrs.get("gene_name", attrs.get("Parent"))

                    if key not in res.keys():
                        res[key] = []

                    res[key].append(int(line[3]) if line[6]
                                    == "+" else int(line[4]))

    return {x: y for x, y in track(res.items()) if len(y) > 1}


def is_all_in(ats, tss, distance=100):
    ats = sorted(ats)
    tss = sorted(tss)

    match = 0
    i, j = 0, 0
    while i < len(ats) and j < len(tss):
        curr_ats = ats[i]
        curr_tss = tss[j]

        if curr_ats < curr_tss - distance:
            i += 1
        elif curr_ats > curr_tss + distance:
            j += 1
        else:
            match += 1
            i += 1

    return match == len(ats)


def load_results(path: str, tss, distance=100):
    res = {}

    header = None
    with open(path) as r:
        for line in r:
            line = line.split()

            if not header:
                header = line
            else:
                temp = {x: y for x, y in zip(header, line)}

                if "gene_name" not in temp.keys():
                    temp["gene_name"] = temp["reference_id"]

                if "," in temp["gene_name"]:
                    continue

                # if "ABRAXAS" not in temp["gene_name"]:
                #     continue

                e = Event.create(temp)
                if e:
                    if temp["gene_name"] in res.keys():
                        # print(res[temp["gene_name"]],  temp["gene_name"])
                        res[temp["gene_name"]] + e
                    else:
                        res[temp["gene_name"]] = e
    return [y for x, y in res.items() if len(y) > 1 and is_all_in(y.sites, tss.get(x, []), distance)]


def call(cmd):
    with open(os.devnull, "w") as w:
        try:
            check_call(cmd, shell=True, stdout=w, stderr=w)
        except CalledProcessError as err:
            print(err)


def exists(path):
    for _, _, files in os.walk(os.path.dirname(path)):
        for f in files:
            if os.path.basename(f) == os.path.basename(path):
                return True
    return False


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def main(
    bam: str,
    output: str,
    res: str,
    species: str = "human",
    distance: int = 100,
    n_jobs: int = 15
):
    gtf = GTF[species]
    tss = load_gtf(gtf)

    res = load_results(res, tss, distance)

    cmds = {}
    events = []

    for sites in res:
        idx = sites.name
        event = sites.region(100)

        o = os.path.join(output, "{}.pdf".format(idx))

        if os.path.exists(o):
            print("skip:", o)
            continue

        cmd = SCRIPT + " plot -S -e {e} --bam {b} --gtf {gtf} -o {o} --color-factor {c} --indicator-lines {site} -p 10 --show-side -t 1000".format(
            e=event, o=o, b=bam, c=3,
            site=sites, gtf=gtf
        )

        if event not in cmds.keys():
            events.append(event)
        cmds[event] = cmd

        # print(cmd)
        # break

    # temp = sorted(events)

    # if len(temp) > 5000:
    #     random.seed(42)
    #     temp = random.choices(temp, k = 5000)

    try:
        with Pool(n_jobs, init_worker) as p:
            list(track(p.imap(call, [cmds[e]
                 for e in events]), total=len(events)))
    except KeyboardInterrupt as err:
        print(err)
        sys.exit(0)


if __name__ == "__main__":
    from fire import Fire
    Fire(main)
