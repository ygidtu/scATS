import os
from subprocess import check_call
from multiprocessing import Pool


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
        return ",".join([str(x) for x in self.sites])

    @classmethod
    def create(cls, data: dict):
        return cls(data["utr"], data["infered_sites"], data["gene_name"])


def format_event(idx: str, span: int = 100):
    chrom, sites, strand = idx.split("_")[-1].split(":")
    sites = sites.split("-")

    chrom, sites, strand = idx.split("_")[0].split(":")
    sites = [int(x) for x in sites.split("-")]

    return "{}:{}-{}:{}".format(chrom, sites[0] - span, sites[-1] + span, strand)


def load_results(path: str):
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
                
                # if "ABCD1" != temp["gene_name"]:
                #     continue

                e = Event.create(temp)
                if e:
                    if temp["gene_name"] in res.keys():
                        # print(res[temp["gene_name"]],  temp["gene_name"])
                        res[temp["gene_name"]] + e
                    else:
                        res[temp["gene_name"]] = e
    return list(res.values())


def call(cmd):
    check_call(cmd, shell=True)


def main(bam: str, output: str, res: str, species: str = "human"):
    res = load_results(res)

    cmds = {}
    events = []

    gtf = GTF[species]
    
    for sites in res:
        idx = sites.name
        event = sites.region(100)  

        o = os.path.join(output, "{}.pdf".format(idx))

        # if os.path.exists(o) or os.path.exists(os.path.join(output, "good", "{}.svg".format(event))) or os.path.exists(os.path.join(output, "not", "{}.svg".format(event))):
        #     print("skip:", o)
        #     continue

        cmd = SCRIPT + " plot -S -e {e} --bam {b} --gtf {gtf} -o {o} --color-factor {c} --indicator-lines {site} -p 10 --show-side -t 1000".format(
            e=event, o=o, b = bam, c=3,
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

    with Pool(15) as p:
        p.map(call, [cmds[e] for e in events])


if __name__ == "__main__":
    from fire import Fire
    Fire(main)
