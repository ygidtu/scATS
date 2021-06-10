import os
from  glob import glob
from subprocess import check_call
import pysam

utr = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes_utr.bed"
output = "/mnt/raid64/ATS/Personal/zhangyiming/infered"

fs = glob("/mnt/raid64/ATS/Personal/zhangyiming/bams/*.bam")
# print(fs)

rbam = "/mnt/raid64/Covid19_Gravida/apamix/bam/"


class Lock:
    def __init__(self, path):
        self.path = path
        self.l = self.path + ".lock"
    
    def lock(self):
        if not os.path.exists(self.l):
            with open(self.l, "w+") as w:
                w.write(str(os.getpid()))

    def unlock(self):
        if os.path.exists(self.l):
            os.remove(self.l)

    def is_locked(self):
        return os.path.exists(self.l)


for f in fs:
    key = os.path.basename(f).replace(".bam", "")

    path = os.path.realpath(f)

    n = os.path.basename(os.path.dirname(os.path.dirname(path)))
    n = n.split("-")[0]


    o = os.path.join(output, f"{key}.txt")

    b = os.path.join(rbam, n + ".bam")

    if os.path.exists(b):
        l = Lock(o)
        if not l.is_locked():
            print(key)
            l.lock()
            cmd = f"julia /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/run.jl -i {utr} -o {o} --using-R 15 --cage-mode {b}"
            print(cmd)
            check_call(cmd, shell=True)
            # l.unlock()

# en_arr = []
# st_arr = []
# with pysam.AlignmentFile("/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam") as r:
#     for rec in r.fetch("1", 6205375, 6206201):

#         if rec.is_qcfail or (rec.has_tag("NH") and rec.get_tag("NH") > 1):
#             continue
        
#         if rec.is_read1:
#             continue
        
#         utr_site = 6206201

#         if rec.mate_is_unmapped:
#             continue
        
#         # R2 needs locate in UTR

#         if 6205375 > rec.reference_start + 1 or 6206201 < rec.reference_end + 1:
#             continue

#         en_arr.append(rec.reference_end)

#         r1 = r.mate(rec)

#         st_arr.append(r1.reference_start)

# print(st_arr)
# print(en_arr)
