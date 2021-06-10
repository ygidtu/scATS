import pysam

with open("/mnt/raid61/Personal_data/zhangyiming/code/afe/0_200.txt") as r:
    first = [line.strip() for line in r]

with open("/mnt/raid61/Personal_data/zhangyiming/code/afe/200_400.txt") as r:
    second = [line.strip() for line in r]

with pysam.AlignmentFile("/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam") as r:
    fw = pysam.AlignmentFile("/mnt/raid61/Personal_data/zhangyiming/code/afe/0_200.bam", mode="wb", template = r)
    sw = pysam.AlignmentFile("/mnt/raid61/Personal_data/zhangyiming/code/afe/200_400.bam", mode="wb", template = r)

    for rec in r.fetch("1", 1212595, 1214738):
        if rec.query_name in first:
            fw.write(rec)
        elif rec.query_name in second:
            sw.write(rec)

fw.close()
sw.close()
pysam.index("/mnt/raid61/Personal_data/zhangyiming/code/afe/0_200.bam")
pysam.index("/mnt/raid61/Personal_data/zhangyiming/code/afe/200_400.bam")