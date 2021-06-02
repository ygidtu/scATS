kernprof -l ./main.py ats \
    --utr data/utr.bed \
    --utr-length 1000 \
    --output data/test.txt \
    --debug \
    data/test.bam

python -m line_profiler main.py.lprof > profiler.txt