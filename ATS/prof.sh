# kernprof -l ./main.py ats \
#     --utr data/utr.bed \
#     --utr-length 1000 \
#     --output data/test.txt \
#     --debug \
#     data/test.bam

kernprof -l test.py
python -m line_profiler test.py.lprof > profiler.txt