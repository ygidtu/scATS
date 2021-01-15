import os, json

with open("data.json") as r:
    data = json.load(r)


for key, url in data.items():
    #if os.path.exists(f"SRR{key}.sra"):
    #    continue
    try:
        check_call(f"wget -c -O SRR{key}.sra {url}", shell=True)
    except Exception a:
        continue