from glob import glob
import requests
from rich import print

url = "https://www.encodeproject.org/search/"

bws  = glob("*/*/*/*.bigWig")

with open("file.csv", "w+") as w:
    for b in bws:
        b = b.replace(".bigWig", "")
        b = b.split("/")

        resp = requests.get(url, params = {"searchTerm": b[-1], "frame": "object", "format": "json"})
        resp = resp.json()
        resp = resp["@graph"]
        resp = resp[0]
        # print(resp)
        b += [resp["assembly"], resp["date_created"], resp["lab"]]
        w.write(f"{','.join(b)}\n")