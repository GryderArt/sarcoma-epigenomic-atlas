#!/usr/bin/env python3
"""Full EGA dataset catalog -> filter to sarcoma. EGA holds controlled-access human data;
its metadata is open, so we can enumerate what EXISTS even where we cannot download it.
That is exactly what a gap map needs: 'this experiment was done, and it is not reusable'."""
import urllib.request, json, time, csv, sys, os, re, collections
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
OUT = DATA
def get(u,tries=5):
    for a in range(tries):
        try:
            r=urllib.request.Request(u,headers={"Accept":"application/json",
                                                "User-Agent":"sarcoma-atlas/1.0"})
            return urllib.request.urlopen(r,timeout=120).read().decode("utf-8","replace")
        except Exception as e:
            if a==tries-1: sys.stderr.write(f"fail {u}: {e}\n"); return None
            time.sleep(2*(a+1))
all_ds=[]; off=0; LIM=200
while True:
    t=get(f"https://metadata.ega-archive.org/datasets?limit={LIM}&offset={off}")
    if t is None: break
    try: d=json.loads(t)
    except Exception: break
    if not d: break
    all_ds.extend(d); off+=LIM
    if off % 2000 == 0: print(f"  {off} datasets...", flush=True)
    if off > 60000: break
print(f"{len(all_ds)} EGA datasets enumerated", flush=True)
json.dump(all_ds, twrite("ega_datasets_raw.json"))
