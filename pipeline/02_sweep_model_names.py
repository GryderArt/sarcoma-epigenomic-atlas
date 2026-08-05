#!/usr/bin/env python3
"""Second sweep keyed on MODEL NAMES (paediatric 427 + adult 156), because many
landmark datasets name only the cell line in their title."""
import json, time, urllib.parse, urllib.request, sys, csv, os, re
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = WORK; T = DATA; O = DATA

GENERIC = {"WT","OS","ES","RT","SS","MRT","CS","NY","HAL","VK","SIM","MIC","POE","CTR","MOS",
           "KAS","GIST","EWS","RMS","SARC","PDX","NA","T","P","M","OSA","COL","2T","RD","HOS",
           "SARG","MHM","CCA","ORS","JJ","CG-1","CS1","CS-1","G292","NONE","UNKNOWN","SIMON",
           "HAMON","DUNN","1273","ES-2","ES2","BIRCH","FUJI","GCT","SKN","MES","LP6","SVR",
           "AS-M","ASM","GOT3","T449","T778","LIS-3","SD-437A"}
LOWPREC = {"SIM/EW27","Dunn osteosarcoma","1273/99","ES2","HAMON","GCT","HT-1080"}

def load_names():
    names = {}
    files = [("RMS_models.tsv",), ("OS_models.tsv",), ("EWS_models.tsv",), ("STS_models.tsv",),
             ("ADULT_models.tsv",)]
    srcs = [os.path.join(T, f[0]) for f in files]
    srcs.append(os.path.join(O, "T2_models.tsv"))
    for p in srcs:
        if not os.path.exists(p): continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            canon = (r.get("model_name") or "").strip()
            if not canon or "NO MODEL" in canon.upper() or canon.upper().startswith("NONE"): continue
            if len(canon) > 45: continue
            if canon in LOWPREC: continue
            cand = [canon] + [a.strip() for a in re.split(r"[;,]", r.get("aliases") or "") if a.strip()]
            for c in cand:
                c = c.strip().strip("()")
                if len(c) < 4 or c.upper() in GENERIC: continue
                if re.fullmatch(r"[A-Za-z]{1,3}", c): continue
                if c.lower() in ("unknown","aliases","not registered","n/a","none"): continue
                names.setdefault(c, canon)
    return names

def eget(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read().decode("utf-8","replace")
        except Exception: time.sleep(2*(i+1))
    return None

def main():
    names = load_names()
    print(f"{len(names)} model-name queries", flush=True)
    hits = {}
    for i, (nm, canon) in enumerate(sorted(names.items())):
        js = eget(f"{EUTILS}/esearch.fcgi?db=gds&term="
                  f"{urllib.parse.quote(chr(34)+nm+chr(34)+'[All Fields] AND GSE[Entry Type]')}"
                  f"&retmax=3000&retmode=json")
        time.sleep(0.34)
        if not js: continue
        try:
            r = json.loads(js)["esearchresult"]
            ids = r.get("idlist", []); n = int(r.get("count", 0))
        except Exception: continue
        if n > 2000:
            print(f"  SKIP too broad ({n}) {nm}", flush=True); continue
        for uid in ids: hits.setdefault(uid, set()).add(canon)
        if n: print(f"  {nm:34s} -> {n:5d}", flush=True)
        if i % 60 == 0: print(f"[{i}/{len(names)}] running UIDs {len(hits)}", flush=True)
    json.dump({k: sorted(v) for k, v in hits.items()}, twrite("model_uid_hits_v2.json"))
    prev = set(json.load(topen("uid_families_v2.json")))
    new = set(hits) - prev
    print(f"\nmodel-sweep UIDs {len(hits)}; NEW beyond disease sweep: {len(new)}", flush=True)

if __name__ == "__main__":
    main()
