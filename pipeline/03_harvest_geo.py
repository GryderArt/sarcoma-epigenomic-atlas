#!/usr/bin/env python3
"""Harvest series summaries + sample-level metadata for the union of the
disease-term sweep and the model-name sweep. Resumable: skips series already
written, so it survives container reclamation."""
import json, time, urllib.parse, urllib.request, sys, csv, os, gzip, threading, queue, re
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = WORK   # raw GEO records: ~450k samples, not tracked in git

SCOLS = ["uid","accession","title","summary","gdstype","taxon","n_samples","pdat",
         "suppfile","pubmed","gpl","entrytype","families","relations"]
GCOLS = ["gsm","gse","families","title","source_name","characteristics","organism",
         "library_strategy","library_source","library_selection","molecule","instrument",
         "platform","sample_type","suppl_files","status_date","contact_institute"]

def eget(url, data=None, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data.encode() if data else None)
            with urllib.request.urlopen(req, timeout=150) as r:
                return r.read().decode("utf-8","replace")
        except Exception:
            time.sleep(2*(i+1))
    return None

def load_uids():
    fam = {k: set(v) for k, v in json.load(topen("uid_families_v2.json")).items()}
    p = f"{OUT}/model_uid_hits_v2.json"
    if os.path.exists(p):
        for u, ms in json.load(open(p)).items():
            fam.setdefault(u, set()).add("MODELSWEEP:" + ";".join(ms)[:120])
    return fam

def stage_summaries(fam):
    path = f"{OUT}/gse_summary_v2.tsv"
    done = set()
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"): done.add(r["uid"])
    todo = [u for u in fam if u not in done]
    print(f"summaries: {len(done)} done, {len(todo)} to fetch", flush=True)
    mode = "a" if done else "w"
    with open(path, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=SCOLS, delimiter="\t", extrasaction="ignore")
        if not done: w.writeheader()
        for i in range(0, len(todo), 300):
            chunk = todo[i:i+300]
            js = eget(f"{EUTILS}/esummary.fcgi", data="db=gds&retmode=json&id="+",".join(chunk))
            time.sleep(0.4)
            if not js: continue
            try: res = json.loads(js).get("result", {})
            except Exception: continue
            for uid in res.get("uids", []):
                d = res[uid]
                w.writerow({"uid":uid,"accession":d.get("accession",""),
                    "title":(d.get("title") or "").replace("\t"," ").replace("\n"," "),
                    "summary":(d.get("summary") or "").replace("\t"," ").replace("\n"," "),
                    "gdstype":d.get("gdstype",""),"taxon":d.get("taxon",""),
                    "n_samples":d.get("n_samples",""),"pdat":d.get("pdat") or d.get("PDAT",""),
                    "suppfile":d.get("suppFile",""),
                    "pubmed":";".join(str(x) for x in (d.get("pubmedids") or [])),
                    "gpl":d.get("gpl",""),"entrytype":d.get("entryType",""),
                    "families":";".join(sorted(fam[uid]))[:300],
                    "relations":""})
            f.flush()
            if i % 3000 == 0: print(f"  summaries {i+len(chunk)}/{len(todo)}", flush=True)
    return path

def parse_gsm(txt, gse, famstr):
    out=[]; cur=None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            if cur: out.append(cur)
            cur={c:"" for c in GCOLS}; cur["gsm"]=line.split("=",1)[1].strip()
            cur["gse"]=gse; cur["families"]=famstr; cur["_ch"]=[]; cur["_sf"]=[]
        elif cur is None: continue
        elif line.startswith("!Sample_title ="): cur["title"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_source_name"): cur["source_name"]=(cur["source_name"]+" | "+line.split("=",1)[1].strip()).strip(" |")
        elif line.startswith("!Sample_characteristics"): cur["_ch"].append(line.split("=",1)[1].strip())
        elif line.startswith("!Sample_organism"): cur["organism"]=(cur["organism"]+";"+line.split("=",1)[1].strip()).strip(";")
        elif line.startswith("!Sample_library_strategy"): cur["library_strategy"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_library_source"): cur["library_source"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_library_selection"): cur["library_selection"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_molecule"): cur["molecule"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_instrument_model"): cur["instrument"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_platform_id"): cur["platform"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_type ="): cur["sample_type"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_status"): cur["status_date"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_contact_institute"): cur["contact_institute"]=line.split("=",1)[1].strip()
        elif line.startswith("!Sample_supplementary_file"):
            cur["_sf"].append(line.split("=",1)[1].strip())
    if cur: out.append(cur)
    for r in out:
        r["characteristics"]=" ;; ".join(r.pop("_ch"))[:2000]
        r["suppl_files"]=" ;; ".join(r.pop("_sf"))[:1500]
    return out

def stage_gsms():
    spath = f"{OUT}/gse_summary_v2.tsv"
    gpath = f"{OUT}/gsm_records_v2.tsv.gz"
    dpath = f"{OUT}/gsm_done_v2.txt"
    fam = {}
    for r in csv.DictReader(open(spath), delimiter="\t"):
        if r["accession"].startswith("GSE"): fam[r["accession"]] = r["families"]
    done = set()
    if os.path.exists(dpath):
        done = set(open(dpath).read().split())
    todo = [g for g in fam if g not in done]
    print(f"GSM harvest: {len(done)} series done, {len(todo)} to go", flush=True)
    fh = gzip.open(gpath, "at" if done else "wt", newline="")
    w = csv.DictWriter(fh, fieldnames=GCOLS, delimiter="\t", extrasaction="ignore")
    if not done: w.writeheader()
    df = open(dpath, "a")
    lock = threading.Lock(); cnt = {"n":0,"g":0,"f":0}
    q = queue.Queue()
    for g in todo: q.put(g)
    def worker():
        while True:
            try: gse = q.get_nowait()
            except queue.Empty: return
            txt = None
            for i in range(4):
                try:
                    with urllib.request.urlopen(
                      f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief",
                      timeout=200) as r:
                        txt = r.read().decode("utf-8","replace"); break
                except Exception: time.sleep(2*(i+1))
            with lock:
                cnt["n"] += 1
                if txt is None: cnt["f"] += 1
                else:
                    rs = parse_gsm(txt, gse, fam.get(gse,""))
                    cnt["g"] += len(rs)
                    for x in rs: w.writerow(x)
                    df.write(gse+"\n")
                if cnt["n"] % 100 == 0:
                    print(f"  {cnt['n']}/{len(todo)} series, {cnt['g']} GSM, {cnt['f']} fail", flush=True)
                    fh.flush(); df.flush()
            time.sleep(0.33)
    ths=[threading.Thread(target=worker) for _ in range(3)]
    [t.start() for t in ths]; [t.join() for t in ths]
    fh.close(); df.close()
    print(f"DONE {cnt['n']} series, {cnt['g']} samples, {cnt['f']} failures", flush=True)

if __name__ == "__main__":
    fam = load_uids()
    print(f"{len(fam)} unique UIDs in union", flush=True)
    stage_summaries(fam)
    stage_gsms()
