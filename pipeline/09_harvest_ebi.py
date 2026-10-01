#!/usr/bin/env python3
"""Harvest the EBI side of the ledger (meeting note #5).

Four resources, all reachable without credentials:
  ArrayExpress / BioStudies  https://www.ebi.ac.uk/biostudies/api/v1/{collection}/search
  ENA portal                 https://www.ebi.ac.uk/ena/portal/api/search
  EGA metadata               https://metadata.ega-archive.org/...      (metadata open, data controlled)
  EpiRR (IHEC registry)      https://www.ebi.ac.uk/vg/epirr/view/all?format=tsv   (returns JSON)

The point is to find what is NOT already in GEO. E-GEOD-* and most ENA/SRA records are
mirrors of GEO; E-MTAB-*, EGAS/EGAD and EpiRR entries are the genuinely additional ones.
"""
import json, urllib.request, urllib.parse, time, csv, os, re, sys, collections
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
OUT = DATA

TERMS = ["sarcoma","rhabdomyosarcoma","Ewing sarcoma","osteosarcoma","synovial sarcoma",
         "liposarcoma","leiomyosarcoma","angiosarcoma","chondrosarcoma","chordoma",
         "MPNST","malignant peripheral nerve sheath tumor","desmoplastic small round cell",
         "rhabdoid tumor","rhabdoid tumor","gastrointestinal stromal tumor",
         "gastrointestinal stromal tumor","epithelioid sarcoma","alveolar soft part sarcoma",
         "clear cell sarcoma","epithelioid hemangioendothelioma","epithelioid hemangioendothelioma",
         "undifferentiated pleomorphic sarcoma","myxofibrosarcoma","solitary fibrous tumor",
         "dermatofibrosarcoma","desmoid","fibrosarcoma","Kaposi sarcoma","PEComa",
         "endometrial stromal sarcoma","giant cell tumor of bone","chondroblastoma"]
EPI_HINT = re.compile(r"chip[-\s]?seq|atac|cut&run|cut&tag|cut ?and ?run|hi-?c|hichip|methylat|"
                      r"bisulfite|h3k\d|histone|chromatin|epigenom|enhancer|dnase|wgbs|rrbs|"
                      r"nucleosome|accessib", re.I)

def get(url, tries=4, accept="application/json"):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept,
                                                       "User-Agent": "sarcoma-atlas/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if a == tries-1: sys.stderr.write(f"  fail {url[:80]}: {e}\n")
            time.sleep(2*(a+1))
    return None

rows = []

# ---------------------------------------------------------------- ArrayExpress
print("ArrayExpress / BioStudies", flush=True)
seen = set()
for t in TERMS:
    off = 0
    while True:
        u = ("https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search?"
             f"query={urllib.parse.quote(t)}&pageSize=100&page={off+1}")
        js = get(u)
        if not js: break
        try: d = json.loads(js)
        except Exception: break
        hits = d.get("hits", [])
        for h in hits:
            acc = h.get("accession","")
            if not acc or acc in seen: continue
            seen.add(acc)
            title = (h.get("title") or "")
            rows.append({"repository":"ArrayExpress/BioStudies","accession":acc,
                "title":title[:300],"query_term":t,
                "is_geo_mirror":"Y" if acc.startswith("E-GEOD") else "",
                "epigenomic_hint":"Y" if EPI_HINT.search(title) else "",
                "release_date":str(h.get("release_date") or ""),
                "n_files":"", "access":"open",
                "url":f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}"})
        off += 1
        if len(hits) < 100 or off > 8: break
    time.sleep(0.25)
print(f"  {len(seen)} unique ArrayExpress accessions", flush=True)

# ---------------------------------------------------------------- ENA
print("ENA portal", flush=True)
ena = set()
for t in TERMS:
    q = urllib.parse.quote(f'study_title="*{t}*"')
    u = (f"https://www.ebi.ac.uk/ena/portal/api/search?result=study&query={q}"
         f"&fields=study_accession,study_title,first_public&limit=1000&format=tsv")
    txt = get(u, accept="text/plain")
    if not txt: continue
    for line in txt.splitlines()[1:]:
        p = line.split("\t")
        if len(p) < 2 or p[0] in ena: continue
        ena.add(p[0])
        rows.append({"repository":"ENA","accession":p[0],"title":p[1][:300],"query_term":t,
            "is_geo_mirror":"", "epigenomic_hint":"Y" if EPI_HINT.search(p[1]) else "",
            "release_date":p[2] if len(p)>2 else "","n_files":"","access":"open",
            "url":f"https://www.ebi.ac.uk/ena/browser/view/{p[0]}"})
    time.sleep(0.25)
print(f"  {len(ena)} unique ENA studies", flush=True)

# ---------------------------------------------------------------- EGA
print("EGA metadata", flush=True)
ega = set()
for t in TERMS:
    u = ("https://metadata.ega-archive.org/studies?limit=100&query=" + urllib.parse.quote(t))
    js = get(u)
    if not js: continue
    try: d = json.loads(js)
    except Exception: continue
    if isinstance(d, dict): d = d.get("results") or d.get("data") or []
    for h in d:
        acc = h.get("accession_id") or h.get("accession") or ""
        if not acc or acc in ega or acc.startswith("DUMMY"): continue
        title = (h.get("title") or "")
        if not re.search(r"sarcom|rhabdo|ewing|osteosarc|gist|chordoma|chondrosarc|desmoid|"
                         r"kaposi|liposarc|leiomyosarc|angiosarc|pecoma|mpnst", title, re.I):
            continue
        ega.add(acc)
        rows.append({"repository":"EGA","accession":acc,"title":title[:300],"query_term":t,
            "is_geo_mirror":"", "epigenomic_hint":"Y" if EPI_HINT.search(title) else "",
            "release_date":str(h.get("released_date") or "")[:10],"n_files":"",
            "access":"CONTROLLED (DAC approval required)",
            "url":f"https://ega-archive.org/studies/{acc}"})
    time.sleep(0.25)
print(f"  {len(ega)} unique EGA studies", flush=True)

# ---------------------------------------------------------------- EpiRR
print("EpiRR (IHEC reference epigenomes)", flush=True)
js = get("https://www.ebi.ac.uk/vg/epirr/view/all?format=tsv")
n_epirr = 0
if js:
    try: d = json.loads(js)
    except Exception: d = []
    print(f"  {len(d)} reference epigenomes in the registry", flush=True)
    for e in d:
        blob = " ".join(str(e.get(k,"")) for k in
                        ("description","local_name","project","full_accession"))
        if not re.search(r"sarcom|rhabdo|ewing|osteosarc|gist|chordoma|chondrosarc|"
                         r"kaposi|liposarc|leiomyosarc|angiosarc|mpnst|desmoid", blob, re.I):
            continue
        n_epirr += 1
        rows.append({"repository":"EpiRR/IHEC","accession":e.get("full_accession",""),
            "title":(e.get("description") or e.get("local_name") or "")[:300],
            "query_term":"registry scan","is_geo_mirror":"",
            "epigenomic_hint":"Y","release_date":"",
            "n_files":str(e.get("status","")),
            "access":"open (points to underlying ENA/GEO/EGA)",
            "url":f"https://www.ebi.ac.uk/vg/epirr/view/{e.get('full_accession','')}"})
    print(f"  {n_epirr} sarcoma-matching reference epigenomes", flush=True)

cols = ["repository","accession","title","query_term","is_geo_mirror","epigenomic_hint",
        "release_date","n_files","access","url"]
with twrite("T12_ebi_all.tsv") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader(); [w.writerow(r) for r in rows]
print(f"\nwrote {len(rows)} rows -> T12_ebi_all.tsv")
c = collections.Counter(r["repository"] for r in rows)
for k,v in c.most_common(): print(f"  {v:6d}  {k}")
native = [r for r in rows if r["repository"]=="ArrayExpress/BioStudies" and not r["is_geo_mirror"]]
epi_native = [r for r in native if r["epigenomic_hint"]=="Y"]
print(f"\nArrayExpress NATIVE (not GEO mirrors): {len(native)}  of which epigenomic-looking: {len(epi_native)}")
for r in epi_native[:12]: print(f"   {r['accession']:16s} {r['title'][:72]}")
