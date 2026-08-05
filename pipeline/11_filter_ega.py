#!/usr/bin/env python3
"""Filter the full EGA dataset catalogue (21,279 datasets) to sarcoma, and resolve each
hit's governing DAC and policy. EGA metadata is open even where the data is not, so this
tells us which experiments EXIST but cannot be reused -- the distinction the white paper
needs to draw between 'never done' and 'done but locked'."""
import json, re, csv, collections, urllib.request, time, sys, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
OUT = DATA

SARC = re.compile(
    r"\bsarcom|rhabdomyosarc|\brms\b|ewing|osteosarc|\bgist\b|gastrointestinal stromal|"
    r"chordoma|chondrosarc|chondroblastom|desmoid|kaposi|liposarc|leiomyosarc|angiosarc|"
    r"haemangiosarc|hemangiosarc|\bmpnst\b|peripheral nerve sheath|dermatofibrosarc|"
    r"\bdsrct\b|desmoplastic small round|rhabdoid|\batrt\b|\bpecoma\b|myxofibrosarc|"
    r"solitary fibrous|haemangioendothelio|hemangioendothelio|giant cell tumou?r of bone|"
    r"\bgctb\b|synovial sarc|epithelioid sarc|alveolar soft part|\bASPS\b|clear cell sarc|"
    r"fibromatosis|myoepithelial carcinom|extraskeletal myxoid|\bCIC-DUX|BCOR-CCNB3|"
    r"fibrosarcom|neurofibrosarcom|\bmalignant fibrous histiocytom\b|"
    r"undifferentiated pleomorphic|epithelioid haemangio|intimal sarcom|"
    r"endometrial stromal sarcom|\bLGFMS\b|low[- ]grade fibromyxoid", re.I)
# terms that make a hit a false positive when they are the ONLY reason it matched
NEG = re.compile(r"sarcoidosis|sarcopenia|sarcolemma|sarcoplasm|sarcomere", re.I)

EPI = re.compile(
    r"chip[- ]?seq|chip[- ]?exo|cut&run|cut ?and ?run|cut&tag|atac|dnase|faire|mnase|"
    r"\bhi-?c\b|hichip|micro-?c|4c-?seq|chia-?pet|repli-?seq|methylat|methylome|bisulfite|"
    r"\bwgbs\b|\brrbs\b|mbd-?seq|medip|\bh3k\d|histone|chromatin|epigenom|epigenetic|"
    r"\b450k\b|\bepic\b|infinium|nucleosome|accessib|enhancer|\bctcf\b|\bbrd4\b|"
    r"bisulphite|oxbs|nanopore.*methyl|methyl.*nanopore", re.I)

def get(u, tries=4):
    for a in range(tries):
        try:
            r = urllib.request.Request(u, headers={"Accept": "application/json",
                                                   "User-Agent": "sarcoma-atlas/1.0"})
            return urllib.request.urlopen(r, timeout=90).read().decode("utf-8", "replace")
        except Exception:
            if a == tries - 1: return None
            time.sleep(1.5 * (a + 1))

def main():
    try:
        ds = json.load(topen("ega_datasets_raw.json"))
    except FileNotFoundError:
        raise SystemExit(
            "ega_datasets_raw.json not found. It is the full EGA catalogue; the repository "
            "ships it gzipped under data/. If it is missing, regenerate it with\n"
            "    python3 10_harvest_ega.py\n"
            "(about 110 paged API calls, a few minutes).")
    print(f"{len(ds):,} EGA datasets in the catalogue")

    hits = []
    for d in ds:
        blob = " ".join(str(d.get(k) or "") for k in ("title", "description"))
        techs = " ".join(d.get("technologies") or [])
        types = " ".join(d.get("dataset_types") or [])
        if not SARC.search(blob): continue
        # drop hits that matched only on a sarco- word that is not a sarcoma
        stripped = NEG.sub(" ", blob)
        if not SARC.search(stripped): continue
        # EGA free text carries newlines and tabs; flatten so the TSV stays one row
        # per dataset and downstream readers do not need a quoted-CSV parser
        flat = lambda v, n: " ".join(str(v or "").split())[:n]
        hits.append({
            "accession": d.get("accession_id",""),
            "title": flat(d.get("title"), 300),
            "description": flat(d.get("description"), 600),
            "n_samples": d.get("num_samples") or 0,
            "technologies": flat(techs, 200),
            "dataset_types": flat(types, 200),
            "access_type": d.get("access_type") or "",
            "released": str(d.get("released_date") or "")[:10],
            "policy": d.get("policy_accession_id") or "",
            "is_epigenomic": "Y" if EPI.search(blob + " " + techs + " " + types) else "",
            "url": f"https://ega-archive.org/datasets/{d.get('accession_id','')}",
        })
    print(f"{len(hits)} sarcoma datasets")

    # resolve each distinct policy to its DAC so the paper can name who controls access
    pols = sorted({h["policy"] for h in hits if h["policy"]})
    print(f"resolving {len(pols)} access policies...", flush=True)
    pmap = {}
    for p in pols:
        t = get(f"https://metadata.ega-archive.org/policies/{p}")
        if not t: continue
        try:
            j = json.loads(t)
            if isinstance(j, list): j = j[0] if j else {}
            pmap[p] = {"policy_title": " ".join(str(j.get("title") or "").split())[:160],
                       "dac": j.get("dac_accession_id") or ""}
        except Exception:
            pass
        time.sleep(0.12)
    dacs = sorted({v["dac"] for v in pmap.values() if v.get("dac")})
    dmap = {}
    for dd in dacs:
        t = get(f"https://metadata.ega-archive.org/dacs/{dd}")
        if not t: continue
        try:
            j = json.loads(t)
            if isinstance(j, list): j = j[0] if j else {}
            dmap[dd] = " ".join(str(j.get("title") or j.get("contact_email") or "").split())[:160]
        except Exception:
            pass
        time.sleep(0.12)
    for h in hits:
        p = pmap.get(h["policy"], {})
        h["policy_title"] = p.get("policy_title", "")
        h["dac"] = p.get("dac", "")
        h["dac_name"] = dmap.get(h["dac"], "")

    cols = ["accession","title","n_samples","is_epigenomic","technologies","dataset_types",
            "access_type","released","policy","policy_title","dac","dac_name","description","url"]
    hits.sort(key=lambda h: (h["is_epigenomic"] != "Y", -int(h["n_samples"] or 0)))
    # EGA free text carries newlines and tabs in several fields; flatten everything on the
    # way out so the TSV is exactly one line per dataset and needs no quoted-CSV parser
    hits = [{k: (" ".join(v.split()) if isinstance(v, str) else v) for k, v in h.items()}
            for h in hits]
    with twrite("T14_ega_sarcoma_datasets.tsv") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader(); [w.writerow(h) for h in hits]

    epi = [h for h in hits if h["is_epigenomic"] == "Y"]
    print(f"\n{len(hits)} sarcoma datasets, {len(epi)} epigenomic")
    print(f"total samples: {sum(int(h['n_samples'] or 0) for h in hits):,} "
          f"(epigenomic: {sum(int(h['n_samples'] or 0) for h in epi):,})")
    print("access:", dict(collections.Counter(h["access_type"] for h in hits)))
    print("\nEPIGENOMIC sarcoma datasets in EGA:")
    for h in epi:
        print(f"  {h['accession']:20s} {h['n_samples']:>5} samp  {h['access_type']:11s} "
              f"{h['title'][:78]}")
    print("\ntop DACs:")
    for k, n in collections.Counter(h["dac_name"] or h["dac"] for h in hits).most_common(10):
        print(f"  {n:4d}  {k[:90]}")

if __name__ == "__main__":
    main()
