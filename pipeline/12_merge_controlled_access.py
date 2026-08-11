#!/usr/bin/env python3
"""Map EGA (controlled-access) and St Jude CSTN holdings onto the atlas entities, and
separate three states that a GEO-only survey conflates:

  PUBLIC      -- open deposit, anyone can download and reuse
  CONTROLLED  -- the experiment was done and the metadata is public, but the data sits
                 behind a DAC/DUA and cannot be reanalysed without per-dataset approval
  ABSENT      -- no record of the experiment anywhere

For a white paper about what the field can actually build on, CONTROLLED is not the same
as PUBLIC and it is not the same as ABSENT.
"""
import csv, re, collections, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

# entity matching, most specific first -- first match wins
ENT = [
 (r"atypical teratoid|\batrt\b|rhabdoid|\bmrt\b|\bsmarcb1\b|\bini1\b", "Rhabdoid tumor/ATRT"),
 (r"tfcp2|rhabdomyosarc|\brms\b|\berms\b|\barms\b", "RMS (any subtype)"),
 (r"ewing|\bES\b tumors|ewing'?s", "Ewing"),
 (r"osteosarc|\bsaos\b|\bSJOS\b", "Osteosarcoma"),
 (r"synovial", "Synovial sarcoma"),
 (r"\bmpnst\b|peripheral nerve sheath", "MPNST"),
 (r"chordoma", "Chordoma"),
 (r"chondrosarc|\bCSA\b|mesenchymal_CSA", "Chondrosarcoma"),
 (r"giant cell tumou?r of bone|\bgctb\b", "GCTB/Chondroblastoma"),
 (r"epithelioid sarcom|\bEpS\b", "Epithelioid sarcoma"),
 (r"solitary fibrous|hemangiopericytoma", "Solitary fibrous tumour"),
 (r"desmoplastic small round|\bdsrct\b", "DSRCT"),
 (r"clear cell sarcom", "Clear cell sarcoma"),
 (r"liposarc", "Liposarcoma (any subtype)"),
 (r"leiomyosarc", "Leiomyosarcoma"),
 (r"angiosarc|haemangiosarc", "Angiosarcoma"),
 (r"\bgist\b|gastrointestinal stromal", "GIST"),
 (r"kaposi", "Kaposi sarcoma"),
 (r"undifferentiated (pleomorphic )?sarcom|\bUPS\b", "UPS / undifferentiated sarcoma"),
 (r"\bCDS\b and \bES\b|CIC-DUX|\bCIC\b", "CIC-DUX4"),
 (r"desmoid|fibromatosis", "Desmoid"),
 (r"sarcom", "Sarcoma, unspecified / mixed cohort"),
]
COMP = [(re.compile(p, re.I), e) for p, e in ENT]
def entity(text):
    for rx, e in COMP:
        if rx.search(text): return e
    return "Sarcoma, unspecified / mixed cohort"

REG = re.compile(r"chip[- ]?seq|chip[- ]?exo|cut&?\s?run|cut&amp;run|cut&?tag|\batac\b|"
                 r"atacseq|atac-seq|dnase|faire|mnase|\bhi-?c\b|hichip|micro-?c|4c-?seq|"
                 r"chia-?pet|\bctcf\b|\bbrd4\b|histone|chromatin|\bh3k\d", re.I)
DNAME = re.compile(r"methylat|methylome|bisulfite|bisulphite|\bwgbs\b|\brrbs\b|mbd-?seq|"
                   r"medip|\b450k\b|\b850k\b|infinium|\bepic\b|oxbs", re.I)

def main():
    ega = list(csv.DictReader(topen("T14_ega_sarcoma_datasets.tsv"), delimiter="\t"))
    rows = []
    for d in ega:
        blob = " ".join((d["title"], d["description"], d["technologies"], d["dataset_types"]))
        reg = bool(REG.search(blob)); dna = bool(DNAME.search(blob))
        if not (reg or dna) and d["is_epigenomic"] != "Y": continue
        rows.append({
            "source": "EGA", "accession": d["accession"], "entity": entity(blob),
            "assay_family": "regulatory" if reg else ("DNA methylation" if dna else "other"),
            "n_samples": int(d["n_samples"] or 0),
            "access": "CONTROLLED", "governing_body": d["dac_name"] or d["dac"],
            "title": d["title"], "url": d["url"], "released": d["released"]})

    # St Jude CSTN -- the epigenomic rows of the inventory the agent assembled
    if os.path.exists(f"{O}/T16_stjude_cstn_inventory.tsv"):
        for d in csv.DictReader(topen("T16_stjude_cstn_inventory.tsv"), delimiter="\t"):
            blob = " ".join(str(v) for v in d.values())
            if not (REG.search(blob) or DNAME.search(blob)): continue
            if d.get("repository","").upper().startswith("EGA"): continue   # already above
            n = re.search(r"\d+", str(d.get("n_samples") or ""))
            rows.append({
                "source": "St Jude CSTN", "accession": d.get("accession",""),
                "entity": entity(str(d.get("diseases","")) + " " + str(d.get("resource",""))),
                "assay_family": "regulatory" if REG.search(str(d.get("assay",""))) else "DNA methylation",
                "n_samples": int(n.group()) if n else 0,
                "access": ("OPEN" if d.get("repository","").upper().startswith("GEO")
                           else "CONTROLLED"),
                "governing_body": d.get("access",""), "title": d.get("resource",""),
                "url": "", "released": ""})

    # ---- CCDI (NCI Childhood Cancer Data Initiative), dbGaP-controlled
    ccdi = os.path.join(DATA, "T24_ccdi_entity_counts.tsv")
    if os.path.exists(ccdi):
        n_ccdi = 0
        for c in csv.DictReader(open(ccdi, newline=""), delimiter="\t"):
            n = int(c["ccdi_participants_with_methylation"] or 0)
            if not n: continue
            rows.append({
                "source": "CCDI", "accession": "dbGaP (multiple studies)",
                "entity": c["atlas_entity"], "assay_family": "DNA methylation",
                "n_samples": n, "access": "CONTROLLED",
                "governing_body": "NCI CCDI / dbGaP Data Access Committee",
                "title": f"CCDI methylation arrays, {c['ccdi_participants']} participants "
                         f"/ {c['ccdi_samples']} samples in cohort",
                "url": "https://ccdi.cancer.gov/explore", "released": ""})
            n_ccdi += n
        print(f"  + CCDI: {n_ccdi:,} participants with methylation arrays across "
              f"{sum(1 for r in rows if r['source']=='CCDI')} entities")
        print("    CCDI regulatory epigenomics for sarcoma: 0 "
              "(no ChIP-seq/ATAC/bisulfite in the sarcoma cohort at all)")

    rows.sort(key=lambda r: (r["entity"], r["assay_family"] != "regulatory", -r["n_samples"]))
    cols = ["source","accession","entity","assay_family","n_samples","access",
            "governing_body","title","released","url"]
    with twrite("T15_controlled_access.tsv") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]

    print(f"{len(rows)} controlled/closed epigenomic datasets mapped to entities\n")
    agg = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        agg[r["entity"]]["ds"] += 1
        agg[r["entity"]][r["assay_family"]] += r["n_samples"]
    print(f"{'entity':38s} {'datasets':>8s} {'regulatory':>11s} {'methylation':>12s}")
    for e in sorted(agg, key=lambda e: -(agg[e]['regulatory'] + agg[e]['DNA methylation'])):
        a = agg[e]
        print(f"{e[:37]:38s} {a['ds']:8d} {a['regulatory']:11d} {a['DNA methylation']:12d}")
    print(f"\ntotal controlled epigenomic samples: "
          f"{sum(r['n_samples'] for r in rows):,} "
          f"(regulatory {sum(r['n_samples'] for r in rows if r['assay_family']=='regulatory'):,})")

if __name__ == "__main__":
    main()
