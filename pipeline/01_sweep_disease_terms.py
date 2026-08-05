#!/usr/bin/env python3
"""Full re-harvest: paediatric + adult sarcoma entities + the bare 'sarcoma' term.

Fixes three failure modes found in the first build:
  1. Pan-sarcoma studies described only as "sarcoma" (e.g. GSE140686, the 1,505-sample
     Koelsche methylation classifier) were never retrieved - the bare term was never queried.
  2. A small number of series were dropped during harvesting (e.g. GSE148724).
  3. Adult sarcoma entities were out of scope entirely.

Writes incrementally so it can be resumed if the container is reclaimed.
"""
import json, time, urllib.parse, urllib.request, sys, csv, os, re
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = WORK   # raw sweep output: large, regenerable, not tracked in git

PAEDIATRIC = {
 "RMS": ['"rhabdomyosarcoma"','"alveolar rhabdomyosarcoma"','"embryonal rhabdomyosarcoma"',
         '"spindle cell rhabdomyosarcoma"','"sclerosing rhabdomyosarcoma"','"PAX3-FOXO1"',
         '"PAX7-FOXO1"','"PAX3-FKHR"','"PAX7-FKHR"'],
 "EWS": ['"Ewing sarcoma"','"Ewings sarcoma"','"EWS-FLI1"','"EWSR1-FLI1"','"EWS-ERG"',
         '"peripheral primitive neuroectodermal tumor"','"Askin tumor"'],
 "EWSlike": ['"CIC-DUX4"','"CIC rearranged sarcoma"','"BCOR-CCNB3"','"BCOR rearranged sarcoma"',
             '"undifferentiated round cell sarcoma"','"Ewing-like sarcoma"','"EWSR1-NFATC2"',
             '"EWSR1-PATZ1"'],
 "OS": ['"osteosarcoma"','"osteogenic sarcoma"'],
 "SS": ['"synovial sarcoma"','"SS18-SSX"','"SYT-SSX"'],
 "MPNST": ['"malignant peripheral nerve sheath tumor"','"MPNST"','"neurofibrosarcoma"'],
 "DSRCT": ['"desmoplastic small round cell tumor"','"EWSR1-WT1"'],
 "ASPS": ['"alveolar soft part sarcoma"','"ASPSCR1-TFE3"'],
 "CCS": ['"clear cell sarcoma"','"EWSR1-ATF1"'],
 "EpiS": ['"epithelioid sarcoma"'],
 "RT": ['"rhabdoid tumor"','"malignant rhabdoid tumor"','"atypical teratoid rhabdoid tumor"',
        '"SMARCB1"AND"tumor"','"INI1 deficient"'],
 "IFS": ['"infantile fibrosarcoma"','"ETV6-NTRK3"','"congenital fibrosarcoma"'],
 "NTRK": ['"NTRK rearranged spindle cell"','"lipofibromatosis-like neural tumor"'],
 "IMT": ['"inflammatory myofibroblastic tumor"','"ALK rearranged"AND"myofibroblastic"'],
 "CHORD": ['"chordoma"','"brachyury"AND"tumor"'],
 "GCTB": ['"giant cell tumor of bone"','"chondroblastoma"','"H3F3A"AND"bone"','"G34W"'],
 "DFSP": ['"dermatofibrosarcoma protuberans"','"COL1A1-PDGFB"'],
 "DESM": ['"desmoid tumor"','"aggressive fibromatosis"'],
 "EMC": ['"extraskeletal myxoid chondrosarcoma"','"EWSR1-NR4A3"'],
 "MYOFIB": ['"infantile myofibromatosis"','"myofibroma"'],
}

ADULT = {
 # --- liposarcoma, the four subtypes explicitly requested ---
 "LPS-DD": ['"dedifferentiated liposarcoma"','"DDLPS"'],
 "LPS-WD": ['"well-differentiated liposarcoma"','"atypical lipomatous tumor"','"WDLPS"'],
 "LPS-MYX": ['"myxoid liposarcoma"','"round cell liposarcoma"','"FUS-DDIT3"','"EWSR1-DDIT3"',
             '"TLS-CHOP"','"myxoid round cell liposarcoma"'],
 "LPS-PLEO": ['"pleomorphic liposarcoma"'],
 "LPS-NOS": ['"liposarcoma"','"MDM2 amplified"AND"sarcoma"'],
 # --- vascular ---
 "ANGIO": ['"angiosarcoma"','"cutaneous angiosarcoma"','"radiation-associated angiosarcoma"',
           '"MYC amplification"AND"angiosarcoma"'],
 "EHE": ['"epithelioid hemangioendothelioma"','"WWTR1-CAMTA1"','"TAZ-CAMTA1"','"YAP1-TFE3"',
         '"kaposiform hemangioendothelioma"','"pseudomyogenic hemangioendothelioma"'],
 "KS": ['"Kaposi sarcoma"'],
 # --- the common adult soft tissue sarcomas ---
 "UPS": ['"undifferentiated pleomorphic sarcoma"','"malignant fibrous histiocytoma"'],
 "MFS": ['"myxofibrosarcoma"'],
 "SFT": ['"solitary fibrous tumor"','"NAB2-STAT6"','"hemangiopericytoma"'],
 "LMS": ['"leiomyosarcoma"','"uterine leiomyosarcoma"'],
 "GIST-ADULT": ['"gastrointestinal stromal tumor"','"KIT mutant"AND"GIST"','"PDGFRA"AND"GIST"',
                '"SDH deficient"AND"GIST"'],
 "FIBRO": ['"fibrosarcoma"','"adult fibrosarcoma"','"myxoinflammatory fibroblastic"'],
 "LGFMS": ['"low grade fibromyxoid sarcoma"','"sclerosing epithelioid fibrosarcoma"',
           '"FUS-CREB3L2"','"EWSR1-CREB3L1"'],
 "PECOMA": ['"PEComa"','"perivascular epithelioid cell tumor"','"malignant PEComa"',
            '"TFE3"AND"PEComa"'],
 "CHS-ADULT": ['"chondrosarcoma"','"dedifferentiated chondrosarcoma"','"mesenchymal chondrosarcoma"',
               '"HEY1-NCOA2"','"IDH1"AND"chondrosarcoma"','"clear cell chondrosarcoma"'],
 "MYOEP": ['"myoepithelial carcinoma"AND"soft tissue"','"EWSR1-POU5F1"'],
 "SYNOV-AD": ['"monophasic synovial sarcoma"','"biphasic synovial sarcoma"'],
 "SPINDLE": ['"spindle cell sarcoma"','"pleomorphic sarcoma"'],
 "ENDOM": ['"endometrial stromal sarcoma"','"JAZF1-SUZ12"','"YWHAE-NUTM2"'],
 "OSTEO-AD": ['"parosteal osteosarcoma"','"periosteal osteosarcoma"','"secondary osteosarcoma"',
              '"Paget"AND"osteosarcoma"'],
 "ALVEOLAR-AD": ['"epithelioid angiosarcoma"'],
 "OTHER-AD": ['"malignant granular cell tumor"','"intimal sarcoma"','"undifferentiated sarcoma"',
              '"NUTM1 rearranged sarcoma"','"CIC-NUTM1"','"sarcoma of unknown primary"'],
 # --- the fix: pan-sarcoma studies that name no entity at all ---
 "PANSARC": ['"sarcoma"','"sarcomas"','"soft tissue sarcoma"','"bone sarcoma"',
             '"sarcoma"AND"methylation"','"mesenchymal tumor"','"mesenchymal neoplasm"'],
}
TERMS = {**PAEDIATRIC, **ADULT}

def eget(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write(f"  retry {i+1} ({e})\n"); time.sleep(3*(i+1))
    return None

def esearch_all(term, retmax=100000):
    ids = []
    q = urllib.parse.quote(term)
    js = eget(f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retmax=0&retmode=json")
    if not js: return ids
    try: cnt = int(json.loads(js)["esearchresult"]["count"])
    except Exception: return ids
    for start in range(0, min(cnt, retmax), 5000):
        j = eget(f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retstart={start}&retmax=5000&retmode=json")
        time.sleep(0.36)
        if not j: continue
        try: ids.extend(json.loads(j)["esearchresult"]["idlist"])
        except Exception: pass
    return ids

def main():
    uid_fam = {}
    counts = {}
    for fam, terms in TERMS.items():
        fu = set()
        for t in terms:
            ids = esearch_all(f"({t}) AND GSE[Entry Type]")
            counts[f"{fam}|{t}"] = len(ids)
            fu.update(ids)
            print(f"{fam:12s} {t:52s} {len(ids):6d}", flush=True)
            time.sleep(0.32)
        for u in fu: uid_fam.setdefault(u, set()).add(fam)
        print(f"  == {fam}: {len(fu)} unique GSE UIDs  (running total {len(uid_fam)})", flush=True)
    json.dump({k: sorted(v) for k, v in uid_fam.items()}, twrite("uid_families_v2.json"))
    json.dump(counts, twrite("term_counts_v2.json"), indent=1)
    print(f"\nTOTAL unique GSE UIDs: {len(uid_fam)}", flush=True)

    # sanity: the series the coauthors named must be present
    for g in ("GSE140686","GSE108028","GSE148724","GSE174376","GSE270506","GSE83726","GSE271176"):
        uid = "200" + g[3:].zfill(6)
        print(f"  check {g}: {'PRESENT' if uid in uid_fam else '*** MISSING ***'}", flush=True)

if __name__ == "__main__":
    main()
