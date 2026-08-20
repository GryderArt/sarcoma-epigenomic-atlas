#!/usr/bin/env python3
"""Enumerate the St Jude CSTN epigenetic browsers at sample x assay resolution.

The CSTN "epigenetic landscape" pages on viz.stjude.cloud are ProteinPaint embeds. The
page ships its whole track manifest inline -- every sample, every assay, every underlying
file -- so the browser contents can be enumerated exactly, which the CSTN portal's own
Shiny interface does not allow.

What this is for: EGA describes these deposits as, for example, "ChIP-Seq files for RMS,
242 samples". That is opaque. The manifest resolves the same material to the atlas's own
grain -- H3K27ac for SJRHB013758_X1 -- so the gap map can say which marks exist on which
model rather than only how many files there are.

ACCESS: the tracks are viewable in the browser; the underlying bigWig and bed files sit
behind Cloudflare bot protection and, for real access, the St Jude / EGA data access
agreement. They are recorded as "formal request" and are NOT added to any sample total,
because the same material is already counted once under EGA.
"""
import urllib.request, html, re, json, csv, os, collections

from _paths import DATA, twrite

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/120 Safari/537.36"}

PAGES = {
    "RHB": ("Epigenetic landscape of rhabdomyosarcoma subtypes",
            "https://viz.stjude.cloud/st-jude-childrens-research-hospital/visualization/"
            "epigenetic-landscape-of-rhabdomyosarcoma-subtypes~67"),
    "SJOS2018": ("Epigenetic landscape of xenografts from osteosarcoma and rare tumors",
                 "https://viz.stjude.cloud/st-jude-childrens-research-hospital/visualization/"
                 "epigenetic-landscape-of-xenografts-from-osteosarcoma-and-rare-tumors~69"),
}

# assay -> atlas assay class / target
REG_ASSAYS = {"BRD4", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3",
              "H3K9-14Ac", "H3K9me3", "RNAPolII", "CTCF", "INPUT"}
DERIVED = {"SE", "SuperEnhancer(SE)", "SE.noK4me3"}      # called peaks, not raw assays
DNAME_ASSAYS = {"WGBS"}

ENTITY = [
    (r"\(ERMS\)|_ERMS", "FN-RMS"), (r"\(ARMS\)|_ARMS", "FP-RMS"),
    (r"\(SCLEROS\)", "RMS-MYOD1"), (r"^SJRHB", "RMS-NOS"),
    (r"^SJOS", "Osteosarcoma"), (r"^SJEWS", "Ewing"), (r"^SJLPS", "Liposarcoma-NOS"),
    (r"^SJHGS", "UPS/MFH"), (r"^SJNORM", "normal/reference"),
]
ENTITY = [(re.compile(p, re.I), e) for p, e in ENTITY]


def entity_of(label):
    for rx, e in ENTITY:
        if rx.search(label):
            return e
    return "Sarcoma NOS"


def fetch(url):
    r = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(r, timeout=120).read().decode("utf-8", "replace")


def extract_config(body):
    i = body.find("runproteinpaint(")
    if i < 0:
        raise RuntimeError("no runproteinpaint( call on the page")
    seg = html.unescape(body[i + len("runproteinpaint("):])
    start = seg.find("{")
    depth, end = 0, None
    for j, ch in enumerate(seg[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    raw = seg[start:end]
    # the holder is a live DOM handle, not JSON
    raw = re.sub(r'"holder"\s*:\s*[A-Za-z0-9_\[\]\(\)\.]+\s*,', '"holder": null,', raw)
    return json.loads(raw)


def main():
    rows = []
    for study, (title, url) in PAGES.items():
        cfg = extract_config(fetch(url))
        sv = cfg["studyview"]
        for assay in sv.get("assays", []):
            blk = sv.get(assay)
            if not isinstance(blk, dict):
                continue
            for sample, inner in blk.items():
                if sample == "config" or not isinstance(inner, dict):
                    continue
                for t in inner.values():
                    if not isinstance(t, dict):
                        continue
                    fam = ("regulatory" if assay in REG_ASSAYS else
                           "DNA methylation" if assay in DNAME_ASSAYS else
                           "derived (called peaks)")
                    rows.append({
                        "study": study, "viz_title": title,
                        "sample_label": sample,
                        "sample_id": t.get("sampleID", ""),
                        "atlas_entity": entity_of(sample),
                        "assay": assay, "assay_family": fam,
                        "is_control": "Y" if assay == "INPUT" else "",
                        "genome": sv.get("genome", ""),
                        "track_type": t.get("type", ""),
                        "file": t.get("file", "") or t.get("url", ""),
                        "access": "FORMAL REQUEST (St Jude CSTN / EGA data access agreement)",
                        "counted_in_totals": "no - same material as the EGA deposit",
                        "viz_url": url,
                    })
    with twrite("T28_stjude_viz_tracks.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); [w.writerow(r) for r in rows]

    print(f"{len(rows)} tracks across {len({r['sample_label'] for r in rows})} samples")
    for study in PAGES:
        g = [r for r in rows if r["study"] == study]
        print(f"  {study:9s} {len(g):4d} tracks · {len({r['sample_label'] for r in g})} samples "
              f"· {len({r['assay'] for r in g})} assays")
    print("\nassay family:", dict(collections.Counter(r["assay_family"] for r in rows)))
    print("\nby atlas entity (samples):")
    ent = collections.defaultdict(set)
    for r in rows: ent[r["atlas_entity"]].add(r["sample_label"])
    for e, s in sorted(ent.items(), key=lambda kv: -len(kv[1])):
        print(f"   {len(s):3d}  {e}")
    k27 = {r["sample_label"] for r in rows if r["assay"] == "H3K27Ac"}
    print(f"\nsamples with H3K27ac: {len(k27)}")


if __name__ == "__main__":
    main()
