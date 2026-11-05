#!/usr/bin/env python3
"""Resolve every organism string the harvest recorded against NCBI Taxonomy.

GEO gives each sample a free-text organism; the atlas carried it through without ever
asking what it was. This stage resolves each distinct string once to a taxid and its
lineage, and writes T34_organism_taxonomy.tsv -- the lookup every later stage uses to
decide whether a record's material could host the sarcoma it was filed under.

`host_class` is the question the atlas actually needs answered:

  human              Homo sapiens
  vertebrate         any other vertebrate -- mouse, rat, dog, zebrafish and the rest can
                     carry a sarcoma, so these are real material, not errors
  invertebrate       fly, worm, sea urchin, butterfly -- no mesenchymal tumours
  fungus / plant / bacterium / archaeon / other
                     cannot host a sarcoma under any reading
  passenger          virus, synthetic construct, metagenome, unparseable -- says nothing
                     about the host, so it never decides a record on its own

Resolution is cached: the table is only re-queried for strings it does not already hold.
"""
import csv, gzip, json, os, re, sys, time, urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import DATA, resolve, dpath

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "sarcoma-epigenomic-atlas/1.0 (research; taxonomy resolution)"}
OUT = "T34_organism_taxonomy.tsv"
COLS = ["organism", "taxid", "scientific_name", "rank", "host_class", "lineage"]

# Strings that describe no organism: they are passengers in a sample, not its material.
PASSENGER = re.compile(r"virus|viroid|phage|synthetic|other\s+sequences|metagenome|"
                       r"unidentified|^\d+$", re.I)
# Lineage divisions, most specific first.
DIVISIONS = [("Vertebrata", "vertebrate"), ("Metazoa", "invertebrate"),
             ("Fungi", "fungus"), ("Viridiplantae", "plant"),
             ("Bacteria", "bacterium"), ("Archaea", "archaeon"),
             ("Viruses", "passenger")]


def eget(url, tries=4):
    for a in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(1.5 * (a + 1))
    return ""


def strings_in_use():
    """Every distinct organism string across the tables that carry one."""
    seen = set()
    for t in ("T4_samples_atomic.tsv", "T30_delta_samples.tsv"):
        try:
            p = resolve(t)
        except FileNotFoundError:
            continue
        op = gzip.open if p.endswith(".gz") else open
        with op(p, "rt", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                for part in re.split(r"[;/]", r.get("organism") or ""):
                    part = part.strip()
                    if part:
                        seen.add(part)
    return seen


def lookup(name):
    """(taxid, scientific_name, rank, host_class, lineage) for one organism string."""
    if PASSENGER.search(name):
        return ("", name, "", "passenger", "")
    # "Mus musculus xenograft" and friends: resolve the species, keep the qualifier out
    q = re.sub(r"\b(xenograft|pdx|strain\s+\S+)\b.*$", "", name, flags=re.I).strip()
    js = eget(f"{EUTILS}/esearch.fcgi?db=taxonomy&retmode=json&term="
              + urllib.parse.quote(f'"{q}"'))
    ids = []
    try:
        ids = json.loads(js)["esearchresult"].get("idlist", [])
    except Exception:
        pass
    if not ids and re.search(r"\s\d+$", q):
        # "Methylobacter tundripaludum 21": a strain number NCBI does not index as a name
        js = eget(f"{EUTILS}/esearch.fcgi?db=taxonomy&retmode=json&term="
                  + urllib.parse.quote('"' + re.sub(r"\s\d+$", "", q) + '"'))
        try:
            ids = json.loads(js)["esearchresult"].get("idlist", [])
        except Exception:
            pass
    if not ids:
        return ("", name, "", "unresolved", "")
    taxid = ids[0]
    xml = eget(f"{EUTILS}/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml")
    sci = re.search(r"<ScientificName>(.*?)</ScientificName>", xml)
    rank = re.search(r"<Rank>(.*?)</Rank>", xml)
    lin = re.search(r"<Lineage>(.*?)</Lineage>", xml)
    lineage = lin.group(1) if lin else ""
    host = "other"
    if (sci.group(1) if sci else "") == "Homo sapiens":
        host = "human"
    else:
        for marker, cls in DIVISIONS:
            if re.search(rf"\b{marker}\b", lineage):
                host = cls
                break
    return (taxid, sci.group(1) if sci else name, rank.group(1) if rank else "",
            host, lineage)


def main():
    path = dpath(OUT)
    have = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
            have[r["organism"]] = r
    todo = sorted(strings_in_use() - set(have))
    print(f"{len(have)} organisms cached; {len(todo)} to resolve")
    for i, name in enumerate(todo, 1):
        taxid, sci, rank, host, lineage = lookup(name)
        have[name] = dict(zip(COLS, [name, taxid, sci, rank, host, lineage]))
        print(f"  [{i}/{len(todo)}] {name!r} -> {host} ({sci}, taxid {taxid or '-'})",
              flush=True)
        time.sleep(0.35)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for k in sorted(have):
            w.writerow(have[k])
    from collections import Counter
    print("\nhost_class:", dict(Counter(r["host_class"] for r in have.values())))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
