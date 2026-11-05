#!/usr/bin/env python3
"""Regression tests for the species gate (correction #21).

The case that prompted it: a collaborator followed the GEO link for ESS-1, an endometrial
stromal sarcoma line, and landed on GSE147927 -- "A high-resolution protein architecture
of the budding yeast genome". "Ess1" is a yeast prolyl isomerase. The organism had been
in the table all along.

Run: python3 test_species_gate.py   (needs data/T34_organism_taxonomy.tsv from stage 40)
"""
import importlib.util, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("c04", os.path.join(HERE, "04_classify.py"))
C4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(C4)


def rec(**kw):
    r = {"title": "", "source_name": "", "characteristics": "", "organism": "",
         "gse_title": "", "model_matched": "", "disease_source": "", "sample_type": ""}
    r.update(kw)
    return r


def rules(r):
    return {x[0] for x in C4.species_material(r)}


CASES = [
    # (name, record, expected rules)
    ("GSE147927: yeast Ess1 read as the line ESS-1",
     rec(title="14462_Ess1_rep1_ChIP-exo", source_name="14462_Ess1_rep1_ChIP-exo",
         characteristics="antibody: TAP-tag: Sigma i5006 ;; growthmedia: YPD",
         organism="Saccharomyces cerevisiae", model_matched="ESS-1",
         gse_title="A high-resolution protein architecture of the budding yeast genome",
         disease_source="model", sample_type="cell_line"), {"21a"}),

    ("Arabidopsis pulled in by 'chromosome arms' matching ARMS",
     rec(title="ROS1 WT rep1", source_name="seedling", organism="Arabidopsis thaliana",
         gse_title="Contrasting epigenetic dynamics in chromosome arms and pericentromeres",
         disease_source="series", sample_type="unspecified"), {"21a"}),

    ("archaeal RNA chemistry filed under UPS/MFH",
     rec(title="TK_WT_1", source_name="cells", organism="Thermococcus kodakarensis",
         disease_source="series", sample_type="unspecified"), {"21a"}),

    ("LPS the endotoxin read as the liposarcoma line LPS1",
     rec(title="HDAC3_KO_LPS_1", source_name="bone marrow-derived macrophages",
         characteristics="strain: C57BL/6 ;; age: 8 weeks ;; agent: LPS",
         organism="Mus musculus", model_matched="LPS1", disease_source="model",
         sample_type="cell_line"), {"21c"}),

    ("Mycoplasma gallisepticum 'Mg63' read as the line MG-63",
     rec(title="Mg63", source_name="tracheal mucosa", organism="Mycoplasmoides gallisepticum",
         model_matched="MG-63", disease_source="model", sample_type="cell_line"), {"21a"}),

    # --- material that must survive ---
    ("a human sarcoma cell line",
     rec(title="RH30 H3K27ac", source_name="RH30 cells", organism="Homo sapiens",
         model_matched="RH30", disease_source="model", sample_type="cell_line"), set()),

    ("a mouse carrying a human line as a xenograft",
     rec(title="RH30 xenograft tumor", source_name="RH30 xenograft",
         characteristics="strain: NSG ;; tissue: tumor from mice injected with RH30 cells",
         organism="Mus musculus", model_matched="RH30", disease_source="model",
         sample_type="xenograft(CDX)"), set()),

    ("MOS-J, a murine osteosarcoma line in its own species",
     rec(title="01_N1", source_name="MOSJ-Dkk1 cell",
         characteristics="tissue: MOSJ-Dkk1 cultures", organism="Mus musculus",
         model_matched="MOS-J", disease_source="model", sample_type="cell_line"), set()),

    ("a human tumour whose record names a virus as well",
     rec(title="KS lesion", source_name="Kaposi sarcoma biopsy",
         organism="Homo sapiens;Human gammaherpesvirus 8", disease_source="sample",
         sample_type="primary_tumor"), set()),

    ("an unresolved organism is never held against a record",
     rec(title="patient 12 tumour", source_name="osteosarcoma biopsy", organism="",
         disease_source="sample", sample_type="primary_tumor"), set()),
]

HOSTS = [("Homo sapiens", "human"), ("Mus musculus", "vertebrate"),
         ("Homo sapiens;Mus musculus", "human"), ("Canis lupus familiaris", "vertebrate"),
         ("Saccharomyces cerevisiae", "fungus"), ("Arabidopsis thaliana", "plant"),
         ("Thermococcus kodakarensis", "archaeon"), ("Drosophila melanogaster", "invertebrate"),
         ("Homo sapiens;Human gammaherpesvirus 8", "human"),
         ("synthetic construct", "unknown"), ("", "unknown")]

if __name__ == "__main__":
    bad = 0
    if not C4._taxonomy():
        sys.exit("T34_organism_taxonomy.tsv not found -- run 40_resolve_species.py first")
    for organism, want in HOSTS:
        got = C4.host_class(organism)
        ok = got == want
        bad += not ok
        print(f"  {'ok ' if ok else 'FAIL'}  host_class({organism!r}) = {got}"
              + ("" if ok else f"  expected {want}"))
    print()
    for name, r, want in CASES:
        got = rules(r)
        ok = got == want
        bad += not ok
        print(f"  {'ok ' if ok else 'FAIL'}  {name}"
              + ("" if ok else f"\n          got {got or '{}'}, expected {want or '{}'}"))
    print(f"\n{len(HOSTS) + len(CASES) - bad}/{len(HOSTS) + len(CASES)} passed")
    sys.exit(1 if bad else 0)
