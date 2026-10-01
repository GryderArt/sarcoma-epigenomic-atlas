"""British to American spelling normalisation for the atlas's own text.

Two things in this repository are written in English and they must be treated
differently.

  OUR text  -- entity labels, column values we derive, figure captions, prose. These
               should read in US English, because the atlas is a US-authored resource
               reporting US incidence.

  THEIR text -- the verbatim contents of GEO, EGA, CCDI and DKFZ records, and the
               regular expressions that match them. European depositors write
               "tumour" and "haemangioendothelioma"; rewriting either the records or
               the patterns that find them would corrupt the harvest. These are left
               exactly as deposited.

The separation is enforced two ways: TSV rewriting is column-scoped with an explicit
deny-list of verbatim fields, and Python rewriting masks every raw-string literal
before substituting and restores it afterwards.
"""
import re

# Substring rules: entity labels, and stems that are prefixes of longer words
# ("haemangio" inside "haemangioendothelioma"). Applied longest-first.
SUBSTR = [
    ("Solitary fibrous tumour", "Solitary fibrous tumor"),
    ("Ossifying fibromyxoid tumour", "Ossifying fibromyxoid tumor"),
    ("Inflammatory myofibroblastic tumour", "Inflammatory myofibroblastic tumor"),
    ("Epithelioid haemangioendothelioma", "Epithelioid hemangioendothelioma"),
    ("Rhabdoid tumour", "Rhabdoid tumor"),
    ("Giant cell tumour", "Giant cell tumor"),
    ("haemangio", "hemangio"), ("Haemangio", "Hemangio"),
    ("haemangiosarc", "hemangiosarc"),
    ("paediatric", "pediatric"), ("Paediatric", "Pediatric"),
    ("tumour", "tumor"), ("Tumour", "Tumor"),
    ("per cent", "percent"), ("Per cent", "Percent"),
]
SUBSTR.sort(key=lambda kv: -len(kv[0]))

# Whole-word rules. Written out as complete words so that "analysis" and "analyses",
# which are identical in both varieties, are never touched by an "analyse" stem.
WORDS = {
    "analyse": "analyze", "analysed": "analyzed", "analysing": "analyzing",
    "reanalyse": "reanalyze", "reanalysed": "reanalyzed", "reanalysing": "reanalyzing",
    "normalise": "normalize", "normalised": "normalized", "normalising": "normalizing",
    "normalisation": "normalization",
    "organise": "organize", "organised": "organized", "organising": "organizing",
    "organisation": "organization",
    "recognise": "recognize", "recognised": "recognized",
    "characterise": "characterize", "characterised": "characterized",
    "characterising": "characterizing",
    "summarise": "summarize", "summarised": "summarized",
    "prioritise": "prioritize", "prioritised": "prioritized",
    "prioritising": "prioritizing",
    "utilise": "utilize", "utilised": "utilized",
    "generalise": "generalize", "generalised": "generalized",
    "generalisation": "generalization",
    "standardise": "standardize", "standardised": "standardized",
    "minimise": "minimize", "maximise": "maximize",
    "emphasise": "emphasize", "emphasised": "emphasized",
    "catalogue": "catalog", "catalogues": "catalogs",
    "catalogued": "cataloged", "cataloguing": "cataloging",
    "unlabelled": "unlabeled", "relabelled": "relabeled",
    "labelled": "labeled", "labelling": "labeling", "labels": "labels",
    "modelling": "modeling", "modelled": "modeled",
    "signalling": "signaling", "totalled": "totaled", "cancelled": "canceled",
    "colour": "color", "colours": "colors", "coloured": "colored",
    "behaviour": "behavior", "favour": "favor", "favoured": "favored",
    "centre": "center", "centres": "centers",
    "defence": "defense", "licence": "license",
    "grey": "gray", "ageing": "aging", "judgement": "judgment", "fulfil": "fulfill",
    "whilst": "while", "amongst": "among", "towards": "toward",
    "entity-labelled": "entity-labeled",
}
# add capitalised forms automatically
WORDS.update({k.capitalize(): v.capitalize() for k, v in list(WORDS.items())})
_WORD_RE = re.compile(r"\b(" + "|".join(sorted(map(re.escape, WORDS), key=len,
                                                reverse=True)) + r")\b")

# Columns holding verbatim third-party text. Never rewritten.
VERBATIM_COLS = {
    "title", "source_name", "characteristics", "gse_title", "antibody_raw",
    "contact_institute", "description", "institutional_diagnosis", "diseases",
    "assays_available", "viz_title", "aliases", "key_alterations", "sample_label",
    "viz_url", "file", "portal_diagnosis", "portal_subtype", "admin_label",
    "duplicate_of_gsm", "platform", "library_strategy", "organism",
}
# Files that are wholly third-party inventories.
VERBATIM_FILES = {"T16_stjude_cstn_inventory.tsv", "T17_stjude_opdx_models.tsv",
                  "T28_stjude_viz_tracks.tsv"}

# The (?<![A-Za-z0-9_]) guard matters: without it the trailing "r" of a word such as
# "tumour" followed by a closing quote is read as the start of a raw string, and the
# label before a pattern gets masked along with it.
# The (?<![A-Za-z0-9_]) guard matters: without it the trailing "r" of a word such as
# "tumour" followed by a closing quote reads as the start of a raw string, and the
# label sitting before a pattern gets masked along with it.
_DQ = '(?<![A-Za-z0-9_])r"(?:[^"\\\\]|\\\\.)*"'
_SQ = "(?<![A-Za-z0-9_])r'(?:[^'\\\\]|\\\\.)*'"
_RAWSTR = re.compile(_DQ + "|" + _SQ)


def convert(text):
    """Normalise one string: whole-word rules first, then substring rules."""
    text = _WORD_RE.sub(lambda m: WORDS[m.group(1)], text)
    for a, b in SUBSTR:
        if a in text:
            text = text.replace(a, b)
    return text


def convert_python(src):
    """Convert a Python source file, leaving every raw-string literal untouched.

    Raw strings in this pipeline are regular expressions matched against deposited
    text, which is frequently written in British English. Rewriting them would stop
    them matching.
    """
    held = []

    def stash(m):
        held.append(m.group(0))
        return f"\x00{len(held) - 1}\x00"

    masked = _RAWSTR.sub(stash, src)
    out = convert(masked)
    return re.sub(r"\x00(\d+)\x00", lambda m: held[int(m.group(1))], out)
