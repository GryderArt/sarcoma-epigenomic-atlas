#!/usr/bin/env python3
"""Classifier v2 — paediatric + adult sarcoma ontology.

Carries every correction established during the paediatric build:
  - subtype comes from SAMPLE-level text, never the series title
  - fusion-negative RMS = absence of PAX3/PAX7 fusion and mutant MYOD1
  - PAX3/PAX7 fused to ANY partner is fusion-positive
  - antibody/vendor catalogue text stripped before model matching
  - explicit `enrichment target:` beats blob inference for ChIP inputs
  - H3.3 residue distinguishes bone tumour (G34W/L, K36M) from glioma (K27M, G34R/V)
"""
import csv, gzip, re, json, sys, collections, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
D = WORK; T = DATA; O = DATA
csv.field_size_limit(10**7)

# ---------------------------------------------------------------- disease ontology
DIS = [
 # ---- paediatric core
 ("RMS-MYOD1", r"myod1[^;|]{0,20}(l122r|p\.leu122arg|\bmutant\b|\bmut\b)|(mutant|mut)\s*myod1|"
               r"sclerosing\s+rhabdomyo|spindle\s+cell\s+rhabdomyo"),
 ("FP-RMS",  r"alveolar\s+rhabdomyosarcoma|\barms\b(?!.*adult)|fusion[\-\s]positive\s+rms|"
             r"pax\s*[37]\s*(?:::|--|[-/:])\s*(?!(?:negative|positive|null|low|high|wt|ko|kd)\b)[A-Za-z][A-Za-z0-9]{2,8}|"
             r"pax[37][\-\s/:]*fkhr|\bfp[\-\s]?rms\b|\bp3f\b|\bp7f\b"),
 ("FN-RMS",  r"embryonal\s+rhabdomyosarcoma|\berms\b|fusion[\-\s]negative\s+rms|\bfn[\-\s]?rms\b"),
 ("RMS-NOS", r"rhabdomyosarcoma|\brms\b"),
 ("Ewing",   r"ewing|\bews[\-\s/:]*fli|\bewsr1[\-\s/:]*fli|askin|peripheral\s+primitive\s+neuroecto|\bpnet\b"),
 ("CIC-DUX4",r"cic[\-\s/:]*dux4|cic[\-\s]rearranged|cic[\-\s/:]*nutm1"),
 ("BCOR-sarcoma", r"bcor[\-\s/:]*ccnb3|bcor[\-\s]rearranged|bcor[\-\s]itd"),
 ("EWSR1-NFATC2/PATZ1", r"ewsr1[\-\s/:]*nfatc2|ewsr1[\-\s/:]*patz1|fus[\-\s/:]*nfatc2"),
 ("Osteosarcoma", r"osteosarcoma|osteogenic\s+sarcoma|\bu2os\b|\bu\-2\s?os\b"),
 ("Rhabdoid tumor/ATRT", r"rhabdoid\s+tumou?r|atypical\s+teratoid|\batrt\b|\bmrt\b|"
                         r"smarcb1[\-\s]deficient|ini1[\-\s]deficient"),
 ("DSRCT",   r"desmoplastic\s+small\s+round\s+cell|ewsr1[\-\s/:]*wt1"),
 ("Infantile fibrosarcoma", r"infantile\s+fibrosarcoma|congenital\s+fibrosarcoma|etv6[\-\s/:]*ntrk3"),
 ("GCTB/Chondroblastoma", r"giant\s+cell\s+tumou?r\s+of\s+bone|\bgctb\b|chondroblastoma|g34w|k36m"),
 # ---- shared paediatric/adult
 ("Synovial sarcoma", r"synovial\s+sarcoma|ss18[\-\s/:]*ssx|syt[\-\s/:]*ssx"),
 ("MPNST",   r"malignant\s+peripheral\s+nerve\s+sheath|\bmpnst\b|neurofibrosarcoma"),
 ("ASPS",    r"alveolar\s+soft\s+part|aspscr1[\-\s/:]*tfe3"),
 ("Clear cell sarcoma", r"clear\s+cell\s+sarcoma|ewsr1[\-\s/:]*atf1"),
 ("Epithelioid sarcoma", r"epithelioid\s+sarcoma"),
 ("IMT",     r"inflammatory\s+myofibroblastic"),
 ("Chordoma", r"chordoma|brachyury"),
 ("EMC",     r"extraskeletal\s+myxoid\s+chondrosarcoma|ewsr1[\-\s/:]*nr4a3"),
 ("DFSP",    r"dermatofibrosarcoma|col1a1[\-\s/:]*pdgfb"),
 ("Desmoid", r"desmoid|aggressive\s+fibromatosis"),
 # ---- adult liposarcoma, resolved to subtype
 ("Liposarcoma-myxoid", r"myxoid\s+liposarcoma|round\s+cell\s+liposarcoma|fus[\-\s/:]*ddit3|"
                        r"ewsr1[\-\s/:]*ddit3|tls[\-\s/:]*chop|\bmlps\b|\bmrcls\b"),
 ("Liposarcoma-dediff", r"dedifferentiated\s+liposarcoma|\bddlps\b|\bddls\b"),
 ("Liposarcoma-WD", r"well[\-\s]differentiated\s+liposarcoma|atypical\s+lipomatous|\bwdlps\b|\balt\b(?=.*lipo)"),
 ("Liposarcoma-pleomorphic", r"pleomorphic\s+liposarcoma"),
 ("Liposarcoma-NOS", r"liposarcoma"),
 # ---- adult vascular
 ("Angiosarcoma", r"angiosarcoma|h[ae]mangiosarcoma"),
 ("EHE", r"epithelioid\s+h[ae]mangioendothelioma|wwtr1[\-\s/:]*camta1|taz[\-\s/:]*camta1|"
         r"yap1[\-\s/:]*tfe3|h[ae]mangioendothelioma"),
 ("Kaposi sarcoma", r"kaposi"),
 # ---- other adult soft tissue
 ("UPS/MFH", r"undifferentiated\s+pleomorphic\s+sarcoma|malignant\s+fibrous\s+histiocytoma|\bups\b|\bmfh\b"),
 ("Myxofibrosarcoma", r"myxofibrosarcoma"),
 ("Solitary fibrous tumour", r"solitary\s+fibrous|nab2[\-\s/:]*stat6|h[ae]mangiopericytoma"),
 ("Leiomyosarcoma", r"leiomyosarcoma|\blms\b"),
 ("GIST",    r"gastrointestinal\s+stromal|\bgist\b"),
 ("Endometrial stromal sarcoma", r"endometrial\s+stromal\s+sarcoma|jazf1[\-\s/:]*suz12|ywhae[\-\s/:]*nutm2"),
 ("PEComa",  r"pecoma|perivascular\s+epithelioid"),
 ("LGFMS/SEF", r"low[\-\s]grade\s+fibromyxoid|sclerosing\s+epithelioid\s+fibrosarcoma|"
               r"fus[\-\s/:]*creb3l|ewsr1[\-\s/:]*creb3l"),
 ("Myoepithelial carcinoma", r"myoepithelial\s+(carcinoma|tumou?r).{0,30}soft\s+tissue|ewsr1[\-\s/:]*pou5f1"),
 ("Intimal sarcoma", r"intimal\s+sarcoma"),
 ("Chondrosarcoma", r"chondrosarcoma"),
 ("Fibrosarcoma NOS", r"fibrosarcoma"),
 ("Sarcoma NOS", r"sarcoma|mesenchymal\s+(tumou?r|neoplas)"),
]
DIS_PATS = [(k, re.compile(v, re.I)) for k, v in DIS]

PAED_CORE = {"FP-RMS","FN-RMS","RMS-MYOD1","RMS-NOS","Ewing","CIC-DUX4","BCOR-sarcoma",
             "EWSR1-NFATC2/PATZ1","Osteosarcoma","Rhabdoid tumor/ATRT","DSRCT",
             "Infantile fibrosarcoma","GCTB/Chondroblastoma"}
SHARED = {"Synovial sarcoma","MPNST","ASPS","Clear cell sarcoma","Epithelioid sarcoma","IMT",
          "Chordoma","EMC","DFSP","Desmoid","Liposarcoma-myxoid","Chondrosarcoma","GIST",
          "EHE","Angiosarcoma"}
ADULT_CORE = {"Liposarcoma-dediff","Liposarcoma-WD","Liposarcoma-pleomorphic","Liposarcoma-NOS",
              "UPS/MFH","Myxofibrosarcoma","Solitary fibrous tumour","Leiomyosarcoma",
              "Endometrial stromal sarcoma","PEComa","LGFMS/SEF","Myoepithelial carcinoma",
              "Intimal sarcoma","Fibrosarcoma NOS","Kaposi sarcoma"}
IN_SCOPE_DIS = PAED_CORE | SHARED | ADULT_CORE | {"Sarcoma NOS"}
AGE_CLASS = {**{d: "paediatric" for d in PAED_CORE},
             **{d: "both" for d in SHARED},
             **{d: "adult" for d in ADULT_CORE}, "Sarcoma NOS": "both"}

# ---------------------------------------------------------------- assay + target rules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ASSAY_RULES = [
 ("HiChIP", r"hichip|hi-chip|plac[\-\s]?seq"), ("ChIA-PET", r"chia[\-\s]?pet"),
 ("Micro-C", r"micro[\-\s]?c\b"), ("Capture-HiC", r"capture\s*hi[\-\s]?c|promoter\s*capture"),
 ("Hi-C", r"\bhi[\-\s]?c\b|in\s?situ\s?hic"), ("4C-seq", r"\b4c[\-\s]?seq\b"),
 ("Repli-seq", r"repli[\-\s]?seq"),
 ("CUT&Tag", r"cut&tag|cut\s*and\s*tag|cut[\-\s]?tag|cuttag"),
 ("CUT&RUN", r"cut&run|cut\s*and\s*run|cut[\-\s]?run|cutrun"),
 ("ChIP-exo", r"chip[\-\s]?exo"),
 ("scATAC-seq", r"\bsc[\-\s]?atac|single[\-\s]cell\s+atac|sn[\-\s]?atac|single[\-\s]nucle\w+\s+atac|multiome"),
 ("ATAC-seq", r"\batac\b|assay\s+for\s+transposase"),
 ("DNase-seq", r"dnase[\-\s]?seq|dnase\s?i?\s?hypersensit|dhs\b"),
 ("FAIRE-seq", r"faire"), ("MNase-seq", r"mnase[\-\s]?seq|nucleosome\s+occupancy"),
 ("WGBS", r"wgbs|whole[\-\s]genome\s+bisulfite|em[\-\s]?seq|methyl[\-\s]?seq"),
 ("RRBS", r"\brrbs\b|reduced\s+representation\s+bisulfite"),
 ("Methyl-array", r"450k|850k|epic\s*(array|beadchip|v2)|methylation\s*(bead)?(chip|array)|"
                  r"humanmethylation|infinium.*methyl"),
 ("MeDIP/hMeDIP", r"me[\-]?dip|hmedip|hmc[\-\s]?seq|5hmc|tab[\-\s]?seq"),
 ("Bisulfite-PCR", r"bisulfite\s+(pcr|amplicon|pyroseq)|pyrosequencing"),
 ("scRNA-seq", r"\bsc[\-\s]?rna|single[\-\s]cell\s+rna|sn[\-\s]?rna[\-\s]?seq|single[\-\s]nucle\w+\s+rna|"
               r"10x\s+genomics|smart[\-\s]?seq|drop[\-\s]?seq|cite[\-\s]?seq"),
 ("Spatial", r"visium|spatial\s+transcript|geomx|merfish|xenium|slide[\-\s]?seq|spatial\s+atac"),
 ("GRO/PRO-seq", r"gro[\-\s]?seq|pro[\-\s]?seq|nascent|net[\-\s]?seq|start[\-\s]?seq"),
 ("Ribo-seq", r"ribo[\-\s]?seq|ribosome\s+profiling"),
 ("CLIP-seq", r"clip[\-\s]?seq|par[\-\s]?clip|iclip"),
 ("CRISPR-screen", r"crispr\s+screen|sgrna\s+library|genome[\-\s]wide\s+screen|shrna\s+screen"),
 ("ChIP-seq", r"chip[\-\s]?seq|chip[\-\s]?sequencing|chromatin\s+immunoprecip"),
 ("ChIP-chip", r"chip[\-\s]?chip|chip[\-\s]?on[\-\s]?chip"),
 ("WGS", r"whole[\-\s]genome\s+sequencing|\bwgs\b"),
 ("WES", r"whole[\-\s]exome|\bwes\b|exome[\-\s]?seq"),
 ("SNP-array", r"snp\s*(6|array|genotyping)|cytoscan|omniexpress"),
 ("aCGH", r"\bacgh\b|array\s*cgh|comparative\s+genomic"),
 ("miRNA", r"mirna|micro[\-\s]?rna|small\s*rna"),
 ("Proteomics", r"proteom|mass\s+spec|lc[\-\s]?ms|tmt\b|silac"),
 ("RNA-seq", r"rna[\-\s]?seq|transcriptome\s+sequencing|total\s+rna|mrna[\-\s]?seq|polya"),
 ("Expr-array", r"affymetrix|agilent|illumina\s+bead|expression\s+array|microarray|u133|gene\s*chip"),
]
ASSAY_PATS = [(k, re.compile(v, re.I)) for k, v in ASSAY_RULES]
LIBSTRAT = {"RNA-Seq":"RNA-seq","ChIP-Seq":"ChIP-seq","ATAC-seq":"ATAC-seq","Bisulfite-Seq":"WGBS",
 "DNase-Hypersensitivity":"DNase-seq","Hi-C":"Hi-C","WGS":"WGS","WXS":"WES","miRNA-Seq":"miRNA",
 "ncRNA-Seq":"ncRNA-seq","MeDIP-Seq":"MeDIP/hMeDIP","MRE-Seq":"MeDIP/hMeDIP","MBD-Seq":"MeDIP/hMeDIP",
 "FAIRE-seq":"FAIRE-seq","CUT&RUN":"CUT&RUN","CUT&Tag":"CUT&Tag","RIP-Seq":"RIP-seq",
 "CLIP-Seq":"CLIP-seq","Targeted-Capture":"Targeted-seq","AMPLICON":"Targeted-seq","MNase-Seq":"MNase-seq"}
MARK_RULES = [("H2BNTac",r"h2bnt[\-\s]?ac|h2bk5ac|h2bk12ac|h2bk15ac|h2bk20ac"),
 ("H3K27ac",r"h3k27[\-\s]?ac"),("H3K27me3",r"h3k27[\-\s]?me3"),("H3K27me2",r"h3k27[\-\s]?me2"),
 ("H3K4me1",r"h3k4[\-\s]?me1"),("H3K4me2",r"h3k4[\-\s]?me2"),("H3K4me3",r"h3k4[\-\s]?me3"),
 ("H3K9me3",r"h3k9[\-\s]?me3"),("H3K9me2",r"h3k9[\-\s]?me2"),("H3K9ac",r"h3k9[\-\s]?ac"),
 ("H3K36me3",r"h3k36[\-\s]?me3"),("H3K36me2",r"h3k36[\-\s]?me2"),("H3K79me2",r"h3k79[\-\s]?me"),
 ("H3K18la",r"h3k18[\-\s]?la|lactyl"),("H4K20me",r"h4k20[\-\s]?me"),
 ("H4ac",r"h4k(5|8|12|16)ac|h4[\-\s]?ac\b|pan[\-\s]?h4"),("H2AK119ub",r"h2ak119|h2aub|ubh2a"),
 ("H3K122ac",r"h3k122"),("H2AZ",r"h2a\.?z"),("5hmC",r"5hmc|hydroxymethyl"),
 ("5mC",r"\b5mc\b|dna\s+methylation|cpg\s+methylation")]
MARK_PATS = [(k, re.compile(v, re.I)) for k, v in MARK_RULES]
TF_RULES = [("PAX3-FOXO1",r"pax3[\-\s/:]*foxo1|pax3[\-\s/:]*fkhr|\bp3f\b"),
 ("PAX7-FOXO1",r"pax7[\-\s/:]*foxo1|\bp7f\b"),("EWS-FLI1",r"ews(r1)?[\-\s/:]*fli1?"),
 ("EWS-ERG",r"ews(r1)?[\-\s/:]*erg"),("SS18-SSX",r"ss18[\-\s/:]*ssx|syt[\-\s/:]*ssx"),
 ("EWSR1-WT1",r"ews(r1)?[\-\s/:]*wt1"),("CIC-DUX4",r"cic[\-\s/:]*dux4"),
 ("BCOR-CCNB3",r"bcor[\-\s/:]*ccnb3"),("TAZ-CAMTA1",r"(wwtr1|taz)[\-\s/:]*camta1"),
 ("FUS-DDIT3",r"fus[\-\s/:]*ddit3|tls[\-\s/:]*chop"),("NAB2-STAT6",r"nab2[\-\s/:]*stat6"),
 ("MDM2",r"\bmdm2\b"),("CDK4",r"\bcdk4\b"),("MYOD1",r"\bmyod1?\b"),("MYOG",r"\bmyog(enin)?\b"),
 ("MYCN",r"\bmycn\b|n[\-\s]?myc"),("MYC",r"\bc?[\-\s]?myc\b"),("RUNX2",r"\brunx2\b"),
 ("FOSL1",r"\bfosl1\b|\bfra[\-\s]?1\b"),("SOX2",r"\bsox2\b"),("BRD4",r"\bbrd4\b|\bbrd[23]\b"),
 ("EP300/CBP",r"\bep300\b|\bp300\b|\bcrebbp\b|\bcbp\b"),
 ("POLR2A",r"pol\s?(ii|2)|polr2a|rnapii|rpb1"),("CTCF",r"\bctcf\b"),
 ("Cohesin",r"\brad21\b|\bsmc1|\bsmc3|\bstag[12]\b|\bnipbl\b"),
 ("SMARCA4",r"smarca4|\bbrg1\b"),("SMARCA2",r"smarca2|\bbrm\b"),
 ("SMARCB1",r"smarcb1|\bini1\b|\bbaf47\b|\bsnf5\b"),("ARID1A",r"arid1a|\bbaf250"),
 ("EZH2/PRC2",r"\bezh2\b|\bsuz12\b|\beed\b|\bprc2\b"),
 ("RING1B/PRC1",r"ring1b|\brnf2\b|\bcbx\d|\bprc1\b"),("HDAC",r"\bhdac\d?\b"),
 ("LSD1/KDM1A",r"lsd1|kdm1a"),("NKX2-2",r"nkx2[\-\s]?2"),("TEAD",r"\btead\d?\b|\byap1?\b"),
 ("TP53",r"\btp53\b|\bp53\b"),("input",r"\binput\b|\bwce\b|whole\s+cell\s+extract"),
 ("IgG",r"\bigg\b|\bmock\b|no\s+antibody")]
TF_PATS = [(k, re.compile(v, re.I)) for k, v in TF_RULES]
EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq","DNase-seq",
       "FAIRE-seq","MNase-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR",
       "Hi-C","HiChIP","ChIA-PET","Micro-C","Capture-HiC","4C-seq","Repli-seq"}

# ---------------------------------------------------------------- model index
GENERIC_ALIAS = {"WT","OS","ES","RT","SS","MRT","CS","NY","HAL","VK","SIM","MIC","POE","CTR","MOS",
 "KAS","GIST","EWS","RMS","SARC","PDX","N/A","NA","-","UNKNOWN","NONE","T","P","M","OSA","COL",
 "HOS 58","2T","SIMON","HAMON","DUNN","1273","ES-2","ES2","BIRCH","CG-1","JJ","SARG","GCT",
 "SKN","MES","LP6","SVR","AS-M","ASM","GOT3","T449","T778","LIS-3"}
_BAD_KEY = re.compile(r"^(os|es|cs|mos|ny|rt|ms|sj)\d{1,3}$|^(nb|ab|sc|cst|pa|ma)\d{3,}$")
_AB_STRIP = re.compile(r"(chip\s*)?antibod(y|ies)\s*[:=][^;|]*", re.I)

def load_models():
    idx = {}
    for fn in ("RMS_models.tsv","OS_models.tsv","EWS_models.tsv","STS_models.tsv","ADULT_models.tsv"):
        p = os.path.join(T, fn)
        if not os.path.exists(p): continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            canon = (r.get("model_name") or "").strip()
            if not canon or "NO MODEL" in canon.upper() or canon.upper().startswith("NONE"): continue
            dis = (r.get("disease") or "").strip()
            sub = (r.get("subtype_or_fusion") or r.get("subtype_code") or r.get("subtype_if_known")
                   or r.get("subtype/fusion") or "").strip()
            mt = (r.get("model_type") or r.get("model_type(cell_line/PDX/organoid)") or "").strip()
            rrid = (r.get("RRID_or_Cellosaurus") or "").strip()
            prob = (r.get("problematic_flag") or "").strip()
            for n in [canon] + [a.strip() for a in re.split(r"[;,]", r.get("aliases") or "") if a.strip()]:
                n = n.strip().strip("()")
                if not n or len(n) < 3 or n.upper() in GENERIC_ALIAS: continue
                if n.lower() in ("unknown","aliases","not registered","none","n/a"): continue
                idx.setdefault(n, (canon, dis, sub, mt, rrid, prob))
    p2 = os.path.join(O, "T2_models.tsv")
    if os.path.exists(p2):
        for r in csv.DictReader(open(p2), delimiter="\t"):
            canon = (r.get("model_name") or "").strip()
            if not canon: continue
            for n in [canon] + [a.strip() for a in re.split(r"[;,]", r.get("aliases") or "") if a.strip()]:
                n = n.strip().strip("()")
                if not n or len(n) < 3 or n.upper() in GENERIC_ALIAS: continue
                idx.setdefault(n, (canon, r.get("disease",""), r.get("subtype_or_fusion",""),
                                   r.get("model_type",""), r.get("RRID_or_Cellosaurus",""),
                                   r.get("problematic_flag","")))
    return idx
MODELS = load_models()
def _norm(x): return re.sub(r"[^a-z0-9]", "", x.lower())
MODEL_KEYS = {}
_MAXKEY = 3
for _n, _v in MODELS.items():
    _k = _norm(_n)
    if len(_k) < 3 or _BAD_KEY.match(_k): continue
    if _k.isdigit() and len(_k) < 4: continue
    _MAXKEY = max(_MAXKEY, len(re.findall(r"[A-Za-z0-9]+", _n)))
    MODEL_KEYS.setdefault(_k, _v)
_MAXKEY = min(_MAXKEY, 6)
_TOK = re.compile(r"[a-z0-9]+")
_RD_CUE = re.compile(r"rhabdomyosarc|\brms\b|myoblast|myogenic|\bmyod", re.I)
_RD = next((v for n, v in MODELS.items() if v[0] == "RD"), None)

def match_model(text):
    text = _AB_STRIP.sub(" ", text)
    toks = _TOK.findall(text.lower())
    L = len(toks)
    if not L: return None
    get = MODEL_KEYS.get; best = None; bestlen = 0
    for i in range(L):
        acc = toks[i]
        if len(acc) <= 48:
            h = get(acc)
            if h is not None and len(acc) > bestlen: best, bestlen = h, len(acc)
        for j in range(i+1, min(i+_MAXKEY, L)):
            acc += toks[j]
            if len(acc) > 48: break
            h = get(acc)
            if h is not None and len(acc) > bestlen: best, bestlen = h, len(acc)
    if best is None and _RD and "rd" in toks and _RD_CUE.search(text): best = _RD
    return best

AB_RE = re.compile(r"(?:chip[\s_]*)?anti(?:body|gen)?\s*[:=]\s*([^;|]{2,60})", re.I)
ENRICH = re.compile(r"enrichment\s*target\s*[:=]\s*([^;|]{1,60})", re.I)
INPUTY = re.compile(r"\binput\b|\bwce\b|whole[\-\s]cell\s+extract|\bigg\b|\bmock\b", re.I)
NONSARC = re.compile(r"\b(melanoma|carcinoma|adenocarcinoma|leukemi|lymphoma|myeloma|gliom|"
 r"glioblastoma|medulloblastoma|neuroblastoma|retinoblastoma|meningioma|wilms|nephroblastoma|"
 r"hepatoblastoma|craniopharyngioma|ependymoma|breast\s+cancer|lung\s+cancer|colorectal|"
 r"prostate\s+cancer|pancreatic\s+cancer|ovarian\s+cancer|gastric\s+cancer|hepatocellular|"
 r"cholangiocarcinoma|mesothelioma|\bdipg\b|diffuse\s+midline)\b", re.I)
TISSUEY = re.compile(r"\b(tissue|specimen|biopsy|resection|surgical|fresh[\-\s]frozen|snap[\-\s]frozen|"
 r"ffpe|formalin|patient|case|donor|autopsy|primary\s+tumou?r|tumou?r\s+sample)\b", re.I)
NORMALY = re.compile(r"\b(normal|healthy|non[\-\s]?(tumou?r|malignant)|human\s+muscle\s+tissue|"
 r"skeletal\s+muscle|myoblast|myotube|\bmsc\b|mesenchymal\s+stem|preadipocyte|adipose[\-\s]derived|"
 r"primary\s+(fibroblast|osteoblast|chondrocyte)|\bhfob\b|\bimr[\-\s]?90\b|huvec)\b", re.I)
TUMORY = re.compile(r"tumou?r|sarcom|carcinom|blastom|malignan|cancer|metasta|neoplas|biopsy|"
 r"chemo(sensitive|resistant)|patient|\bpdx\b|xenograft", re.I)
K27M_GLIOMA = re.compile(r"k27m|g34r|g34v", re.I)
GLIOMA_CTX = re.compile(r"gliom|\bdipg\b|diffuse\s+midline|astrocytom|\bgbm\b|brain\s*stem|pons|thalam", re.I)

def first(pats, blob): return [k for k, p in pats if p.search(blob)]

def sample_type(blob, organism):
    b = blob.lower()
    if re.search(r"\bpdx\b|patient[\-\s]derived\s+xenograft|o-pdx", b): return "PDX"
    if re.search(r"xenograft|\bcdx\b|orthotopic", b): return "xenograft(CDX)"
    if re.search(r"organoid|spheroid|tumoroid", b): return "organoid"
    if re.search(r"primary\s+tumou?r|patient\s+(tumou?r|sample|biopsy)|biopsy|resection|ffpe|"
                 r"tumou?r\s+(tissue|specimen)|surgical|autopsy", b): return "primary_tumor"
    if re.search(r"metasta", b): return "metastasis"
    if re.search(r"\bcell\s+line\b|\bcells?\b.*\b(cultur|passag)|\bcell\s+cultur", b): return "cell_line"
    if NORMALY.search(b): return "normal/reference"
    if organism and "Homo sapiens" not in organism and re.search(r"mouse|murine|\bmus\b|transgenic", b):
        return "mouse_model"
    return "unspecified"

def main():
    ctx = {}
    for r in csv.DictReader(open(f"{D}/gse_summary_v2.tsv"), delimiter="\t"):
        ctx[r["accession"]] = (r["title"] + " || " + r["summary"], r["pubmed"], r["pdat"],
                               r["taxon"], r["families"])
    out = f"{D}/gsm_annotated_v2.tsv.gz"
    n = 0; w = None
    with gzip.open(out, "wt", newline="") as fo:
        seen = set()
        def reader():
            try:
                for x in csv.DictReader(gzip.open(f"{D}/gsm_records_v2.tsv.gz","rt"), delimiter="\t"):
                    g = x.get("gsm")
                    if not g or g in seen: continue
                    seen.add(g); yield x
            except (EOFError, OSError, csv.Error) as e:
                sys.stderr.write(f"[tolerated: {e}]\n")
        for r in reader():
            gse = r.get("gse","")
            gctx, pmid, pdat, taxon, fams = ctx.get(gse, ("","","","",""))
            sblob = " ;; ".join([r.get("title",""), r.get("source_name",""),
                                 r.get("characteristics",""), r.get("suppl_files","")[:300]])
            mblob = " ;; ".join([r.get("title","")[:200], r.get("source_name","")[:200],
                                 r.get("characteristics","")[:500]])
            blob = sblob + " ;; " + gctx
            ls = r.get("library_strategy","").strip()
            hits = first(ASSAY_PATS, sblob) or first(ASSAY_PATS, blob)
            assay = hits[0] if hits else ""
            if ls in ("CUT&RUN","CUT&Tag"): assay = ls
            elif ls == "ATAC-seq": assay = "ATAC-seq"
            elif ls == "Hi-C" and assay not in ("HiChIP","Micro-C","ChIA-PET","Capture-HiC"): assay = "Hi-C"
            elif not assay and ls in LIBSTRAT: assay = LIBSTRAT[ls]
            if not assay:
                st = r.get("sample_type",""); lsrc = r.get("library_source","")
                assay = ("RNA-seq" if "TRANSCRIPTOMIC" in lsrc.upper() else
                         "Expr-array" if st == "RNA" else
                         "genomic-array" if st == "genomic" else "other/unknown")
            if assay in ("RNA-seq","Expr-array") and re.search(
                    r"single[\-\s]cell|single[\-\s]nucle|\bscrna|\bsnrna|10x\s+genomics|smart[\-\s]?seq2?", blob, re.I):
                assay = "scRNA-seq"
            if assay == "ATAC-seq" and re.search(r"single[\-\s]cell|single[\-\s]nucle|\bscatac|\bsnatac|multiome", blob, re.I):
                assay = "scATAC-seq"

            marks = first(MARK_PATS, sblob); tfs = first(TF_PATS, sblob)
            target = ""; ab = ""
            if assay in ("ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","HiChIP","ChIA-PET"):
                m = ENRICH.search(r.get("characteristics",""))
                if m and INPUTY.search(m.group(1)): target = "input/none"
                elif INPUTY.search(r.get("title","")): target = "input/none"
                else:
                    a = AB_RE.search(r.get("characteristics",""))
                    ab = (a.group(1).strip() if a else "")[:40]
                    combo = ab + " ;; " + sblob
                    mk = first(MARK_PATS, combo); tf = first(TF_PATS, combo)
                    target = mk[0] if mk else (tf[0] if tf else (ab or "unspecified"))
                    if target in ("input","IgG"): target = "input/none"
            elif assay in ("WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"):
                target = "5hmC" if "5hmC" in marks else "5mC"

            hit = match_model(mblob)
            model = hit[0] if hit else ""
            mdis  = hit[1] if hit else ""
            msub  = hit[2] if hit else ""
            rrid  = hit[4] if hit else ""
            mprob = hit[5] if hit else ""

            dh_s = [k for k, p in DIS_PATS if p.search(sblob)]
            dh_c = [k for k, p in DIS_PATS if p.search(gctx)]
            sp_s = [d for d in dh_s if d != "Sarcoma NOS"]
            sp_c = [d for d in dh_c if d != "Sarcoma NOS"]
            disease = (sp_s or sp_c or dh_s or dh_c or [""])[0]
            dsrc = ("sample" if sp_s else "series" if sp_c else
                    "sample" if dh_s else "series" if dh_c else "none")
            if not disease and mdis: disease, dsrc = mdis, "model"
            # H3.3 residue guard: keep glioma out of the bone-tumour bin
            if disease == "GCTB/Chondroblastoma" and (K27M_GLIOMA.search(blob) or GLIOMA_CTX.search(blob)) \
               and not re.search(r"g34w|g34l|k36m", blob, re.I):
                disease = ""

            stype = sample_type(sblob, r.get("organism",""))
            if stype == "unspecified" and model:
                mt = (MODELS.get(model, ("","","","","",""))[3] or "").lower()
                stype = "PDX" if "pdx" in mt else "cell_line" if "cell" in mt else stype
            if stype == "unspecified" and re.search(r"cell\s*line|\bcells\b", sblob, re.I): stype = "cell_line"
            if stype == "unspecified" and not model and TISSUEY.search(sblob): stype = "primary_tumor"
            if NORMALY.search(r.get("source_name","")) and not TUMORY.search(r.get("source_name","")):
                stype = "normal/reference"

            if model: relevance = "model_confirmed"
            elif sp_s: relevance = "sample_explicit"
            elif sp_c: relevance = "series_explicit"
            elif dh_s or dh_c: relevance = "sarcoma_nos"
            else: relevance = "off_target"
            in_scope = "Y" if (relevance != "off_target" and
                               (disease in IN_SCOPE_DIS or model)) else "N"
            if NONSARC.search(sblob) and not model and not sp_s and disease in ("", "Sarcoma NOS"):
                relevance, in_scope = "off_target", "N"

            rec = {**r, "assay_class": assay, "epi_target_norm": target, "antibody_raw": ab,
                   "sample_type": stype, "model_matched": model, "model_rrid": rrid,
                   "model_subtype": msub, "model_problematic": mprob,
                   "disease": disease, "disease_source": dsrc,
                   "age_class": AGE_CLASS.get(disease, ""),
                   "relevance": relevance, "in_scope": in_scope,
                   "is_epigenomic": "Y" if assay in EPI else "N",
                   "gse_title": gctx.split(" || ")[0][:300], "gse_pubmed": pmid,
                   "gse_date": pdat, "gse_taxon": taxon, "gse_families": fams}
            if w is None:
                w = csv.DictWriter(fo, fieldnames=list(rec.keys()), delimiter="\t", extrasaction="ignore")
                w.writeheader()
            w.writerow(rec); n += 1
    print(f"annotated {n} -> {out}")

if __name__ == "__main__":
    main()
