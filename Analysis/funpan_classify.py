#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_classify.py
==================
Phenotype classification engine for FunPan notebooks.

Implements a weighted rules engine with:
1. Weighted scoring (hard/medium/soft signals)
2. Two-tier classification: Provenance x Phenotype
3. NCBI BioProject enrichment, geocoding, visualization
4. Audit trail with matched patterns and confidence scores

Extracted from Theme1_Analysis_Functions.py.
"""

from __future__ import annotations
import os
import re
import json
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import requests


# =============================================================================
# NCBI EUTILS API FUNCTIONS (for BioProject enrichment)
# =============================================================================

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_SESSION = requests.Session()  # reuse TCP connections


def _params(base: dict, email: Optional[str], api_key: Optional[str]) -> dict:
    """Add email and api_key to NCBI API params."""
    p = dict(base)
    if email:   p["email"] = email
    if api_key: p["api_key"] = api_key
    p.setdefault("tool", "strain_classifier")
    return p


def _get(endpoint: str, params: dict, timeout: int = 30,
         max_retries: int = 6, backoff_base: float = 0.5) -> requests.Response:
    """GET with exponential backoff on 429/5xx and respect Retry-After."""
    url = f"{EUTILS}/{endpoint}"
    for attempt in range(max_retries):
        r = _SESSION.get(url, params=params, timeout=timeout)
        if r.status_code == 200:
            return r
        if r.status_code in (429, 500, 502, 503, 504):
            retry_after = r.headers.get("Retry-After")
            try:
                delay = float(retry_after) if retry_after else backoff_base * (2 ** attempt)
            except Exception:
                delay = backoff_base * (2 ** attempt)
            time.sleep(min(delay, 30))
            continue
        r.raise_for_status()
    r.raise_for_status()


def esearch_ids(db: str, term: str, email: Optional[str], api_key: Optional[str], max_retries: int = 3) -> List[str]:
    """Search NCBI database and return IDs with retry on JSON errors."""
    for attempt in range(max_retries):
        try:
            resp = _get("esearch.fcgi", _params({"db": db, "term": term, "retmode": "json"}, email, api_key))
            js = resp.json()
            return js.get("esearchresult", {}).get("idlist", []) or []
        except json.JSONDecodeError as e:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            print(f"Warning: esearch failed for {term} after {max_retries} retries: {e}")
            return []
        except Exception as e:
            print(f"Warning: esearch failed for {term}: {e}")
            return []
    return []


def elink_ids(dbfrom: str, db: str, ids: List[str], email: Optional[str], api_key: Optional[str], max_retries: int = 3) -> List[str]:
    """Link IDs from one database to another with retry on JSON errors."""
    if not ids: return []
    for attempt in range(max_retries):
        try:
            resp = _get("elink.fcgi", _params({"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"}, email, api_key))
            js = resp.json()
            break  # Success
        except json.JSONDecodeError as e:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))  # Exponential backoff: 2, 4, 6 seconds
                continue
            print(f"Warning: elink failed for {dbfrom}->{db} after {max_retries} retries: {e}")
            return []
        except Exception as e:
            print(f"Warning: elink failed for {dbfrom}->{db}: {e}")
            return []
    out = []
    for ls in js.get("linksets", []):
        for ldb in ls.get("linksetdbs", []):
            if ldb.get("dbto") == db:
                out += [str(x) for x in ldb.get("links", [])]
    # de-dup preserve order
    seen = set()
    uniq = []
    for u in out:
        if u not in seen:
            seen.add(u)
            uniq.append(u)
    return uniq


def _strip_ns(root: ET.Element) -> ET.Element:
    """Strip XML namespaces."""
    for el in root.iter():
        if isinstance(el.tag, str) and "}" in el.tag:
            el.tag = el.tag.split("}", 1)[1]
    return root


def _collect_text(node: Optional[ET.Element]) -> Optional[str]:
    """Recursively collect all text from an XML element."""
    if node is None: return None
    parts = []
    def rec(x: ET.Element):
        if x.text and x.text.strip(): parts.append(x.text.strip())
        for c in list(x): rec(c)
        if x.tail and x.tail.strip(): parts.append(x.tail.strip())
    rec(node)
    s = " ".join(parts).strip()
    return s or None


def bioproject_summary(uid_or_acc: str, email: Optional[str], api_key: Optional[str]) -> Dict[str, Optional[str]]:
    """
    XML-only: one efetch call to get accession, uid, title, description, grants.
    Accepts either numeric UID or PRJ* accession as id=; parses <ArchiveID> for both.
    """
    out = {"accession": None, "uid": None, "title": None, "description": None, "grants": None}
    try:
        xml = _get("efetch.fcgi", _params({"db": "bioproject", "id": uid_or_acc, "retmode": "xml"}, email, api_key)).text
        root = _strip_ns(ET.fromstring(xml))

        arch = root.find(".//Project/ProjectID/ArchiveID")
        acc = arch.get("accession") if arch is not None else None
        uid = arch.get("id") if arch is not None else None

        title = root.findtext(".//Project/ProjectDescr/Title")
        desc = _collect_text(root.find(".//Project/ProjectDescr/Description"))

        grants_txts = []
        for g in root.findall(".//Project/ProjectDescr/Grant"):
            gid = g.get("GrantId") or (g.findtext("GrantId") or "")
            gtitle = g.findtext("Title") or ""
            agency_el = g.find("Agency")
            agency_text = (agency_el.text or "").strip() if (agency_el is not None and agency_el.text) else ""
            abbr = (agency_el.get("abbr") if agency_el is not None else None) or ""
            meta = []
            if gid: meta.append(f"Grant ID {gid}")
            if agency_text and abbr and abbr != agency_text:
                meta.append(f"{agency_text} [{abbr}]")
            elif agency_text:
                meta.append(agency_text)
            bits = []
            if gtitle: bits.append(gtitle)
            if meta:   bits.append("(" + ", ".join(meta) + ")")
            s = " ".join(bits).strip()
            if s:
                grants_txts.append(s)
        grants = " | ".join(grants_txts) if grants_txts else None

        out.update({"accession": acc or uid_or_acc, "uid": uid, "title": title, "description": desc, "grants": grants})
    except Exception:
        pass
    return out


def fetch_bioprojects_for_assembly(assembly_acc: str, email: Optional[str] = None,
                                    api_key: Optional[str] = None, pause: float = 0.35) -> List[Dict[str, Optional[str]]]:
    """
    assembly accession -> BioProject(s) via esearch/elink -> ONE XML efetch per BioProject
    """
    acc = (assembly_acc or "").strip()
    if not acc: return []
    asm_ids = esearch_ids("assembly", f"{acc}[Assembly Accession]", email, api_key) or esearch_ids("assembly", acc, email, api_key)
    if not asm_ids: return []
    results: List[Dict[str, Optional[str]]] = []

    bp_uids = elink_ids("assembly", "bioproject", asm_ids, email, api_key)
    if not bp_uids:
        bs = elink_ids("assembly", "biosample", asm_ids, email, api_key)
        if bs:
            bp_uids = elink_ids("biosample", "bioproject", bs, email, api_key)

    seen = set()
    for uid in bp_uids:
        if uid in seen: continue
        seen.add(uid)
        results.append(bioproject_summary(uid, email, api_key))
        time.sleep(pause)

    if not results:
        try:
            js = _get("esummary.fcgi", _params({"db": "assembly", "id": ",".join(asm_ids), "retmode": "json"}, email, api_key)).json()
            text = json.dumps(js).upper()
            accs = set()
            for tok in text.replace(",", " ").split():
                if tok.startswith(("PRJNA", "PRJEB", "PRJDB")):
                    accs.add(tok.strip('";,[] '))
            for a in sorted(accs):
                results.append(bioproject_summary(a, email, api_key))
                time.sleep(pause)
        except Exception:
            pass

    def richness(x):
        return int(bool(x.get("description"))) + int(bool(x.get("title"))) + int(bool(x.get("grants")))
    dedup = {}
    for r in results:
        key = r.get("accession") or r.get("uid")
        if not key: continue
        cur = dedup.get(key)
        if cur is None or richness(r) > richness(cur):
            dedup[key] = r
    return list(dedup.values())


def enrich_df_with_bioproject(df: pd.DataFrame,
                              assembly_col: str = "Assembly Accession",
                              email: Optional[str] = None,
                              api_key: Optional[str] = None,
                              pause: float = 0.35) -> pd.DataFrame:
    """
    Adds BioProject Accession/Title/Description/Grants columns (XML-only retrieval).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with assembly metadata.
    assembly_col : str
        Column name containing assembly accessions.
    email : str, optional
        Email for NCBI API.
    api_key : str, optional
        NCBI API key for higher rate limits.
    pause : float
        Pause between API calls in seconds.

    Returns
    -------
    pd.DataFrame
        DataFrame with added BioProject columns.
    """
    out = df.copy()
    cache: Dict[str, List[Dict[str, Optional[str]]]] = {}
    bp_accs = []
    bp_titles = []
    bp_descs = []
    bp_grants = []

    for acc in out[assembly_col].astype(str).fillna(""):
        key = acc.strip()
        if not key:
            bp_accs.append(None)
            bp_titles.append(None)
            bp_descs.append(None)
            bp_grants.append(None)
            continue
        if key not in cache:
            cache[key] = fetch_bioprojects_for_assembly(key, email=email, api_key=api_key, pause=pause)
        items = cache[key]
        bp_accs.append(" | ".join(filter(None, [x.get("accession") for x in items])) or None)
        bp_titles.append(" | ".join(filter(None, [x.get("title") for x in items])) or None)
        bp_descs.append(" | ".join(filter(None, [x.get("description") for x in items])) or None)
        bp_grants.append(" | ".join(filter(None, [x.get("grants") for x in items])) or None)

    out["BioProject Accession"] = bp_accs
    out["BioProject Title"] = bp_titles
    out["BioProject Description"] = bp_descs
    out["BioProject Grants"] = bp_grants
    return out


def fetch_wgs_title(wgs_acc: str, email: Optional[str] = None,
                    api_key: Optional[str] = None) -> Dict[str, Optional[str]]:
    """
    Fetch the title/definition for a WGS project accession from NCBI.

    Returns both the short definition and the full reference title (if available).

    Parameters
    ----------
    wgs_acc : str
        WGS project accession (e.g., "LKBF01", "ABDB01").
    email : str, optional
        Email for NCBI API.
    api_key : str, optional
        NCBI API key for higher rate limits.

    Returns
    -------
    dict
        Dictionary with "definition" (short) and "reference_title" (full project title).
    """
    result = {"definition": None, "reference_title": None}
    acc = (wgs_acc or "").strip()
    if not acc:
        return result

    # Normalize accession: LKBF01 -> LKBF00000000 (master record)
    # WGS accessions are typically 4-6 letters + 2 digits for version
    master_acc = None
    if re.match(r'^[A-Z]{4,6}\d{2}$', acc.upper()):
        # Convert LKBF01 to LKBF00000000
        prefix = re.match(r'^([A-Z]{4,6})', acc.upper()).group(1)
        master_acc = f"{prefix}00000000"

    try:
        # Search nuccore for the master accession
        search_term = master_acc or acc
        ids = esearch_ids("nuccore", f"{search_term}[Accession]", email, api_key)
        if not ids:
            ids = esearch_ids("nuccore", search_term, email, api_key)
        if not ids:
            return result

        # Fetch the full GenBank record to get the reference title
        gb_text = _get("efetch.fcgi", _params({
            "db": "nuccore",
            "id": ids[0],
            "rettype": "gb",
            "retmode": "text"
        }, email, api_key)).text

        # Parse DEFINITION
        def_match = re.search(r'DEFINITION\s+(.+?)(?=\nACCESSION)', gb_text, re.DOTALL)
        if def_match:
            definition = re.sub(r'\s+', ' ', def_match.group(1)).strip()
            result["definition"] = definition

        # Parse first REFERENCE TITLE (usually the main publication/project title)
        # Look for TITLE in REFERENCE section
        title_match = re.search(r'REFERENCE\s+1[^\n]*\n.*?TITLE\s+(.+?)(?=\n\s{0,2}[A-Z])', gb_text, re.DOTALL)
        if title_match:
            ref_title = re.sub(r'\s+', ' ', title_match.group(1)).strip()
            result["reference_title"] = ref_title

    except Exception:
        pass

    return result


def enrich_df_with_wgs_title(df: pd.DataFrame,
                              wgs_col: str = "WGS project accession",
                              email: Optional[str] = None,
                              api_key: Optional[str] = None,
                              pause: float = 0.35) -> pd.DataFrame:
    """
    Add WGS Project Title columns by fetching from NCBI GenBank records.

    Adds two columns:
    - "WGS Definition": Short definition (e.g., "Aspergillus niger strain L2, whole genome shotgun sequencing project")
    - "WGS Reference Title": Full reference title (e.g., "Comparative genomics of Aspergillus niger strains...")

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with WGS project accessions.
    wgs_col : str
        Column name containing WGS project accessions.
    email : str, optional
        Email for NCBI API.
    api_key : str, optional
        NCBI API key for higher rate limits.
    pause : float
        Pause between API calls in seconds.

    Returns
    -------
    pd.DataFrame
        DataFrame with added WGS title columns.
    """
    out = df.copy()

    if wgs_col not in out.columns:
        print(f"Warning: Column '{wgs_col}' not found in DataFrame")
        return out

    cache: Dict[str, Dict[str, Optional[str]]] = {}
    definitions = []
    ref_titles = []

    total = len(out)
    for i, acc in enumerate(out[wgs_col].astype(str).fillna("")):
        key = acc.strip()
        if not key or key.lower() in ("nan", "none", ""):
            definitions.append(None)
            ref_titles.append(None)
            continue

        if key not in cache:
            if (i + 1) % 20 == 0:
                print(f"  Fetching WGS titles: {i + 1}/{total}...")
            cache[key] = fetch_wgs_title(key, email=email, api_key=api_key)
            time.sleep(pause)

        info = cache[key]
        definitions.append(info.get("definition"))
        ref_titles.append(info.get("reference_title"))

    out["WGS Definition"] = definitions
    out["WGS Reference Title"] = ref_titles
    return out


# =============================================================================
# CLASSIFICATION RESULT DATACLASS (Two-tier: Provenance + Phenotype)
# =============================================================================

@dataclass
class ClassificationResult:
    """
    Result of two-tier classification:

    A) Provenance (where physically isolated from):
       - Clinical: host-derived (human/animal patient samples)
       - Environmental: soil/air/water/plant debris/built environment
       - Industrial-origin: factory/fermentor/production site/fermentation starters
       - Culture-derived: ATCC/NRRL/CBS/"from culture"/maintained lab stocks
       - Unknown: insufficient provenance information

    B) Phenotype (biological association - for pan-GWAS):
       - Human-pathogenic: isolated from human infection OR validated virulence
       - Animal-pathogenic: isolated from animal infection/disease
       - Plant-pathogenic: mycotoxigenic, crop contamination, plant disease
       - Industrial-trait: citric acid producer, enzyme producer, fermentation
       - Environmental: general environmental isolate, no specific association
       - Lab: reference strain, lab workhorse, model organism
       - Unknown: insufficient phenotype information
    """
    # Provenance classification
    provenance: str                                 # Where isolated from
    provenance_confidence: str                      # "high", "medium", "low"
    provenance_scores: Dict[str, float]             # Scores per provenance category
    provenance_matched: Dict[str, List[str]]        # Matched patterns for provenance

    # Phenotype classification (for pan-GWAS)
    phenotype: str                                  # Biological association
    phenotype_confidence: str                       # "high", "medium", "low"
    phenotype_scores: Dict[str, float]              # Scores per phenotype category
    phenotype_matched: Dict[str, List[str]]         # Matched patterns for phenotype

    # Metadata
    field_provenance: Dict[str, List[str]]          # Which fields contributed
    negations_applied: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame storage."""
        return {
            # Provenance (where isolated from)
            "Provenance": self.provenance,
            "ProvenanceConfidence": self.provenance_confidence,
            "ProvenanceScores": json.dumps(self.provenance_scores),
            "ProvenanceMatched": json.dumps(self.provenance_matched),
            # Phenotype (for pan-GWAS)
            "Phenotype": self.phenotype,
            "PhenotypeConfidence": self.phenotype_confidence,
            "PhenotypeScores": json.dumps(self.phenotype_scores),
            "PhenotypeMatched": json.dumps(self.phenotype_matched),
            # Legacy columns for backward compatibility
            "IsolationClass": self.phenotype,  # Phenotype is the main analysis class
            "UsageClass": self.phenotype,
            # Metadata
            "FieldProvenance": json.dumps(self.field_provenance),
            "NegationsApplied": json.dumps(self.negations_applied),
        }


# =============================================================================
# WEIGHT CONSTANTS
# =============================================================================

WEIGHT_HARD = 3.0      # Strong evidence
WEIGHT_MEDIUM = 1.5    # Moderate evidence
WEIGHT_SOFT = 0.5      # Weak evidence (common/ambiguous terms)
WEIGHT_OVERRIDE = 10.0 # Known strain override


# =============================================================================
# TWO-TIER PATTERN DEFINITIONS
# =============================================================================
# A) Provenance: WHERE the strain was physically isolated from
# B) Phenotype: WHAT biological association/trait the strain has

# ===========================================================================
# A) PROVENANCE PATTERNS (physical source location)
# ===========================================================================

# ---------------------------------------------------------------------------
# PROVENANCE: Clinical (any host-derived sample)
# ---------------------------------------------------------------------------
PROV_CLINICAL = [
    # Human clinical samples
    (r"\bclinical(?:ly)?\s+isolat(?:e|ed)\b", "clinical isolate"),
    (r"\bclinical\s+(?:laboratory\s+)?sample\b", "clinical sample"),
    (r"\bpatient\b", "patient"),
    (r"\bhospital(?:ized)?\b", "hospital"),
    (r"\bICU\b", "ICU"),
    (r"\bBALF?\b|\bbronchoalveolar\s+lavage\b", "BAL"),
    (r"\bsputum\b", "sputum"),
    (r"\btracheal\b|\bbronchial\b", "respiratory specimen"),
    (r"\blung\s+tissue\b|\bpulmonary\b", "lung tissue"),
    (r"\bcornea\b|\bcontact\s+lens\b|\beye\b", "eye sample"),
    (r"\bsinus\b|\bparanasal\b|\bnasal\b", "sinus/nasal"),
    (r"\bwound\b|\bulcer\b|\bpus\b|\babscess\b", "wound"),
    (r"\bblood\s+(?:culture|sample)\b|\bblood\b", "blood"),
    (r"\burine\b", "urine"),
    (r"\bswab\b|\bbiopsy\b", "clinical specimen"),
    (r"\bear\b(?!\s*of\s*corn)", "ear"),
    (r"\boral\b|\bgingiva\b|\bdental\b", "oral"),
    (r"\bdrain(?:age)?\b|\bperitoneal\b|\babdominal\b", "abdominal"),
    # Fecal samples (human/animal)
    (r"\bhuman\s+fec(?:es|al)\b|\bhu[mn]an\s+feces\b", "human feces"),
    (r"\bfec(?:es|al)\b|\bstool\b", "fecal"),
    # Host-derived samples
    (r"\bhomo\s+sapiens\b|\bhuman\b", "human"),
    # Animal clinical samples
    (r"\bvet(?:erinary)?\s+(?:clinic|hospital|sample)\b", "veterinary clinical"),
    (r"\bdiseased\s+(?:animal|bird|chicken)\b", "diseased animal"),
    (r"\binfected\s+(?:animal|bird|tissue|host)\b", "infected host"),
]

# ---------------------------------------------------------------------------
# PROVENANCE: Environmental (soil/air/water/built environment)
# ---------------------------------------------------------------------------
PROV_ENVIRONMENTAL = [
    # Soil/Earth
    (r"\bsoil\b", "soil"),
    (r"\bearth\b|\bground\b", "ground"),
    (r"\brhizosphere\b", "rhizosphere"),
    (r"\bcompost\b|\bdecaying\b", "compost"),
    (r"\bpermafrost\b", "permafrost"),
    # Air
    (r"\bair\b(?!\s*(?:conditioning|filter))", "air"),
    (r"\batmospher(?:e|ic)\b", "atmospheric"),
    (r"\baerosol\b|\bairborne\b", "airborne"),
    (r"\bdust\b", "dust"),
    # Water
    (r"\bwater\b(?!\s*(?:treatment|plant))", "water"),
    (r"\baquatic\b|\bfreshwater\b|\bmarine\b", "aquatic"),
    (r"\bsediment\b", "sediment"),
    # Plant debris/agriculture (non-disease context)
    (r"\bplant\s+debris\b|\blitter\b|\bdetritus\b", "plant debris"),
    (r"\bforest\b|\bwoodland\b|\brain[-\s]?forest\b", "forest"),
    (r"\bagricultural\s+(?:soil|field|land)\b", "agricultural field"),
    (r"\bcrop\s+field\b|\bfarm(?:land)?\b", "farmland"),
    # Food/plant products (environmental, not industrial production)
    (r"\bcashew\s*nuts?\b|\bpeanuts?\b|\bnuts?\b", "nuts"),
    (r"\bfood\b(?!\s*(?:factory|processing|production))", "food"),
    (r"\bendophyte\b|\bepiphyte\b", "endophyte"),
    (r"\binsect\b|\bbee\b|\bhoney\s*bee\b", "insect-associated"),
    (r"\broot\b|\bleaf\b|\bplant\b", "plant"),
    # Built environment
    (r"\bindoor\b|\bbuilding\b|\bhouse\b", "indoor/building"),
    (r"\boutdoor\b|\benvironment(?:al)?\b", "outdoor/environmental"),
    (r"\bsurvey\b|\bcollection\b|\bsampling\b", "environmental survey"),
    (r"\bISS\b|\bInternational\s+Space\s+Station\b", "ISS"),
    # ENVO ontology terms (common environmental descriptors)
    (r"\bENVO:\d+\b", "ENVO ontology term"),
]

# ---------------------------------------------------------------------------
# PROVENANCE: Industrial-origin (factory/production site)
# ---------------------------------------------------------------------------
PROV_INDUSTRIAL = [
    # Fermentation sites
    (r"\bkoji\b|\bnuru?k\b|\bdaqu\b|\bqu\b|\bmeju\b", "fermentation starter"),
    (r"\bferment(?:ation|er|or)?\s+(?:plant|facility|site)\b", "fermentation facility"),
    (r"\bferment(?:ed|ing)?\s+foods?\b", "fermented food"),
    (r"\bsoy\s*sauce\s+(?:factory|plant|production)\b", "soy sauce production"),
    (r"\bbrew(?:ery|ing)\b|\bwinery\b", "brewery/winery"),
    (r"\bdistillery\b", "distillery"),
    # Traditional Asian fermented foods/beverages (Aspergillus-based)
    (r"\bsake\b|\bnihonshu\b|\brice\s+wine\b", "sake/rice wine"),
    (r"\bmiso\b|\bdoenjang\b", "miso/fermented soybean paste"),
    (r"\bsoy\s*sauce\b|\bshoyu\b|\bganjang\b", "soy sauce"),
    (r"\bamazake\b|\bshiokoji\b", "fermented rice product"),
    (r"\bvinegar\b.*\bproduction\b|\brice\s+vinegar\b", "vinegar production"),
    (r"\bbaijiu\b|\bhuangjiu\b|\bchinese\s+(?:wine|liquor)\b", "Chinese fermented alcohol"),
    (r"\bmakgeolli\b|\btakju\b|\bkorean\s+rice\s+wine\b", "Korean rice wine"),
    (r"\btempeh\b|\boncom\b", "tempeh/fermented soybean"),
    # Fermented beverages/foods
    (r"\bpu[-\s]?er\s+tea\b|\bkombucha\b|\btea\s+fung(?:us|i)\b", "fermented tea"),
    # Production facilities
    (r"\bfactory\b|\bplant\b(?=.*production)", "factory"),
    (r"\bbioreactor\b|\bfermentor\b|\bfermenter\b", "bioreactor"),
    (r"\bproduction\s+(?:facility|site|plant|line)\b", "production facility"),
    (r"\bindustrial\s+(?:isolate|strain|facility|site)\b", "industrial site"),
    # Food processing
    (r"\bcheese\s+(?:factory|production|facility)\b", "cheese production"),
    (r"\bdairy\s+(?:plant|factory|production)\b", "dairy production"),
    (r"\bfood\s+(?:factory|processing|production)\b", "food processing"),
]

# ---------------------------------------------------------------------------
# PROVENANCE: Culture-derived (lab stocks/collections)
# ---------------------------------------------------------------------------
PROV_CULTURE = [
    # Culture collections
    (r"\bATCC\b", "ATCC"),
    (r"\bCBS\b", "CBS"),
    (r"\bNRRL\b", "NRRL"),
    (r"\bNBRC\b", "NBRC"),
    (r"\bJCM\b", "JCM"),
    (r"\bDSM\b|\bDSMZ\b", "DSMZ"),
    (r"\bFGSC\b", "FGSC"),
    (r"\bIMI\b", "IMI"),
    # Lab/culture context
    (r"\bfrom\s+culture\b|\bculture\s+collection\b", "culture collection"),
    (r"\btype\s+strain\b", "type strain"),
    (r"\breference\s+(?:strain|isolate)\b", "reference strain"),
    (r"\blab(?:oratory)?\s+(?:strain|stock|isolate)\b", "lab strain"),
    (r"\bmaintained\s+(?:strain|culture)\b", "maintained culture"),
    (r"\bderived\s+from\b.*\bstrain\b", "derived strain"),
    (r"\bparent(?:al)?\s+strain\b", "parental strain"),
    # Standalone lab/laboratory (as isolation source)
    (r"\blab(?:oratory)?\b", "lab"),
]

# ===========================================================================
# B) PHENOTYPE PATTERNS (biological association - for pan-GWAS)
# ===========================================================================

# ---------------------------------------------------------------------------
# PHENOTYPE: Human-pathogenic association
# ---------------------------------------------------------------------------
PHENO_HUMAN_PATHOGEN = [
    (r"\baspergillosis\b", "aspergillosis"),
    (r"\baspergilloma\b", "aspergilloma"),
    (r"\bCAPA\b", "CAPA"),
    (r"\bIPA\b|\binvasive\s+(?:pulmonary\s+)?aspergillosis\b", "invasive aspergillosis"),
    (r"\bhuman[-\s]?pathogen(?:ic)?\b", "human pathogen"),
    (r"\bpathogenic\s+to\s+humans?\b", "pathogenic to humans"),
    (r"\bmajor\s+etiological\s+agent\b", "etiological agent"),
    (r"\bkeratitis\b|\bfungal\s+keratitis\b", "keratitis"),
    (r"\binvasive\s+(?:fungal\s+)?(?:infection|disease)\b", "invasive infection"),
    (r"\bcystic\s+fibrosis\b", "CF-associated"),
    (r"\bNRZ[-\s]?\d+", "NRZ reference (pathogen)"),
    # Note: FDA-ARGOS removed - database includes environmental/food samples
    (r"(?:\banti[-\s]?fungal|azole|triazole)\b.*\bresistan\w+\b", "antifungal resistance"),
    (r"\bvirulen(?:t|ce)\b", "virulence"),
    # "pathogen" requires context - avoid matching labels like "Pathogen: environmental sample"
    (r"\b(?:known|opportunistic|major|primary|human)\s+pathogen\b", "pathogen (with context)"),
    (r"\bclinical\s+(?:significance|relevance|importance)\b", "clinical significance"),
]

# ---------------------------------------------------------------------------
# PHENOTYPE: Animal-pathogenic association
# ---------------------------------------------------------------------------
PHENO_ANIMAL_PATHOGEN = [
    (r"\banimal\s+(?:pathogen|disease|infection)\b", "animal pathogen"),
    (r"\bavian\s+aspergillosis\b", "avian aspergillosis"),
    (r"\bvet(?:erinary)?\s+pathogen\b", "veterinary pathogen"),
    (r"\bkakapo\b", "kakapo pathogen"),
    (r"\btick\s+(?:egg|pathogen)\b|\bRhipicephalus\b", "tick pathogen"),
    (r"\bpoultry\s+(?:disease|pathogen)\b", "poultry pathogen"),
    (r"\bzoonotic\b|\bzoonosis\b", "zoonotic"),
    (r"\blivestock\s+(?:disease|pathogen)\b", "livestock pathogen"),
    (r"\bwildlife\s+(?:disease|pathogen)\b", "wildlife pathogen"),
]

# ---------------------------------------------------------------------------
# PHENOTYPE: Plant-pathogenic association (mycotoxigenic, crop contamination)
# ---------------------------------------------------------------------------
PHENO_PLANT_PATHOGEN = [
    # Mycotoxins
    (r"\baflatoxin(?:s|ogenic)?\b", "aflatoxin"),
    (r"\bmycotoxin(?:s|ogenic)?\b", "mycotoxin"),
    (r"\bochratoxin\b", "ochratoxin"),
    (r"\bfumonisin\b", "fumonisin"),
    (r"\bcyclopiazonic\s+acid\b", "cyclopiazonic acid"),
    # Plant disease
    (r"\bplant\s+(?:pathogen|disease)\b", "plant pathogen"),
    (r"\bcrop\s+(?:contamination|disease|pathogen)\b", "crop contamination"),
    (r"\bgrain\s+(?:contamination|spoilage)\b", "grain contamination"),
    (r"\bpost[-\s]?harvest\s+(?:rot|decay|disease)\b", "postharvest disease"),
    (r"\bseed(?:borne)?\s+pathogen\b", "seedborne pathogen"),
    (r"\brot\b.*\b(?:corn|maize|peanut|grain)\b", "crop rot"),
    # Contamination context
    (r"\bcontaminated\s+(?:grain|corn|maize|peanut|food)\b", "contaminated crop"),
    (r"\btoxigenic\b", "toxigenic"),
    # Crop/nut sources (A. flavus on these = aflatoxin contamination context)
    (r"\bcashew\s*nut", "cashew (crop contamination)"),
    (r"\bpeanut", "peanut (crop contamination)"),
    (r"\bpistachio", "pistachio (crop contamination)"),
    (r"\b(?:moldy|mouldy)\s+(?:corn|maize|peanut|nut|grain)", "moldy crop"),
    (r"\bcorn\b(?!\s*(?:field|farm))", "corn (crop contamination)"),
    (r"\bmaize\b(?!\s*(?:field|farm))", "maize (crop contamination)"),
    (r"\btree\s*nut", "tree nut (crop contamination)"),
    (r"\balmond", "almond (crop contamination)"),
    (r"\bwalnut", "walnut (crop contamination)"),
    (r"\bhazelnut", "hazelnut (crop contamination)"),
    (r"\btiger\s*nut", "tiger nut (crop contamination)"),
]

# ---------------------------------------------------------------------------
# PHENOTYPE: Industrial trait (production capabilities)
# ---------------------------------------------------------------------------
PHENO_INDUSTRIAL_TRAIT = [
    # Acid production
    (r"\bcitric\s+acid\b.*\b(?:produc|high|yield)\b", "citric acid producer"),
    (r"\bitaconic\s+acid\b", "itaconic acid"),
    (r"\bgluconic\s+acid\b", "gluconic acid"),
    (r"\bkojic\s+acid\b", "kojic acid"),
    # Enzyme production
    (r"\benzyme\s+(?:produc|hyper[-\s]?produc)\b", "enzyme producer"),
    (r"\b(?:gluco)?amylase\b.*produc\b", "amylase producer"),
    (r"\bcellulase\b.*produc\b", "cellulase producer"),
    (r"\bprotease\b.*produc\b", "protease producer"),
    (r"\blipase\b.*produc\b", "lipase producer"),
    # Fermentation performance
    (r"\bfermentation\s+(?:performance|trait|capability)\b", "fermentation trait"),
    (r"\bsake\b|\bshochu\b|\bmiso\b|\bsoy\s*sauce\b", "traditional fermentation"),
    (r"\bindustrial\s+(?:production|biotechnology|strain)\b", "industrial production"),
    (r"\bhigh\s+(?:productivity|yield|production)\b", "high producer"),
    (r"\bGRAS\b", "GRAS status"),
    # Biocontrol (industrial application)
    (r"\bbiocontrol\b|\bAfla[-\s]?(?:Safe|Guard)\b", "biocontrol"),
    (r"\bnon[-\s]?(?:aflatoxigenic|toxigenic)\b", "non-toxigenic (biocontrol)"),
]

# ---------------------------------------------------------------------------
# PHENOTYPE: Lab (reference/model strains)
# ---------------------------------------------------------------------------
PHENO_LAB = [
    (r"\breference\s+(?:strain|genome|isolate)\b", "reference strain"),
    (r"\btype\s+strain\b", "type strain"),
    (r"\bmodel\s+(?:organism|strain)\b", "model organism"),
    (r"\blab(?:oratory)?\s+(?:strain|horse|workhorse)\b", "lab strain"),
    (r"\bmutagenesis\b|\bmutagenized\b|\bmutant\b", "mutant"),
    (r"\bUV[-\s]?mutat(?:ed|ion)\b", "UV mutant"),
    (r"\bgene\s+(?:deletion|knockout|disruption)\b", "gene knockout"),
    (r"\bexperimental\s+evolution\b", "experimental evolution"),
    (r"\blab(?:oratory)?[-\s]?adapted\b", "lab-adapted"),
    (r"\bgenome\s+(?:sequenc|reference|resource)\b", "genome resource"),
]

# ---------------------------------------------------------------------------
# PHENOTYPE: Environmental (default, no specific association)
# ---------------------------------------------------------------------------
PHENO_ENVIRONMENTAL = [
    (r"\benvironmental\s+(?:isolate|strain|sample)\b", "environmental isolate"),
    (r"^environmental$", "environmental (explicit)"),  # Just "environmental" alone
    (r"\bsoil\s+isolate\b", "soil isolate"),
    (r"\bsaprob(?:e|ic)\b|\bsaprophyt(?:e|ic)\b", "saprophyte"),
    (r"\bbiodiversity\b|\bsurvey\b|\bpopulation\s+study\b", "biodiversity study"),
    (r"\becolog(?:y|ical)\b", "ecological"),
    # Agricultural research contexts (field studies, not lab)
    (r"\bexperimental\s+(?:farm|orchard|field|plot)\b", "experimental farm/orchard"),
    (r"\bagricultural\s+field\b", "agricultural field"),
    (r"\bfield\s+soil\b", "field soil"),
    (r"\bcitrus\s+field\b", "citrus field"),
    (r"\bgrant\s+farm\b|\bwolfskill\b", "research farm"),
]

# ---------------------------------------------------------------------------
# NEGATION PATTERNS (reduce scores)
# ---------------------------------------------------------------------------
NEGATION_PATTERNS = [
    (r"\bnon[-\s]?pathogen(?:ic)?\b", "non-pathogenic"),
    (r"\bnon[-\s]?(?:aflatoxigenic|toxigenic|mycotoxigenic)\b", "non-toxigenic"),
    (r"\bavirulent\b", "avirulent"),
    (r"\battenuated\b", "attenuated"),
    (r"\bnot\s+(?:pathogenic|virulent|toxigenic)\b", "not pathogenic"),
]

# ---------------------------------------------------------------------------
# CULTURE MEDIUM PATTERNS (should not be interpreted as pathogenic context)
# These are laboratory culture media, not sources of infection
# ---------------------------------------------------------------------------
CULTURE_MEDIUM_PATTERNS = [
    (r"\bsheep\s+blood(?:\s+agar)?\b", "sheep blood agar (culture medium)"),
    (r"\bMTB\s+medium\b", "MTB medium (culture medium)"),
    (r"\bPDA\b|\bpotato\s+dextrose\b", "PDA (culture medium)"),
    (r"\bagar\b", "agar (culture medium)"),
    (r"\bRose\s+bengal\b", "Rose bengal (culture medium)"),
    (r"\bczapek\b|\bCzapek\b", "Czapek medium (culture medium)"),
    (r"\bBrain\s+heart\s+infusion\b|\bBHI\b", "BHI (culture medium)"),
    (r"\bYeast\s+extract\s+peptone\b|\bYEP\b", "YEP (culture medium)"),
    (r"\bChromogenic\s+Medium\b", "Chromogenic medium (culture medium)"),
]

# ---------------------------------------------------------------------------
# AGRICULTURAL RESEARCH PATTERNS (not Lab phenotype - these are field studies)
# ---------------------------------------------------------------------------
AGRICULTURAL_RESEARCH_PATTERNS = [
    (r"\bexperimental\s+(?:farm|orchard|field|plot)\b", "experimental farm/orchard"),
    (r"\bagricultural\s+field\b", "agricultural field"),
    (r"\bfield\s+soil\b", "field soil"),
    (r"\bcitrus\s+field\b", "citrus field"),
    (r"\bgrant\s+farm\b|\bwolfskill\b", "research farm"),
]


# =============================================================================
# LEGACY PATTERN DEFINITIONS (for backward compatibility)
# =============================================================================

# ---------------------------------------------------------------------------
# CLINICAL-ASSOCIATED PATTERNS
# ---------------------------------------------------------------------------

CLINICAL_HARD = [
    (r"\bclinical(?:ly)?\s+isolat(?:e|ed)\b", "clinical isolate"),
    (r"\bclinical\s+(?:laboratory\s+)?sample\b", "clinical sample"),
    (r"\bNRZ[-\s]?\d+", "NRZ strain (German National Reference Centre)"),
    (r"\bpatient\b", "patient"),
    (r"\bhospital(?:ized)?\b", "hospital"),
    (r"\bICU\b", "ICU"),
    (r"\baspergillosis\b", "aspergillosis"),
    (r"\baspergilloma\b", "aspergilloma"),
    (r"\bCAPA\b", "CAPA"),
    (r"\bIPA\b", "IPA"),
    (r"\bhuman[-\s]?pathogen\b", "human pathogen"),
    (r"\bpathogenic\s+to\s+humans?\b", "pathogenic to humans"),
    (r"\bmajor\s+etiological\s+agent\b", "major etiological agent"),
    (r"\betiological\s+agent\s+of\b", "etiological agent"),
    (r"\bkeratitis\b", "keratitis"),
    (r"\binvasive\s+(?:fungal\s+)?(?:infection|disease)\b", "invasive infection"),
    (r"\bCOVID-19\b|\bSARS[-\s]?CoV[-\s]?2\b", "COVID-19 associated"),
    (r"\bcystic\s+fibrosis\b", "cystic fibrosis"),
    (r"(?:\b(?:anti[-\s]?fungal|azole|triazole|voriconazole|itraconazole|posaconazole)\b.*\bresistan\w+\b)|(?:\bresistan\w+\b.*\b(?:anti[-\s]?fungal|azole|triazole|voriconazole|itraconazole|posaconazole)\b)", "antifungal resistance"),
]

CLINICAL_MEDIUM = [
    (r"\bBALF?\b|\bbronchoalveolar\s+lavage(?:\s+fluid)?\b", "BAL/BALF"),
    (r"\bsputum\b|\bexpectorated\s+sputum\b", "sputum"),
    (r"\btracheal\s+(?:aspiration|secretion)\b", "tracheal specimen"),
    (r"\bbronchial\s+secretion\b|\bbronchoscopy\b", "bronchial specimen"),
    (r"\blung\s+tissue\b|\bpulmonary\s+lesion\b", "lung tissue"),
    (r"\bcornea\b|\bcontact\s+lens\b", "corneal/eye"),
    (r"\bsinus\b|\bparanasal\b", "sinus"),
    (r"\bwound\b|\bulcer\b|\bpus\b", "wound/ulcer"),
    (r"\bnecrotic\s+tissue\b", "necrotic tissue"),
    (r"\bblood\s+(?:culture|sample|isolate)\b", "blood culture"),
    (r"\brespiratory\s+specimen\b", "respiratory specimen"),
    (r"\bmedical\s+mycolog(?:y|ical)\b", "medical mycology"),
    (r"\bMedical\s+Mycology\s+Research\s+Center\b|\bMMRC\b", "MMRC"),
    (r"\bInstitute\s+of\s+Tropical\s+Medicine\b|\bNekken\b", "ITM/Nekken"),
    # Oral/dental
    (r"\bgingiva\b|\boral\b|\bmouth\b|\bdental\b", "oral/dental"),
    # Abdominal/peritoneal
    (r"\bperitoneal\b|\bintraperitoneal\b|\babdominal\b", "peritoneal/abdominal"),
    (r"\bdrain(?:age)?\b|\babscess\b", "drainage/abscess"),
    (r"\bcavity\b", "body cavity"),
    # Discharge
    (r"\bdischarge\b|\bexudate\b|\bsecretion\b", "discharge/exudate"),
]

CLINICAL_SOFT = [
    (r"\bblood\b", "blood"),
    (r"\burine\b", "urine"),
    (r"\bswab\b", "swab"),
    (r"\bbiopsy\b", "biopsy"),
    (r"\brespiratory\b", "respiratory"),
    (r"\bear\b(?!\s*(?:canal|of\s+corn))", "ear"),
    (r"\beye\b", "eye"),
    (r"\bnose\b|\bnasal\b", "nasal"),
    (r"\bfec(?:es|al)\b|\bstool\b", "fecal"),
]

# ---------------------------------------------------------------------------
# ANIMAL PATHOGEN PATTERNS
# ---------------------------------------------------------------------------

ANIMAL_CONTEXT = [
    (r"\banimal\b", "animal"),
    (r"\bavian\b|\bbird(?:s)?\b", "avian/bird"),
    (r"\bvet(?:erinary)?\b", "veterinary"),
    (r"\bkakapo\b", "kakapo"),
    (r"\btick\b|\bRhipicephalus\b", "tick"),
    (r"\bhen(?:s)?\b|\bchicken(?:s)?\b|\bpoultry\b", "poultry"),
    (r"\bwild(?:life)?\s+rehabilitation\b", "wildlife rehab"),
    (r"\bdog\b|\bcat\b|\bcanine\b|\bfeline\b", "companion animal"),
    (r"\bhorse\b|\bequine\b", "equine"),
    (r"\bcattle\b|\bbovine\b|\bswine\b|\bporcine\b", "livestock"),
]

ANIMAL_DISEASE = [
    (r"\binfect(?:ion|ed|s)\b", "infection"),
    (r"\bdisease(?:d)?\b", "disease"),
    (r"\blesion\b", "lesion"),
    (r"\bnecropsy\b", "necropsy"),
    (r"\bmycosis\b|\bmycotic\b", "mycosis"),
    (r"\bpathogen(?:ic)?\b", "pathogen"),
    # Arthropod egg isolation often indicates opportunistic infection
    (r"\b(?:tick|mite|insect|arthropod)\b.*\beggs?\b", "arthropod eggs"),
    (r"\beggs?\b.*\b(?:tick|mite|insect|arthropod|Rhipicephalus)\b", "arthropod eggs"),
]

ANIMAL_NONDISEASE = [
    (r"\bfecal\s+sample\b|\bfeces\b", "fecal sample (survey)"),
    (r"\bgut\s+(?:fungi|mycobiome|microbiome)\b", "gut microbiome"),
    (r"\bmycobiome\b|\bmicrobiome\b", "microbiome study"),
]

# ---------------------------------------------------------------------------
# PLANT DISEASE-ASSOCIATED PATTERNS
# ---------------------------------------------------------------------------

PLANT_HOSTS = [
    (r"\bplant(?:s)?\b", "plant"),
    (r"\bcrop(?:s)?\b", "crop"),
    (r"\borchard(?:s)?\b", "orchard"),
    (r"\brhizosphere\b", "rhizosphere"),
    (r"\broot(?:s)?\b", "root"),
    (r"\bleaf(?:ves)?\b", "leaf"),
    (r"\bstem(?:s)?\b", "stem"),
    (r"\bfruit(?:s)?\b", "fruit"),
    (r"\bseed(?:s)?\b", "seed"),
    (r"\bkernel(?:s)?\b", "kernel"),
    (r"\bear(?:s)?\b(?=.*(?:corn|maize))", "ear (corn)"),
    (r"\bspike(?:s)?\b", "spike"),
    (r"\bbud(?:s)?\b", "bud"),
    (r"\bpeanut(?:s)?\b", "peanut"),
    (r"\bcashew(?:s)?(?:\s*nuts?)?\b", "cashew"),
    (r"\bpistachio(?:s)?\b", "pistachio"),
    (r"\bcorn(?:s)?\b|\bmaize\b", "corn/maize"),
    (r"\bwheat\b", "wheat"),
    (r"\bcotton(?:seed)?\b", "cotton"),
    (r"\balmond(?:s)?\b", "almond"),
    (r"\bgrape(?:s)?\b", "grape"),
    (r"\bcoffee\b", "coffee"),
    (r"\brice\b", "rice"),
    (r"\bsorghum\b", "sorghum"),
    (r"\bdate(?:s)?(?:\s*fruit)?\b", "date"),
    (r"\btiger\s*nuts?\b", "tiger nut"),
]

PLANT_DISEASE_HARD = [
    (r"\bphytopath(?:ogen|ogenic)\b", "phytopathogen"),
    (r"\bplant[-\s]?pathogen\b", "plant pathogen"),
    (r"\bdiseased\b", "diseased"),
    (r"\bsymptomatic\b", "symptomatic"),
    (r"\b(?:blight|wilt|rot|smut|mildew|canker|rust|anthracnose|scab|damping[-\s]?off)\b", "plant disease"),
    (r"\bmoldy\b|\bmould(?:y|ed)\b", "moldy"),
    (r"\baflatoxin\s+contaminat(?:ion|ed)\b", "aflatoxin contamination"),
]

PLANT_DISEASE_MEDIUM = [
    (r"\binfect(?:ion|ed|s)\b", "infection"),
    (r"\b(?:opportunistic|facultative)\s+pathogen\b", "opportunistic pathogen"),
    (r"\bpost[-\s]?harvest\b", "post-harvest"),
]

PLANT_NONDISEASE = [
    (r"\bendophyte(?:s)?\b|\bendophytic\b", "endophyte"),
    (r"\bepiphyte(?:s)?\b|\bepiphytic\b", "epiphyte"),
    (r"\bphyllosphere\b", "phyllosphere"),
    (r"\bplant\s+microbiome\b", "plant microbiome"),
]

# ---------------------------------------------------------------------------
# INDUSTRIAL-ASSOCIATED PATTERNS
# ---------------------------------------------------------------------------

INDUSTRIAL_FERMENT_HARD = [
    (r"\bferment(?:ed|ation|er|or|ing)?\b", "fermentation"),
    (r"\bkoji\b", "koji"),
    (r"\bmiso\b", "miso"),
    (r"\bsoy\s*sauce\b|\bshoyu\b|\bganjang\b", "soy sauce"),
    # Asian fermentation starters - qu (曲) is used in many forms
    (r"\bnuru?k\b|\bdaqu\b|\bjiuqu\b|\bmeju\b|\bdoenjang\b", "Asian fermented"),
    (r"\b(?:wheat|rice|barley|sorghum|grain)\s+qu\b", "qu fermentation starter"),
    (r"\bqu\s*(?:starter|culture|sample|block|cake)\b", "qu fermentation"),
    (r"\bhuangjiu\b|\bbaijiu\b|\bmakgeolli\b|\bsake\b|\bsoju\b", "Asian liquor"),
    (r"\bqu\b.*\b(?:brew|ferment|isolat)", "qu brewing/fermentation"),
    # Catch "Wheat Qu" or similar patterns where qu is at word boundary
    (r"\bqu\b", "qu (fermentation starter)"),
    (r"\bstarter\s+culture\b", "starter culture"),
    (r"\bbioreactor\b|\bfermento?r\b", "bioreactor/fermentor"),
    (r"\bproduction\s+strain\b", "production strain"),
    (r"\bGRAS\b", "GRAS"),
    (r"\bfood[-\s]?grade\b", "food-grade"),
    (r"\bdomesticat(?:ed|ion)\b", "domesticated"),
    (r"\bcheese\b|\bcamembert\b|\bbrie\b|\broquefort\b|\bgorgonzola\b", "cheese"),
    (r"\bmold[-\s]?ripened\b|\bsurface[-\s]?ripened\b", "ripened cheese"),
    (r"\bdairy\b|\bcreamery\b|\bfromager(?:ie|y)\b", "dairy"),
    (r"\bbrew(?:ery|ing|house)?\b", "brewing"),
    (r"\bwinery\b", "winery"),
    (r"\bvinegar\b", "vinegar"),
    (r"\bliquor\b", "liquor"),
    # Industrial strain indicators
    (r"\bindustrial\s+strain\b", "industrial strain"),
    (r"\bwidely\s+used\b.*\b(?:brew|ferment|produc)", "widely used industrial"),
]

INDUSTRIAL_ENZYME_HARD = [
    (r"\benzyme\s+production\b", "enzyme production"),
    (r"\bglucose\s+oxidase\b|\bGOX\b", "glucose oxidase"),
    (r"\bcellulase\b", "cellulase"),
    (r"\b(?:gluco)?amylase\b|\bamylase\b", "amylase"),
    (r"\bprotease\b", "protease"),
    (r"\blipase\b", "lipase"),
    (r"\bxylanase\b", "xylanase"),
    (r"\bpectinase\b", "pectinase"),
    (r"\blaccase\b", "laccase"),
    (r"\bphytase\b", "phytase"),
    (r"\bglucanase\b", "glucanase"),
    (r"\benzyme\s+(?:hyper|over)[-\s]?produc(?:tion|er)\b", "enzyme hyperproducer"),
    (r"\blignocellulos(?:e|ic)\b", "lignocellulosic"),
    (r"\bbioconversion\b", "bioconversion"),
    (r"\bbiorefin(?:e|ery)\b", "biorefinery"),
]

INDUSTRIAL_CHEM_HARD = [
    (r"\bchemical\s+production\b", "chemical production"),
    (r"\bmetabolite\s+production\b", "metabolite production"),
    (r"\bpenicillin\b", "penicillin"),
    (r"\bcephalosporin\b", "cephalosporin"),
    (r"\blactam\b", "beta-lactam"),
    (r"\bantibiotic\s+production\b", "antibiotic production"),
    (r"\bcitric\s+acid\b", "citric acid"),
    (r"\bitaconic\s+acid\b", "itaconic acid"),
    (r"\bgluconic\s+acid\b", "gluconic acid"),
    (r"\bkojic\s+acid\b", "kojic acid"),
    (r"\bindustrial\s+production\b", "industrial production"),
    (r"\bsecondary\s+metabolite\s+production\b", "secondary metabolite production"),
]

INDUSTRIAL_BIOCONTROL_HARD = [
    (r"\bbiocontrol\b|\bbio[-\s]?control\b", "biocontrol"),
    (r"\bAflasafe\b", "Aflasafe"),
    (r"\bAfla[-\s]?Guard\b", "Afla-Guard"),
    (r"\bnon[-\s]?(?:aflatoxigenic|toxigenic)\b", "non-toxigenic"),
    (r"\bbiostimulant\b", "biostimulant"),
    (r"\bbiofertiliz(?:er)?\b", "biofertilizer"),
    (r"\bplant\s+growth\s+promot(?:ion|ing)\b|\bPGPR\b", "plant growth promotion"),
]

INDUSTRIAL_SOFT = [
    (r"\bfactory\b", "factory"),
    (r"\bproduction\s+room\b", "production room"),
    (r"\bbioprocess\b", "bioprocess"),
]

# ---------------------------------------------------------------------------
# LAB/REFERENCE PATTERNS
# ---------------------------------------------------------------------------

LAB_REFERENCE_HARD = [
    # FDA-ARGOS is a genome project, not necessarily lab strains - removed
    (r"\breference\s+(?:strain|isolate)\b", "reference strain"),
    (r"\btype\s+strain\b", "type strain"),
    (r"\bdiagnostic\s+(?:panel|reference|standard)\b", "diagnostic reference"),
    (r"\bquality[-\s]?controlled\s+reference\b", "QC reference"),
    (r"\blab(?:oratory)?\s+(?:evolution|adaptation|strain)\b", "lab strain"),
    (r"\bexperimental\s+evolution\b", "experimental evolution"),
    (r"\bmutagenesis\b|\bmutagenized\b", "mutagenesis"),
    (r"\bUV[-\s]?mutat(?:ed|ion|ant)\b", "UV mutant"),
    (r"\bgene\s+(?:deletion|knockout|disruption)\b", "gene knockout"),
    (r"\btransform(?:ed|ant)\b", "transformant"),
]

LAB_REFERENCE_MEDIUM = [
    # Actual lab strain indicators (NOT culture collection IDs - those don't indicate lab origin)
    (r"\bparental\s+strain\b", "parental strain"),
    (r"\blaboratory[-\s]?adapted\b", "laboratory adapted"),
    (r"\bpassage\s+\d+\b", "passage number"),
    (r"\bderived\s+from\b.*\b(?:strain|isolate)\b", "derived strain"),
    (r"\bstrain\s+collection\b|\bculture\s+collection\b", "culture collection"),
    (r"\bbiological\s+resource\s+center\b|\bBRC\b", "BRC"),
    (r"\brepository\b", "repository"),
]

LAB_REFERENCE_SOFT = [
    # Removed generic "annotation" patterns - all genomes have annotations
    (r"\blab(?:oratory)?\s+(?:horse|workhorse)\b", "lab horse"),
    (r"\bmodel\s+(?:organism|strain)\b", "model organism"),
]

# ---------------------------------------------------------------------------
# ENVIRONMENTAL PATTERNS
# ---------------------------------------------------------------------------

ENVIRONMENTAL_MEDIUM = [
    (r"\benvironmental\s+(?:sample|isolate|source)\b", "environmental isolate"),
    (r"\bforest\b|\brain[-\s]?forest\b", "forest"),
    (r"\bmountain\b", "mountain"),
    (r"\bbiome\b", "biome"),
    (r"\bsavann?ah?\b", "savanna"),
    (r"\bcoastal\b", "coastal"),
    (r"\bmangrove\b", "mangrove"),
    (r"\bpermafrost\b", "permafrost"),
    (r"\bwild[-\s]?type\b", "wild-type"),
    (r"\binsect[-\s]?associated\b", "insect-associated"),
    (r"\bbee\b", "bee"),
    (r"\bcompost\b", "compost"),
    # Clear environmental sample types (moved from SOFT for proper scoring)
    (r"\bsoil\b", "soil"),
    (r"\bair\b(?!\s*(?:way|borne))", "air"),  # air but not airway
    (r"\bdust\b", "dust"),
    (r"\bindoor\b|\bbuilt\s+environment\b", "indoor/built environment"),
    (r"\bInternational\s+Space\s+Station\b|\bISS\b", "ISS"),
]

ENVIRONMENTAL_SOFT = [
    (r"\bwater\b", "water"),  # water is ambiguous (could be lab water)
    (r"\blake\b|\bpond\b|\bocean\b|\bsea\b", "aquatic"),
    (r"\bgroundwater\b", "groundwater"),
]

# ---------------------------------------------------------------------------
# NEGATION PATTERNS
# ---------------------------------------------------------------------------

PATHOGENIC_NEGATIONS = [
    (r"\bnon[-\s]?pathogenic\b", "non-pathogenic"),
    (r"\bavirulen\w*\b", "avirulent"),
    (r"\battenuated\b", "attenuated"),
    (r"\bsafe\s+strain\b", "safe strain"),
    (r"\bnon[-\s]?(?:tox(?:ic|igenic)|aflatoxigenic)\b", "non-toxigenic"),
]

MYCOTOXIN_NEGATIONS = [
    (r"\bnon[-\s]?aflatoxigenic\b", "non-aflatoxigenic"),
    (r"\bnon[-\s]?toxigenic\b", "non-toxigenic"),
    (r"\batoxigenic\b", "atoxigenic"),
]

# ---------------------------------------------------------------------------
# TAXONOMY/SURVEY CONTEXT
# ---------------------------------------------------------------------------

TAXONOMY_CONTEXT = [
    (r"\btaxonom(?:y|ic)\b", "taxonomy"),
    (r"\bspecies\s+limit\b", "species delimitation"),
    (r"\bphylogenet(?:ic|ics)\b", "phylogenetics"),
    (r"\bpopulation\s+(?:genomics|structure|genetics)\b", "population study"),
    (r"\bmating[-\s]?type\b", "mating type"),
    (r"\bsexual\s+reproduction\b", "sexual reproduction"),
    (r"\bheterokaryon\b", "heterokaryon"),
    (r"\bdiploid\b|\bmeiosis\b", "ploidy study"),
]

SURVEY_CONTEXT = [
    (r"\bcommunity\s+(?:composition|structure)\b", "community study"),
    (r"\bdiversity\b", "diversity study"),
    (r"\bsurvey\b", "survey"),
    (r"\bmetagenom\w+\b", "metagenomics"),
    (r"\bamplicon\b", "amplicon"),
    (r"\bITS1?\b|\bITS2\b", "ITS sequencing"),
    (r"\bmetabarcoding\b", "metabarcoding"),
]

# ---------------------------------------------------------------------------
# MYCOTOXIN PATTERNS
# ---------------------------------------------------------------------------

MYCOTOXIN_POSITIVE = [
    (r"\baflatoxin(?:s)?\b", "aflatoxin"),
    (r"\baflatoxigen(?:ic|esis)\b", "aflatoxigenic"),
    (r"\bmycotoxin(?:s)?\b", "mycotoxin"),
    (r"\bochratoxin\b", "ochratoxin"),
    (r"\bfumonisin\b", "fumonisin"),
    (r"\bcyclopiazonic\s+acid\b|\bCPA\b", "cyclopiazonic acid"),
    (r"\bsterigmatocystin\b", "sterigmatocystin"),
    (r"\bgliotoxin\b", "gliotoxin"),
]

# ---------------------------------------------------------------------------
# ESTABLISHED STRAIN OVERRIDES
# ---------------------------------------------------------------------------

ESTABLISHED_STRAINS = [
    # Industrial strains (enzyme/chemical production) - patterns match strain names in metadata
    {"name": "Aspergillus niger CBS 513.88", "patterns": [r"\bCBS\s*513\.?88\b", r"\bFGSC\s*A1513\b", r"\bASM285v2\b"], "class": "Industrial", "subtype": "Enzyme production"},
    {"name": "Aspergillus niger ATCC 1015", "patterns": [r"\bATCC\s*1015(?:\b|[D-])", r"\bUSDA\s*3528\.7\b", r"\bFGSC\s*A1144\b"], "class": "Industrial", "subtype": "Chemical production"},
    {"name": "Aspergillus terreus ATCC 20542", "patterns": [r"\bATCC\s*20542\b"], "class": "Industrial", "subtype": "Chemical production"},
    # Industrial fermentation strains
    {"name": "Aspergillus oryzae RIB40", "patterns": [r"\brib[-\s]?40\b", r"\bNBRC\s*100959\b", r"\bATCC\s*42149\b"], "class": "Industrial", "subtype": "Fermentation"},
    # Environmental strains (biocontrol)
    {"name": "Aspergillus flavus NRRL 21882", "patterns": [r"\bNRRL\s*21882\b", r"\bAfla[-\s]?Guard\b"], "class": "Environmental", "subtype": "Biocontrol"},
    {"name": "Aspergillus flavus AF36", "patterns": [r"\bAF[-\s]?36\b", r"\bNRRL\s*18543\b"], "class": "Environmental", "subtype": "Biocontrol"},
    # Note: FDA-ARGOS removed - it's a database with samples from various sources;
    # phenotype should be determined by actual isolation source, not project name
    # Reference clinical strains (reclassified from Lab → original phenotype)
    {"name": "Aspergillus fumigatus Af293", "patterns": [r"\bAf\s*293\b|\bAf293\b", r"\bATCC\s*MYA[-\s]?4609\b", r"\bCBS\s*101355\b"], "class": "Human", "subtype": "Clinical reference strain"},
    {"name": "Aspergillus fumigatus A1163", "patterns": [r"\bA1163\b", r"\bCEA\s*10\b", r"\bCBS\s*144\.?89\b"], "class": "Human", "subtype": "Clinical reference strain"},
    {"name": "Aspergillus fumigatus CEA17", "patterns": [r"\bCEA[-\s]?17\b"], "class": "Human", "subtype": "Clinical reference strain"},
    # Reference plant-pathogenic strains
    {"name": "Aspergillus flavus NRRL 3357", "patterns": [r"\bNRRL\s*3357\b"], "class": "Plant", "subtype": "Aflatoxigenic reference"},
    {"name": "Aspergillus flavus CA14", "patterns": [r"\bCA[-\s]?14\b", r"\bKuPG"], "class": "Plant", "subtype": "Plant-associated reference"},
    {"name": "Aspergillus flavus AF13", "patterns": [r"\bAF[-\s]?13\b"], "class": "Plant", "subtype": "Aflatoxigenic reference"},
    {"name": "Aspergillus flavus ATCC 22546", "patterns": [r"\bATCC\s*22546\b"], "class": "Plant", "subtype": "Plant-associated reference"},
    # Reference industrial strains
    {"name": "Aspergillus niger N402", "patterns": [r"\baspergillus\s+niger\b.*\bN402\b", r"\bATCC\s*64974\b"], "class": "Industrial", "subtype": "Industrial reference (ATCC 1015 derivative)"},
    # Other reference strains (non-Aspergillus or different section)
    {"name": "Aspergillus nidulans FGSC A4", "patterns": [r"\bFGSC\s*A4\b", r"\baspergillus\s+nidulans\b.*\bA4\b"], "class": "Environmental", "subtype": "Reference strain"},

    # ── Strain-name reclassifications for previously Unknown strains ──
    # These strains have no/insufficient isolation source metadata in NCBI.
    # Classifications assigned based on published literature.

    # Fumigatus: clinical strains with no isolation source metadata
    # Ref: Fedorova et al. (2008) PLoS Genet 4:e1000046 — CEA10 lineage
    # Ref: Rybak et al. (2019) mBio 10:e02213-18 — A1160 is CEA10 pyrG1 derivative
    {"name": "Aspergillus fumigatus A1160", "patterns": [r"\bA1160\b"], "class": "Human", "subtype": "Clinical reference (CEA10 pyrG1 derivative)"},
    # Ref: Valsecchi et al. (2022) Microbiol Spectr — pksP deletion in clinical isolates
    {"name": "Aspergillus fumigatus IP_23", "patterns": [r"\bIP[-_\s]?23\b"], "class": "Human", "subtype": "Clinical pksP variant, France"},
    {"name": "Aspergillus fumigatus IP_24", "patterns": [r"\bIP[-_\s]?24\b"], "class": "Human", "subtype": "Clinical pksP variant, France"},
    # Ref: Zheng et al. (2007) Appl Microbiol Biotechnol 74:233-239 — dye decolorization from rice straw
    {"name": "Aspergillus fumigatus XC6", "patterns": [r"\bXC[-\s]?6\b"], "class": "Environmental", "subtype": "Isolated from rice straw"},

    # Niger: known strains with ambiguous metadata
    # Ref: ATCC product page (https://www.atcc.org/products/10864) — type strain from soil
    # Ref: Braunsdorf et al. (2017) Stand Genomic Sci 12:53 — biofilm-forming cellulase producer
    {"name": "Aspergillus niger ATCC 10864", "patterns": [r"\bATCC\s*10864\b"], "class": "Environmental", "subtype": "Soil isolate, cellulase producer"},
    # Ref: Samson et al. (2007) Stud Mycol 59:1-10 — A. niger neotype, gallic acid production
    {"name": "Aspergillus niger CBS 554.65", "patterns": [r"\bCBS\s*554\.?65\b"], "class": "Industrial", "subtype": "Neotype strain, tannase/gallic acid production"},
    # Ref: Seekles et al. (2022) G3 12:jkac124 — 24 A. niger sensu stricto strains; leather isolate
    {"name": "Aspergillus niger CBS 113.50", "patterns": [r"\bCBS\s*113\.?50\b"], "class": "Environmental", "subtype": "Leather isolate"},
    # Ref: Seekles et al. (2022) G3 12:jkac124 — 24 A. niger sensu stricto strains; raisin isolate
    {"name": "Aspergillus niger CBS 147323", "patterns": [r"\bCBS\s*147323\b"], "class": "Plant", "subtype": "Raisin isolate"},
    # Ref: Takahashi et al. (2023) Microbiol Spectr — taxonomy of Aspergillus series Nigri
    {"name": "Aspergillus niger IFM 59636", "patterns": [r"\bIFM\s*59636\b"], "class": "Human", "subtype": "Abdominal drain isolate"},
    # CCTCC = China Center of Industrial Culture Collection — industrial deposit
    {"name": "Aspergillus niger CCTCC 206047", "patterns": [r"\bCCTCC\s*206047\b"], "class": "Industrial", "subtype": "Industrial culture collection"},

    # Oryzae: fermentation/koji strains with no metadata
    # Ref: Watarai et al. (2019) Microbiol Spectr — soy sauce koji comparison
    # Ref: Zhao et al. (2023) Microbiol Spectr 11:e00836-22 — industrial A. oryzae comparison
    {"name": "Aspergillus oryzae RIB326", "patterns": [r"\bRIB[-\s]?326\b"], "class": "Industrial", "subtype": "Soy sauce koji strain"},
    # Ref: Umemura et al. (2012) Genome Biol Evol 4:1169-1186 — Tokyo Tech koji collection
    # Ref: Machida et al. (2008) Adv Appl Microbiol 64:159-167 — history of koji mold genomics
    {"name": "Aspergillus oryzae TK-13", "patterns": [r"\bTK[-\s]?13\b"], "class": "Industrial", "subtype": "Tokyo Tech fermentation strain"},
    {"name": "Aspergillus oryzae TK-17", "patterns": [r"\bTK[-\s]?17\b"], "class": "Industrial", "subtype": "Tokyo Tech fermentation strain"},
    {"name": "Aspergillus oryzae TK-18", "patterns": [r"\bTK[-\s]?18\b"], "class": "Industrial", "subtype": "Tokyo Tech fermentation strain"},
    {"name": "Aspergillus oryzae TK-21", "patterns": [r"\bTK[-\s]?21\b"], "class": "Industrial", "subtype": "Tokyo Tech fermentation strain"},
    {"name": "Aspergillus oryzae TK-27", "patterns": [r"\bTK[-\s]?27\b"], "class": "Industrial", "subtype": "Tokyo Tech fermentation strain"},
    {"name": "Aspergillus oryzae TK-29", "patterns": [r"\bTK[-\s]?29\b"], "class": "Industrial", "subtype": "Tokyo Tech koji strain"},
    # Originally missubmitted as Neurospora intermedia (corrected by Aalborg University)
    {"name": "Aspergillus oryzae CBS 466.91", "patterns": [r"\bCBS\s*466\.?91\b"], "class": "Industrial", "subtype": "Culture collection koji strain"},
]


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _compile_patterns(pattern_list: List[Tuple[str, str]]) -> List[Tuple[re.Pattern, str]]:
    """Compile regex patterns with their labels."""
    return [(re.compile(p, re.IGNORECASE), label) for p, label in pattern_list]


def _match_patterns(text: str, compiled_patterns: List[Tuple[re.Pattern, str]]) -> List[str]:
    """Return list of matched pattern labels."""
    return [label for pat, label in compiled_patterns if pat.search(text)]


def _score_patterns(text: str, compiled_patterns: List[Tuple[re.Pattern, str]], weight: float) -> Tuple[float, List[str]]:
    """Score text against patterns, return (score, matched_labels)."""
    matched = _match_patterns(text, compiled_patterns)
    return (len(matched) * weight, matched)


# Pre-compile all patterns at module load
_R_CLINICAL_HARD = _compile_patterns(CLINICAL_HARD)
_R_CLINICAL_MEDIUM = _compile_patterns(CLINICAL_MEDIUM)
_R_CLINICAL_SOFT = _compile_patterns(CLINICAL_SOFT)
_R_ANIMAL_CONTEXT = _compile_patterns(ANIMAL_CONTEXT)
_R_ANIMAL_DISEASE = _compile_patterns(ANIMAL_DISEASE)
_R_ANIMAL_NONDISEASE = _compile_patterns(ANIMAL_NONDISEASE)
_R_PLANT_HOSTS = _compile_patterns(PLANT_HOSTS)
_R_PLANT_DISEASE_HARD = _compile_patterns(PLANT_DISEASE_HARD)
_R_PLANT_DISEASE_MEDIUM = _compile_patterns(PLANT_DISEASE_MEDIUM)
_R_PLANT_NONDISEASE = _compile_patterns(PLANT_NONDISEASE)
_R_IND_FERMENT = _compile_patterns(INDUSTRIAL_FERMENT_HARD)
_R_IND_ENZYME = _compile_patterns(INDUSTRIAL_ENZYME_HARD)
_R_IND_CHEM = _compile_patterns(INDUSTRIAL_CHEM_HARD)
_R_IND_BIOCONTROL = _compile_patterns(INDUSTRIAL_BIOCONTROL_HARD)
_R_IND_SOFT = _compile_patterns(INDUSTRIAL_SOFT)
_R_LAB_HARD = _compile_patterns(LAB_REFERENCE_HARD)
_R_LAB_MEDIUM = _compile_patterns(LAB_REFERENCE_MEDIUM)
_R_LAB_SOFT = _compile_patterns(LAB_REFERENCE_SOFT)
_R_ENV_MEDIUM = _compile_patterns(ENVIRONMENTAL_MEDIUM)
_R_ENV_SOFT = _compile_patterns(ENVIRONMENTAL_SOFT)
_R_PATH_NEG = _compile_patterns(PATHOGENIC_NEGATIONS)
_R_MYCOTOX_NEG = _compile_patterns(MYCOTOXIN_NEGATIONS)
_R_TAXONOMY = _compile_patterns(TAXONOMY_CONTEXT)
_R_SURVEY = _compile_patterns(SURVEY_CONTEXT)
_R_MYCOTOXIN_POS = _compile_patterns(MYCOTOXIN_POSITIVE)

# Compile established strain patterns
_ESTABLISHED_COMPILED = []
for _strain_data in ESTABLISHED_STRAINS:
    _compiled_pats = [re.compile(p, re.IGNORECASE) for p in _strain_data["patterns"]]
    _ESTABLISHED_COMPILED.append({
        "name": _strain_data["name"],
        "patterns": _compiled_pats,
        "class": _strain_data["class"],
        "subtype": _strain_data["subtype"]
    })

# Pre-compile TWO-TIER classification patterns
# A) Provenance patterns
_R_PROV_CLINICAL = _compile_patterns(PROV_CLINICAL)
_R_PROV_ENVIRONMENTAL = _compile_patterns(PROV_ENVIRONMENTAL)
_R_PROV_INDUSTRIAL = _compile_patterns(PROV_INDUSTRIAL)
_R_PROV_CULTURE = _compile_patterns(PROV_CULTURE)

# B) Phenotype patterns
_R_PHENO_HUMAN_PATHOGEN = _compile_patterns(PHENO_HUMAN_PATHOGEN)
_R_PHENO_ANIMAL_PATHOGEN = _compile_patterns(PHENO_ANIMAL_PATHOGEN)
_R_PHENO_PLANT_PATHOGEN = _compile_patterns(PHENO_PLANT_PATHOGEN)
_R_PHENO_INDUSTRIAL_TRAIT = _compile_patterns(PHENO_INDUSTRIAL_TRAIT)
_R_PHENO_LAB = _compile_patterns(PHENO_LAB)
_R_PHENO_ENVIRONMENTAL = _compile_patterns(PHENO_ENVIRONMENTAL)
_R_NEGATION = _compile_patterns(NEGATION_PATTERNS)


def _check_established_strain(text: str) -> Optional[Dict]:
    """Check if text matches an established strain."""
    for strain_info in _ESTABLISHED_COMPILED:
        if any(p.search(text) for p in strain_info["patterns"]):
            return strain_info
    return None


# =============================================================================
# MAIN CLASSIFICATION FUNCTIONS
# =============================================================================

def classify_isolation_source(
    fields: Dict[str, Optional[str]],
    debug: bool = False
) -> ClassificationResult:
    """
    Classify isolation source using TWO-TIER classification system.

    A) Provenance (where physically isolated from):
       - Clinical: host-derived (human/animal patient samples)
       - Environmental: soil/air/water/plant debris/built environment
       - Industrial-origin: factory/fermentor/production site
       - Culture-derived: ATCC/NRRL/CBS/lab stocks
       - Unknown: insufficient provenance information

    B) Phenotype (biological association - for pan-GWAS):
       - Human-pathogenic: human infection OR validated virulence
       - Animal-pathogenic: animal infection/disease
       - Plant-pathogenic: mycotoxigenic, crop contamination, plant disease
       - Industrial-trait: citric acid producer, enzyme producer, fermentation
       - Environmental: general environmental isolate
       - Lab: reference strain, lab workhorse
       - Unknown: insufficient phenotype information

    Parameters
    ----------
    fields : Dict[str, Optional[str]]
        Dictionary of field names to their values.

    debug : bool
        If True, print detailed scoring information.

    Returns
    -------
    ClassificationResult
        Contains provenance, phenotype, confidence scores, and matched patterns.
    """
    # Normalize and combine text
    field_texts = {}
    for key, val in fields.items():
        if val is None:
            continue
        if isinstance(val, float) and np.isnan(val):
            continue
        text = str(val).strip()
        if not text or text.lower() in {"missing", "not applicable", "na", "none", "unknown", "nan", ""}:
            continue
        normalized = re.sub(r'[\s\r\n\t]+', ' ', text).strip()
        if normalized:
            field_texts[key] = normalized

    combined = " | ".join(field_texts.values()) if field_texts else ""

    # Empty text -> Unknown for both
    if not combined:
        return ClassificationResult(
            provenance="Unknown",
            provenance_confidence="low",
            provenance_scores={"Unknown": 0},
            provenance_matched={},
            phenotype="Unknown",
            phenotype_confidence="low",
            phenotype_scores={"Unknown": 0},
            phenotype_matched={},
            field_provenance={},
            negations_applied=[]
        )

    # Check negations
    negations = _match_patterns(combined, _R_NEGATION)

    # =========================================================================
    # A) PROVENANCE SCORING (where physically isolated from)
    # =========================================================================
    prov_scores = {
        "Clinical": 0.0,
        "Environmental": 0.0,
        "Industrial-origin": 0.0,
        "Culture-derived": 0.0,
    }
    prov_matched = {k: [] for k in prov_scores}

    # Clinical provenance (host-derived samples)
    clin_score, clin_matched = _score_patterns(combined, _R_PROV_CLINICAL, WEIGHT_HARD)
    prov_scores["Clinical"] = clin_score
    prov_matched["Clinical"] = clin_matched

    # Environmental provenance (soil/air/water/built env)
    env_score, env_matched = _score_patterns(combined, _R_PROV_ENVIRONMENTAL, WEIGHT_MEDIUM)
    prov_scores["Environmental"] = env_score
    prov_matched["Environmental"] = env_matched

    # Industrial-origin provenance (factory/production site)
    ind_score, ind_matched = _score_patterns(combined, _R_PROV_INDUSTRIAL, WEIGHT_HARD)
    prov_scores["Industrial-origin"] = ind_score
    prov_matched["Industrial-origin"] = ind_matched

    # Culture-derived provenance (lab stocks/collections)
    cult_score, cult_matched = _score_patterns(combined, _R_PROV_CULTURE, WEIGHT_MEDIUM)
    prov_scores["Culture-derived"] = cult_score
    prov_matched["Culture-derived"] = cult_matched

    # Determine provenance winner
    prov_sorted = sorted(prov_scores.items(), key=lambda x: x[1], reverse=True)
    prov_top, prov_top_score = prov_sorted[0]
    prov_second_score = prov_sorted[1][1] if len(prov_sorted) > 1 else 0
    prov_margin = prov_top_score - prov_second_score

    if prov_top_score == 0:
        final_prov = "Unknown"
        prov_confidence = "low"
    elif prov_margin >= WEIGHT_HARD * 2:
        final_prov = prov_top
        prov_confidence = "high"
    elif prov_margin >= WEIGHT_HARD or prov_top_score >= WEIGHT_HARD:
        final_prov = prov_top
        prov_confidence = "medium"
    else:
        final_prov = prov_top
        prov_confidence = "low"

    # =========================================================================
    # B) PHENOTYPE SCORING (biological association - for pan-GWAS)
    # =========================================================================
    pheno_scores = {
        "Human-pathogenic": 0.0,
        "Animal-pathogenic": 0.0,
        "Plant-pathogenic": 0.0,
        "Industrial-trait": 0.0,
        "Environmental": 0.0,
    }
    pheno_matched = {k: [] for k in pheno_scores}

    # Human-pathogenic phenotype
    if not negations:
        hp_score, hp_matched = _score_patterns(combined, _R_PHENO_HUMAN_PATHOGEN, WEIGHT_HARD)
        pheno_scores["Human-pathogenic"] = hp_score
        pheno_matched["Human-pathogenic"] = hp_matched
    else:
        # Negations reduce pathogenic scores
        hp_score, hp_matched = _score_patterns(combined, _R_PHENO_HUMAN_PATHOGEN, WEIGHT_HARD)
        pheno_scores["Human-pathogenic"] = max(0, hp_score - len(negations) * WEIGHT_HARD)
        pheno_matched["Human-pathogenic"] = hp_matched

    # Animal-pathogenic phenotype
    if not negations:
        ap_score, ap_matched = _score_patterns(combined, _R_PHENO_ANIMAL_PATHOGEN, WEIGHT_HARD)
        pheno_scores["Animal-pathogenic"] = ap_score
        pheno_matched["Animal-pathogenic"] = ap_matched

    # Plant-pathogenic phenotype (mycotoxigenic, crop contamination)
    if not negations:
        pp_score, pp_matched = _score_patterns(combined, _R_PHENO_PLANT_PATHOGEN, WEIGHT_HARD)
        pheno_scores["Plant-pathogenic"] = pp_score
        pheno_matched["Plant-pathogenic"] = pp_matched
    else:
        # Non-toxigenic negation specifically reduces plant-pathogenic
        pp_score, pp_matched = _score_patterns(combined, _R_PHENO_PLANT_PATHOGEN, WEIGHT_HARD)
        pheno_scores["Plant-pathogenic"] = max(0, pp_score - len(negations) * WEIGHT_HARD)
        pheno_matched["Plant-pathogenic"] = pp_matched

    # Industrial-trait phenotype (production capabilities)
    it_score, it_matched = _score_patterns(combined, _R_PHENO_INDUSTRIAL_TRAIT, WEIGHT_HARD)
    pheno_scores["Industrial-trait"] = it_score
    pheno_matched["Industrial-trait"] = it_matched

    # Lab phenotype removed — reference/model strains are classified by their
    # original isolation phenotype (Human, Plant, Industrial, etc.)
    # Lab scoring patterns still run but contribute to provenance (Culture-derived),
    # not to phenotype.

    # Environmental phenotype (default, no specific association)
    env_pheno_score, env_pheno_matched = _score_patterns(combined, _R_PHENO_ENVIRONMENTAL, WEIGHT_MEDIUM)
    # Non-pathogenic/non-toxigenic strains boost Environmental
    if negations:
        env_pheno_score += len(negations) * WEIGHT_MEDIUM
        env_pheno_matched.extend([f"Negation: {n}" for n in negations])
    pheno_scores["Environmental"] = env_pheno_score
    pheno_matched["Environmental"] = env_pheno_matched

    # Determine phenotype winner
    pheno_sorted = sorted(pheno_scores.items(), key=lambda x: x[1], reverse=True)
    pheno_top, pheno_top_score = pheno_sorted[0]
    pheno_second_score = pheno_sorted[1][1] if len(pheno_sorted) > 1 else 0
    pheno_margin = pheno_top_score - pheno_second_score

    if pheno_top_score == 0:
        final_pheno = "Unknown"
        pheno_confidence = "low"
    elif pheno_margin >= WEIGHT_HARD * 2:
        final_pheno = pheno_top
        pheno_confidence = "high"
    elif pheno_margin >= WEIGHT_HARD or pheno_top_score >= WEIGHT_HARD:
        final_pheno = pheno_top
        pheno_confidence = "medium"
    else:
        final_pheno = pheno_top
        pheno_confidence = "low"

    # =========================================================================
    # ENVIRONMENTAL PROVENANCE GUARD
    # =========================================================================
    # When the isolation source is clearly environmental (soil, air, water, etc.)
    # and provenance was classified as Environmental, the phenotype should not be
    # overridden by ambiguous keywords from project descriptions (e.g., "antifungal
    # resistance" in a project title doesn't make a soil isolate Human-pathogenic,
    # "non-toxigenic" doesn't make a soil isolate Industrial-trait, and "ochratoxin"
    # in a taxonomy study doesn't make a fermentation isolate Plant-pathogenic).
    if final_prov == "Environmental":
        iso_text = field_texts.get("isolation_source", "").lower()
        is_clearly_environmental = any(w in iso_text for w in [
            "soil", "compost", "air ", "water", "sediment", "leaf litter",
            "dust", "cave", "marine", "mangrove", "permafrost",
        ])

        if is_clearly_environmental and final_pheno in ("Human-pathogenic", "Industrial-trait"):
            # Only override if the phenotype signal came from project-level keywords
            # (not from direct clinical/industrial evidence in the isolation source
            # or strong evidence elsewhere in the metadata)
            has_direct_evidence = False
            if final_pheno == "Human-pathogenic":
                has_direct_evidence = any(w in iso_text for w in [
                    "patient", "sputum", "lung", "blood", "clinical", "bronch",
                    "infection", "aspergillosis", "cornea", "wound",
                ])
            elif final_pheno == "Industrial-trait":
                # Check isolation source AND full combined text for industrial evidence
                # (e.g., "citric acid producing strain" in WGS title)
                all_text = combined.lower()
                has_direct_evidence = any(w in iso_text for w in [
                    "ferment", "koji", "sake", "miso", "soy", "industrial",
                    "enzyme", "citric", "daqu", "nuruk",
                ]) or any(w in all_text for w in [
                    "citric acid produc", "enzyme produc", "industrial produc",
                    "fermentation strain",
                ])

            if not has_direct_evidence:
                final_pheno = "Environmental"
                pheno_confidence = "medium"
                pheno_scores["Environmental"] = WEIGHT_OVERRIDE
                pheno_matched["Environmental"].append(
                    f"Environmental provenance guard: isolation source is '{iso_text[:50]}', "
                    f"overriding ambiguous {final_pheno} signal from project metadata"
                )

    # =========================================================================
    # MYCOTOXIN IN TAXONOMY/SURVEY STUDIES GUARD
    # =========================================================================
    # When the phenotype is Plant-pathogenic solely from mycotoxin keywords
    # (ochratoxin, fumonisin) but the study is a taxonomy/survey and the
    # isolation source suggests fermentation/industrial, correct to Industrial.
    if final_pheno == "Plant-pathogenic":
        iso_text = field_texts.get("isolation_source", "").lower()
        submitter_text = field_texts.get("submitter", "").lower()
        is_industrial_source = any(w in iso_text for w in [
            "liquor", "daqu", "qu ", "ferment", "vinegar", "baijiu",
        ])
        is_industrial_submitter = "industrial" in submitter_text

        if is_industrial_source or is_industrial_submitter:
            # Check if the Plant-pathogenic score came only from mycotoxin keywords
            pp_matched = pheno_matched.get("Plant-pathogenic", [])
            mycotoxin_only = all(any(m in match for m in [
                "ochratoxin", "fumonisin", "mycotoxin", "aflatoxin"
            ]) for match in pp_matched) if pp_matched else False

            if mycotoxin_only:
                final_pheno = "Industrial-trait"
                pheno_confidence = "high"
                pheno_scores["Industrial-trait"] = WEIGHT_OVERRIDE
                pheno_matched["Industrial-trait"].append(
                    f"Industrial source guard: isolation='{iso_text[:50]}', "
                    f"mycotoxin keywords from taxonomy study, not crop contamination"
                )

    # =========================================================================
    # INDUSTRIAL OVERRIDE: Industrial provenance overrides all other phenotypes
    # =========================================================================
    # Industrial strains (fermentation, enzyme/chemical production) are selected
    # to be non-toxigenic. If a strain is used industrially, it means no mycotoxin
    # production - this overrides any Plant-pathogenic signals from project
    # descriptions mentioning mycotoxins (which may be for comparison/taxonomy).
    if final_prov == "Industrial-origin":
        final_pheno = "Industrial-trait"
        pheno_confidence = "high"
        pheno_scores["Industrial-trait"] = WEIGHT_OVERRIDE
        pheno_matched["Industrial-trait"] = ["Industrial provenance override (non-toxigenic industrial strain)"]

    # =========================================================================
    # ESTABLISHED STRAIN OVERRIDES
    # =========================================================================
    strain_match = _check_established_strain(combined)
    if strain_match:
        # Map established strain class to phenotype
        class_to_pheno = {
            "Human": "Human-pathogenic",
            "Animal": "Animal-pathogenic",
            "Plant": "Plant-pathogenic",
            "Industrial": "Industrial-trait",
            "Environmental": "Environmental",
        }
        # Established strains get high confidence overrides
        established_pheno = class_to_pheno.get(strain_match["class"], "Unknown")
        final_pheno = established_pheno
        pheno_confidence = "high"
        pheno_scores[established_pheno] = WEIGHT_OVERRIDE
        pheno_matched[established_pheno] = [f"Established strain: {strain_match['name']}"]

    # =========================================================================
    # FALLBACK: Use provenance to infer phenotype when Unknown
    # =========================================================================
    if final_pheno == "Unknown" and final_prov != "Unknown":
        # If we know WHERE it came from but not WHAT association, infer from provenance
        if final_prov == "Clinical":
            # Check what specific markers matched for Clinical provenance
            clinical_markers = prov_matched.get("Clinical", [])
            clinical_markers_str = str(clinical_markers).lower()

            # Also check original combined text for culture medium exclusion
            combined_lower = combined.lower()

            # Pathogenic context markers (actual disease/infection context)
            pathogenic_markers = ["patient", "hospital", "icu", "bal", "sputum",
                                 "lung", "respiratory", "wound", "abscess", "pus",
                                 "infection", "aspergillosis", "clinical isolate",
                                 "clinical sample", "biopsy", "swab",
                                 "cornea", "eye", "sinus", "nasal", "ear", "oral",
                                 "diseased", "infected", "veterinary"]

            # "blood" is only pathogenic if NOT preceded by "sheep" (culture medium)
            has_real_blood = "blood" in clinical_markers_str and "sheep" not in combined_lower

            # Commensal/survey markers (NOT pathogenic - gut fungi from surveys)
            commensal_markers = ["fec", "stool", "microbiome", "survey"]

            # Culture medium markers (should NOT count as pathogenic)
            culture_medium_markers = ["sheep blood", "agar", "pda", "medium", "broth"]

            has_pathogenic_context = any(m in clinical_markers_str for m in pathogenic_markers) or has_real_blood
            has_commensal = any(m in clinical_markers_str for m in commensal_markers)
            has_culture_medium = any(m in combined_lower for m in culture_medium_markers)

            # If fecal sample with culture medium, it's commensal regardless of medium type
            if has_commensal and has_culture_medium:
                has_pathogenic_context = False  # Culture medium doesn't make it pathogenic

            if has_pathogenic_context and not has_commensal:
                # Actual pathogenic context -> Human-pathogenic
                final_pheno = "Human-pathogenic"
                pheno_confidence = "low"
                pheno_matched["Human-pathogenic"] = ["Inferred from Clinical provenance (pathogenic context)"]
            elif has_commensal:
                # Fecal/microbiome samples are commensal, NOT pathogenic
                final_pheno = "Environmental"
                pheno_confidence = "low"
                pheno_matched["Environmental"] = ["Inferred from Clinical provenance (commensal sample)"]
            elif has_pathogenic_context and has_commensal:
                # Both present but no culture medium - prefer pathogenic context
                final_pheno = "Human-pathogenic"
                pheno_confidence = "low"
                pheno_matched["Human-pathogenic"] = ["Inferred from Clinical provenance (mixed context)"]
            else:
                # Just "human" marker without specific context -> Unknown, let other fallbacks handle
                final_pheno = "Unknown"
                pheno_confidence = "low"
        elif final_prov == "Environmental":
            # Environmental samples without specific association -> Environmental phenotype
            final_pheno = "Environmental"
            pheno_confidence = "low"
            pheno_matched["Environmental"] = [f"Inferred from Environmental provenance"]
        elif final_prov == "Industrial-origin":
            # Industrial site samples -> Industrial-trait
            final_pheno = "Industrial-trait"
            pheno_confidence = "low"
            pheno_matched["Industrial-trait"] = [f"Inferred from Industrial-origin provenance"]
        elif final_prov == "Culture-derived":
            # Culture collection without specific phenotype -> Unknown
            # (Lab category removed; culture-derived strains should be classified
            # by their original isolation source, not their culture status)
            final_pheno = "Unknown"
            pheno_confidence = "low"
            pheno_matched["Unknown"] = [f"Culture-derived with no specific phenotype information"]

    # Build field provenance tracking
    field_prov_tracking = {}
    if final_prov != "Unknown":
        field_prov_tracking[f"Provenance:{final_prov}"] = list(field_texts.keys())
    if final_pheno != "Unknown":
        field_prov_tracking[f"Phenotype:{final_pheno}"] = list(field_texts.keys())

    if debug:
        print(f"\n=== Two-Tier Classification Debug ===")
        print(f"Combined text: {combined[:200]}...")
        print(f"\nPROVENANCE scores: {prov_scores}")
        print(f"PROVENANCE matched: {prov_matched}")
        print(f"Final PROVENANCE: {final_prov} (confidence: {prov_confidence})")
        print(f"\nPHENOTYPE scores: {pheno_scores}")
        print(f"PHENOTYPE matched: {pheno_matched}")
        print(f"Final PHENOTYPE: {final_pheno} (confidence: {pheno_confidence})")

    return ClassificationResult(
        provenance=final_prov,
        provenance_confidence=prov_confidence,
        provenance_scores=prov_scores,
        provenance_matched={k: v for k, v in prov_matched.items() if v},
        phenotype=final_pheno,
        phenotype_confidence=pheno_confidence,
        phenotype_scores=pheno_scores,
        phenotype_matched={k: v for k, v in pheno_matched.items() if v},
        field_provenance=field_prov_tracking,
        negations_applied=negations
    )


def classify_strain_from_row(
    row: pd.Series,
    field_mapping: Optional[Dict[str, List[str]]] = None,
    debug: bool = False
) -> ClassificationResult:
    """
    Classify a DataFrame row.

    Parameters
    ----------
    row : pd.Series
        A row from a metadata DataFrame.
    field_mapping : Dict[str, List[str]], optional
        Mapping from expected field names to list of actual column names to try.
    debug : bool
        If True, print debug information.

    Returns
    -------
    ClassificationResult
    """
    default_mapping = {
        "isolation_source": ["isolation_source", "Isolation Source", "isolationSource"],
        "host": ["host", "Host"],
        "env_biome": ["env_broad_scale", "env_biome", "env_local_scale", "env_medium"],
        "geo_loc_name": ["geo_loc_name", "geoLocName", "geographic_location", "Geo Location"],
        "description": ["BioSample Description Title", "biosample_title", "title", "description"],
        "comment": ["Comment", "comment", "Comments"],
        "strain": ["Organism Infraspecific Names Strain", "Strain", "strain", "Organism Infraspecific Names Isolate", "isolate", "Isolate"],
        "submitter": ["Assembly Submitter", "submitter"],
        "bioproject_title": ["BioProject Title", "bioproject_title"],
        "bioproject_description": ["BioProject Description", "bioproject_description"],
        "annotation_text": ["annotation_text", "CombinedText"],
        "species": ["Species", "species", "Organism Name"],
    }

    if field_mapping:
        for key, candidates in default_mapping.items():
            if key not in field_mapping:
                field_mapping[key] = candidates
    else:
        field_mapping = default_mapping

    fields = {}
    for field_name, candidates in field_mapping.items():
        if isinstance(candidates, str):
            candidates = [candidates]
        for col_name in candidates:
            if col_name in row.index:
                val = row[col_name]
                if pd.notna(val) and str(val).strip():
                    fields[field_name] = str(val)
                    break

    return classify_isolation_source(fields, debug=debug)


def classify_dataframe(
    df: pd.DataFrame,
    field_mapping: Optional[Dict[str, List[str]]] = None,
    add_legacy_columns: bool = True,
    debug: bool = False
) -> pd.DataFrame:
    """
    Classify all rows in a DataFrame using TWO-TIER classification.

    Adds both Provenance (where isolated from) and Phenotype (biological
    association for pan-GWAS) columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with strain metadata.
    field_mapping : Dict[str, List[str]], optional
        Mapping from expected field names to column names.
    add_legacy_columns : bool
        If True, also add IsolationClass/UsageClass columns for
        backward compatibility (maps Phenotype to IsolationClass).
    debug : bool
        If True, print debug for first row.

    Returns
    -------
    pd.DataFrame
        Copy of input DataFrame with classification columns added:
        - Provenance, ProvenanceConfidence
        - Phenotype, PhenotypeConfidence
        - Plus legacy columns if requested
    """
    result_df = df.copy()
    results = [
        classify_strain_from_row(row, field_mapping, debug=(debug and i == 0))
        for i, (_, row) in enumerate(df.iterrows())
    ]

    # =========================================================================
    # A) PROVENANCE COLUMNS (where isolated from)
    # =========================================================================
    result_df["Provenance"] = [r.provenance for r in results]
    result_df["ProvenanceConfidence"] = [r.provenance_confidence for r in results]
    result_df["ProvenanceScores"] = [json.dumps(r.provenance_scores) for r in results]
    result_df["ProvenanceMatched"] = [json.dumps(r.provenance_matched) for r in results]

    # =========================================================================
    # B) PHENOTYPE COLUMNS (biological association - for pan-GWAS)
    # =========================================================================
    result_df["Phenotype"] = [r.phenotype for r in results]
    result_df["PhenotypeConfidence"] = [r.phenotype_confidence for r in results]
    result_df["PhenotypeScores"] = [json.dumps(r.phenotype_scores) for r in results]
    result_df["PhenotypeMatched"] = [json.dumps(r.phenotype_matched) for r in results]

    # Metadata columns
    result_df["FieldProvenance"] = [json.dumps(r.field_provenance) for r in results]
    result_df["NegationsApplied"] = [json.dumps(r.negations_applied) for r in results]

    # =========================================================================
    # SPECIES-AWARE POST-PROCESSING OVERRIDES
    # =========================================================================
    # These rules apply biological knowledge about specific species

    # Check if Species column exists
    species_col = None
    for col in ['Species', 'species', 'Organism Name']:
        if col in result_df.columns:
            species_col = col
            break

    if species_col:
        # Get isolation source column
        iso_col = None
        for col in ['Isolation Source', 'isolation_source', 'IsolationSource']:
            if col in result_df.columns:
                iso_col = col
                break

        if iso_col:
            # ANY species + "fermented" in isolation source → Industrial-trait
            # Fermented foods are industrial products (sake, miso, soy sauce, etc.)
            fermented_mask = (
                result_df[iso_col].str.contains(r'ferment', case=False, na=False) &
                (result_df["Phenotype"] == "Environmental")
            )
            if fermented_mask.any():
                result_df.loc[fermented_mask, "Phenotype"] = "Industrial-trait"
                result_df.loc[fermented_mask, "PhenotypeConfidence"] = "medium"

            # A. oryzae + "food"/"fermented food" → Industrial-trait
            # A. oryzae is a domesticated fermentation species; food context = fermentation
            oryzae_food_mask = (
                result_df[species_col].str.contains('oryzae', case=False, na=False) &
                result_df[iso_col].str.contains(r'food', case=False, na=False, regex=True) &
                (result_df["Phenotype"] == "Environmental")
            )
            if oryzae_food_mask.any():
                result_df.loc[oryzae_food_mask, "Phenotype"] = "Industrial-trait"
                result_df.loc[oryzae_food_mask, "PhenotypeConfidence"] = "medium"

            # A. oryzae soil in Industrial-trait → Environmental (soil is not industrial)
            # BUT only if there's no explicit industrial context in the matched patterns OR description
            oryzae_soil_industrial = (
                result_df[species_col].str.contains('oryzae', case=False, na=False) &
                result_df[iso_col].str.contains(r'^soil$', case=False, na=False, regex=True) &
                (result_df["Phenotype"] == "Industrial-trait")
            )
            if oryzae_soil_industrial.any():
                for idx in result_df[oryzae_soil_industrial].index:
                    # Check multiple sources for explicit industrial context
                    pheno_matched_str = str(result_df.loc[idx, "PhenotypeMatched"]).lower()
                    prov_matched_str = str(result_df.loc[idx, "ProvenanceMatched"]).lower()
                    bp_desc = str(result_df.loc[idx, 'BioProject Description']).lower() if 'BioProject Description' in result_df.columns else ''

                    industrial_terms = [
                        "citric acid", "enzyme", "ferment", "producer", "production",
                        "industrial", "biotechnology", "bioprocess", "sake", "miso",
                        "soy sauce", "shoyu", "koji"
                    ]
                    has_explicit_industrial = (
                        any(term in pheno_matched_str for term in industrial_terms) or
                        any(term in prov_matched_str for term in industrial_terms) or
                        any(term in bp_desc for term in industrial_terms)
                    )
                    if not has_explicit_industrial:
                        result_df.loc[idx, "Phenotype"] = "Environmental"
                        result_df.loc[idx, "PhenotypeConfidence"] = "medium"

            # A. niger soil in Industrial-trait → Environmental
            # BUT only if there's no explicit industrial context (e.g., citric acid producer)
            niger_soil_industrial = (
                result_df[species_col].str.contains('niger', case=False, na=False) &
                result_df[iso_col].str.contains(r'^soil$', case=False, na=False, regex=True) &
                (result_df["Phenotype"] == "Industrial-trait")
            )
            if niger_soil_industrial.any():
                for idx in result_df[niger_soil_industrial].index:
                    # Check multiple sources for explicit industrial context
                    pheno_matched_str = str(result_df.loc[idx, "PhenotypeMatched"]).lower()
                    prov_matched_str = str(result_df.loc[idx, "ProvenanceMatched"]).lower()
                    bp_desc = str(result_df.loc[idx, 'BioProject Description']).lower() if 'BioProject Description' in result_df.columns else ''

                    industrial_terms = [
                        "citric acid", "enzyme", "ferment", "producer", "production",
                        "industrial", "biotechnology", "bioprocess"
                    ]
                    has_explicit_industrial = (
                        any(term in pheno_matched_str for term in industrial_terms) or
                        any(term in prov_matched_str for term in industrial_terms) or
                        any(term in bp_desc for term in industrial_terms)
                    )
                    if not has_explicit_industrial:
                        result_df.loc[idx, "Phenotype"] = "Environmental"
                        result_df.loc[idx, "PhenotypeConfidence"] = "medium"

            # Isolation source = "laboratory" alone needs context check
            # Check BioProject description for actual context
            lab_only_mask = (
                result_df[iso_col].str.match(r'^lab(?:oratory)?$', case=False, na=False) &
                (result_df["Phenotype"].isin(["Unknown"]))
            )
            if lab_only_mask.any() and 'BioProject Description' in result_df.columns:
                for idx in result_df[lab_only_mask].index:
                    bp_desc = str(result_df.loc[idx, 'BioProject Description']).lower()
                    if 'industrial' in bp_desc:
                        result_df.loc[idx, "Phenotype"] = "Industrial-trait"
                        result_df.loc[idx, "PhenotypeConfidence"] = "medium"
                    elif 'marine' in bp_desc or 'sponge' in bp_desc:
                        result_df.loc[idx, "Phenotype"] = "Environmental"
                        result_df.loc[idx, "PhenotypeConfidence"] = "medium"
                    else:
                        result_df.loc[idx, "Phenotype"] = "Unknown"
                        result_df.loc[idx, "PhenotypeConfidence"] = "low"

    # =========================================================================
    # BIOPROJECT DESCRIPTION-BASED FIXES FOR UNKNOWN SAMPLES
    # =========================================================================
    if 'BioProject Description' in result_df.columns:
        unknown_mask = result_df["Phenotype"] == "Unknown"
        if unknown_mask.any():
            for idx in result_df[unknown_mask].index:
                bp_desc = str(result_df.loc[idx, 'BioProject Description']).lower()
                species = str(result_df.loc[idx, species_col]).lower() if species_col else ''

                # A. oryzae + domestication context → Industrial-trait (fermentation evolution)
                if 'oryzae' in species and 'domestication' in bp_desc:
                    result_df.loc[idx, "Phenotype"] = "Industrial-trait"
                    result_df.loc[idx, "PhenotypeConfidence"] = "medium"
                # Lignocellulosic degradation studies → Environmental
                elif 'lignocellulosic' in bp_desc or 'degradation potential' in bp_desc:
                    result_df.loc[idx, "Phenotype"] = "Environmental"
                    result_df.loc[idx, "PhenotypeConfidence"] = "low"

    # =========================================================================
    # SUBMITTER-BASED OVERRIDES
    # =========================================================================
    # FDA submitter → Human-pathogenic (diagnostic reference genomes for pathogen identification)
    submitter_col = None
    for col in ['Assembly Submitter', 'submitter']:
        if col in result_df.columns:
            submitter_col = col
            break

    if submitter_col:
        fda_mask = (
            result_df[submitter_col].str.contains('FDA|Food and Drug Administration',
                                                   case=False, na=False, regex=True) &
            (result_df["Phenotype"] != "Human-pathogenic")
        )
        if fda_mask.any():
            result_df.loc[fda_mask, "Phenotype"] = "Human-pathogenic"
            result_df.loc[fda_mask, "PhenotypeConfidence"] = "medium"

    # =========================================================================
    # LEGACY COLUMNS (for backward compatibility)
    # =========================================================================
    if add_legacy_columns:
        # Map Phenotype to old IsolationClass naming
        pheno_to_iso = {
            "Human-pathogenic": "Human",
            "Animal-pathogenic": "Animal",
            "Plant-pathogenic": "Plant",
            "Industrial-trait": "Industrial",
            "Environmental": "Environmental",
            "Unknown": "Unknown",
        }
        result_df["IsolationClass"] = result_df["Phenotype"].map(pheno_to_iso)
        result_df["IsolationSubclass"] = result_df["Phenotype"]  # Use phenotype as subclass
        result_df["IsolationSubtype"] = result_df["Phenotype"]

        # Use phenotype confidence for classification confidence
        result_df["ClassConfidence"] = result_df["PhenotypeConfidence"]
        result_df["ClassificationConfidence"] = result_df["PhenotypeConfidence"]

        # Legacy UsageClass mapping
        usage_mapping = {
            "Human-pathogenic": "Pathogenic",
            "Animal-pathogenic": "Pathogenic",
            "Plant-pathogenic": "Pathogenic",
            "Industrial-trait": "Industrial",
            "Environmental": "Wild-type",
            "Unknown": "Unknown",
        }
        result_df["UsageClass"] = result_df["Phenotype"].map(usage_mapping)

    return result_df


def summarize_classification(df: pd.DataFrame, class_col: str = "Phenotype") -> pd.DataFrame:
    """
    Summarize TWO-TIER classification results.

    Parameters
    ----------
    df : pd.DataFrame
        Classified DataFrame with Provenance and Phenotype columns.
    class_col : str
        Column name for primary classification (default: "Phenotype").

    Returns
    -------
    pd.DataFrame
        Summary statistics.
    """
    print("\n" + "=" * 70)
    print("TWO-TIER ISOLATION SOURCE CLASSIFICATION SUMMARY")
    print("=" * 70)
    print(f"\nTotal samples: {len(df)}")

    # =========================================================================
    # A) PROVENANCE SUMMARY (where isolated from)
    # =========================================================================
    if "Provenance" in df.columns:
        prov_counts = df["Provenance"].value_counts()
        print("\n" + "-" * 70)
        print("A) PROVENANCE (where physically isolated from)")
        print("-" * 70)
        for prov, count in prov_counts.items():
            print(f"  {prov}: {count} ({100 * count / len(df):.1f}%)")

        if "ProvenanceConfidence" in df.columns:
            print(f"\n  Confidence by provenance:")
            conf_table = df.groupby(["Provenance", "ProvenanceConfidence"]).size().unstack(fill_value=0)
            print("  " + conf_table.to_string().replace("\n", "\n  "))

    # =========================================================================
    # B) PHENOTYPE SUMMARY (biological association - for pan-GWAS)
    # =========================================================================
    if "Phenotype" in df.columns:
        pheno_counts = df["Phenotype"].value_counts()
        print("\n" + "-" * 70)
        print("B) PHENOTYPE (biological association - for pan-GWAS)")
        print("-" * 70)
        for pheno, count in pheno_counts.items():
            print(f"  {pheno}: {count} ({100 * count / len(df):.1f}%)")

        if "PhenotypeConfidence" in df.columns:
            print(f"\n  Confidence by phenotype:")
            conf_table = df.groupby(["Phenotype", "PhenotypeConfidence"]).size().unstack(fill_value=0)
            print("  " + conf_table.to_string().replace("\n", "\n  "))

    # =========================================================================
    # CROSS-TABULATION: Provenance x Phenotype
    # =========================================================================
    if "Provenance" in df.columns and "Phenotype" in df.columns:
        print("\n" + "-" * 70)
        print("CROSS-TABULATION: Provenance x Phenotype")
        print("-" * 70)
        cross_tab = pd.crosstab(df["Provenance"], df["Phenotype"], margins=True)
        print(cross_tab.to_string())

    # Legacy support
    if class_col not in df.columns and "Phenotype" in df.columns:
        class_col = "Phenotype"
    if class_col not in df.columns:
        raise ValueError(f"Column '{class_col}' not found")

    class_counts = df[class_col].value_counts()
    return class_counts.to_frame("count")


def audit_classification(
    df: pd.DataFrame,
    sample_n: int = 5,
    focus_phenotype: Optional[str] = None,
    focus_provenance: Optional[str] = None,
    focus_confidence: Optional[str] = None
) -> None:
    """
    Print detailed audit information for sample classifications.

    Parameters
    ----------
    df : pd.DataFrame
        Classified DataFrame.
    sample_n : int
        Number of samples to audit.
    focus_phenotype : str, optional
        If specified, only audit this phenotype.
    focus_provenance : str, optional
        If specified, only audit this provenance.
    focus_confidence : str, optional
        If specified, only audit this confidence level.
    """
    subset = df.copy()

    if focus_phenotype and "Phenotype" in subset.columns:
        subset = subset[subset["Phenotype"] == focus_phenotype]

    if focus_provenance and "Provenance" in subset.columns:
        subset = subset[subset["Provenance"] == focus_provenance]

    if focus_confidence:
        if "PhenotypeConfidence" in subset.columns:
            subset = subset[subset["PhenotypeConfidence"] == focus_confidence]
        elif "ClassificationConfidence" in subset.columns:
            subset = subset[subset["ClassificationConfidence"] == focus_confidence]

    if len(subset) == 0:
        print("No samples match the audit criteria.")
        return

    sample = subset.sample(min(sample_n, len(subset)))

    print("\n" + "=" * 80)
    print("TWO-TIER CLASSIFICATION AUDIT")
    print("=" * 80)

    for idx, row in sample.iterrows():
        print(f"\n--- {row.get('Assembly Accession', idx)} ---")
        print(f"Species: {row.get('Species', 'N/A')}")

        # Two-tier classification
        prov = row.get('Provenance', 'N/A')
        prov_conf = row.get('ProvenanceConfidence', 'N/A')
        pheno = row.get('Phenotype', 'N/A')
        pheno_conf = row.get('PhenotypeConfidence', 'N/A')

        print(f"PROVENANCE: {prov} (confidence: {prov_conf})")
        print(f"PHENOTYPE:  {pheno} (confidence: {pheno_conf})")

        for fld in ["Isolation Source", "Strain", "Comment", "BioProject Title"]:
            if fld in row and pd.notna(row[fld]):
                val = str(row[fld])[:80]
                print(f"{fld}: {val}{'...' if len(str(row[fld])) > 80 else ''}")

        # Show matched patterns for both tiers
        if "ProvenanceMatched" in row and pd.notna(row["ProvenanceMatched"]):
            try:
                prov_patterns = json.loads(row["ProvenanceMatched"])
                if prov_patterns:
                    print(f"Provenance matched: {prov_patterns}")
            except Exception:
                pass

        if "PhenotypeMatched" in row and pd.notna(row["PhenotypeMatched"]):
            try:
                pheno_patterns = json.loads(row["PhenotypeMatched"])
                if pheno_patterns:
                    print(f"Phenotype matched: {pheno_patterns}")
            except Exception:
                pass


# =============================================================================
# GEOGRAPHIC VISUALIZATION FUNCTIONS  
# =============================================================================

# Try importing optional geographic packages
try:
    from geopy.geocoders import Nominatim
    from geopy.extra.rate_limiter import RateLimiter
    HAS_GEOPY = True
except ImportError:
    HAS_GEOPY = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False

# Default paths for geographic data (stored in NB0_Results)
RESULTS_BASE = '/datadrive/Analysis/NB0_Results'
GEOCODE_CACHE_PATH = os.path.join(RESULTS_BASE, 'geocode_cache.json')
MANUAL_LATLON_PATH = os.path.join(RESULTS_BASE, 'manual_latlon_additions.txt')


def load_geocode_cache(cache_path: str = None) -> dict:
    """Load geocoding cache from JSON file."""
    path = cache_path or GEOCODE_CACHE_PATH
    if os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}


def save_geocode_cache(cache: dict, cache_path: str = None):
    """Save geocoding cache to JSON file."""
    path = cache_path or GEOCODE_CACHE_PATH
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(cache, f)


def normalize_location_key(loc: str) -> str:
    """Normalize location string for matching (lowercase, collapse whitespace)."""
    if not loc:
        return ""
    s = str(loc).lower().strip()
    s = re.sub(r'\s+', ' ', s)  # Collapse whitespace
    s = re.sub(r'\s*:\s*', ':', s)  # Remove spaces around colons
    s = re.sub(r'\s*,\s*', ', ', s)  # Normalize comma spacing
    return s


def load_manual_latlon(latlon_path: str = None) -> dict:
    """Load manually curated lat/lon coordinates from txt file.
    
    Returns a dict mapping normalized location strings to (lat, lon) tuples.
    """
    path = latlon_path or MANUAL_LATLON_PATH
    manual_coords = {}
    if not os.path.exists(path):
        return manual_coords
    
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
        
        current_location = None
        current_lat = None
        current_lon = None
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if line.startswith('Location:'):
                current_location = line.replace('Location:', '').strip()
            elif line.startswith('Latitude:'):
                try:
                    current_lat = float(line.replace('Latitude:', '').strip())
                except ValueError:
                    current_lat = None
            elif line.startswith('Longitude:'):
                try:
                    current_lon = float(line.replace('Longitude:', '').strip())
                except ValueError:
                    current_lon = None
            elif line.startswith('Notes:'):
                if current_location and current_lat is not None and current_lon is not None:
                    key = normalize_location_key(current_location)
                    manual_coords[key] = (current_lat, current_lon)
                current_location = None
                current_lat = None
                current_lon = None
        
        if current_location and current_lat is not None and current_lon is not None:
            key = normalize_location_key(current_location)
            manual_coords[key] = (current_lat, current_lon)
        
        if manual_coords:
            print(f"Loaded {len(manual_coords)} manual lat/lon entries from {path}")
    except Exception as e:
        print(f"Warning: Could not load manual lat/lon file: {e}")
    
    return manual_coords


def parse_lat_lon(lat_lon_str: str) -> Tuple[Optional[float], Optional[float]]:
    """Parse Lat/Lon string like '35.6762 N 139.6503 E' to (lat, lon)."""
    if pd.isna(lat_lon_str) or not lat_lon_str:
        return None, None
    if str(lat_lon_str).lower() in ['not collected', 'nan', 'missing']:
        return None, None
    try:
        s = str(lat_lon_str).strip()
        parts = s.replace(',', ' ').split()
        if len(parts) >= 4:
            lat = float(parts[0])
            if parts[1].upper() == 'S':
                lat = -lat
            lon = float(parts[2])
            if parts[3].upper() == 'W':
                lon = -lon
            return lat, lon
        elif len(parts) == 2:
            return float(parts[0]), float(parts[1])
    except:
        pass
    return None, None


def geocode_location(location_str: str, geolocator, cache: dict) -> Tuple[Optional[float], Optional[float]]:
    """Geocode a location string to (lat, lon) using cache."""
    if pd.isna(location_str) or not location_str:
        return None, None
    if str(location_str).lower() in ['not collected', 'missing', 'unknown', 'nan']:
        return None, None
    
    loc_key = str(location_str).strip().lower()
    if loc_key in cache:
        return tuple(cache[loc_key]) if cache[loc_key][0] is not None else (None, None)
    
    try:
        loc_clean = location_str.replace(':', ',').strip()
        location = geolocator.geocode(loc_clean, timeout=10)
        if location:
            result = [location.latitude, location.longitude]
            cache[loc_key] = result
            return tuple(result)
    except Exception:
        pass
    
    cache[loc_key] = [None, None]
    return None, None


def extract_coordinates(
    df: pd.DataFrame,
    lat_lon_col: str = 'Lat/Lon',
    geo_loc_col: str = 'Geo Location',
    species_col: str = 'Species',
    accession_col: str = 'Assembly Accession',
    isolation_col: str = 'Isolation Source',
    strain_col: str = 'Organism Infraspecific Names Strain',
    use_geocoding: bool = True,
    cache_path: str = None,
    manual_latlon_path: str = None
) -> Tuple[pd.DataFrame, dict]:
    """Extract geographic coordinates from a DataFrame.
    
    Tries multiple sources in order:
    1. Lat/Lon column (direct coordinates)
    2. Geo Location (geocoded via API/cache)
    3. Manual additions (from txt file)
    
    Returns:
        Tuple of (coords_df, stats_dict)
    """
    # Initialize geocoder
    geolocator = None
    if use_geocoding and HAS_GEOPY:
        geolocator = Nominatim(user_agent="aspergillus_strain_analysis_v1")
    
    geo_cache = load_geocode_cache(cache_path)
    initial_cache_size = len(geo_cache)
    manual_latlon = load_manual_latlon(manual_latlon_path)
    
    # Check column existence
    lat_lon_col = lat_lon_col if lat_lon_col in df.columns else None
    geo_loc_col = geo_loc_col if geo_loc_col in df.columns else None
    
    coords_data = []
    n_from_latlon = 0
    n_from_geocode = 0
    n_from_cache = 0
    n_from_manual = 0
    
    for idx, row in df.iterrows():
        lat, lon = None, None
        source = None
        
        # Try Lat/Lon first
        if lat_lon_col and pd.notna(row.get(lat_lon_col)):
            lat, lon = parse_lat_lon(row[lat_lon_col])
            if lat is not None:
                source = 'Lat/Lon'
                n_from_latlon += 1
        
        # Fall back to Geo Location (geocoding)
        if lat is None and geo_loc_col and pd.notna(row.get(geo_loc_col)) and geolocator:
            loc_key = str(row[geo_loc_col]).strip().lower()
            was_cached = loc_key in geo_cache
            lat, lon = geocode_location(row[geo_loc_col], geolocator, geo_cache)
            if lat is not None:
                source = 'Geo Location'
                if was_cached:
                    n_from_cache += 1
                else:
                    n_from_geocode += 1
        
        # Fall back to manual lat/lon entries
        if lat is None and manual_latlon:
            if geo_loc_col and pd.notna(row.get(geo_loc_col)):
                loc_key = normalize_location_key(row[geo_loc_col])
                if loc_key in manual_latlon:
                    lat, lon = manual_latlon[loc_key]
                    source = 'Manual'
                    n_from_manual += 1
            
            # Try matching strain name
            if lat is None and strain_col in df.columns:
                strain = str(row.get(strain_col, '')).strip()
                for manual_key, (m_lat, m_lon) in manual_latlon.items():
                    if strain and strain.lower() in manual_key.lower():
                        lat, lon = m_lat, m_lon
                        source = 'Manual (strain)'
                        n_from_manual += 1
                        break
        
        if lat is not None and lon is not None:
            coords_data.append({
                'lat': lat,
                'lon': lon,
                'species': row.get(species_col, 'Unknown'),
                'accession': row.get(accession_col, ''),
                'isolation_source': row.get(isolation_col, ''),
                'source': source
            })
    
    # Save geocode cache if it grew
    if len(geo_cache) > initial_cache_size:
        save_geocode_cache(geo_cache, cache_path)
    
    df_coords = pd.DataFrame(coords_data)
    
    stats = {
        'total': len(df),
        'mapped': len(df_coords),
        'from_latlon': n_from_latlon,
        'from_geocode': n_from_geocode,
        'from_cache': n_from_cache,
        'from_manual': n_from_manual
    }
    
    return df_coords, stats


def plot_strain_world_map(
    df_coords: pd.DataFrame,
    species_colors: dict = None,
    output_path: str = None,
    figsize: Tuple[int, int] = (16, 10),
    title: str = 'Aspergillus Strain Isolation Locations'
):
    """Plot strain locations on a world map."""
    import matplotlib.pyplot as plt
    
    if species_colors is None:
        species_colors = {
            'A. fumigatus': '#e74c3c',
            'A. flavus': '#ffd500',
            'A. niger': '#1a1aff',
            'A. oryzae': '#27ae60',
        }
    
    if HAS_CARTOPY:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.Robinson()})
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='wheat', edgecolor='none')
        ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':')
        
        for species in df_coords['species'].unique():
            subset = df_coords[df_coords['species'] == species]
            color = species_colors.get(species, '#888888')
            ax.scatter(
                subset['lon'], subset['lat'],
                c=color, s=50, alpha=0.7, label=species,
                transform=ccrs.PlateCarree(), edgecolor='white', linewidth=0.5
            )
        
        ax.legend(loc='lower left', fontsize=10)
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        for species in df_coords['species'].unique():
            subset = df_coords[df_coords['species'] == species]
            color = species_colors.get(species, '#888888')
            ax.scatter(subset['lon'], subset['lat'], c=color, s=50, alpha=0.7, label=species)
        ax.legend()
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title(title)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Map saved to: {output_path}")
    


# =============================================================================
# SAMPLE ID NORMALIZATION FUNCTIONS
# =============================================================================

def normalize_sample_id(sample_id: str) -> str:
    """Normalize sample ID by removing version suffix."""
    if pd.isna(sample_id):
        return ""
    s = str(sample_id).strip()
    # Remove version suffix like .1, .2
    if re.match(r'^GC[AF]_\d+\.\d+$', s):
        return s.rsplit('.', 1)[0]
    return s


def create_sample_id_mappings(df: pd.DataFrame, id_col: str = 'Assembly Accession') -> dict:
    """Create mappings between full and normalized sample IDs."""
    mappings = {}
    for idx, row in df.iterrows():
        full_id = row.get(id_col, '')
        norm_id = normalize_sample_id(full_id)
        if norm_id:
            mappings[norm_id] = full_id
            mappings[full_id] = full_id
    return mappings


def create_cross_sample_mappings(
    pheno_sample_ids: List[str],
    pca_sample_ids: List[str],
    feature_sample_ids: List[str]
) -> dict:
    """Create cross-mappings between phenotype, PCA, and feature matrix sample IDs.

    Handles format mismatches like:
    - GCA_000002855.2 (metadata format)
    - GCF_000002855.4_ASM285v2 (PCA/kinship format)
    """
    # Normalize all IDs
    pheno_norm = {normalize_sample_id(s): s for s in pheno_sample_ids}
    pca_norm = {normalize_sample_id(s): s for s in pca_sample_ids}
    feature_norm = {normalize_sample_id(s): s for s in feature_sample_ids}

    # Create mappings
    pheno_to_pca = {}
    pheno_to_feature = {}

    for norm_id, pheno_id in pheno_norm.items():
        if norm_id in pca_norm:
            pheno_to_pca[pheno_id] = pca_norm[norm_id]
        if norm_id in feature_norm:
            pheno_to_feature[pheno_id] = feature_norm[norm_id]

    # Find common samples (in all three)
    common_norm = set(pheno_norm.keys()) & set(pca_norm.keys()) & set(feature_norm.keys())
    common_samples = [pheno_norm[n] for n in common_norm]

    return {
        'pheno_to_pca': pheno_to_pca,
        'pheno_to_feature': pheno_to_feature,
        'common_samples': common_samples,
        'pheno_norm': pheno_norm,
        'pca_norm': pca_norm,
        'feature_norm': feature_norm
    }


# =============================================================================
# PHENOTYPE FILE CREATION FUNCTIONS
# =============================================================================

def create_phenotype_files(samples: pd.DataFrame, contrasts: dict, output_dir: str):
    """
    Create binary phenotype files for each contrast.
    
    strict: high-confidence only
    broad: high + medium confidence
    lenient: all confidence levels (high + medium + low)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    phenotype_files = {}
    
    for contrast_name, contrast_config in contrasts.items():
        case_labels = contrast_config['case']
        control_labels = contrast_config['control']
        
        for stringency in ['strict', 'broad', 'lenient']:
            if stringency == 'strict':
                valid_conf = ['high']
            elif stringency == 'broad':
                valid_conf = ['high', 'medium']
            else:
                valid_conf = ['high', 'medium', 'low']
            
            # Filter samples
            case_mask = (samples['label_short'].isin(case_labels)) & (samples['confidence'].isin(valid_conf))
            ctrl_mask = (samples['label_short'].isin(control_labels)) & (samples['confidence'].isin(valid_conf))
            
            case_samples = samples[case_mask]['sample_id'].tolist()
            ctrl_samples = samples[ctrl_mask]['sample_id'].tolist()
            
            # Create phenotype DataFrame
            pheno_data = []
            for s in case_samples:
                pheno_data.append({'sample_id': s, 'phenotype': 1})
            for s in ctrl_samples:
                pheno_data.append({'sample_id': s, 'phenotype': 0})
            
            if len(pheno_data) > 0:
                pheno_df = pd.DataFrame(pheno_data)
                filename = f"{contrast_name}_{stringency}.tsv"
                filepath = os.path.join(output_dir, filename)
                pheno_df.to_csv(filepath, sep='\t', index=False)
                phenotype_files[f"{contrast_name}_{stringency}"] = {
                    'path': filepath,
                    'n_case': len(case_samples),
                    'n_control': len(ctrl_samples)
                }
    
    return phenotype_files


# =============================================================================
# ASSOCIATION TESTING FUNCTIONS
# =============================================================================

def run_association_test(
    feature_matrix: pd.DataFrame,
    phenotype: pd.Series,
    kinship: pd.DataFrame = None,
    pcs: pd.DataFrame = None,
    n_pcs: int = 3,
    min_maf: float = 0.05
) -> pd.DataFrame:
    """
    Run association test for each feature against phenotype.
    
    Uses logistic regression with optional PC covariates.
    Mixed model with kinship can be approximated by including PCs.
    """
    from scipy import stats
    import warnings
    warnings.filterwarnings('ignore')
    
    # Align samples
    common_samples = list(set(feature_matrix.index) & set(phenotype.index))
    if pcs is not None:
        common_samples = list(set(common_samples) & set(pcs.index))
    
    X = feature_matrix.loc[common_samples]
    y = phenotype.loc[common_samples]
    
    if pcs is not None:
        pc_cols = [c for c in pcs.columns if c.startswith('PC')][:n_pcs]
        covariates = pcs.loc[common_samples, pc_cols]
    else:
        covariates = None
    
    results = []
    
    for feature in X.columns:
        feat_values = X[feature].values
        
        # Calculate MAF
        maf = np.mean(feat_values)
        if maf > 0.5:
            maf = 1 - maf
        
        if maf < min_maf:
            continue
        
        try:
            # Simple chi-square test for binary features
            if len(np.unique(feat_values)) == 2:
                contingency = pd.crosstab(feat_values, y.values)
                if contingency.shape == (2, 2):
                    chi2, pval, dof, expected = stats.chi2_contingency(contingency)
                else:
                    pval = 1.0
            else:
                # Mann-Whitney U for continuous
                group0 = feat_values[y == 0]
                group1 = feat_values[y == 1]
                if len(group0) > 0 and len(group1) > 0:
                    stat, pval = stats.mannwhitneyu(group0, group1, alternative='two-sided')
                else:
                    pval = 1.0
            
            # Effect size (case vs control mean difference)
            case_mean = np.mean(feat_values[y == 1])
            ctrl_mean = np.mean(feat_values[y == 0])
            effect = case_mean - ctrl_mean
            
            results.append({
                'feature': feature,
                'pvalue': pval,
                'effect': effect,
                'case_mean': case_mean,
                'ctrl_mean': ctrl_mean,
                'maf': maf,
                'n_samples': len(common_samples)
            })
        except Exception:
            continue
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        _, results_df['pvalue_adj'], _, _ = multipletests(
            results_df['pvalue'], method='fdr_bh'
        )
    
    return results_df


# =============================================================================
# DIAGNOSTIC PLOTTING FUNCTIONS
# =============================================================================

def calculate_genomic_inflation(pvalues: np.ndarray) -> float:
    """Calculate genomic inflation factor (lambda)."""
    from scipy import stats
    pvalues = np.array(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues > 0]
    
    chi2_observed = stats.chi2.ppf(1 - pvalues, df=1)
    lambda_gc = np.median(chi2_observed) / stats.chi2.ppf(0.5, df=1)
    return lambda_gc


def apply_genomic_control(results_df, lambda_gc=None):
    """Apply genomic control correction to pan-GWAS results.

    Divides chi-squared statistics by lambda_gc to correct for
    population structure inflation, then recomputes p-values and FDR.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    df = results_df.copy()
    pvals = df['pvalue'].values.astype(float)

    if lambda_gc is None:
        lambda_gc = calculate_genomic_inflation(pvals)

    if lambda_gc <= 1.0:
        df['pvalue_gc'] = pvals
    else:
        valid = (pvals > 0) & (pvals < 1) & ~np.isnan(pvals)
        pvals_gc = np.full_like(pvals, np.nan)
        chi2 = stats.chi2.ppf(1 - pvals[valid], df=1)
        chi2_corrected = chi2 / lambda_gc
        pvals_gc[valid] = 1 - stats.chi2.cdf(chi2_corrected, df=1)
        df['pvalue_gc'] = pvals_gc

    mask = ~np.isnan(df['pvalue_gc'])
    df['fdr_gc'] = np.nan
    if mask.sum() > 0:
        df.loc[mask, 'fdr_gc'] = multipletests(
            df.loc[mask, 'pvalue_gc'], method='fdr_bh'
        )[1]

    return df


def plot_qq(pvalues: np.ndarray, ax=None, title: str = 'QQ Plot'):
    """Plot QQ plot of p-values."""
    import matplotlib.pyplot as plt
    from scipy import stats
    
    pvalues = np.array(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues > 0]
    
    observed = -np.log10(np.sort(pvalues))
    expected = -np.log10(np.linspace(1/len(pvalues), 1, len(pvalues)))
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    
    ax.scatter(expected, observed, alpha=0.5, s=10)
    max_val = max(max(expected), max(observed))
    ax.plot([0, max_val], [0, max_val], 'r--', label='Expected')
    
    lambda_gc = calculate_genomic_inflation(pvalues)
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')
    ax.set_title(f'{title}\nλ = {lambda_gc:.3f}')
    ax.legend()
    
    return ax


def plot_manhattan(
    results: pd.DataFrame,
    pval_col: str = 'pvalue',
    feature_col: str = 'feature',
    ax=None,
    title: str = 'Manhattan Plot',
    significance_threshold: float = 0.05
):
    """Plot Manhattan-style plot of association results."""
    import matplotlib.pyplot as plt
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    
    results = results.sort_values(feature_col)
    results['idx'] = range(len(results))
    
    log_p = -np.log10(results[pval_col].values)
    
    ax.scatter(results['idx'], log_p, alpha=0.5, s=10)
    
    # Significance line
    if significance_threshold:
        ax.axhline(-np.log10(significance_threshold), color='red', linestyle='--', 
                   label=f'p = {significance_threshold}')
        
        # Bonferroni line
        bonf = significance_threshold / len(results)
        ax.axhline(-np.log10(bonf), color='blue', linestyle=':', 
                   label=f'Bonferroni ({bonf:.2e})')
    
    ax.set_xlabel('Feature Index')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(title)
    ax.legend()
    
    return ax


def create_diagnostic_plots(
    results: pd.DataFrame,
    output_path: str = None,
    title_prefix: str = ''
):
    """Create QQ and Manhattan diagnostic plots."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    plot_qq(results['pvalue'].values, ax=axes[0], 
            title=f'{title_prefix} QQ Plot')
    plot_manhattan(results, ax=axes[1], 
                   title=f'{title_prefix} Manhattan Plot')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Diagnostic plots saved to: {output_path}")
    


# --- Functions extracted from NB0_DataPrep notebook ---

# =============================================================================
# Strain world map
# =============================================================================

try:
    from geopy.geocoders import Nominatim
    HAS_GEOPY = True
except ImportError:
    print('WARNING: geopy not installed. Geocoding disabled.')
    HAS_GEOPY = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    print('Note: cartopy not installed. Using simple scatter plot.')
    HAS_CARTOPY = False

# ---- Geocoding helpers ----
def load_geocode_cache():
    # Try NB0 cache first, use NB0 cache
    nb0_cache = os.path.join(RESULTS_BASE, 'geocode_cache.json')
    for path in [nb0_cache, GEOCODE_CACHE_PATH]:
        if os.path.exists(path):
            with open(path, 'r') as f:
                return json.load(f)
    return {}

def save_geocode_cache(cache):
    out = os.path.join(RESULTS_BASE, 'geocode_cache.json')
    with open(out, 'w') as f:
        json.dump(cache, f)

def load_manual_latlon():
    manual_coords = {'by_location': {}, 'by_strain': {}}
    if not os.path.exists(MANUAL_LATLON_PATH):
        return manual_coords
    current_location = current_lat = current_lon = current_notes = None
    with open(MANUAL_LATLON_PATH, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('Location:'):
            if current_location and current_lat is not None and current_lon is not None:
                strain_match = re.search(r'\(([A-Z]+\d+)\)', current_location)
                if strain_match:
                    manual_coords['by_strain'][strain_match.group(1).upper()] = (current_lat, current_lon, current_notes)
                else:
                    manual_coords['by_location'][current_location.lower().strip()] = (current_lat, current_lon, current_notes)
            current_location = line.replace('Location:', '').strip()
            current_lat = current_lon = current_notes = None
        elif line.startswith('Latitude:'):
            try: current_lat = float(line.replace('Latitude:', '').strip())
            except ValueError: pass
        elif line.startswith('Longitude:'):
            try: current_lon = float(line.replace('Longitude:', '').strip())
            except ValueError: pass
        elif line.startswith('Notes:'):
            current_notes = line.replace('Notes:', '').strip()
    if current_location and current_lat is not None and current_lon is not None:
        strain_match = re.search(r'\(([A-Z]+\d+)\)', current_location)
        if strain_match:
            manual_coords['by_strain'][strain_match.group(1).upper()] = (current_lat, current_lon, current_notes)
        else:
            manual_coords['by_location'][current_location.lower().strip()] = (current_lat, current_lon, current_notes)
    print(f'Loaded manual coordinates: {len(manual_coords["by_location"])} by location, {len(manual_coords["by_strain"])} by strain')
    return manual_coords

def parse_lat_lon(lat_lon_str):
    if pd.isna(lat_lon_str) or not lat_lon_str or str(lat_lon_str).lower() in ['not collected', 'nan', 'missing']:
        return None, None
    try:
        s = str(lat_lon_str).strip()
        parts = s.replace(',', ' ').split()
        if len(parts) >= 4:
            lat = float(parts[0])
            if parts[1].upper() == 'S': lat = -lat
            lon = float(parts[2])
            if parts[3].upper() == 'W': lon = -lon
            return lat, lon
        elif len(parts) == 2:
            return float(parts[0]), float(parts[1])
    except:
        pass
    return None, None

def lookup_manual_coords(row, manual_coords):
    strain_cols = [
        'Organism Infraspecific Names Strain', 'Strain', 'Infraspecific Name',
        'Assembly Accession', 'sample_id', 'BioSample Description Title', 'Assembly Name',
    ]
    for col in strain_cols:
        if col in row.index and pd.notna(row.get(col)):
            val = str(row[col]).upper()
            for pattern in [r'(NRRL\d+)', r'(CBS[\d\.]+)', r'(ATCC\d+)', r'(IFO\d+)', r'(FGSC\d+)']:
                match = re.search(pattern, val)
                if match:
                    sid = match.group(1).replace('.', '')
                    if sid in manual_coords['by_strain']:
                        return manual_coords['by_strain'][sid][:2]
    geo_loc = row.get('Geo Location', '')
    if pd.notna(geo_loc) and geo_loc:
        loc_key = str(geo_loc).lower().strip()
        if loc_key in manual_coords['by_location']:
            return manual_coords['by_location'][loc_key][:2]
        for mk, coords in manual_coords['by_location'].items():
            if mk in loc_key or loc_key in mk:
                return coords[:2]
    return None, None

def geocode_location(location_str, geolocator, cache):
    if pd.isna(location_str) or not location_str:
        return None, None
    if str(location_str).lower() in ['not collected', 'missing', 'unknown', 'nan']:
        return None, None
    loc_key = str(location_str).strip().lower()
    if loc_key in cache:
        return tuple(cache[loc_key]) if cache[loc_key][0] is not None else (None, None)
    try:
        loc_clean = location_str.replace(':', ',').strip()
        location = geolocator.geocode(loc_clean, timeout=10)
        if location:
            result = [location.latitude, location.longitude]
            cache[loc_key] = result
            return tuple(result)
    except:
        pass
    cache[loc_key] = [None, None]
    return None, None

def build_strain_coordinates(df_analysis, species_colors=None):
    """Extract geographic coordinates for all strains in df_analysis.

    Returns
    -------
    df_coords : pd.DataFrame
        DataFrame with lat, lon, species, accession, isolation_source, source columns.
    """
    if species_colors is None:
        from funpan_utils import SPECIES_COLORS
        species_colors = SPECIES_COLORS

    print('Extracting coordinates for strain locations...')
    manual_coords = load_manual_latlon()
    geolocator = Nominatim(user_agent='aspergillus_strain_analysis_v1') if HAS_GEOPY else None
    geo_cache = load_geocode_cache()
    initial_cache_size = len(geo_cache)

    coords_data = []
    lat_lon_col = 'Lat/Lon' if 'Lat/Lon' in df_analysis.columns else None
    geo_loc_col = 'Geo Location' if 'Geo Location' in df_analysis.columns else None
    n_from_latlon = n_from_manual = n_from_geocode = n_from_cache = 0

    for idx, row in df_analysis.iterrows():
        lat, lon = None, None
        source = None

        if lat_lon_col and pd.notna(row.get(lat_lon_col)):
            lat, lon = parse_lat_lon(row[lat_lon_col])
            if lat is not None:
                source = 'Lat/Lon'; n_from_latlon += 1

        if lat is None:
            lat, lon = lookup_manual_coords(row, manual_coords)
            if lat is not None:
                source = 'Manual'; n_from_manual += 1

        if lat is None and geo_loc_col and pd.notna(row.get(geo_loc_col)) and geolocator:
            loc_key = str(row[geo_loc_col]).strip().lower()
            was_cached = loc_key in geo_cache
            lat, lon = geocode_location(row[geo_loc_col], geolocator, geo_cache)
            if lat is not None:
                source = 'Geo Location'
                if was_cached: n_from_cache += 1
                else: n_from_geocode += 1

        if lat is not None and lon is not None:
            coords_data.append({
                'lat': lat, 'lon': lon,
                'species': row['Species'],
                'accession': row['Assembly Accession'],
                'isolation_source': row.get('Isolation Source', ''),
                'source': source,
            })

    if len(geo_cache) > initial_cache_size:
        save_geocode_cache(geo_cache)
        print(f'Geocode cache updated: {len(geo_cache)} locations')

    df_coords = pd.DataFrame(coords_data)
    print(f'\nCoordinates found for {len(df_coords)}/{len(df_analysis)} samples '
          f'({100 * len(df_coords) / len(df_analysis):.1f}%)')
    print(f'  From Lat/Lon column: {n_from_latlon}')
    print(f'  From Manual additions: {n_from_manual}')
    print(f'  From Geo Location (geocoded): {n_from_geocode}')
    print(f'  From Geo Location (cached): {n_from_cache}')

    return df_coords


def plot_strain_world_map_from_coords(df_coords, results_base, species_colors=None,
                                       fontsize=15, dpi=400,
                                       marker_size=90, figsize=(18, 11)):
    """Plot a world map of strain locations and save to results_base.

    Parameters
    ----------
    df_coords : pd.DataFrame
        Output of ``build_strain_coordinates``.
    results_base : str
        Directory for saving the plot.
    species_colors : dict, optional
        Species-to-color mapping.
    fontsize : int
        Base font size (title is +3, legend title is +1).
    dpi : int
        Resolution of saved PNG.
    marker_size : int
        Scatter marker area in points^2 (matplotlib ``s`` argument). Default 90.
        Try 150-250 for poster figures.
    figsize : tuple
        Figure size in inches.
    """
    import matplotlib.pyplot as plt
    if species_colors is None:
        from funpan_utils import SPECIES_COLORS_DISPLAY
        species_colors = SPECIES_COLORS_DISPLAY

    if len(df_coords) == 0:
        print('\nNo coordinates available for mapping.')
        return

    color_lookup = dict(species_colors)

    FS = fontsize
    if HAS_CARTOPY:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
        ax.add_feature(cfeature.OCEAN, facecolor='#e6f3ff')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':', alpha=0.5)
        for species in sorted(df_coords['species'].unique()):
            sp_data = df_coords[df_coords['species'] == species]
            color = color_lookup.get(species, '#333333')
            ax.scatter(sp_data['lon'], sp_data['lat'], c=color, s=marker_size, alpha=0.75,
                       transform=ccrs.PlateCarree(),
                       label=f'{species} (n={len(sp_data)})',
                       edgecolors='white', linewidth=0.7, zorder=5)
        ax.legend(loc='lower left', fontsize=FS, title='Species',
                  title_fontsize=FS + 1, framealpha=0.95)
        ax.set_title(f'Geographic Distribution of Aspergillus Strains (n={len(df_coords)})',
                     fontsize=FS + 3, fontweight='bold')
    else:
        fig, ax = plt.subplots(figsize=figsize)
        for species in sorted(df_coords['species'].unique()):
            sp_data = df_coords[df_coords['species'] == species]
            color = color_lookup.get(species, '#333333')
            ax.scatter(sp_data['lon'], sp_data['lat'], c=color, s=marker_size, alpha=0.75,
                       label=f'{species} (n={len(sp_data)})',
                       edgecolors='white', linewidth=0.7)
        ax.set_xlabel('Longitude', fontsize=FS); ax.set_ylabel('Latitude', fontsize=FS)
        ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
        ax.tick_params(axis='both', labelsize=FS)
        ax.legend(loc='lower left', fontsize=FS, title='Species',
                  title_fontsize=FS + 1, framealpha=0.95)
        ax.set_title(f'Geographic Distribution of Aspergillus Strains (n={len(df_coords)})',
                     fontsize=FS + 3, fontweight='bold')
        ax.grid(True, alpha=0.3); ax.set_facecolor('#f5f5f5')

    plt.tight_layout()
    map_path = os.path.join(results_base, 'strain_world_map.png')
    fig.savefig(map_path, dpi=dpi, bbox_inches='tight')
    print(f'\nMap saved to: {map_path}')


# =============================================================================
# Per-species phenotype file generation
# =============================================================================

LABEL_MAP = {
    'Human-pathogenic': 'Human',
    'Animal-pathogenic': 'Animal',
    'Plant-pathogenic': 'Plant',
    'Industrial-trait': 'Industrial',
    'Environmental': 'Environmental',
    'Lab': 'Lab',
    'Unknown': 'Unknown',
}

MIN_SAMPLES = 3  # Minimum samples per group for GWAS


def create_phenotype_files(samples, contrasts, output_dir):
    """Create binary phenotype files for each contrast (strict/broad/lenient)."""
    phenotype_files = {}
    for case_label, control_label, contrast_name in contrasts:
        for mode in ['strict', 'broad', 'lenient']:
            if mode == 'strict':
                conf_mask = samples['confidence'] == 'high'
            elif mode == 'broad':
                conf_mask = samples['confidence'].isin(['high', 'medium'])
            else:
                conf_mask = pd.Series([True] * len(samples), index=samples.index)

            label_mask = samples['label_short'].isin([case_label, control_label])
            subset = samples[conf_mask & label_mask].copy()
            subset['phenotype'] = (subset['label_short'] == case_label).astype(int)
            pheno_df = subset[['sample_id', 'phenotype']].copy()

            filename = f'pheno_{contrast_name}_{mode}.tsv'
            filepath = os.path.join(output_dir, filename)
            pheno_df.to_csv(filepath, sep='\t', index=False)

            n_case = (pheno_df['phenotype'] == 1).sum()
            n_control = (pheno_df['phenotype'] == 0).sum()
            phenotype_files[f'{contrast_name}_{mode}'] = {
                'path': filepath, 'n_case': n_case, 'n_control': n_control,
            }
            print(f'  {filename}: {case_label}={n_case}, {control_label}={n_control}')
    return phenotype_files


def generate_per_species_phenotype_files(df_classified, species_list, species_display,
                                          results_base, label_map=None,
                                          min_samples=None):
    """Generate per-species phenotype master sheets and binary contrast files.

    Parameters
    ----------
    df_classified : pd.DataFrame
        Classified metadata (output of ``classify_dataframe``).
    species_list : list
        List of species keys (e.g. ``['fumigatus', 'flavus', ...]``).
    species_display : dict
        ``{species_key: display_name}`` mapping.
    results_base : str
        Base output directory.
    label_map : dict, optional
        Phenotype label shortening map. Defaults to ``LABEL_MAP``.
    min_samples : int, optional
        Minimum samples per group for GWAS contrasts. Defaults to ``MIN_SAMPLES``.
    """
    if label_map is None:
        label_map = LABEL_MAP
    if min_samples is None:
        min_samples = MIN_SAMPLES

    for species in species_list:
        sp_display = species_display[species]
        df_sp = df_classified[df_classified['Species'] == sp_display].copy()
        if len(df_sp) == 0:
            print(f'\n{sp_display}: no samples, skipping')
            continue

        sp_results = os.path.join(results_base, species)
        os.makedirs(sp_results, exist_ok=True)

        print(f'\n{"=" * 60}')
        print(f'{sp_display} ({len(df_sp)} samples)')
        print('=' * 60)

        # Master sample sheet
        class_col = 'Phenotype'
        conf_col = 'PhenotypeConfidence'
        cols_to_use = ['Assembly Accession', class_col, 'Provenance', conf_col,
                       'Isolation Source', 'Geo Location', 'Assembly BioProject Accession']
        cols_to_use = [c for c in cols_to_use if c in df_sp.columns]

        samples_df = df_sp[cols_to_use].copy()
        rename_map = {
            'Assembly Accession': 'sample_id',
            class_col: 'label',
            'Provenance': 'subclass',
            conf_col: 'confidence',
            'Isolation Source': 'isolation_text',
            'Geo Location': 'country',
            'Assembly BioProject Accession': 'bioproject',
        }
        samples_df = samples_df.rename(columns={k: v for k, v in rename_map.items() if k in samples_df.columns})
        samples_df['label_short'] = samples_df['label'].map(label_map).fillna(samples_df['label'])

        samples_path = os.path.join(sp_results, f'{species}_samples.tsv')
        samples_df.to_csv(samples_path, sep='\t', index=False)
        print(f'  Saved master sample sheet: {samples_path}')

        phenotype_df = samples_df.copy()
        phenotype_df = phenotype_df.rename(columns={'label': 'isolation_class', 'isolation_text': 'isolation_source'})
        phenotype_path = os.path.join(sp_results, 'phenotype_data.tsv')
        phenotype_df.to_csv(phenotype_path, sep='\t', index=False)
        print(f'  Saved phenotype data: {phenotype_path}')

        print(f'\n  Label distribution:')
        print(samples_df['label_short'].value_counts().to_string())

        # Binary contrast files
        gwas_samples = samples_df[~samples_df['label_short'].isin(['Lab', 'Unknown'])].copy()
        label_counts = gwas_samples['label_short'].value_counts()

        contrasts = []
        if label_counts.get('Human', 0) >= min_samples:
            if label_counts.get('Environmental', 0) >= min_samples:
                contrasts.append(('Human', 'Environmental', 'human_vs_env'))
            if label_counts.get('Industrial', 0) >= min_samples:
                contrasts.append(('Human', 'Industrial', 'human_vs_ind'))
        if label_counts.get('Industrial', 0) >= min_samples and label_counts.get('Environmental', 0) >= min_samples:
            contrasts.append(('Industrial', 'Environmental', 'ind_vs_env'))
        if label_counts.get('Plant', 0) >= min_samples:
            if label_counts.get('Environmental', 0) >= min_samples:
                contrasts.append(('Plant', 'Environmental', 'plant_vs_env'))
            if label_counts.get('Human', 0) >= min_samples:
                contrasts.append(('Human', 'Plant', 'human_vs_plant'))
        if label_counts.get('Animal', 0) >= min_samples:
            if label_counts.get('Environmental', 0) >= min_samples:
                contrasts.append(('Animal', 'Environmental', 'animal_vs_env'))

        if contrasts:
            pheno_dir = os.path.join(sp_results, 'phenotypes')
            os.makedirs(pheno_dir, exist_ok=True)
            print(f'\n  Creating phenotype contrast files...')
            print(f'    strict = high confidence only')
            print(f'    broad  = high + medium confidence')
            print(f'    lenient = all confidence levels')
            pf = create_phenotype_files(gwas_samples, contrasts, pheno_dir)
            print(f'  Created {len(pf)} phenotype files in {pheno_dir}')
        else:
            print(f'  No valid contrasts (need >= {min_samples} samples per group)')

