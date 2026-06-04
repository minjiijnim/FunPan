#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_utils.py
===============
Shared utilities for FunPan analysis notebooks.
Contains: constants, metadata loaders, annotation parsers, QC visualization.

Extracted for FunPan analysis.
"""

from __future__ import annotations
import os
import re
import sys
import glob
import json
import random
import itertools
from pathlib import Path
from typing import Dict, List, Iterable, Optional
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

from scipy.signal import find_peaks
from scipy.stats import fisher_exact
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# ---- Biopython for GenBank parsing ----
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None
    print("[WARN] biopython not available - BiG-SCAPE functions may not work")


# =============================================================================
# CONSTANTS
# =============================================================================

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")

# EggNOG annotation fields to keep
EGGNOG_KEEP = [
    "Description_EggNog", "COG_category", "GOs_EggNog", "EC_EggNog", "KEGG_ko", "KEGG_Pathway",
    "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy_EggNog",
    "BiGG_Reaction", "PFAMs_EggNog"
]

# dbCAN fields (overview.tsv only)
DBCAN_KEEP = [
    "CAZy_dbCAN", "dbcan_Recommend", "EC_dbCAN", "dbcan_NumTools", "dbcan_Substrate",
]

# InterProScan fields
INTERPRO_KEEP = [
    "interpro_IPR", "GOs_Interpro", "interpro_Pathways", "PFAMs_Interpro",
    "EC_Interpro", "Description_Interpro",
]

# SignalP fields
SIGNALP_KEEP = [
    "signalp_Prediction", "signalp_SP_Probability", "signalp_Cleavage_Site",
]

_CAZY_FAM_RE = re.compile(r"\b(?:GH|GT|PL|CE|CBM|AA)\d+(?:_\d+)?\b")
_MULTI_SEP_RE = re.compile(r"[;,|]")

# COG Category Descriptions (NCBI reference)
COG_DESCRIPTIONS = {
    # INFORMATION STORAGE AND PROCESSING
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    # CELLULAR PROCESSES AND SIGNALING
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    # METABOLISM
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    # POORLY CHARACTERIZED
    "R": "General function prediction only",
    "S": "Function unknown",
    "X": "Mobilome: prophages, transposons",
}


# =============================================================================
# UNINFORMATIVE-TERM FILTER FOR FUNCTIONAL ENRICHMENT
# =============================================================================
# Standard practice in fungal pangenome / functional-enrichment work is to
# exclude annotation terms that are not interpretable functional categories
# (e.g. COG-S "Function unknown", Pfam DUF families, GO root nodes). Including
# them inflates the multiple-testing burden and can produce trivial "findings"
# that say nothing about biology -- e.g. a rare compartment "enriched" for
# "Function unknown" simply means it has more uncharacterised proteins, which
# is tautological for a strain-private compartment.

# COG categories explicitly flagged as uninformative.
UNINFORMATIVE_COG_TERMS = {
    "S",   # Function unknown
    "R",   # General function prediction only (deprecated in modern COG schemes)
}

# GO root terms -- always present, never informative.
UNINFORMATIVE_GO_TERMS = {
    "GO:0003674",  # molecular_function (root)
    "GO:0008150",  # biological_process (root)
    "GO:0005575",  # cellular_component (root)
}

# EC placeholder term used by some annotators for "no enzyme classification".
UNINFORMATIVE_EC_TERMS = {"-.-.-.-", "-", ""}

# Description substrings (lowercased) that mark a term as uninformative.
# Applied to InterPro / KEGG_ko / Pfam descriptions when available.
UNINFORMATIVE_DESC_PATTERNS = (
    "uncharacterised",
    "uncharacterized",
    "hypothetical protein",
    "protein of unknown function",
    "domain of unknown function",
)


def is_informative_term(layer, term, description=None):
    """Return False for annotation terms that should be excluded from
    functional enrichment / convergence tests.

    Filters applied per layer:

    - ``COG_category``: drop ``S`` (Function unknown) and ``R`` (General
      function prediction only).
    - ``PFAMs``: drop any term name beginning with ``DUF`` (case-insensitive),
      and any term whose description matches an uninformative pattern.
    - ``interpro_IPR``: drop terms whose description matches an uninformative
      pattern (uncharacterised / DUF / hypothetical).
    - ``KEGG_ko``: drop KOs whose description matches an uninformative pattern.
    - ``GOs``: drop the three root terms (molecular_function,
      biological_process, cellular_component).
    - ``EC``: drop placeholder ``-.-.-.-`` and empty terms.

    All other layers (``CAZy``, ``KEGG_Pathway``, ``KEGG_TC``,
    ``dbcan_Substrate``, ``KEGG_Module``, ``KEGG_Reaction``, ``BRITE``,
    ``BiGG_Reaction``) pass through; they are already curated functional
    catalogues without "function unknown" sentinels.

    Parameters
    ----------
    layer : str
        Annotation layer name as used in the OG consensus column headers.
    term : str
        Term identifier (e.g. ``"S"``, ``"PF11111"``, ``"GO:0003824"``).
    description : str, optional
        Human-readable description; only consulted when the term ID alone is
        not enough to judge (e.g. KEGG_ko / interpro_IPR / Pfam descriptions).

    Returns
    -------
    bool
        True if the term should be tested, False if it should be skipped.
    """
    if term is None:
        return False
    try:
        if pd.isna(term):
            return False
    except (TypeError, ValueError):
        pass
    term_s = str(term).strip()
    if term_s == "" or term_s.lower() == "nan":
        return False

    desc_s = ""
    if description is not None:
        try:
            if not pd.isna(description):
                desc_s = str(description).lower()
        except (TypeError, ValueError):
            desc_s = str(description).lower()

    if layer == "COG_category":
        return term_s.upper() not in UNINFORMATIVE_COG_TERMS
    if layer == "PFAMs":
        if term_s.upper().startswith("DUF"):
            return False
        if desc_s and any(p in desc_s for p in UNINFORMATIVE_DESC_PATTERNS):
            return False
        return True
    if layer == "interpro_IPR":
        if desc_s and any(p in desc_s for p in UNINFORMATIVE_DESC_PATTERNS):
            return False
        # Also catch the case where the term name is a DUF passed through
        if term_s.upper().startswith("DUF"):
            return False
        return True
    if layer == "KEGG_ko":
        if desc_s and any(p in desc_s for p in UNINFORMATIVE_DESC_PATTERNS):
            return False
        return True
    if layer == "GOs":
        return term_s not in UNINFORMATIVE_GO_TERMS
    if layer == "EC":
        return term_s not in UNINFORMATIVE_EC_TERMS

    # Default: keep
    return True


# Protease PFAM families (peptidases)
PROTEASE_PFAMS = {
    "PF00082": "Subtilase (S8 serine protease)",
    "PF00026": "Eukaryotic aspartyl protease",
    "PF00089": "Trypsin-like serine protease",
    "PF01546": "Peptidase family M20",
    "PF00246": "Metallopeptidase family M12",
    "PF01435": "Peptidase family M48",
    "PF00814": "Peptidase family M1 (aminopeptidase)",
    "PF00675": "Peptidase family C14 (caspase)",
    "PF00112": "Papain family cysteine protease",
    "PF01650": "Peptidase family C13",
    "PF00560": "Leucine-rich repeat (often in proteases)",
    "PF03572": "Peptidase family S41",
    "PF00557": "Metallopeptidase family M24",
    "PF01432": "Peptidase family M3",
}

# Transporter PFAM families
TRANSPORTER_PFAMS = {
    "PF00005": "ABC transporter",
    "PF00664": "ABC transporter transmembrane region",
    "PF00083": "Sugar transporter (MFS)",
    "PF07690": "Major Facilitator Superfamily (MFS)",
    "PF00324": "Amino acid permease",
    "PF01490": "Transmembrane amino acid transporter",
    "PF01061": "ABC-2 type transporter",
    "PF00854": "POT family proton oligopeptide transporter",
    "PF01699": "Na+/solute symporter",
    "PF00520": "Ion channel",
    "PF01566": "Natural resistance-associated macrophage protein",
    "PF03092": "BET1 SNARE superfamily",
}

# Secretion-related GO terms
SECRETION_GO_TERMS = {
    "GO:0005576": "extracellular region",
    "GO:0005615": "extracellular space",
    "GO:0031012": "extracellular matrix",
    "GO:0005618": "cell wall",
    "GO:0009986": "cell surface",
    "GO:0046930": "pore complex",
    "GO:0006887": "exocytosis",
    "GO:0032940": "secretion by cell",
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _split_multi_terms(val: str):
    """Split multi-valued annotation fields into a list of terms."""
    if not isinstance(val, str):
        return []
    val = val.strip()
    if not val or val == "-":
        return []
    return [p.strip() for p in _MULTI_SEP_RE.split(val) if p.strip()]


def unique_join(vals):
    """Order-preserving, de-duplicating join for lists/iterables of strings."""
    seen, out = set(), []
    for v in vals:
        if not isinstance(v, str):
            continue
        s = v.strip()
        if not s or s in seen:
            continue
        seen.add(s)
        out.append(s)
    return ", ".join(out)


def combine_annotations(*cols):
    """Combine multiple annotation columns into a single column with unique values."""
    all_terms = []
    for col_val in cols:
        if pd.isna(col_val) or col_val == "" or col_val == "-":
            continue
        terms = _split_multi_terms(str(col_val))
        all_terms.extend(terms)
    return unique_join(all_terms)


# =============================================================================
# METADATA AND QC FUNCTIONS (from notebook Cell 2, 3)
# =============================================================================

def load_all_metadata(genus, species_list, json_basename="busco_updated.jsonl"):
    """Load metadata for ALL downloaded genomes (before filtering)."""
    frames = []
    for species in species_list:
        species_dir = os.path.join(f"/datadrive/Species/{genus}", species)
        json_path = os.path.join(species_dir, json_basename)

        if not os.path.isfile(json_path):
            print(f"[warn] No metadata found for {species}")
            continue

        records = []
        with open(json_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    records.append(json.loads(line))

        df = pd.DataFrame(records)
        df["Species"] = species

        col_map = {}
        for c in df.columns:
            cl = c.lower().replace(" ", "").replace("_", "")
            if "numberofcontigs" in cl:
                col_map[c] = "contigs"
            elif "contign50" in cl:
                col_map[c] = "n50"
            elif "buscocomplete" in cl:
                col_map[c] = "busco"
            elif "buscototal" in cl:
                col_map[c] = "busco_total"
            elif "totalsequencelength" in cl:
                col_map[c] = "total_len"
            elif "assemblyaccession" in cl:
                col_map[c] = "Assembly Accession"

        df = df.rename(columns=col_map)

        for col in ["contigs", "n50", "busco", "busco_total", "total_len"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        if "busco" in df.columns and "busco_total" in df.columns:
            df["busco"] = (df["busco"] / df["busco_total"] * 100)

        frames.append(df)

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def get_filtered_accessions(genus, species_list):
    """Get set of accessions that passed filtering (exist in funannotate_output)."""
    filtered = {}
    acc_re = re.compile(r"(GC[AF]_\d+\.\d+)")

    for species in species_list:
        funannotate_dir = os.path.join(f"/datadrive/Species/{genus}", species, "funannotate_output")
        accs = set()
        if os.path.isdir(funannotate_dir):
            for folder in os.listdir(funannotate_dir):
                m = acc_re.search(folder)
                if m:
                    accs.add(m.group(1))
        filtered[species] = accs

    return filtered


def get_proteome_accessions(genus, species_list):
    """Get set of accessions that have proteomes in filtered_protein directory."""
    proteomes = {}
    acc_re = re.compile(r"(GC[AF]_\d+\.\d+)")

    for species in species_list:
        protein_dir = os.path.join(f"/datadrive/Species/{genus}", species, "filtered_protein")
        accs = set()
        if os.path.isdir(protein_dir):
            for f in os.listdir(protein_dir):
                if f.endswith(".faa"):
                    m = acc_re.search(f)
                    if m:
                        accs.add(m.group(1))
        proteomes[species] = accs

    return proteomes


def calc_quality_score(row):
    """Calculate quality score: busco_pct - (contigs * 0.01) + (n50 * 0.000001)"""
    busco = row.get("busco", np.nan)
    contigs = row.get("contigs", np.nan)
    n50 = row.get("n50", np.nan)

    if pd.isna(busco) or pd.isna(contigs) or pd.isna(n50):
        return np.nan

    busco_pct = busco if busco > 1 else busco * 100
    return busco_pct - (contigs * 0.01) + (n50 * 0.000001)


def find_references(df, species_list, filtered_accs, proteome_accs):
    """Find reference genome and proteome per species."""
    references = {}

    for species in species_list:
        sp_data = df[df["Species"] == species].copy()

        if sp_data.empty:
            references[species] = {"ref_genome": None, "ref_proteome": None}
            continue

        sp_data["quality_score"] = sp_data.apply(calc_quality_score, axis=1)
        sp_data["is_gcf"] = sp_data["Assembly Accession"].str.startswith("GCF_")

        # Reference GENOME
        passed_accs = filtered_accs.get(species, set())
        passed_data = sp_data[sp_data["Assembly Accession"].isin(passed_accs)]

        ref_genome = None
        ref_genome_score = None
        if not passed_data.empty:
            gcf_passed = passed_data[passed_data["is_gcf"] == True]
            if not gcf_passed.empty and gcf_passed["quality_score"].notna().any():
                best_idx = gcf_passed["quality_score"].idxmax()
                ref_genome = gcf_passed.loc[best_idx, "Assembly Accession"]
                ref_genome_score = gcf_passed.loc[best_idx, "quality_score"]
            elif passed_data["quality_score"].notna().any():
                best_idx = passed_data["quality_score"].idxmax()
                ref_genome = passed_data.loc[best_idx, "Assembly Accession"]
                ref_genome_score = passed_data.loc[best_idx, "quality_score"]

        # Reference PROTEOME
        prot_accs = proteome_accs.get(species, set())
        prot_data = sp_data[sp_data["Assembly Accession"].isin(prot_accs)]

        ref_proteome = None
        ref_proteome_score = None
        if not prot_data.empty:
            gcf_prot = prot_data[prot_data["is_gcf"] == True]
            if not gcf_prot.empty and gcf_prot["quality_score"].notna().any():
                best_idx = gcf_prot["quality_score"].idxmax()
                ref_proteome = gcf_prot.loc[best_idx, "Assembly Accession"]
                ref_proteome_score = gcf_prot.loc[best_idx, "quality_score"]
            elif prot_data["quality_score"].notna().any():
                best_idx = prot_data["quality_score"].idxmax()
                ref_proteome = prot_data.loc[best_idx, "Assembly Accession"]
                ref_proteome_score = prot_data.loc[best_idx, "quality_score"]

        references[species] = {
            "ref_genome": ref_genome,
            "ref_genome_score": ref_genome_score,
            "ref_proteome": ref_proteome,
            "ref_proteome_score": ref_proteome_score,
            "same_reference": ref_genome == ref_proteome if ref_genome and ref_proteome else None
        }

        if ref_genome:
            g_type = "RefSeq" if "GCF" in ref_genome else "GenBank"
            print(f"  {species}: Ref Genome = {ref_genome} ({g_type})")
        if ref_proteome:
            p_type = "RefSeq" if "GCF" in ref_proteome else "GenBank"
            same_note = " [SAME]" if ref_genome == ref_proteome else " [DIFFERENT]"
            print(f"  {species}: Ref Proteome = {ref_proteome} ({p_type}){same_note}")

    return references


def plot_qc_metrics_per_species(df, species_list, references, metrics=None, figsize=(16, 12), ani_excluded=None):
    """Plot QC metrics (contigs, BUSCO, N50) for all downloaded genomes per species."""
    if metrics is None:
        metrics = [
            ("contigs", "Number of Contigs", True),
            ("busco", "BUSCO Completeness (%)", False),
            ("n50", "Contig N50 (bp)", True),
        ]

    n_species = len(species_list)
    n_metrics = len(metrics)

    species_colors = dict(zip(species_list, sns.color_palette("Set2", n_species)))
    grey_color = "#AAAAAA"
    ref_color = "#E53935"

    fig, axes = plt.subplots(n_metrics, n_species, figsize=figsize, sharey='row')

    if n_species == 1:
        axes = axes.reshape(-1, 1)
    if n_metrics == 1:
        axes = axes.reshape(1, -1)

    for col_idx, species in enumerate(species_list):
        sp_data = df[df["Species"] == species].copy()

        # Separate ANI-excluded samples for distinct coloring
        ani_excluded_accs = set()
        if ani_excluded is not None and isinstance(ani_excluded, dict):
            ani_excluded_accs = ani_excluded.get(species, set())

        ani_excluded_df = sp_data[sp_data["Assembly Accession"].isin(ani_excluded_accs)]

        passed_non_ref = sp_data[(sp_data["passed_filter"] == True) & (sp_data["reference_type"].isna())]
        filtered_non_ref = sp_data[(sp_data["passed_filter"] == False) & (sp_data["reference_type"].isna()) & (~sp_data["Assembly Accession"].isin(ani_excluded_accs))]
        ref_genome_only = sp_data[sp_data["reference_type"] == "genome"]
        ref_proteome_only = sp_data[sp_data["reference_type"] == "proteome"]
        ref_both = sp_data[sp_data["reference_type"] == "both"]

        for row_idx, (metric_col, metric_name, use_log) in enumerate(metrics):
            ax = axes[row_idx, col_idx]

            passed_vals = passed_non_ref[metric_col].dropna()
            filtered_vals = filtered_non_ref[metric_col].dropna()

            jitter_width = 0.3

            if len(filtered_vals) > 0:
                x_filtered = np.random.uniform(-jitter_width, jitter_width, len(filtered_vals))
                ax.scatter(x_filtered, filtered_vals, c=grey_color, alpha=0.5, s=40,
                          edgecolors='none', label='Filtered out', zorder=1)

            if len(passed_vals) > 0:
                x_passed = np.random.uniform(-jitter_width, jitter_width, len(passed_vals))
                ax.scatter(x_passed, passed_vals, c=[species_colors[species]], alpha=0.8, s=50,
                          edgecolors='white', linewidths=0.5, label='Passed', zorder=2)

            # ANI-excluded: distinct red X marker
            if not ani_excluded_df.empty:
                ani_vals = ani_excluded_df[metric_col].dropna()
                if len(ani_vals) > 0:
                    x_ani = np.random.uniform(-jitter_width, jitter_width, len(ani_vals))
                    ax.scatter(x_ani, ani_vals, c='#FF1744', marker='X', s=100,
                              edgecolors='black', linewidths=1, label='ANI excluded', zorder=6)
                    for idx_a, (_, arow) in enumerate(ani_excluded_df.dropna(subset=[metric_col]).iterrows()):
                        ax.annotate(arow["Assembly Accession"], (x_ani[idx_a], arow[metric_col]),
                                   fontsize=7, color='#FF1744', ha='left', va='bottom',
                                   xytext=(4, 2), textcoords='offset points', zorder=7)

            if not ref_genome_only.empty:
                ref_g = ref_genome_only[[metric_col, "Assembly Accession"]].dropna(subset=[metric_col])
                if len(ref_g) > 0:
                    ax.scatter([-0.15] * len(ref_g), ref_g[metric_col].values, c=ref_color, marker='*', s=100,
                              edgecolors='white', linewidths=0.5, label='Ref Genome', zorder=4)
                    for _, rrow in ref_g.iterrows():
                        ax.annotate(rrow["Assembly Accession"], (-0.15, rrow[metric_col]),
                                   fontsize=8, color=ref_color, ha='left', va='bottom',
                                   xytext=(4, 2), textcoords='offset points', zorder=6)

            if not ref_proteome_only.empty:
                ref_p = ref_proteome_only[[metric_col, "Assembly Accession"]].dropna(subset=[metric_col])
                if len(ref_p) > 0:
                    ax.scatter([0.15] * len(ref_p), ref_p[metric_col].values, c=ref_color, marker='D', s=100,
                              edgecolors='white', linewidths=0.5, label='Ref Proteome', zorder=4)
                    for _, rrow in ref_p.iterrows():
                        ax.annotate(rrow["Assembly Accession"], (0.15, rrow[metric_col]),
                                   fontsize=8, color=ref_color, ha='left', va='bottom',
                                   xytext=(4, 2), textcoords='offset points', zorder=6)

            if not ref_both.empty:
                ref_b = ref_both[[metric_col, "Assembly Accession"]].dropna(subset=[metric_col])
                if len(ref_b) > 0:
                    ax.scatter([0] * len(ref_b), ref_b[metric_col].values, c=ref_color, marker='h', s=100,
                              edgecolors='white', linewidths=0.5, label='Ref Both', zorder=5)
                    for _, rrow in ref_b.iterrows():
                        ax.annotate(rrow["Assembly Accession"], (0, rrow[metric_col]),
                                   fontsize=8, color=ref_color, ha='left', va='bottom',
                                   xytext=(4, 2), textcoords='offset points', zorder=6)

            if use_log:
                ax.set_yscale('log')
                ax.yaxis.set_major_formatter(ScalarFormatter())

            ax.set_xlim(-0.6, 0.6)
            ax.set_xticks([0])
            if row_idx == n_metrics - 1:
                ax.set_xticklabels(["Strains"], fontsize=9)
            else:
                ax.set_xticklabels([])

            n_passed = len(sp_data[sp_data["passed_filter"] == True])
            n_filtered = len(sp_data[sp_data["passed_filter"] == False])
            n_total = n_passed + n_filtered
            ax.text(0.02, 0.98, f"n={n_total}\n({n_passed} passed)",
                   transform=ax.transAxes, fontsize=9, va='top', ha='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            if row_idx == 0:
                ax.set_title(f"A. {species}", fontsize=12, fontweight='bold')

            if col_idx == 0:
                ax.set_ylabel(metric_name, fontsize=10)

            ax.grid(True, axis='y', linestyle=':', alpha=0.5)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#66BB6A',
               markersize=10, label='Passed filter'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=grey_color,
               markersize=10, alpha=0.5, label='Filtered out'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor=ref_color,
               markeredgecolor='white', markersize=14, label='Ref Genome (parsnp)'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor=ref_color,
               markeredgecolor='white', markersize=7, label='Ref Proteome (funannotate)'),
        Line2D([0], [0], marker='h', color='w', markerfacecolor=ref_color,
               markeredgecolor='white', markersize=10, label='Both (same strain)'),
        Line2D([0], [0], marker='X', color='w', markerfacecolor='#FF1744',
               markeredgecolor='black', markersize=10, label='ANI excluded (wrong species)')
    ]
    fig.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.99, 0.99),
               fontsize=9, frameon=True)

    fig.suptitle(f"Genome QC Metrics for Aspergillus Species",
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()



# =============================================================================
# EGGNOG I/O FUNCTIONS (from functions_for_analysis_1_2.py)
# =============================================================================

def read_emapper_annotations(path: str) -> pd.DataFrame:
    """Robust reader for eggNOG-mapper .emapper.annotations files."""
    hdr_idx = None
    queries_scanned = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if line.startswith("#query") or line.startswith("query"):
                hdr_idx = i
                break
            if "queries scanned" in line:
                parts = line.split()
                if len(parts) >= 3:
                    queries_scanned = parts[1]

    if hdr_idx is None:
        if queries_scanned == "0":
            raise ValueError(f"EggNOG file appears empty (0 queries scanned): {path}")
        else:
            raise ValueError(f"Could not find header (query/#query) in {path}")

    df = pd.read_csv(path, sep="\t", header=0, skiprows=hdr_idx, dtype=str, na_filter=False)

    colmap = {}
    for c in df.columns:
        cl = c.lstrip("#").strip().lower()
        if cl in ("query", "query_name"):
            colmap[c] = "Protein_ID"
        elif cl == "description":
            colmap[c] = "Description_EggNog"
        elif cl in ("cog_category", "cog", "cogs"):
            colmap[c] = "COG_category"
        elif cl in ("go", "gos", "go_terms"):
            colmap[c] = "GOs_EggNog"
        elif cl in ("ec", "ec_number", "ec_numbers"):
            colmap[c] = "EC_EggNog"
        elif cl in ("kegg_ko", "ko"):
            colmap[c] = "KEGG_ko"
        elif cl in ("kegg_pathway", "pathway"):
            colmap[c] = "KEGG_Pathway"
        elif cl in ("kegg_module", "module"):
            colmap[c] = "KEGG_Module"
        elif cl in ("kegg_reaction", "reaction"):
            colmap[c] = "KEGG_Reaction"
        elif cl in ("kegg_rclass", "rclass"):
            colmap[c] = "KEGG_rclass"
        elif cl in ("brite", "kegg_brite"):
            colmap[c] = "BRITE"
        elif cl in ("kegg_tc", "tcdb", "tc"):
            colmap[c] = "KEGG_TC"
        elif cl in ("cazy", "cazymes"):
            colmap[c] = "CAZy_EggNog"
        elif cl in ("bigg_reaction", "bigg"):
            colmap[c] = "BiGG_Reaction"
        elif cl in ("pfams", "pfam"):
            colmap[c] = "PFAMs_EggNog"
    if colmap:
        df = df.rename(columns=colmap)

    for need in ["Protein_ID"] + EGGNOG_KEEP:
        if need not in df.columns:
            df[need] = ""

    df["Protein_ID"] = df["Protein_ID"].astype(str).str.strip()
    return df[["Protein_ID"] + EGGNOG_KEEP]


def load_emapper_for_accessions(eggnog_dir: str, accessions) -> pd.DataFrame:
    """Load/stack eggNOG annotations for a set of assembly accessions."""
    frames = []
    for acc in sorted(pd.Series(accessions).dropna().unique()):
        ann = os.path.join(eggnog_dir, f"{acc}.emapper.annotations")
        if not os.path.exists(ann):
            hits = glob.glob(os.path.join(eggnog_dir, f"{acc}*.emapper.annotations"))
            if not hits:
                print(f"[eggNOG] Missing for {acc} — skipping")
                continue
            ann = hits[0]
        try:
            df = read_emapper_annotations(ann)
        except Exception as e:
            print(f"[eggNOG] Failed to parse {os.path.basename(ann)}: {e}")
            continue
        df["Assembly Accession"] = acc
        frames.append(df)

    if frames:
        out = pd.concat(frames, ignore_index=True)
    else:
        out = pd.DataFrame(columns=["Protein_ID", "Assembly Accession"] + EGGNOG_KEEP)
    return out


# =============================================================================
# DBCAN I/O FUNCTIONS (from functions_for_analysis_1_2.py)
# =============================================================================

def _extract_cazy_families(*fields) -> str:
    """Pull CAZy families like GH5, GT2, AA9, GH16_1 from given strings."""
    fams = []
    for f in fields:
        if not isinstance(f, str) or f.strip() in ("", "-"):
            continue
        fams.extend(_CAZY_FAM_RE.findall(f))
    return unique_join(fams)


def read_dbcan_overview(path: str) -> pd.DataFrame:
    """Parse dbCAN overview.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str, na_filter=False)

    rename_exact = {
        "Gene ID": "Protein_ID",
        "EC#": "dbcan_ECs",
        "dbCAN_hmm": "dbcan_hmm_raw",
        "dbCAN_sub": "dbcan_sub_raw",
        "DIAMOND": "dbcan_diamond_raw",
        "#ofTools": "dbcan_NumTools",
        "Recommend Results": "dbcan_Recommend",
        "Substrate": "dbcan_Substrate",
    }
    colmap = {}
    for c in df.columns:
        k = c.strip()
        colmap[c] = rename_exact.get(k, c)
    df = df.rename(columns=colmap)

    for need in ["Protein_ID", "dbcan_Recommend", "dbcan_NumTools", "dbcan_ECs",
                 "dbcan_Substrate", "dbcan_hmm_raw", "dbcan_sub_raw", "dbcan_diamond_raw"]:
        if need not in df.columns:
            df[need] = ""

    df["Protein_ID"] = df["Protein_ID"].astype(str).str.strip()

    fams = []
    for rec, hmm, sub, dia in zip(df["dbcan_Recommend"], df["dbcan_hmm_raw"],
                                  df["dbcan_sub_raw"], df["dbcan_diamond_raw"]):
        fam_rec = _extract_cazy_families(rec)
        if fam_rec:
            fams.append(fam_rec)
        else:
            fams.append(_extract_cazy_families(hmm, sub, dia))

    out = pd.DataFrame({
        "Protein_ID": df["Protein_ID"],
        "CAZy_dbCAN": fams,
        "dbcan_Recommend": df["dbcan_Recommend"].astype(str).str.strip(),
        "EC_dbCAN": df["dbcan_ECs"].astype(str).str.replace(r"\s+", "", regex=True),
        "dbcan_NumTools": df["dbcan_NumTools"].astype(str).str.strip(),
        "dbcan_Substrate": df["dbcan_Substrate"].astype(str).str.strip(),
    })
    return out


def load_dbcan_for_accessions(dbcan_dir: str, accessions) -> pd.DataFrame:
    """overview.tsv-only loader for dbCAN."""
    frames = []
    for acc in sorted(pd.Series(accessions).dropna().unique()):
        acc_dir = os.path.join(dbcan_dir, acc)
        if not os.path.isdir(acc_dir):
            hits = glob.glob(os.path.join(dbcan_dir, f"{acc}*"))
            if hits:
                acc_dir = hits[0]
            else:
                print(f"[dbCAN] Missing dir for {acc} — skipping")
                continue

        overview = os.path.join(acc_dir, "overview.tsv")
        if not os.path.exists(overview):
            print(f"[dbCAN] No overview.tsv for {acc} — skipping")
            continue

        try:
            ov_df = read_dbcan_overview(overview)
        except Exception as e:
            print(f"[dbCAN] Failed to parse overview for {acc}: {e}")
            continue

        ov_df["Assembly Accession"] = acc
        frames.append(ov_df[["Protein_ID", "Assembly Accession"] + DBCAN_KEEP])

    if frames:
        return pd.concat(frames, ignore_index=True)
    else:
        return pd.DataFrame(columns=["Protein_ID", "Assembly Accession"] + DBCAN_KEEP)


# =============================================================================
# INTERPROSCAN I/O FUNCTIONS (from functions_for_analysis_1_2.py)
# =============================================================================

def read_interproscan_tsv(path: str) -> pd.DataFrame:
    """Parse InterProScan TSV output."""
    df = pd.read_csv(
        path, sep="\t", header=None, dtype=str, na_filter=False,
        names=["Protein_ID", "MD5", "Length", "Analysis", "Signature_Acc", "Signature_Desc",
               "Start", "Stop", "Score", "Status", "Date", "IPR_Acc", "IPR_Desc", "GO", "Pathways"]
    )

    grouped = df.groupby("Protein_ID", as_index=False).agg({
        "IPR_Acc": lambda x: unique_join([v for v in x if v and v != "-"]),
        "GO": lambda x: unique_join([v for v in x if v and v != "-"]),
        "Pathways": lambda x: unique_join([v for v in x if v and v != "-"]),
        "IPR_Desc": lambda x: unique_join([v for v in x if v and v != "-"]),
    })

    pfam_df = (
        df[df["Analysis"] == "Pfam"]
        .groupby("Protein_ID", as_index=False)
        .agg({"Signature_Acc": lambda x: unique_join([v for v in x if v and v != "-" and v.startswith("PF")])})
        .rename(columns={"Signature_Acc": "PFAMs_Interpro"})
    )

    def extract_ec(pathways_str):
        if not pathways_str or pathways_str == "-":
            return ""
        ec_pattern = re.compile(r'\b(EC:)?(\d+\.\d+\.\d+\.\d+)\b')
        matches = ec_pattern.findall(pathways_str)
        ecs = [m[1] for m in matches if m[1]]
        return unique_join(ecs)

    ec_df = grouped[["Protein_ID", "Pathways"]].copy()
    ec_df["EC_Interpro"] = ec_df["Pathways"].apply(extract_ec)

    out = pd.DataFrame({
        "Protein_ID": grouped["Protein_ID"],
        "interpro_IPR": grouped["IPR_Acc"],
        "GOs_Interpro": grouped["GO"],
        "interpro_Pathways": grouped["Pathways"],
        "Description_Interpro": grouped["IPR_Desc"],
    })

    out = out.merge(pfam_df, on="Protein_ID", how="left")
    out["PFAMs_Interpro"] = out["PFAMs_Interpro"].fillna("")

    out = out.merge(ec_df[["Protein_ID", "EC_Interpro"]], on="Protein_ID", how="left")
    out["EC_Interpro"] = out["EC_Interpro"].fillna("")

    return out


def load_interproscan_for_accessions(interproscan_dir: str, accessions) -> pd.DataFrame:
    """Load InterProScan TSV files."""
    frames = []
    for acc in sorted(pd.Series(accessions).dropna().unique()):
        acc_pattern = os.path.join(interproscan_dir, f"{acc}*")
        hits = [h for h in glob.glob(acc_pattern) if os.path.isdir(h)]
        if not hits:
            print(f"[InterProScan] Missing dir for {acc} — skipping")
            continue
        acc_dir = hits[0]

        tsv_pattern = os.path.join(acc_dir, f"{acc}*.tsv")
        tsv_files = glob.glob(tsv_pattern)
        if not tsv_files:
            print(f"[InterProScan] No TSV for {acc} in {acc_dir} — skipping")
            continue

        tsv_file = tsv_files[0]

        try:
            df = read_interproscan_tsv(tsv_file)
        except Exception as e:
            print(f"[InterProScan] Failed to parse {os.path.basename(tsv_file)}: {e}")
            continue

        df["Assembly Accession"] = acc
        frames.append(df)

    if frames:
        out = pd.concat(frames, ignore_index=True)
    else:
        out = pd.DataFrame(columns=["Protein_ID", "Assembly Accession"] + INTERPRO_KEEP)

    for col in INTERPRO_KEEP:
        if col not in out.columns:
            out[col] = ""

    return out


# =============================================================================
# SIGNALP I/O FUNCTIONS (from functions_for_analysis_1_2.py)
# =============================================================================

def read_signalp6_txt(path: str) -> pd.DataFrame:
    """Parse SignalP 6.x prediction_results.txt output."""
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            raw_id = parts[0].strip()
            prot_id = raw_id.split()[0] if " " in raw_id else raw_id

            prediction = parts[1].strip() if len(parts) > 1 else ""
            other_prob = parts[2].strip() if len(parts) > 2 else ""
            sp_prob = parts[3].strip() if len(parts) > 3 else ""
            cs_position = parts[4].strip() if len(parts) > 4 else ""

            rows.append({
                "Protein_ID": prot_id,
                "signalp_Prediction": prediction,
                "signalp_SP_Probability": sp_prob,
                "signalp_OTHER_Probability": other_prob,
                "signalp_Cleavage_Site": cs_position,
            })

    if not rows:
        return pd.DataFrame(columns=["Protein_ID", "signalp_Prediction",
                                      "signalp_SP_Probability", "signalp_Cleavage_Site"])

    return pd.DataFrame(rows)


def load_signalp_for_accessions(signalp_dir: str, accessions) -> pd.DataFrame:
    """Load/stack SignalP predictions for a set of assembly accessions."""
    frames = []
    for acc in sorted(pd.Series(accessions).dropna().unique()):
        acc_dir = os.path.join(signalp_dir, acc)
        if not os.path.isdir(acc_dir):
            hits = glob.glob(os.path.join(signalp_dir, f"{acc}*"))
            if hits and os.path.isdir(hits[0]):
                acc_dir = hits[0]
            else:
                print(f"[SignalP] Missing dir for {acc} — skipping")
                continue

        sp6_file = os.path.join(acc_dir, "prediction_results.txt")
        if os.path.exists(sp6_file):
            try:
                df = read_signalp6_txt(sp6_file)
                df["Assembly Accession"] = acc
                frames.append(df)
                continue
            except Exception as e:
                print(f"[SignalP] Failed to parse SignalP6 for {acc}: {e}")

        print(f"[SignalP] No output files for {acc} — skipping")

    if frames:
        out = pd.concat(frames, ignore_index=True)
    else:
        out = pd.DataFrame(columns=["Protein_ID", "Assembly Accession"] + SIGNALP_KEEP)
    return out




# =============================================================================
# SPECIES CONFIGURATION (shared across notebooks)
# =============================================================================

GENUS = 'Aspergillus'

SPECIES_LIST = ['fumigatus', 'flavus', 'niger', 'oryzae']

# ANI-failed accessions to exclude from per-species analyses at runtime.
# These were below the 95% ANI threshold to the GCF reference and likely
# represent misidentified species (kept on disk so the audit trail is preserved,
# but filtered when loaded by the analysis functions).
ANI_EXCLUDED = {
    'flavus':    {'GCA_023653635'},   # A. parasiticus, ANI ~94.08% to A. flavus reference
    'fumigatus': set(),
    'niger':     set(),               # A. tubingensis already removed at directory level
    'oryzae':    set(),
}


def is_ani_excluded(accession_or_label, species):
    """Return True if the accession matches an ANI-excluded entry for the species."""
    excluded = ANI_EXCLUDED.get(species, set())
    if not excluded:
        return False
    s = str(accession_or_label)
    return any(acc in s for acc in excluded)


def filter_ani_excluded(df, species, axis='auto'):
    """Drop rows/columns matching ANI-excluded accessions for ``species``.

    ``axis`` may be ``'rows'``, ``'cols'``, ``'both'``, or ``'auto'`` which
    drops whichever axis carries the matching accession.
    """
    excluded = ANI_EXCLUDED.get(species, set())
    if not excluded or df is None:
        return df
    df_out = df
    if axis in ('cols', 'both', 'auto'):
        keep_cols = [c for c in df_out.columns
                     if not any(a in str(c) for a in excluded)]
        if len(keep_cols) < df_out.shape[1]:
            df_out = df_out[keep_cols]
    if axis in ('rows', 'both', 'auto'):
        # Index has .map (not .apply); cast to Series so boolean indexing works
        keep = ~pd.Series(df_out.index.astype(str), index=df_out.index).map(
            lambda s: any(a in s for a in excluded))
        if keep.sum() < len(df_out):
            df_out = df_out[keep]
    return df_out

SPECIES_DISPLAY = {
    'fumigatus': 'A. fumigatus',
    'flavus':    'A. flavus',
    'niger':     'A. niger',
    'oryzae':    'A. oryzae',
}

SPECIES_COLORS = {
    'fumigatus': '#e74c3c',
    'flavus':    '#ffd500',
    'niger':     '#1a1aff',
    'oryzae':    '#27ae60',
}

# Display-name keyed version (used by world map and other plots)
SPECIES_COLORS_DISPLAY = {
    'A. fumigatus': '#e74c3c',
    'A. flavus':    '#ffd500',
    'A. niger':     '#1a1aff',
    'A. oryzae':    '#27ae60',
}

ISOLATION_CLASSES = [
    'Human-pathogenic', 'Animal-pathogenic', 'Plant-pathogenic',
    'Industrial-trait', 'Environmental', 'Lab',
]

ISOLATION_CLASS_COLORS = {
    'Human-pathogenic':  '#E64B35',
    'Animal-pathogenic': '#F39B7F',
    'Plant-pathogenic':  '#91D1C2',
    'Industrial-trait':  '#3C5488',
    'Environmental':     '#00A087',
    'Lab':               '#B09C85',
    'Unknown':           '#CCCCCC',
}


def load_and_prepare_metadata(genus=None, species_list=None):
    """Load all genome metadata, apply QC filter flags, and identify references.

    Returns
    -------
    df_all : pd.DataFrame
        All genomes with ``passed_filter``, ``reference_type``, and ``is_reference`` columns.
    filtered_accs : dict
        ``{species: set_of_accessions}`` that passed QC.
    proteome_accs : dict
        ``{species: set_of_accessions}`` with proteomes available.
    references : dict
        Per-species reference genome/proteome info.
    """
    if genus is None:
        genus = GENUS
    if species_list is None:
        species_list = SPECIES_LIST

    print('Loading metadata for all downloaded genomes...')
    df_all = load_all_metadata(genus, species_list)
    print(f'Total genomes downloaded: {len(df_all)}')

    # Get filtered accessions (passed QC -> funannotate_output)
    filtered_accs = get_filtered_accessions(genus, species_list)
    for sp, accs in filtered_accs.items():
        print(f'  {sp}: {len(accs)} passed filtering')

    # Get proteome accessions
    proteome_accs = get_proteome_accessions(genus, species_list)
    print('\nProteomes available (filtered_protein):')
    for sp, accs in proteome_accs.items():
        print(f'  {sp}: {len(accs)} proteomes')

    # Mark which assemblies passed filtering
    def _mark_filtered(row):
        species = row['Species']
        acc = row.get('Assembly Accession', '')
        return acc in filtered_accs.get(species, set())

    df_all['passed_filter'] = df_all.apply(_mark_filtered, axis=1)
    passed = df_all['passed_filter'].sum()
    failed = len(df_all) - passed
    print(f'\nFiltering summary: {passed} passed, {failed} filtered out')

    # Identify reference genomes and proteomes per species
    print('Identifying reference genomes (parsnp) and proteomes (funannotate)...')
    references = find_references(df_all, species_list, filtered_accs, proteome_accs)

    def _get_reference_type(row):
        species = row['Species']
        acc = row.get('Assembly Accession', '')
        ref_info = references.get(species, {})
        is_ref_genome = acc == ref_info.get('ref_genome')
        is_ref_proteome = acc == ref_info.get('ref_proteome')
        if is_ref_genome and is_ref_proteome:
            return 'both'
        elif is_ref_genome:
            return 'genome'
        elif is_ref_proteome:
            return 'proteome'
        return None

    df_all['reference_type'] = df_all.apply(_get_reference_type, axis=1)
    df_all['is_reference'] = df_all['reference_type'].notna()

    print(f'\nReference summary:')
    print(f'  Ref genomes: {(df_all["reference_type"] == "genome").sum() + (df_all["reference_type"] == "both").sum()}')
    print(f'  Ref proteomes: {(df_all["reference_type"] == "proteome").sum() + (df_all["reference_type"] == "both").sum()}')

    return df_all, filtered_accs, proteome_accs, references


def apply_ani_exclusions(df_all, filtered_accs, ani_excluded):
    """Remove ANI-excluded accessions from filtered_accs and update df_all.

    Parameters
    ----------
    df_all : pd.DataFrame
        Full metadata with ``passed_filter`` column.
    filtered_accs : dict
        ``{species: set_of_accessions}`` -- modified in-place.
    ani_excluded : dict
        ``{species: set_of_accessions_to_exclude}``.

    Returns
    -------
    all_passed_accs : set
        Union of all passed accessions after ANI exclusion.
    """
    for species, excluded_accs in ani_excluded.items():
        if species in filtered_accs:
            before = len(filtered_accs[species])
            filtered_accs[species] -= excluded_accs
            removed = before - len(filtered_accs[species])
            if removed > 0:
                print(f'  ANI exclusion: removed {removed} from {species}')
                for acc in excluded_accs:
                    mask = df_all['Assembly Accession'].str.startswith(acc.split('.')[0])
                    df_all.loc[mask, 'passed_filter'] = False

    passed = df_all['passed_filter'].sum()
    print(f'  After ANI exclusion: {passed} passed')

    all_passed_accs = set()
    for sp, accs in filtered_accs.items():
        print(f'  {sp}: {len(accs)} passed QC')
        all_passed_accs.update(accs)
    print(f'\nTotal QC-passed assemblies: {len(all_passed_accs)}')

    return all_passed_accs


def _run_fastani_for_species(species, genus='Aspergillus', threads=8):
    """Run fastANI for one species: all filtered genomes vs GCF reference.

    Returns list of dicts with 'accession' and 'ani', or empty list on failure.
    """
    import subprocess, tempfile
    genome_dir = f'/datadrive/Species/{genus}/{species}/filtered_genome'
    if not os.path.isdir(genome_dir):
        print(f'  {species}: No filtered_genome directory')
        return []

    # Find GCF reference
    ref = None
    for f in os.listdir(genome_dir):
        if f.startswith('GCF_') and f.endswith('.fna') and '_renamed' not in f:
            ref = os.path.join(genome_dir, f)
            break
    if ref is None:
        print(f'  {species}: No GCF reference genome found')
        return []

    # Build query list (all .fna except renamed)
    queries = [os.path.join(genome_dir, f)
               for f in os.listdir(genome_dir)
               if f.endswith('.fna') and '_renamed' not in f]
    if not queries:
        return []

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as ql:
        ql.write('\n'.join(queries))
        ql_path = ql.name

    out_path = ql_path + '.ani'
    try:
        subprocess.run(
            ['fastANI', '--ql', ql_path, '--ref', ref, '-o', out_path, '-t', str(threads)],
            capture_output=True, text=True, timeout=600
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f'  {species}: fastANI failed ({e})')
        os.unlink(ql_path)
        return []

    records = []
    if os.path.exists(out_path):
        with open(out_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    m = re.search(r'(GC[AF]_\d+\.\d+)', parts[0])
                    if m:
                        records.append({'accession': m.group(1), 'ani': float(parts[2])})
        os.unlink(out_path)
    os.unlink(ql_path)
    return records


def load_ani_results(species_list, results_dir, genus='Aspergillus', threads=8):
    """Load ANI results from cached TSV files, computing with fastANI if missing.

    On first run, fastANI is called for each species and results are saved to
    ``results_dir/{species}/{species}_ani_results.tsv`` so subsequent runs are instant.

    Returns
    -------
    ani_results_all : dict
        ``{species: pd.DataFrame}`` with columns ``accession`` and ``ani``.
    """
    ani_results_all = {}
    for species in species_list:
        sp_dir = os.path.join(results_dir, species)
        os.makedirs(sp_dir, exist_ok=True)
        ani_path = os.path.join(sp_dir, f'{species}_ani_results.tsv')

        if os.path.exists(ani_path):
            # Load cached results
            records = []
            with open(ani_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        m = re.search(r'(GC[AF]_\d+\.\d+)', parts[0])
                        ani_val = float(parts[2])
                        if m:
                            records.append({'accession': m.group(1), 'ani': ani_val})
        else:
            # Compute with fastANI
            print(f'  {species}: No cached ANI results — running fastANI...')
            records = _run_fastani_for_species(species, genus, threads)
            if records:
                # Save for next time
                with open(ani_path, 'w') as f:
                    for r in records:
                        f.write(f"{r['accession']}\tref\t{r['ani']}\n")
                print(f'  {species}: Saved {len(records)} ANI results to {ani_path}')

        if records:
            ani_results_all[species] = pd.DataFrame(records)
            n_below = sum(1 for r in records if r['ani'] < 95)
            print(f'  {species}: {len(records)} genomes, '
                  f'ANI range: {min(r["ani"] for r in records):.2f} - '
                  f'{max(r["ani"] for r in records):.2f}%, '
                  f'{n_below} below 95%')
    return ani_results_all


def plot_ani_verification(ani_results_all, species_list, results_base):
    """Plot ANI species verification jitter plots.

    Returns the matplotlib Figure or None if no data.
    """
    if not ani_results_all:
        print('\nNo ANI results found. Run fastANI first (see ani_and_filter_genome_QC.sh)')
        return None

    n_species = len(ani_results_all)
    fig, axes = plt.subplots(1, n_species, figsize=(4 * n_species, 5), sharey=True)
    if n_species == 1:
        axes = [axes]
    sp_colors = dict(zip(species_list, sns.color_palette('Set2', len(species_list))))

    for idx, species in enumerate(species_list):
        ax = axes[idx]
        if species not in ani_results_all:
            ax.set_title(f'A. {species}\n(no data)')
            continue
        df_ani = ani_results_all[species]
        passed_ani = df_ani[df_ani['ani'] >= 95]
        failed_ani = df_ani[df_ani['ani'] < 95]

        jitter = np.random.uniform(-0.3, 0.3, len(passed_ani))
        ax.scatter(jitter, passed_ani['ani'], c=[sp_colors.get(species, 'steelblue')],
                   s=50, alpha=0.7, edgecolors='white', linewidths=0.5, zorder=2)
        if len(failed_ani) > 0:
            jitter_f = np.random.uniform(-0.3, 0.3, len(failed_ani))
            ax.scatter(jitter_f, failed_ani['ani'], c='#FF1744', marker='X', s=120,
                       edgecolors='black', linewidths=1.5, zorder=5)
            for _, row in failed_ani.iterrows():
                ax.annotate(row['accession'], (0, row['ani']),
                            fontsize=7, color='#FF1744', ha='center', va='top',
                            xytext=(0, -8), textcoords='offset points')

        ax.axhline(y=95, color='red', linestyle='--', alpha=0.7, linewidth=1.5)
        ax.text(0.4, 95.05, 'ANI 95%', color='red', fontsize=8, va='bottom')
        ax.set_xlim(-0.6, 0.6)
        ax.set_xticks([0])
        ax.set_xticklabels(['Strains'])
        ax.set_title(f'A. {species}\n(n={len(df_ani)}, {len(failed_ani)} excluded)', fontweight='bold')
        if idx == 0:
            ax.set_ylabel('ANI to species reference (%)')
        ax.grid(axis='y', linestyle=':', alpha=0.5)
        ax.text(0.02, 0.02, f'pass: {len(passed_ani)}\nfail: {len(failed_ani)}',
                transform=ax.transAxes, fontsize=9, va='bottom',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('ANI Species Verification (fastANI vs GCF reference)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(results_base, 'ani_verification.png'), dpi=150, bbox_inches='tight')
    print('\nThreshold: 95% ANI -- samples below are likely different species')


def build_qc_summary_tables(df_all, species_list, references):
    """Build QC filtering summary and QC metrics summary DataFrames.

    Returns
    -------
    summary_df : pd.DataFrame
        Genome counts and reference info per species.
    qc_stats_df : pd.DataFrame
        QC metric statistics for passed genomes.
    """
    summary_stats = []
    for species in species_list:
        sp_data = df_all[df_all['Species'] == species]
        passed_sp = sp_data[sp_data['passed_filter'] == True]
        filtered_sp = sp_data[sp_data['passed_filter'] == False]
        ref_info = references.get(species, {})
        ref_genome = ref_info.get('ref_genome', None)
        ref_proteome = ref_info.get('ref_proteome', None)
        g_type = 'RefSeq' if ref_genome and 'GCF' in str(ref_genome) else 'GenBank' if ref_genome else '-'
        p_type = 'RefSeq' if ref_proteome and 'GCF' in str(ref_proteome) else 'GenBank' if ref_proteome else '-'
        same_strain = 'Yes' if ref_genome and ref_genome == ref_proteome else 'No' if ref_genome and ref_proteome else '-'
        summary_stats.append({
            'Species': f'A. {species}',
            'Downloaded': len(sp_data),
            'Passed': len(passed_sp),
            'Filtered': len(filtered_sp),
            'Pass %': f'{100 * len(passed_sp) / len(sp_data):.1f}' if len(sp_data) > 0 else '-',
            'Ref Genome': ref_genome if ref_genome else '-',
            'G Type': g_type,
            'Ref Proteome': ref_proteome if ref_proteome else '-',
            'P Type': p_type,
            'Same?': same_strain,
        })
    summary_df = pd.DataFrame(summary_stats)

    qc_stats = []
    for species in species_list:
        passed_sp = df_all[(df_all['Species'] == species) & (df_all['passed_filter'] == True)]
        if len(passed_sp) > 0:
            qc_stats.append({
                'Species': f'A. {species}',
                'Median Contigs': f'{passed_sp["contigs"].median():.0f}',
                'Median BUSCO': f'{passed_sp["busco"].median():.1f}%',
                'Median N50': f'{passed_sp["n50"].median() / 1e6:.2f} Mb',
                'Min BUSCO': f'{passed_sp["busco"].min():.1f}%',
                'Max Contigs': f'{passed_sp["contigs"].max():.0f}',
            })
    qc_stats_df = pd.DataFrame(qc_stats)

    return summary_df, qc_stats_df


def load_or_build_enriched_metadata(all_passed_accs, enriched_cache, results_base,
                                     genus=None, species_list=None,
                                     species_display=None):
    """Load pre-enriched metadata from cache or build from scratch.

    Returns
    -------
    df_meta : pd.DataFrame
        Enriched metadata with ``CombinedText`` column ready for classification.
    """
    import funpan_classify as classify

    if genus is None:
        genus = GENUS
    if species_list is None:
        species_list = SPECIES_LIST
    if species_display is None:
        species_display = SPECIES_DISPLAY

    if os.path.exists(enriched_cache):
        print(f'Loading pre-enriched metadata from {enriched_cache}')
        df_meta = pd.read_csv(enriched_cache)

        existing_accs = set(
            df_meta['Assembly Accession']
            .str.extract(r'(GC[AF]_\d+\.\d+)')[0]
            .dropna()
        )
        missing = all_passed_accs - existing_accs
        extra = existing_accs - all_passed_accs
        if missing:
            print(f'  WARNING: {len(missing)} new QC-passed assemblies not in saved file')
            print(f'  Delete {enriched_cache} to re-enrich all assemblies')
        if extra:
            print(f'  Note: {len(extra)} assemblies in saved file no longer pass QC')
            df_meta = df_meta[
                df_meta['Assembly Accession']
                .str.extract(r'(GC[AF]_\d+\.\d+)')[0]
                .isin(all_passed_accs)
            ]
        print(f'  Loaded {len(df_meta)} assemblies')

        if 'WGS Reference Title' not in df_meta.columns:
            print('\nEnriching with WGS Reference Titles from NCBI...')
            df_meta = classify.enrich_df_with_wgs_title(
                df_meta, wgs_col='WGS project accession', pause=0.35
            )
            print('  Done.')

    else:
        print('No pre-enriched metadata found. Loading fresh...')
        df_meta_raw = load_all_metadata(genus, species_list)
        print(f'Loaded {len(df_meta_raw)} total assemblies from metadata')

        df_meta_raw['acc_short'] = (
            df_meta_raw['Assembly Accession']
            .str.extract(r'(GC[AF]_\d+\.\d+)')[0]
        )
        df_meta = df_meta_raw[df_meta_raw['acc_short'].isin(all_passed_accs)].copy()
        df_meta = df_meta.drop(columns=['acc_short'])
        print(f'Filtered to {len(df_meta)} QC-passed assemblies')

        print(f'\nEnriching {len(df_meta)} assemblies with BioProject info from NCBI...')
        df_meta = classify.enrich_df_with_bioproject(
            df_meta, assembly_col='Assembly Accession', pause=0.35
        )
        print('BioProject enrichment complete!')

        print('\nEnriching with WGS Reference Titles from NCBI...')
        df_meta = classify.enrich_df_with_wgs_title(
            df_meta, wgs_col='WGS project accession', pause=0.35
        )
        print('WGS title enrichment complete!')

    # Standardize Species column
    species_map = {sp: species_display[sp] for sp in species_list}
    if 'Species' in df_meta.columns and df_meta['Species'].iloc[0] in species_map:
        df_meta['Species'] = df_meta['Species'].map(species_map)

    # Save enriched metadata
    enriched_out = os.path.join(results_base, 'qc_passed_metadata_enriched.csv')
    df_meta.to_csv(enriched_out, index=False)
    print(f'\nSaved enriched metadata to {enriched_out}')

    print(f'\nSpecies distribution:')
    print(df_meta['Species'].value_counts())

    for col in ['BioProject Description', 'WGS Reference Title']:
        if col in df_meta.columns:
            n = (df_meta[col].notna() & (df_meta[col] != '')).sum()
            print(f'\n{col}: {n}/{len(df_meta)} non-empty')

    # Build CombinedText for classification
    text_cols = [
        'Isolation Source', 'Assembly Submitter', 'Comment',
        'BioSample Description Title', 'Organism Infraspecific Names Strain',
        'BioProject Title', 'BioProject Description', 'BioProject Grants',
        'WGS Reference Title',
    ]
    text_cols = [c for c in text_cols if c in df_meta.columns]
    df_meta['CombinedText'] = (
        df_meta[text_cols].fillna('').astype(str).agg(' | '.join, axis=1)
    )

    return df_meta


def save_classified_outputs(df_classified, results_base):
    """Save classified metadata to CSV files and return analysis subset.

    Returns
    -------
    df_analysis : pd.DataFrame
        Subset excluding Unknown and Lab phenotypes.
    """
    # All classified (including Unknown)
    classified_all_path = os.path.join(results_base, 'qc_passed_classified_all.csv')
    df_classified.to_csv(classified_all_path, index=False)
    print(f'Saved FULL classified data (including Unknown) to:')
    print(f'  {classified_all_path}')

    # Known only (excluding Unknown and Lab)
    df_analysis = df_classified[
        ~df_classified['Phenotype'].isin(['Unknown', 'Lab'])
    ].copy()

    n_unknown = (df_classified['Phenotype'] == 'Unknown').sum()
    n_lab = (df_classified['Phenotype'] == 'Lab').sum()
    print(f'\nAnalysis dataset (excluding {n_unknown} Unknown + {n_lab} Lab phenotypes):')
    print(df_analysis['Phenotype'].value_counts().to_string())

    classified_known_path = os.path.join(results_base, 'qc_passed_classified_known.csv')
    df_analysis.to_csv(classified_known_path, index=False)
    print(f'\nSaved analysis data (excluding Unknown/Lab) to:')
    print(f'  {classified_known_path}')

    gwas_path = os.path.join(results_base, 'phenotype_classified_for_gwas.csv')
    df_analysis.to_csv(gwas_path, index=False)
    print(f'\nSaved pan-GWAS phenotype data to:')
    print(f'  {gwas_path}')

    print(f'\n{"=" * 60}')
    print('Available dataframes:')
    print(f'  df_classified: {len(df_classified)} assemblies (ALL, including Unknown/Lab)')
    print(f'  df_analysis:   {len(df_analysis)} assemblies (for pan-GWAS, no Unknown/Lab)')
    print('=' * 60)

    return df_analysis


def plot_phenotype_distribution(df_analysis, results_base, class_col='Phenotype',
                                fontsize=15, dpi=400,
                                pie_figsize=(10, 11),
                                bar_figsize=(13, 9),
                                pie_legend_ncol=1,
                                bar_ylim_factor=1.45):
    """Plot phenotype distribution as TWO SEPARATE figures (pie and per-species bar).

    Parameters
    ----------
    fontsize : int
        Base font size (title is +2, legend labels are -1).
    dpi : int
        Resolution of saved PNGs.
    pie_figsize, bar_figsize : tuple
        Figure sizes in inches for the pie and bar figures.
    pie_legend_ncol : int
        Number of columns in the pie legend (default 1 = tall single column).
    bar_ylim_factor : float
        Multiplier applied to the tallest bar to size the y-axis. Higher values
        leave more headroom for the bar-plot legend (default 1.45).

    Returns
    -------
    (fig_pie, fig_bar) : matplotlib Figure pair
    """
    UNIFIED_COLORS = {
        'Clinical':          '#c0392b',
        'Industrial-origin': '#2980b9',
        'Culture-derived':   '#8e44ad',
        'Human-pathogenic':  '#e74c3c',
        'Animal-pathogenic': '#7d3c98',
        'Plant-pathogenic':  '#e67e22',
        'Industrial-trait':  '#3498db',
        'Lab':               '#95a5a6',
        'Environmental':     '#27ae60',
        'Unknown':           '#bdc3c7',
    }

    FS = fontsize
    pheno_counts = df_analysis[class_col].value_counts()
    colors = [UNIFIED_COLORS.get(c, '#7f8c8d') for c in pheno_counts.index]

    # ============================== FIGURE 1: PIE ==============================
    fig_pie, ax_pie = plt.subplots(figsize=pie_figsize)

    n_pheno = len(pheno_counts)
    legend_rows = -(-n_pheno // pie_legend_ncol)  # ceiling division
    legend_height = 0.10 + 0.055 * legend_rows * (FS / 15.0)
    legend_height = min(legend_height, 1.4)
    pie_radius = max(0.35, 1.0 - legend_height / 2 - 0.08)
    pie_cy = 1.0 - pie_radius - 0.05
    pie_center = (0.0, pie_cy)
    legend_y = pie_cy - pie_radius - 0.05

    wedges, _ = ax_pie.pie(
        pheno_counts.values, labels=None, autopct=None,
        colors=colors, startangle=90,
        wedgeprops=dict(linewidth=0.5, edgecolor='white'),
        radius=pie_radius, center=pie_center,
    )
    ax_pie.set_xlim(-1.2, 1.2)
    ax_pie.set_ylim(legend_y - legend_height - 0.05, 1.1)
    ax_pie.set_aspect('equal')
    total = pheno_counts.values.sum()
    cx, cy = pie_center
    for wedge, val in zip(wedges, pheno_counts.values):
        pct = 100 * val / total
        ang = (wedge.theta1 + wedge.theta2) / 2
        rad = np.radians(ang)
        r = 0.55 * pie_radius if pct >= 5.0 else 1.15 * pie_radius
        ax_pie.text(cx + r * np.cos(rad), cy + r * np.sin(rad),
                    f'{pct:.1f}%', ha='center', va='center',
                    fontsize=FS, fontweight='bold',
                    color='white' if pct >= 5.0 else 'black')

    ax_pie.legend(wedges, pheno_counts.index,
                  title='Phenotype', title_fontsize=FS, fontsize=FS - 1,
                  loc='upper center', bbox_to_anchor=(0.0, legend_y),
                  bbox_transform=ax_pie.transData,
                  ncol=pie_legend_ncol, framealpha=0.95)
    ax_pie.set_title(f'Phenotype Distribution for Pan-GWAS\n(n={len(df_analysis)}, excl. Unknown/Lab)',
                     fontsize=FS + 2, fontweight='bold')

    fig_pie.tight_layout()
    pie_path = os.path.join(results_base, 'phenotype_distribution_pie.png')
    fig_pie.savefig(pie_path, dpi=dpi, bbox_inches='tight')
    fig_pie.savefig(pie_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f'Saved to {pie_path}  (+ .pdf)')

    # ============================== FIGURE 2: BAR ==============================
    fig_bar, ax_bar = plt.subplots(figsize=bar_figsize)
    species_pheno = df_analysis.groupby(['Species', class_col]).size().unstack(fill_value=0)
    species_pheno = species_pheno[[c for c in pheno_counts.index if c in species_pheno.columns]]
    colors_bar = [UNIFIED_COLORS.get(c, '#7f8c8d') for c in species_pheno.columns]
    species_pheno.plot(kind='bar', stacked=True, ax=ax_bar, color=colors_bar, edgecolor='white')
    ax_bar.set_title('Phenotype Distribution by Species', fontsize=FS + 2, fontweight='bold')
    ax_bar.set_xlabel('Species', fontsize=FS)
    ax_bar.set_ylabel('Number of Strains', fontsize=FS)
    ax_bar.legend(title='Phenotype', loc='upper right',
                  fontsize=FS - 1, title_fontsize=FS, framealpha=0.95)
    ax_bar.tick_params(axis='x', rotation=45, labelsize=FS)
    ax_bar.tick_params(axis='y', labelsize=FS)
    max_total = species_pheno.sum(axis=1).max()
    ax_bar.set_ylim(0, max_total * bar_ylim_factor)
    for x, total_val in enumerate(species_pheno.sum(axis=1).values):
        ax_bar.text(x, total_val + max_total * 0.02, f'n={int(total_val)}',
                    ha='center', va='bottom', fontsize=FS, fontweight='bold')

    fig_bar.tight_layout()
    bar_path = os.path.join(results_base, 'phenotype_distribution_bar.png')
    fig_bar.savefig(bar_path, dpi=dpi, bbox_inches='tight')
    fig_bar.savefig(bar_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f'Saved to {bar_path}  (+ .pdf)')

    return fig_pie, fig_bar


def print_nb0_summary(df_classified, df_analysis, species_list, results_base):
    """Print final NB0 cross-tabulation summary tables."""
    print('Genome counts per species per phenotype:')
    print('=' * 70)
    crosstab = pd.crosstab(
        df_classified['Species'],
        df_classified['Phenotype'],
        margins=True,
    )

    print('\nGenome counts per species per provenance:')
    print('=' * 70)
    crosstab_prov = pd.crosstab(
        df_classified['Species'],
        df_classified['Provenance'],
        margins=True,
    )

    print(f'\n{"=" * 70}')
    print('NB0 Data Preparation complete.')
    print(f'All outputs saved to: {results_base}')
    print(f'  - qc_passed_metadata_enriched.csv')
    print(f'  - qc_passed_classified_all.csv  ({len(df_classified)} genomes)')
    print(f'  - qc_passed_classified_known.csv ({len(df_analysis)} genomes)')
    print(f'  - phenotype_classified_for_gwas.csv')
    print(f'  - phenotype_distribution.png')
    print(f'  - strain_world_map.png')
    for sp in species_list:
        print(f'  - {sp}/phenotype_data.tsv')
        print(f'  - {sp}/phenotypes/pheno_*.tsv')
    print('=' * 70)

    return crosstab, crosstab_prov

