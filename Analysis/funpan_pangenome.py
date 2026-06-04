#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_pangenome.py
===================
Pangenome construction, classification, enrichment analysis, Heap's law,
SNP PCA/GRM, and BGC/GCF matrix analysis for FunPan notebooks.

Extracted from Main_Analysis_Functions.py (lines 881+).
"""

from funpan_utils import *
from funpan_utils import _split_multi_terms, _CAZY_FAM_RE, _MULTI_SEP_RE, _extract_cazy_families

# =============================================================================
# ORTHOFINDER I/O FUNCTIONS (from functions_for_analysis_1_2.py)
# =============================================================================

def load_orthogroups_long(orthogroups_tsv: str, species_name: str = None) -> pd.DataFrame:
    """Read OrthoFinder Orthogroups.tsv and reshape to long format."""
    og = pd.read_csv(orthogroups_tsv, sep="\t")

    rename = {}
    for c in og.columns:
        if c == "Orthogroup":
            continue
        m = ACC_RE.search(c)
        if m:
            rename[c] = m.group(1)
    og = og.rename(columns=rename)

    long = (
        og.melt(id_vars="Orthogroup", var_name="Assembly Accession", value_name="Protein_IDs")
          .dropna(subset=["Protein_IDs"])
    )
    long["Protein_ID"] = long["Protein_IDs"].astype(str).str.split(r",\s*")
    long = long.explode("Protein_ID", ignore_index=True)[["Orthogroup", "Assembly Accession", "Protein_ID"]]
    long["Protein_ID"] = long["Protein_ID"].str.strip()
    long = long[long["Protein_ID"] != ""].reset_index(drop=True)

    if species_name:
        long["Species"] = species_name

    return long


def consensus_terms_for_og(df_og: pd.DataFrame, column: str,
                           min_fraction: float = 0.5, min_support: int = 1,
                           max_terms: int = 1) -> str:
    """Single most-common consensus term for an annotation column within one
    Orthogroup. Matches the Theme convention (Counter.most_common(1)[0][0])
    so that an OG appears in exactly one term-bucket per layer downstream.

    Behaviour:
      - if any term passes (count >= min_support, fraction >= min_fraction)
        return the highest-count one (tie-break by fraction).
      - else return the most-common term as fallback (count >= 1).
    """
    if column not in df_og.columns:
        return ""

    vals = df_og[column].dropna().astype(str)
    if vals.empty:
        return ""

    per_protein_terms = [set(_split_multi_terms(v)) for v in vals]
    all_terms = set().union(*per_protein_terms) if per_protein_terms else set()
    if not all_terms:
        return ""

    counts = {}
    for term in all_terms:
        c = sum(1 for tset in per_protein_terms if term in tset)
        counts[term] = c

    n_proteins = len(per_protein_terms)
    term_stats = [(t, c, c / n_proteins) for t, c in counts.items()]

    filtered = [
        (t, c, f) for (t, c, f) in term_stats
        if c >= min_support and f >= min_fraction
    ]

    if not filtered:
        t_best, _, _ = max(term_stats, key=lambda x: x[1])
        return t_best

    filtered.sort(key=lambda x: (-x[1], -x[2]))
    return filtered[0][0]


def build_og_consensus_table(long_annot: pd.DataFrame, unknown_label: str = "Unannotated") -> pd.DataFrame:
    """Build one consensus-annotation row per Orthogroup."""
    og_groups = long_annot.groupby("Orthogroup", group_keys=False)

    consensus_cols = []
    combined_cols = ["Description", "GOs", "PFAMs", "CAZy", "EC"]
    for c in combined_cols:
        if c in long_annot.columns:
            consensus_cols.append(c)

    for c in EGGNOG_KEEP:
        if c not in ["Description_EggNog", "GOs_EggNog", "PFAMs_EggNog", "CAZy_EggNog", "EC_EggNog"]:
            if c in long_annot.columns:
                consensus_cols.append(c)

    for c in DBCAN_KEEP:
        if c not in ["CAZy_dbCAN", "EC_dbCAN"]:
            if c in long_annot.columns:
                consensus_cols.append(c)

    for c in INTERPRO_KEEP:
        if c not in ["GOs_Interpro", "PFAMs_Interpro", "EC_Interpro", "Description_Interpro"]:
            if c in long_annot.columns:
                consensus_cols.append(c)

    for c in SIGNALP_KEEP:
        if c in long_annot.columns:
            consensus_cols.append(c)

    rows = []
    for og, df_og in og_groups:
        row = {"Orthogroup": og}
        any_annotation = False

        for col in consensus_cols:
            cons = consensus_terms_for_og(df_og, col)
            if cons:
                any_annotation = True
            row[f"Consensus_{col}"] = cons

        if not any_annotation:
            row["Consensus_Description"] = unknown_label

        rows.append(row)

    og_consensus = pd.DataFrame(rows)
    return og_consensus


def build_ortho_annot_table(base_dir: str, species_name: str, out_path: str) -> pd.DataFrame:
    """Build the *long* per-protein OrthoFinder + annotation table."""
    if os.path.isfile(out_path):
        print(f"[i] Found existing {out_path} — loading long per-protein table.")
        df = pd.read_csv(out_path, sep="\t", dtype=str).fillna("")
        return df

    eggnog_dir = os.path.join(base_dir, "eggnog_output")
    dbcan_dir = os.path.join(base_dir, "dbcan_output")
    interpro_dir = os.path.join(base_dir, "interproscan_output")
    signalp_dir = os.path.join(base_dir, "signalp_output")

    og_paths = glob.glob(os.path.join(base_dir, "orthofinder_output", "*", "Orthogroups", "Orthogroups.tsv"))
    if not og_paths:
        raise FileNotFoundError("Could not find Orthogroups.tsv under orthofinder_output/")
    og_path = og_paths[0]
    print(f"[i] Using Orthogroups: {og_path}")

    long = load_orthogroups_long(og_path, species_name=species_name)

    eggnog = load_emapper_for_accessions(eggnog_dir, long["Assembly Accession"].unique()) \
             if os.path.isdir(eggnog_dir) else pd.DataFrame()
    dbcan = load_dbcan_for_accessions(dbcan_dir, long["Assembly Accession"].unique()) \
            if os.path.isdir(dbcan_dir) else pd.DataFrame()
    interpro = load_interproscan_for_accessions(interpro_dir, long["Assembly Accession"].unique()) \
               if os.path.isdir(interpro_dir) else pd.DataFrame()
    signalp = load_signalp_for_accessions(signalp_dir, long["Assembly Accession"].unique()) \
              if os.path.isdir(signalp_dir) else pd.DataFrame()

    merged = long.copy()

    if not eggnog.empty:
        merged = merged.merge(eggnog, on=["Protein_ID", "Assembly Accession"], how="left")
    else:
        for c in EGGNOG_KEEP:
            if c not in merged.columns:
                merged[c] = ""

    if not dbcan.empty:
        merged = merged.merge(dbcan, on=["Protein_ID", "Assembly Accession"], how="left")
    else:
        for c in DBCAN_KEEP:
            if c not in merged.columns:
                merged[c] = ""

    if not interpro.empty:
        merged = merged.merge(interpro, on=["Protein_ID", "Assembly Accession"], how="left")
    else:
        for c in INTERPRO_KEEP:
            if c not in merged.columns:
                merged[c] = ""

    if not signalp.empty:
        merged = merged.merge(signalp, on=["Protein_ID", "Assembly Accession"], how="left")
    else:
        for c in SIGNALP_KEEP:
            if c not in merged.columns:
                merged[c] = ""

    for col in ["GOs_EggNog", "GOs_Interpro", "PFAMs_EggNog", "PFAMs_Interpro",
                "CAZy_EggNog", "CAZy_dbCAN", "EC_EggNog", "EC_dbCAN", "EC_Interpro",
                "Description_EggNog", "Description_Interpro"]:
        if col not in merged.columns:
            merged[col] = ""

    merged["GOs"] = merged.apply(
        lambda row: combine_annotations(row.get("GOs_EggNog", ""), row.get("GOs_Interpro", "")), axis=1
    )
    merged["PFAMs"] = merged.apply(
        lambda row: combine_annotations(row.get("PFAMs_EggNog", ""), row.get("PFAMs_Interpro", "")), axis=1
    )
    merged["CAZy"] = merged["CAZy_dbCAN"].fillna("").astype(str)
    merged["EC"] = merged.apply(
        lambda row: combine_annotations(row.get("EC_EggNog", ""), row.get("EC_dbCAN", ""), row.get("EC_Interpro", "")),
        axis=1
    )

    def merge_description(row):
        eggnog_desc = str(row.get("Description_EggNog", "")).strip()
        interpro_desc = str(row.get("Description_Interpro", "")).strip()
        if eggnog_desc and eggnog_desc not in ["", "-", "nan", "Unannotated"]:
            return eggnog_desc
        elif interpro_desc and interpro_desc not in ["", "-", "nan"]:
            return interpro_desc
        return ""

    merged["Description"] = merged.apply(merge_description, axis=1)

    merged.to_csv(out_path, sep="\t", index=False)
    print(f"[✓] Saved long per-protein OG annotation table: {out_path}")

    return merged


# =============================================================================
# PAV & PANGENOME CLASSIFICATION (from functions_for_analysis_1_2.py)
# =============================================================================

def build_pav(ortho_annot_df: pd.DataFrame) -> pd.DataFrame:
    """Build Presence/Absence matrix (Orthogroup × Assembly Accession)."""
    df = ortho_annot_df.copy()
    if df["Protein_ID"].astype(str).str.contains(",").any():
        df["Protein_ID"] = df["Protein_ID"].astype(str).str.split(r",\s*")
        df = df.explode("Protein_ID", ignore_index=True)

    pav = df.pivot_table(
        index='Orthogroup',
        columns='Assembly Accession',
        values='Protein_ID',
        aggfunc=lambda x: 1 if pd.notnull(x).any() else 0
    ).fillna(0).astype(int)
    return pav


def determine_core_and_rare_thresholds(pav_matrix: pd.DataFrame,
                                       travel_fraction=0.90):
    """S-curve (Hyun/Palsson) inflection point method for pangenome thresholds.

    Based on Hyun et al. 2022 (mSphere) and Chauhan et al. 2024 (E. coli pangenome).

    Fits the gene frequency histogram with two power functions whose sum is
    U-shaped (rare peak at x=1, core peak at x=N):

        P(x) = c1 * x^(-a1) + c2 * (N+1-x)^(-a2)

    The inflection point x* is the minimum of this U-curve. Thresholds are
    placed travel_fraction of the way from x* to each endpoint:

        Core  >=  x* + travel_fraction * (N - x*)
        Rare  <=  x* - travel_fraction * (x* - 1)

    To prevent the core peak (which is typically orders of magnitude taller
    than the rare peak) from dominating the fit, the histogram is fit in
    log-space.

    Parameters
    ----------
    pav_matrix : pd.DataFrame
        Presence/absence matrix (orthogroups x genomes).
    travel_fraction : float
        Fraction of distance from inflection point to endpoint (default 0.90,
        matching Hyun/Chauhan convention).

    Returns
    -------
    core_rate, core_n, rare_rate, rare_n : float, int, float, int
    """
    from scipy.optimize import curve_fit

    num_strains = pav_matrix.shape[1]
    gene_presence_counts = pav_matrix.sum(axis=1)

    # Histogram of gene presence counts
    bins = np.arange(1, num_strains + 1, dtype=float)
    freq_counts = np.array([int((gene_presence_counts == t).sum()) for t in bins],
                           dtype=float)
    # Avoid log(0) by adding 1
    log_freq = np.log(freq_counts + 1.0)

    def two_power(x, c1, a1, c2, a2):
        return c1 * np.power(x, -a1) + c2 * np.power(num_strains + 1 - x, -a2)

    def log_two_power(x, c1, a1, c2, a2):
        return np.log(two_power(x, c1, a1, c2, a2) + 1.0)

    inflection_x = None
    try:
        # Fit in log space so core and rare peaks contribute equally
        p0 = [max(freq_counts[0], 1.0), 1.0, max(freq_counts[-1], 1.0), 1.0]
        bounds = ([0, 0.01, 0, 0.01], [np.inf, 10, np.inf, 10])
        popt, _ = curve_fit(log_two_power, bins, log_freq,
                            p0=p0, bounds=bounds, maxfev=20000)

        # Find minimum of the fitted U-curve in the interior
        x_fine = np.linspace(1.0, float(num_strains), 2000)
        y_fit = two_power(x_fine, *popt)
        inflection_x = float(x_fine[np.argmin(y_fit)])
    except (RuntimeError, ValueError, TypeError):
        pass

    if inflection_x is None:
        # Fallback: median of the histogram (ignoring core peak)
        interior = slice(max(1, int(0.05 * num_strains)),
                         max(2, int(0.95 * num_strains)))
        inflection_x = float(bins[interior][int(np.argmin(freq_counts[interior]))])

    # Travel fraction from inflection point toward each endpoint
    core_n = int(np.round(inflection_x + travel_fraction * (num_strains - inflection_x)))
    rare_n = int(np.round(inflection_x - travel_fraction * (inflection_x - 1)))
    core_n = min(max(core_n, 1), num_strains)
    rare_n = max(1, min(rare_n, num_strains))

    core_rate = core_n / num_strains
    rare_rate = rare_n / num_strains
    return core_rate, core_n, rare_rate, rare_n


def classify_pangenome(ortho_annot_df: pd.DataFrame, pav: pd.DataFrame,
                       core_n: int, rare_n: int) -> pd.DataFrame:
    """Add Pangenome_Class (Core/Accessory/Rare/Other) to the long table."""
    sums = pav.sum(axis=1)
    core_idx = set(sums[sums >= core_n].index)
    accessory_idx = set(sums[(sums < core_n) & (sums > rare_n)].index)
    rare_idx = set(sums[sums <= rare_n].index)

    def _cls(og):
        if og in core_idx:
            return "Core"
        if og in accessory_idx:
            return "Accessory"
        if og in rare_idx:
            return "Rare"
        return "Other"

    df = ortho_annot_df.copy()
    df["Pangenome_Class"] = df["Orthogroup"].map(_cls)
    return df


def plot_gene_freq_hist(pav: pd.DataFrame, core_n: int, rare_n: int, ax=None, annotate=True):
    """Plot gene frequency histogram with core/rare thresholds."""
    gene_freq = pav.sum(axis=1).values
    num_strains = pav.shape[1]
    bins = np.arange(1, num_strains + 2)
    hist, bin_edges = np.histogram(gene_freq, bins=bins)

    core_cnt = int((gene_freq >= core_n).sum())
    rare_cnt = int((gene_freq <= rare_n).sum())
    accessory_cnt = int(((gene_freq < core_n) & (gene_freq > rare_n)).sum())

    core_rate = core_n / num_strains
    rare_rate = rare_n / num_strains

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    ax.bar(bin_edges[:-1], hist, width=0.75, alpha=0.6, label="Gene frequency histogram")
    ax.axvline(core_n, color="green", ls="--", lw=2)
    ax.axvline(rare_n, color="purple", ls="--", lw=2)
    ax.axvspan(core_n, num_strains, color="green", alpha=0.08)
    ax.axvspan(1, rare_n, color="purple", alpha=0.08)

    ax.set_xlabel("Number of genomes gene is present in")
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene presence frequency across genomes")
    ax.set_xticks(np.arange(1, num_strains + 1, max(1, num_strains // 15)))

    if annotate:
        ymax = hist.max() if len(hist) else 1
        ytop = ymax * 1.1

        ax.text(core_n, ytop, f"Core ≥ {core_n}\n({core_rate:.0%})",
                ha="center", va="bottom", color="green", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="green", alpha=0.6))
        ax.text(rare_n, ytop, f"Rare ≤ {rare_n}\n({rare_rate:.0%})",
                ha="center", va="bottom", color="purple", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="purple", alpha=0.6))

        summary_txt = (f"Core genes: {core_cnt}\n"
                       f"Accessory:  {accessory_cnt}\n"
                       f"Rare genes: {rare_cnt}")
        ax.text(0.02, 0.95, summary_txt, transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.5", alpha=0.9))

    return ax


# =============================================================================
# ENRICHMENT FUNCTIONS (from notebook cells 15, 16)
# =============================================================================

def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (nan-safe)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    p_clean = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(p_clean)
    q = np.empty_like(p_clean, dtype=float)
    prev = 1.0
    for i, idx in enumerate(order[::-1]):
        rank = n - i
        val = min(prev, p_clean[idx] * n / rank)
        q[idx] = val
        prev = val
    q[np.isnan(p)] = np.nan
    return q


def build_og_consensus_with_pangenome(long_annot, pav, core_n, rare_n):
    """
    Build orthogroup-level consensus annotations with pangenome classification.

    Uses ONLY the combined annotation columns (GOs, PFAMs, CAZy, EC, Description)
    for consensus calculation. Source-specific columns are excluded.

    Also extracts:
    - Proteases (from PFAM peptidase families)
    - Transporters (from PFAM transporter families + KEGG_TC)
    - Secretome markers (from SignalP + secretion GO terms)
    """
    # First, get pangenome classification for each orthogroup
    sums = pav.sum(axis=1)
    core_idx = set(sums[sums >= core_n].index)
    accessory_idx = set(sums[(sums < core_n) & (sums > rare_n)].index)
    rare_idx = set(sums[sums <= rare_n].index)

    def _cls(og):
        if og in core_idx: return "Core"
        if og in accessory_idx: return "Accessory"
        if og in rare_idx: return "Rare"
        return "Other"

    # Build consensus per orthogroup using COMBINED columns only
    consensus_cols = [
        "GOs",           # Combined: GOs_EggNog + GOs_Interpro
        "PFAMs",         # Combined: PFAMs_EggNog + PFAMs_Interpro
        "CAZy",          # Only from CAZy_dbCAN
        "EC",            # Combined: EC_EggNog + EC_dbCAN + EC_Interpro
        "Description",   # EggNog priority, Interpro fallback
        "COG_category",  # Single source (EggNog only)
        "KEGG_ko",       # Single source (EggNog only)
        "KEGG_Pathway",  # Single source (EggNog only)
        "KEGG_TC",       # Transporter Classification (EggNog)
        "interpro_IPR",  # Single source (Interpro only)
        "dbcan_Substrate",  # CAZyme substrate specificity
    ]

    # Filter to columns that exist
    available_cols = [c for c in consensus_cols if c in long_annot.columns]

    og_groups = long_annot.groupby("Orthogroup")

    rows = []
    for og, df_og in og_groups:
        row = {
            "Orthogroup": og,
            "Pangenome_Class": _cls(og),
            "n_genomes": pav.loc[og].sum() if og in pav.index else 0,
            "n_proteins": len(df_og),
        }

        for col in available_cols:
            cons = consensus_terms_for_og(df_og, col, min_fraction=0.5, min_support=1, max_terms=5)
            row[col] = cons

        # Derive additional annotations
        pfams_str = row.get("PFAMs", "")
        gos_str = row.get("GOs", "")
        kegg_tc = row.get("KEGG_TC", "")

        # Proteases: check if any protease PFAMs present
        protease_hits = [pf for pf in PROTEASE_PFAMS.keys() if pf in pfams_str]
        row["is_protease"] = len(protease_hits) > 0
        row["protease_families"] = "; ".join(protease_hits) if protease_hits else ""

        # Transporters: check PFAM transporters + KEGG_TC
        transporter_hits = [pf for pf in TRANSPORTER_PFAMS.keys() if pf in pfams_str]
        has_kegg_tc = bool(kegg_tc and kegg_tc.strip() and kegg_tc != "-")
        row["is_transporter"] = len(transporter_hits) > 0 or has_kegg_tc
        row["transporter_families"] = "; ".join(transporter_hits) if transporter_hits else ""

        # Secretion GO terms
        secretion_hits = [go for go in SECRETION_GO_TERMS.keys() if go in gos_str]
        row["has_secretion_GO"] = len(secretion_hits) > 0

        rows.append(row)

    og_consensus = pd.DataFrame(rows)

    # Add SignalP consensus at OG level (majority vote)
    if "signalp_Prediction" in long_annot.columns:
        signalp_og = long_annot.groupby("Orthogroup").apply(
            lambda x: (x["signalp_Prediction"].fillna("").str.upper().str.contains("SP")).mean()
        ).to_dict()
        og_consensus["signalp_fraction"] = og_consensus["Orthogroup"].map(signalp_og).fillna(0)
        og_consensus["is_secreted_signalp"] = og_consensus["signalp_fraction"] >= 0.5
    else:
        og_consensus["signalp_fraction"] = 0
        og_consensus["is_secreted_signalp"] = False

    # Combined secretome: SignalP OR secretion GO terms
    og_consensus["is_secreted"] = og_consensus["is_secreted_signalp"] | og_consensus["has_secretion_GO"]

    return og_consensus


def fisher_enrichment_by_pangenome_class_og(og_consensus: pd.DataFrame,
                                             column: str,
                                             classes_order: list = None,
                                             min_count: int = 5,
                                             q_threshold: float = 0.05) -> pd.DataFrame:
    """
    For each term in `column`, test whether it's enriched in each pangenome class
    (Core/Accessory) using Fisher's exact test at the OG level.
    Note: Rare genes are excluded from enrichment analysis by default.
    """
    if classes_order is None:
        classes_order = ["Core", "Accessory"]

    if column not in og_consensus.columns:
        print(f"[WARN] Column '{column}' not in og_consensus. Skipping.")
        return pd.DataFrame()

    df = og_consensus[["Orthogroup", "Pangenome_Class", column]].dropna(subset=[column]).copy()
    df = df[df[column].astype(str).str.strip() != ""]
    df = df[df[column].astype(str).str.strip() != "-"]

    all_terms = []
    for val in df[column]:
        terms = _split_multi_terms(str(val))
        all_terms.extend(terms)

    term_counts = Counter(all_terms)
    common_terms = [t for t, c in term_counts.items() if c >= min_count]

    # Drop uninformative terms (COG-S/R, Pfam DUFs, GO root nodes, EC -.-.-.-)
    # Description-based filters (interpro / KEGG_ko / Pfam) are best-effort:
    # we only have term IDs here, not descriptions, so the description-based
    # cleanup happens at the per-term level downstream where descriptions
    # are joined in.
    from funpan_utils import is_informative_term
    common_terms = [t for t in common_terms if is_informative_term(column, t)]

    if not common_terms:
        return pd.DataFrame()

    results = []
    N_total = df["Orthogroup"].nunique()

    for term in common_terms:
        term_mask = df[column].astype(str).apply(lambda x: term in _split_multi_terms(x))
        ogs_with_term = set(df.loc[term_mask, "Orthogroup"])
        n_term = len(ogs_with_term)

        for pangenome_class in classes_order:
            class_mask = df["Pangenome_Class"] == pangenome_class
            ogs_in_class = set(df.loc[class_mask, "Orthogroup"])
            n_class = len(ogs_in_class)

            a = len(ogs_with_term & ogs_in_class)
            b = n_class - a
            c = n_term - a
            d = N_total - (a + b + c)

            if a < 1 or d < 0:
                continue

            try:
                # Use two-sided test to detect both enrichment AND depletion
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            except Exception:
                continue

            expected = (n_term * n_class) / N_total if N_total > 0 else 0
            fold_enrichment = a / expected if expected > 0 else 0

            # Determine direction: enriched (OR > 1) or depleted (OR < 1)
            direction = "enriched" if odds_ratio > 1 else "depleted" if odds_ratio < 1 else "none"

            results.append({
                "Term": term,
                "Pangenome_Class": pangenome_class,
                "Count_in_Class": a,
                "Count_Total": n_term,
                "Class_Size": n_class,
                "Odds_Ratio": odds_ratio,
                "Fold_Enrichment": fold_enrichment,
                "p_value": p_value,
                "Direction": direction,
            })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)
    result_df["q_value"] = multipletests(result_df["p_value"], method="fdr_bh")[1]
    result_df["Significant"] = result_df["q_value"] < q_threshold

    return result_df.sort_values(["Pangenome_Class", "q_value"])


def run_pangenome_enrichment_og(og_consensus: pd.DataFrame,
                                 columns: list = None,
                                 classes_order: list = None,
                                 min_count: int = 5,
                                 q_threshold: float = 0.05) -> dict:
    """Run enrichment analysis for multiple annotation columns.
    Note: Rare genes are excluded from enrichment analysis by default."""
    if columns is None:
        columns = [
            "COG_category",
            "KEGG_TC",
            "CAZy",
            "dbcan_Substrate",
            "GOs",
            "PFAMs",
            "EC",
            "KEGG_ko",
            "KEGG_Pathway",
            "interpro_IPR",
        ]

    if classes_order is None:
        classes_order = ["Core", "Accessory"]

    results = {}

    for col in columns:
        if col not in og_consensus.columns:
            print(f"[SKIP] Column '{col}' not found")
            continue

        enr = fisher_enrichment_by_pangenome_class_og(
            og_consensus, col, classes_order, min_count, q_threshold
        )

        if enr.empty:
            print(f"[SKIP] No enrichment results for '{col}'")
            continue

        n_sig = enr["Significant"].sum()
        print(f"[OK] {col}: {len(enr)} tests, {n_sig} significant (FDR < {q_threshold})")
        results[col] = enr

    return results


def add_cog_descriptions(df: pd.DataFrame) -> pd.DataFrame:
    """Add COG category descriptions to enrichment results."""
    if df.empty:
        return df
    df = df.copy()
    df["Description"] = df["Term"].map(COG_DESCRIPTIONS).fillna("Unknown")
    return df


def plot_enrichment_heatmap_single_species(enr_df: pd.DataFrame, title: str = "",
                                            max_terms: int = 25,
                                            figsize: tuple = (10, 10),
                                            q_threshold: float = 0.05) -> tuple:
    """
    Plot heatmap of enrichment showing log2(odds_ratio) for Core vs Accessory.
    Uses RdBu_r colormap: Red=Enriched, Blue=Depleted.
    Shows significance stars (*, **, ***).
    ONLY shows significant terms (q < threshold).
    Filters out COG S (function unknown).
    """
    if enr_df.empty:
        print("[WARN] No enrichment data to plot")
        return None, None

    df = enr_df.copy()

    # Filter out COG S (function unknown)
    df = df[df["Term"] != "S"]

    # Only work with significant terms
    sig_df = df[df["q_value"] < q_threshold]

    if sig_df.empty:
        print(f"[INFO] No significant terms (q < {q_threshold}) to display in heatmap")
        return None, None

    # Get top terms sorted by maximum absolute log2(OR) - most extreme effects first.
    # Clip inf (perfect separation) to a finite cap so they don't all tie at the top
    # and push moderately-strong-but-finite effects off the chart.
    LOG2_SORT_CLIP = 8.0
    def max_abs_log2or(term_df):
        ors = term_df["Odds_Ratio"].values
        log2ors = []
        for x in ors:
            if not np.isfinite(x) or x <= 0:
                log2ors.append(LOG2_SORT_CLIP if (x == np.inf or x > 1) else -LOG2_SORT_CLIP)
            else:
                log2ors.append(float(np.clip(np.log2(x), -LOG2_SORT_CLIP, LOG2_SORT_CLIP)))
        return max(abs(v) for v in log2ors) if log2ors else 0

    term_effects = sig_df.groupby("Term").apply(max_abs_log2or)
    top_terms = term_effects.sort_values(ascending=False).head(max_terms).index.tolist()

    if not top_terms:
        print("[INFO] No terms to display")
        return None, None

    # For COG categories, add descriptions
    def add_cog_desc(term):
        if term in COG_DESCRIPTIONS:
            return f"{term}: {COG_DESCRIPTIONS[term]}"
        return term

    top_terms_display = [add_cog_desc(t) for t in top_terms]

    classes = ["Core", "Accessory"]

    # Build matrices for values and p-values
    matrix = np.zeros((len(top_terms), len(classes)))
    pval_matrix = np.ones((len(top_terms), len(classes)))

    LOG2_CLIP = 8.0  # cap for inf/extreme ORs; OR=256 is already biologically extreme
    for i, term in enumerate(top_terms):
        for j, cls in enumerate(classes):
            match = df[(df["Term"] == term) & (df["Pangenome_Class"] == cls)]
            if len(match) > 0:
                or_val = match["Odds_Ratio"].values[0]
                if not np.isfinite(or_val) or or_val <= 0:
                    # inf  -> +cap (perfectly enriched);  0 -> -cap (perfectly depleted)
                    log2_val = LOG2_CLIP if (or_val == np.inf or or_val > 1) else -LOG2_CLIP
                else:
                    log2_val = float(np.clip(np.log2(or_val), -LOG2_CLIP, LOG2_CLIP))
                matrix[i, j] = log2_val
                pval_matrix[i, j] = match["q_value"].values[0]

    fig, ax = plt.subplots(figsize=figsize)

    # Symmetric color scale — use the data's actual range so extremes don't get
    # clipped (the per-cell numeric annotation makes the exact value readable
    # regardless). Floor at 1 so a single weak hit still has visible colour.
    vmax = max(abs(matrix.min()), abs(matrix.max()))
    vmax = max(vmax, 1.0)

    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    # Annotate each cell with its log2OR value + significance stars
    for i in range(len(top_terms)):
        for j in range(len(classes)):
            log2or = matrix[i, j]
            q = pval_matrix[i, j]
            if q < 0.001:
                stars = "***"
            elif q < 0.01:
                stars = "**"
            elif q < 0.05:
                stars = "*"
            else:
                stars = ""
            # Use white text on saturated cells, black on pale cells
            text_color = "white" if abs(log2or) > 0.6 * vmax else "black"
            label = f"{log2or:+.1f}\n{stars}" if stars else (f"{log2or:+.1f}" if log2or != 0 else "")
            if label:
                ax.text(j, i, label, ha="center", va="center",
                        fontsize=8, color=text_color, fontweight="bold")

    ax.set_xticks(np.arange(len(classes)))
    ax.set_xticklabels(classes, fontsize=12, fontweight="bold")
    ax.set_yticks(np.arange(len(top_terms)))
    ax.set_yticklabels(top_terms_display, fontsize=9)

    ax.set_xlabel("Pangenome Class", fontsize=11)
    n_sig = len(sig_df["Term"].unique())
    ax.set_title(f"{title}\n({n_sig} significant terms, * q<0.05, ** q<0.01, *** q<0.001)",
                 fontsize=12, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log2(Odds Ratio)\nRed=Enriched, Blue=Depleted", fontsize=10)

    plt.tight_layout()



def plot_enrichment_bars_single_species(enr_df: pd.DataFrame, title: str = "",
                                         max_terms: int = 12,
                                         q_threshold: float = 0.05,
                                         figsize: tuple = (16, 8)) -> tuple:
    """
    Bar plot showing top enriched AND depleted terms for Core vs Accessory side by side.
    Uses signed -log10(q-value): positive = enriched, negative = depleted.
    Green = enriched, Red = depleted.
    Filters out COG S (function unknown).
    """
    if enr_df.empty:
        print("[WARN] No enrichment data to plot")
        return None, None

    df = enr_df.copy()

    # Filter out COG S (function unknown)
    df = df[df["Term"] != "S"]

    # For COG categories, add descriptions
    def add_cog_desc(term):
        if term in COG_DESCRIPTIONS:
            return f"{term}: {COG_DESCRIPTIONS[term]}"
        return term

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    for idx, pangenome_class in enumerate(["Core", "Accessory"]):
        ax = axes[idx]

        class_df = df[df["Pangenome_Class"] == pangenome_class].copy()
        sig_df = class_df[class_df["q_value"] < q_threshold].copy()

        if len(sig_df) == 0:
            ax.text(0.5, 0.5, f"No significant terms\n(FDR < {q_threshold})",
                   ha="center", va="center", transform=ax.transAxes, fontsize=12)
            ax.set_title(f"{pangenome_class}", fontsize=14, fontweight="bold",
                        color="#2ca02c" if pangenome_class == "Core" else "#e67e22")
            ax.set_xlim(-1, 1)
            continue

        # Get top enriched (OR > 1) and depleted (OR < 1) terms separately
        enriched = sig_df[sig_df["Odds_Ratio"] > 1].nsmallest(max_terms // 2, "q_value")
        depleted = sig_df[sig_df["Odds_Ratio"] < 1].nsmallest(max_terms // 2, "q_value")
        combined = pd.concat([enriched, depleted]).copy()

        if len(combined) == 0:
            ax.text(0.5, 0.5, f"No significant terms\n(FDR < {q_threshold})",
                   ha="center", va="center", transform=ax.transAxes, fontsize=12)
            ax.set_title(f"{pangenome_class}", fontsize=14, fontweight="bold",
                        color="#2ca02c" if pangenome_class == "Core" else "#e67e22")
            continue

        # Calculate signed -log10(q-value): positive for enriched, negative for depleted
        combined["signed_log_q"] = -np.log10(combined["q_value"] + 1e-300)
        combined.loc[combined["Odds_Ratio"] < 1, "signed_log_q"] *= -1
        combined = combined.sort_values("signed_log_q")

        # Labels with COG descriptions
        labels = [add_cog_desc(t) for t in combined["Term"]]

        # Colors: green for enriched (positive), red for depleted (negative)
        colors = ["#d62728" if x < 0 else "#2ca02c" for x in combined["signed_log_q"]]

        ax.barh(range(len(combined)), combined["signed_log_q"], color=colors, alpha=0.8)
        ax.set_yticks(range(len(combined)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("-log10(q-value)\n← Depleted | Enriched →", fontsize=10)

        n_enr = len(enriched)
        n_dep = len(depleted)
        ax.set_title(f"{pangenome_class}\n({n_enr} enriched, {n_dep} depleted)",
                    fontsize=14, fontweight="bold",
                    color="#2ca02c" if pangenome_class == "Core" else "#e67e22")
        ax.axvline(0, color="black", linewidth=1)
        ax.grid(axis="x", alpha=0.3)

    plt.suptitle(f"{title} (FDR < {q_threshold})", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()



def summarize_enrichment_single_species(enr_results: dict, species: str) -> pd.DataFrame:
    """Create summary table of significant enrichments for a single species."""
    summary_rows = []

    for col, enr_df in enr_results.items():
        if enr_df.empty:
            continue

        sig_df = enr_df[enr_df["Significant"] == True]

        for pangenome_class in ["Core", "Accessory"]:
            class_sig = sig_df[sig_df["Pangenome_Class"] == pangenome_class]

            if class_sig.empty:
                continue

            top_terms = class_sig.nsmallest(5, "q_value")["Term"].tolist()

            summary_rows.append({
                "Species": species,
                "Category": col.replace("Consensus_", ""),
                "Pangenome_Class": pangenome_class,
                "N_Significant": len(class_sig),
                "Top_Terms": "; ".join(top_terms[:3]),
                "Min_q_value": class_sig["q_value"].min(),
                "Max_Fold_Enrichment": class_sig["Fold_Enrichment"].max(),
            })

    return pd.DataFrame(summary_rows)


def interpret_enrichment_results(all_enrichment_results: dict, og_consensus_tables: dict,
                                  species_list: list, q_threshold: float = 0.05):
    """
    Provide in-depth biological interpretation of enrichment results.
    Summarizes key findings for Core vs Accessory genome across all species.
    Includes both enriched AND depleted functional categories.
    """
    print("\n" + "="*80)
    print("IN-DEPTH INTERPRETATION OF ENRICHMENT RESULTS")
    print("="*80)

    for species in species_list:
        if species not in all_enrichment_results:
            continue

        print(f"\n{'='*60}")
        print(f"A. {species.upper()}")
        print(f"{'='*60}")

        species_results = all_enrichment_results[species]
        og_consensus = og_consensus_tables.get(species)

        if og_consensus is not None:
            n_core = (og_consensus["Pangenome_Class"] == "Core").sum()
            n_acc = (og_consensus["Pangenome_Class"] == "Accessory").sum()
            print(f"\nPangenome structure: {n_core} Core OGs, {n_acc} Accessory OGs")

        # Track key findings - both enriched AND depleted
        core_enriched = []
        core_depleted = []
        accessory_enriched = []
        accessory_depleted = []

        for col, enr_df in species_results.items():
            if enr_df.empty:
                continue

            sig = enr_df[enr_df["q_value"] < q_threshold]

            for cls in ["Core", "Accessory"]:
                # Enriched (OR > 1)
                cls_enriched = sig[(sig["Pangenome_Class"] == cls) & (sig["Odds_Ratio"] > 1)]
                if not cls_enriched.empty:
                    top = cls_enriched.nsmallest(3, "q_value")
                    terms = [(row["Term"], row["Odds_Ratio"], row["q_value"])
                             for _, row in top.iterrows()]
                    if cls == "Core":
                        core_enriched.append((col, terms))
                    else:
                        accessory_enriched.append((col, terms))

                # Depleted (OR < 1)
                cls_depleted = sig[(sig["Pangenome_Class"] == cls) & (sig["Odds_Ratio"] < 1)]
                if not cls_depleted.empty:
                    top = cls_depleted.nsmallest(3, "q_value")
                    terms = [(row["Term"], row["Odds_Ratio"], row["q_value"])
                             for _, row in top.iterrows()]
                    if cls == "Core":
                        core_depleted.append((col, terms))
                    else:
                        accessory_depleted.append((col, terms))

        # Core genome interpretation - ENRICHED
        print(f"\n📗 CORE GENOME - ENRICHED (over-represented in conserved genes):")
        if core_enriched:
            for col, terms in core_enriched:
                print(f"\n  {col}:")
                for term, odds, q in terms:
                    desc = COG_DESCRIPTIONS.get(term, "")
                    if desc:
                        print(f"    ↑ {term} ({desc}): OR={odds:.1f}, q={q:.2e}")
                    else:
                        print(f"    ↑ {term}: OR={odds:.1f}, q={q:.2e}")
        else:
            print("  No significant enrichments detected in core genome")

        # Core genome interpretation - DEPLETED
        print(f"\n📕 CORE GENOME - DEPLETED (under-represented in conserved genes):")
        if core_depleted:
            for col, terms in core_depleted:
                print(f"\n  {col}:")
                for term, odds, q in terms:
                    desc = COG_DESCRIPTIONS.get(term, "")
                    if desc:
                        print(f"    ↓ {term} ({desc}): OR={odds:.2f}, q={q:.2e}")
                    else:
                        print(f"    ↓ {term}: OR={odds:.2f}, q={q:.2e}")
        else:
            print("  No significant depletions detected in core genome")

        # Accessory genome interpretation - ENRICHED
        print(f"\n📙 ACCESSORY GENOME - ENRICHED (over-represented in variable genes):")
        if accessory_enriched:
            for col, terms in accessory_enriched:
                print(f"\n  {col}:")
                for term, odds, q in terms:
                    desc = COG_DESCRIPTIONS.get(term, "")
                    if desc:
                        print(f"    ↑ {term} ({desc}): OR={odds:.1f}, q={q:.2e}")
                    else:
                        print(f"    ↑ {term}: OR={odds:.1f}, q={q:.2e}")
        else:
            print("  No significant enrichments detected in accessory genome")

        # Accessory genome interpretation - DEPLETED
        print(f"\n📘 ACCESSORY GENOME - DEPLETED (under-represented in variable genes):")
        if accessory_depleted:
            for col, terms in accessory_depleted:
                print(f"\n  {col}:")
                for term, odds, q in terms:
                    desc = COG_DESCRIPTIONS.get(term, "")
                    if desc:
                        print(f"    ↓ {term} ({desc}): OR={odds:.2f}, q={q:.2e}")
                    else:
                        print(f"    ↓ {term}: OR={odds:.2f}, q={q:.2e}")
        else:
            print("  No significant depletions detected in accessory genome")

        # Biological insights
        print(f"\n🔬 BIOLOGICAL INSIGHTS:")

        # Check for secreted proteins
        if og_consensus is not None and "is_secreted" in og_consensus.columns:
            core_sec = og_consensus[(og_consensus["Pangenome_Class"] == "Core") &
                                    (og_consensus["is_secreted"] == True)]
            acc_sec = og_consensus[(og_consensus["Pangenome_Class"] == "Accessory") &
                                   (og_consensus["is_secreted"] == True)]
            core_total = (og_consensus["Pangenome_Class"] == "Core").sum()
            acc_total = (og_consensus["Pangenome_Class"] == "Accessory").sum()

            core_pct = 100 * len(core_sec) / core_total if core_total > 0 else 0
            acc_pct = 100 * len(acc_sec) / acc_total if acc_total > 0 else 0

            print(f"  • Secreted proteins: {core_pct:.1f}% of Core vs {acc_pct:.1f}% of Accessory")
            if acc_pct > core_pct * 1.2:
                print(f"    → Accessory genome enriched for secreted proteins (potential niche adaptation)")

        # Check for CAZymes
        if og_consensus is not None and "CAZy" in og_consensus.columns:
            core_cazy = og_consensus[(og_consensus["Pangenome_Class"] == "Core") &
                                     (og_consensus["CAZy"].fillna("").str.len() > 0)]
            acc_cazy = og_consensus[(og_consensus["Pangenome_Class"] == "Accessory") &
                                    (og_consensus["CAZy"].fillna("").str.len() > 0)]

            core_pct = 100 * len(core_cazy) / core_total if core_total > 0 else 0
            acc_pct = 100 * len(acc_cazy) / acc_total if acc_total > 0 else 0

            print(f"  • CAZymes: {core_pct:.1f}% of Core vs {acc_pct:.1f}% of Accessory")

        # Check for transporters
        if og_consensus is not None and "is_transporter" in og_consensus.columns:
            core_trans = og_consensus[(og_consensus["Pangenome_Class"] == "Core") &
                                      (og_consensus["is_transporter"] == True)]
            acc_trans = og_consensus[(og_consensus["Pangenome_Class"] == "Accessory") &
                                     (og_consensus["is_transporter"] == True)]

            core_pct = 100 * len(core_trans) / core_total if core_total > 0 else 0
            acc_pct = 100 * len(acc_trans) / acc_total if acc_total > 0 else 0

            print(f"  • Transporters: {core_pct:.1f}% of Core vs {acc_pct:.1f}% of Accessory")

    # Cross-species summary
    print("\n" + "="*80)
    print("CROSS-SPECIES SUMMARY")
    print("="*80)
    print("""
INTERPRETING ENRICHMENT AND DEPLETION:

↑ ENRICHED (OR > 1): Category is over-represented in that genome class
↓ DEPLETED (OR < 1): Category is under-represented in that genome class

Key patterns to look for:

CORE GENOME:
  • Enriched → Essential, conserved functions (housekeeping, primary metabolism)
  • Depleted → Functions that vary between strains (not conserved)

ACCESSORY GENOME:
  • Enriched → Strain-specific adaptations, niche specialization
  • Depleted → Functions that are consistently present (moved to core)

Note: Core-enriched ≈ Accessory-depleted (and vice versa) - they are complementary views

Common findings in fungal pangenomes:
  • Core enriched: Translation (J), Energy production (C), Amino acid metabolism (E)
  • Core depleted: Secondary metabolism (Q), Defense mechanisms (V)
  • Accessory enriched: Secondary metabolism (Q), Carbohydrate transport (G), Defense (V)
  • Accessory depleted: Translation (J), Replication/repair (L)

Species-specific considerations:
  • A. fumigatus: Clinical isolates may have unique virulence factors in accessory
  • A. flavus: Aflatoxin genes may vary between strains
  • A. niger: Industrial enzyme diversity in accessory genome
  • A. oryzae: Domestication may have selected specific accessory genes
""")


# =============================================================================
# BIGSCAPE/ANTISMASH FUNCTIONS (from functions_for_analysis_1_6.py)
# =============================================================================

_CANON_USAGE = {
    "environmental": "environmental",
    "environment": "environmental",
    "env": "environmental",
    "industrial": "industrial",
    "industry": "industrial",
    "clinical": "clinical",
    "clinical/lab": "clinical",
    "lab": "lab",
    "laboratory": "lab",
    "unknown": "unknown",
    "na": "unknown",
    "": "unknown",
    None: "unknown",
}


def canon_usage(x: Optional[str]) -> str:
    """Canonicalize usage class strings."""
    if x is None:
        return "unknown"
    s = str(x).strip().lower()
    return _CANON_USAGE.get(s, s) or "unknown"


def list_region_files_for_assembly(antismash_dir: Path, assembly: str) -> List[Path]:
    """Find antiSMASH region GBKs for an assembly."""
    antismash_dir = Path(antismash_dir)
    assembly = str(assembly)
    hits: List[Path] = []

    uniq = antismash_dir / "antismash_gbk_unique"
    if uniq.exists():
        hits.extend(sorted(
            p for p in uniq.rglob("*.region*.gbk")
            if p.name.startswith(assembly) or assembly in p.name
        ))

    if not hits:
        sub = antismash_dir / assembly
        if sub.exists():
            hits.extend(sorted(sub.glob("*.region*.gbk")))
            if not hits:
                hits.extend(sorted(sub.glob("*.gbk")))

    if not hits:
        candidates = list(antismash_dir.rglob("*.region*.gbk")) or list(antismash_dir.rglob("*.gbk"))
        hits = [p for p in candidates if (assembly in p.name) or (assembly in str(p))]

    seen = set()
    dedup = []
    for p in hits:
        if p not in seen:
            seen.add(p)
            dedup.append(p)
    return dedup


def parse_antismash_regions_from_files(gbk_paths: Iterable[Path]) -> pd.DataFrame:
    """Parse antiSMASH .gbk files."""
    if SeqIO is None:
        print("[ERROR] biopython not available")
        return pd.DataFrame()

    rows = []
    for gbk in gbk_paths:
        try:
            for rec in SeqIO.parse(str(gbk), "genbank"):
                region_id = rec.name or rec.id
                sc = rec.annotations.get("structured_comment", {}).get("antiSMASH-Data", {})
                cluster_number = sc.get("cluster-number", "")
                product = sc.get("product", "")

                if not product:
                    predicted: List[str] = []
                    for feat in rec.features:
                        if "product" in feat.qualifiers:
                            predicted.extend(feat.qualifiers.get("product", []))
                    if predicted:
                        product = ";".join(sorted(set(predicted)))

                rows.append({
                    "bgc_path": str(gbk),
                    "region_id": region_id,
                    "cluster_number": cluster_number,
                    "raw_product": product or "",
                })
        except Exception as e:
            print(f"[WARN] Failed to parse {gbk}: {e}", file=sys.stderr)
            continue
    return pd.DataFrame(rows)


def build_regions_from_df(samples: pd.DataFrame,
                          species_root_map: Dict[str, Path],
                          antismash_subdir: str = "antismash_output") -> pd.DataFrame:
    """Build regions dataframe from samples with antiSMASH data."""
    req = {"Species", "Assembly Accession", "UsageClass"}
    if not req.issubset(samples.columns):
        raise ValueError(f"samples must contain columns: {req}")

    parts = []
    for _, row in samples.iterrows():
        species = str(row["Species"])
        assembly = str(row["Assembly Accession"])
        usage = canon_usage(row.get("UsageClass"))

        root = species_root_map.get(species)
        if not root or not Path(root).exists():
            continue

        antismash_dir = Path(root) / antismash_subdir
        if not antismash_dir.exists():
            continue

        files = list_region_files_for_assembly(antismash_dir, assembly)
        if not files:
            continue

        reg = parse_antismash_regions_from_files(files)
        if reg.empty:
            continue

        reg["species"] = species
        reg["strain"] = assembly
        reg["genome_id"] = f"{species}-{assembly}"
        reg["usage_class"] = usage
        parts.append(reg)

    if not parts:
        return pd.DataFrame(columns=[
            "bgc_path", "region_id", "cluster_number", "raw_product",
            "species", "strain", "genome_id", "usage_class"
        ])
    return pd.concat(parts, ignore_index=True)


_ANTISMASH_NOISE = {"hypothetical protein", "hypothetical", "protein", "domain", "like protein", "unknown"}
_TRNA_RE = re.compile(r"\btrna-[a-z]+\b", re.I)


def _tokens_from_antismash_product(raw_product: str) -> List[str]:
    """Split antiSMASH 'product' into tokens."""
    s = (raw_product or "").replace("_", "-")
    s = _TRNA_RE.sub(" ", s)
    toks = [t.strip() for t in re.split(r"[;|,/]+", s) if t.strip()]
    toks = [t for t in toks if t.lower() not in _ANTISMASH_NOISE]
    return toks


def add_antismash_direct_labels(regions: pd.DataFrame,
                                raw_col: str = "raw_product",
                                direct_col: str = "antismash_primary_direct",
                                collapsed_col: str = "antismash_primary_collapsed",
                                tokens_col: str = "antismash_tokens",
                                raw_label_col: str = "antismash_label_raw") -> pd.DataFrame:
    """Add antiSMASH-derived label columns to regions dataframe."""
    df = regions.copy()
    df[raw_label_col] = df[raw_col].fillna("")
    toks = df[raw_col].apply(_tokens_from_antismash_product)
    df[tokens_col] = toks.apply(lambda xs: ";".join(xs))
    df[direct_col] = toks.apply(lambda xs: xs[0] if len(xs) else "Unassigned")

    def _collapse(tok: str) -> str:
        t = (tok or "").lower()
        if t in {"t1pks", "t2pks", "t3pks", "polyketide"}:
            return "PKS"
        if t.startswith("nrps"):
            return "NRPS"
        if t.startswith("terpene"):
            return "Terpene"
        if t in {"fungal-ripp"}:
            return "RiPP"
        if t in {"ripp-like", "napaa"}:
            return "RiPP-like"
        if t in {"indole", "dmats"}:
            return "Indole"
        if t in {"ni-siderophore", "siderophore"}:
            return "Siderophore"
        if t in {"beta-lactone", "betalactone"}:
            return "Beta-lactone"
        if t == "phosphonate":
            return "Phosphonate"
        if t == "isocyanide":
            return "Isocyanide"
        if t == "saccharide":
            return "Saccharide"
        if t == "melanin":
            return "Melanin"
        if t == "bacteriocin":
            return "Bacteriocin"
        return tok

    df[collapsed_col] = toks.apply(lambda xs: "Unassigned" if not xs else _collapse(xs[0]))
    return df


def vivid_colors(n: int):
    """Vivid categorical colors built from tab20 / tab20b / tab20c / tab10."""
    banks = [
        list(plt.get_cmap("tab20").colors),
        list(plt.get_cmap("tab20b").colors),
        list(plt.get_cmap("tab20c").colors),
        list(plt.get_cmap("tab10").colors),
    ]
    base = list(itertools.chain.from_iterable(banks))
    if n <= len(base):
        return base[:n]
    reps = (n + len(base) - 1) // len(base)
    return (base * reps)[:n]


def compute_tag_enrichment(
    regions: pd.DataFrame,
    tag_col: str,
    usage_col: str = "usage_class",
    exclude_tags: tuple = ("Unassigned",),
    min_tag_total: int = 3,
    alternative: str = "greater",
) -> pd.DataFrame:
    """For each tag and usage class, test enrichment using Fisher's exact test."""
    if regions is None or regions.empty:
        return pd.DataFrame(columns=["usage_class", "tag", "a", "b", "c", "d", "or", "log2_or",
                                      "pval", "qval", "-log10(qval)", "frac_in_class", "frac_overall",
                                      "total_in_class", "total_tag", "total_all"])

    def _haldane_or(a, b, c, d) -> float:
        a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        return (a_ * d_) / (b_ * c_)

    df = regions[[usage_col, tag_col]].dropna().copy()
    N = len(df)
    classes = df[usage_col].astype(str).unique()
    tags = [t for t in df[tag_col].astype(str).unique() if t not in exclude_tags]

    out = []
    for tag in tags:
        tag_mask = (df[tag_col].astype(str) == tag)
        total_tag = int(tag_mask.sum())
        if total_tag < min_tag_total:
            continue
        for cls in classes:
            cls_mask = (df[usage_col].astype(str) == cls)
            total_in_class = int(cls_mask.sum())
            a = int((tag_mask & cls_mask).sum())
            b = int(total_in_class - a)
            c = int(total_tag - a)
            d = int(N - (a + b + c))
            OR = _haldane_or(a, b, c, d)
            try:
                _, p = fisher_exact([[a, b], [c, d]], alternative=alternative)
            except Exception:
                p = np.nan
            out.append({
                "usage_class": cls, "tag": tag,
                "a": a, "b": b, "c": c, "d": d,
                "or": float(OR), "log2_or": float(np.log2(OR)),
                "pval": float(p),
                "frac_in_class": (a / total_in_class) if total_in_class else 0.0,
                "frac_overall": (total_tag / N) if N else 0.0,
                "total_in_class": total_in_class, "total_tag": total_tag, "total_all": N
            })
    res = pd.DataFrame(out)
    if res.empty:
        return res
    if res["pval"].notna().any():
        res["qval"] = _bh_fdr(res["pval"].to_numpy())
        res["-log10(qval)"] = res["qval"].apply(lambda q: -np.log10(q) if (pd.notna(q) and q > 0) else 0.0)
    else:
        res["qval"] = np.nan
        res["-log10(qval)"] = 0.0
    return res.sort_values(["usage_class", "qval", "or"], ascending=[True, True, False]).reset_index(drop=True)


def load_bigscape_full(base_dir: Path, species_list: List[str], similarity_cutoff: str = "c0.3") -> pd.DataFrame:
    """Load BiG-SCAPE results for all species."""
    rows = []
    for species in species_list:
        species_dir = Path(base_dir) / species / "bigscape_output" / "output_files"
        if not species_dir.exists():
            print(f"[SKIP] {species}: no output_files dir at {species_dir}")
            continue

        ann_candidates = glob.glob(str(species_dir / "**" / "record_annotations.tsv"), recursive=True)
        if not ann_candidates:
            continue
        ann_path = max(ann_candidates, key=os.path.getmtime)
        ann_df = pd.read_csv(ann_path, sep="\t", dtype=str)

        def _ci_cols(df):
            return {c.lower(): c for c in df.columns}

        ann_cols = _ci_cols(ann_df)
        record_col_ann = ann_cols.get("record") or ann_cols.get("bgc") or list(ann_df.columns)[0]
        ann_df = ann_df.rename(columns={record_col_ann: "Record"})

        cutoff_dirs = glob.glob(str(species_dir / f"*_{similarity_cutoff}"))
        if not cutoff_dirs:
            continue
        cutoff_dir = max(cutoff_dirs, key=os.path.getmtime)

        class_dirs = [d for d in os.listdir(cutoff_dir) if os.path.isdir(os.path.join(cutoff_dir, d))]
        if not class_dirs:
            continue

        for bgc_class in class_dirs:
            clustering_file = os.path.join(cutoff_dir, bgc_class, f"{bgc_class}_clustering_{similarity_cutoff}.tsv")
            if not os.path.exists(clustering_file):
                continue
            clustering_df = pd.read_csv(clustering_file, sep="\t", dtype=str)
            cc = _ci_cols(clustering_df)
            record_col = cc.get("record") or cc.get("bgc") or cc.get("cluster") or list(clustering_df.columns)[0]
            fam_col = cc.get("family") or (list(clustering_df.columns)[1] if len(clustering_df.columns) > 1 else None)
            cc_col = cc.get("cc") or cc.get("clusterclass")
            clustering_df = clustering_df.rename(columns={
                record_col: "Record",
                **({fam_col: "Family"} if fam_col else {}),
                **({cc_col: "CC"} if cc_col else {})
            })
            merged = clustering_df.merge(ann_df, on="Record", how="left")
            merged["BGC_Class"] = bgc_class
            if "Family" not in merged.columns:
                merged["Family"] = pd.NA
            if "CC" not in merged.columns:
                merged["CC"] = pd.NA

            merged["Species"] = species
            merged["Family_ID"] = merged["Species"].astype(str) + "_" + merged["Family"].astype(str)
            merged["CC_ID"] = merged["Species"].astype(str) + "_" + merged["CC"].astype(str)
            rows.append(merged)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def derive_bgc_key_from_path(p: str) -> str:
    """Unique key based on file stem."""
    s = str(p).strip().strip('"').strip("'")
    try:
        stem = Path(s).stem
    except Exception:
        stem = s.rsplit(".", 1)[0]
    return stem.lower()


def bigscape_to_gcf_map(bigscape_full_df: pd.DataFrame) -> pd.DataFrame:
    """Create mapping from BGC key to GCF."""
    if bigscape_full_df is None or bigscape_full_df.empty:
        return pd.DataFrame(columns=["bgc_key", "GCF"])
    df = bigscape_full_df.copy()
    rec_col = next((c for c in ["Record", "record", "BGC", "cluster"] if c in df.columns), None)
    if not rec_col:
        return pd.DataFrame(columns=["bgc_key", "GCF"])
    df["bgc_key"] = df[rec_col].astype(str).apply(derive_bgc_key_from_path)
    if "Family_ID" in df.columns and df["Family_ID"].notna().any():
        df["GCF"] = df["Family_ID"].astype(str)
    elif "Family" in df.columns:
        df["GCF"] = df["Family"].astype(str)
    else:
        df["GCF"] = pd.NA
    return df[["bgc_key", "GCF"]].dropna().drop_duplicates()


def join_regions_with_bigscape_full(regions_df: pd.DataFrame, bigscape_full_df: pd.DataFrame) -> pd.DataFrame:
    """Join regions with BiG-SCAPE GCF assignments."""
    mapping = bigscape_to_gcf_map(bigscape_full_df)
    out = regions_df.copy()
    out["bgc_key"] = out["bgc_path"].astype(str).apply(derive_bgc_key_from_path)
    if mapping.empty:
        out["GCF"] = pd.NA
        return out
    return out.merge(mapping, how="left", on="bgc_key")


def compute_gcf_uniqueness_and_enrichment(regions_with_gcf: pd.DataFrame,
                                          usage_col: str = "usage_class",
                                          gcf_col: str = "GCF",
                                          min_members: int = 3,
                                          top_frac_threshold: float = 0.8) -> pd.DataFrame:
    """Compute GCF uniqueness and enrichment statistics."""
    if regions_with_gcf is None or regions_with_gcf.empty or gcf_col not in regions_with_gcf.columns:
        return pd.DataFrame(columns=["GCF", "n_members", "top_class", "top_frac", "is_unique",
                                      "is_enriched", "classes_breakdown"])
    df = regions_with_gcf.dropna(subset=[gcf_col])
    rows = []
    for gcf, g in df.groupby(gcf_col):
        counts = g[usage_col].value_counts()
        n = int(counts.sum())
        top_class = counts.idxmax()
        top_frac = float(counts.max()) / n if n else 0.0
        is_unique = (len(counts) == 1)
        is_enriched = (counts.max() >= min_members) and (top_frac >= top_frac_threshold)
        rows.append({
            "GCF": gcf,
            "n_members": n,
            "top_class": top_class,
            "top_frac": round(top_frac, 3),
            "is_unique": bool(is_unique),
            "is_enriched": bool(is_enriched),
            "classes_breakdown": ";".join(f"{k}:{v}" for k, v in counts.items())
        })
    return pd.DataFrame(rows).sort_values(["is_unique", "is_enriched", "top_frac", "n_members"],
                                          ascending=[False, False, False, False])


def annotate_gcf_with_class(regions_with_gcf: pd.DataFrame,
                            class_col: str = "antismash_primary_collapsed",
                            gcf_col: str = "GCF") -> pd.DataFrame:
    """Assign a representative class to each GCF."""
    if regions_with_gcf is None or regions_with_gcf.empty:
        return pd.DataFrame(columns=[gcf_col, "GCF_class"])
    if gcf_col not in regions_with_gcf.columns or class_col not in regions_with_gcf.columns:
        return pd.DataFrame(columns=[gcf_col, "GCF_class"])
    df = regions_with_gcf.dropna(subset=[gcf_col]).copy()
    rep = (df.groupby(gcf_col)[class_col]
             .agg(lambda s: s.mode().iat[0] if not s.mode().empty else "Unassigned")
             .reset_index()
             .rename(columns={class_col: "GCF_class"}))
    return rep


def plot_enrichment_heatmap(enr_df: pd.DataFrame,
                            value: str = "log2_or",
                            mask_by_q: float = 0.1,
                            title: str = None,
                            savepath: str = None,
                            cmap_name: str = "viridis"):
    """Heatmap of tag enrichment across usage classes."""
    if enr_df is None or enr_df.empty:
        print("[plot] enrichment table empty; skipping heatmap.")
        return None, None

    df = enr_df.copy()
    if (mask_by_q is not None) and ("qval" in df.columns):
        df.loc[df["qval"].isna() | (df["qval"] > mask_by_q), value] = 0.0

    mat = (df.pivot(index="usage_class", columns="tag", values=value)
             .fillna(0.0)
             .sort_index())
    data = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots()
    im = ax.imshow(data, aspect="auto", cmap=plt.get_cmap(cmap_name))
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(mat.index.tolist())
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(mat.columns.tolist(), rotation=45, ha="right")
    ax.set_title(title or f"Enrichment heatmap: {value}")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(value)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")


# =============================================================================
# ADDITIONAL NOTEBOOK-SPECIFIC FUNCTIONS
# =============================================================================

def plot_orthofinder_summary(orthofinder_data, species_list):
    """
    Create summary visualization of OrthoFinder results per species.
    - Bar chart of orthogroup counts (total, single-copy, species-specific)
    - Pie chart of gene distribution (in OGs vs unassigned)
    """
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.3, wspace=0.25)
    
    # Prepare data
    species_names = []
    total_ogs = []
    single_copy = []
    species_specific = []
    all_species_ogs = []
    genes_in_og = []
    unassigned = []
    pct_in_og = []
    n_genomes_list = []
    
    for sp in species_list:
        data = orthofinder_data.get(sp, {})
        stats = data.get("stats", {})
        
        species_names.append(f"A. {sp}")
        total_ogs.append(stats.get("Number of orthogroups", 0))
        single_copy.append(stats.get("Number of single-copy orthogroups", 0))
        species_specific.append(stats.get("Number of species-specific orthogroups", 0))
        all_species_ogs.append(stats.get("Number of orthogroups with all species present", 0))
        
        total_genes = stats.get("Number of genes", 0)
        in_og = stats.get("Number of genes in orthogroups", 0)
        genes_in_og.append(in_og)
        unassigned.append(total_genes - in_og)
        pct_in_og.append(stats.get("Percentage of genes in orthogroups", 0))
        # Use actual PAV matrix column count (respects ANI exclusions) if available,
        # otherwise fall back to OrthoFinder stats (may include excluded samples)
        pav = data.get("pav")
        if pav is not None:
            n_genomes_list.append(pav.shape[1])
        else:
            n_genomes_list.append(stats.get("Number of species", 0))
    
    # 1. Orthogroup counts bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    x = np.arange(len(species_names))
    width = 0.2
    
    ax1.bar(x - 1.5*width, total_ogs, width, label='Total OGs', color='steelblue')
    ax1.bar(x - 0.5*width, all_species_ogs, width, label='Core OGs (all genomes)', color='forestgreen')
    ax1.bar(x + 0.5*width, single_copy, width, label='Single-copy OGs', color='goldenrod')
    ax1.bar(x + 1.5*width, species_specific, width, label='Species-specific OGs', color='coral')
    
    ax1.set_xlabel('Species')
    ax1.set_ylabel('Number of Orthogroups')
    ax1.set_title('OrthoFinder: Orthogroup Statistics per Species')
    ax1.set_xticks(x)
    ax1.set_xticklabels(species_names, rotation=15)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(axis='y', linestyle=':', alpha=0.5)
    
    # Add genome count annotations
    for i, n in enumerate(n_genomes_list):
        ax1.annotate(f'n={n}', (i, total_ogs[i]), ha='center', va='bottom', fontsize=9)
    
    # 2. Gene assignment percentages
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.bar(x - 0.2, genes_in_og, 0.4, label='Genes in OGs', color='forestgreen', alpha=0.8)
    ax2.bar(x + 0.2, unassigned, 0.4, label='Unassigned genes', color='lightcoral', alpha=0.8)
    
    ax2.set_xlabel('Species')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Gene Assignment to Orthogroups')
    ax2.set_xticks(x)
    ax2.set_xticklabels(species_names, rotation=15)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(axis='y', linestyle=':', alpha=0.5)
    
    # Add percentage annotations
    for i, pct in enumerate(pct_in_og):
        ax2.annotate(f'{pct:.1f}%', (i-0.2, genes_in_og[i]), ha='center', va='bottom', fontsize=9)
    
    # 3. Core genome size comparison (OGs present in all genomes)
    ax3 = fig.add_subplot(gs[1, 0])
    colors = sns.color_palette("Set2", len(species_names))
    bars = ax3.bar(species_names, all_species_ogs, color=colors)
    
    for i, (bar, n, total) in enumerate(zip(bars, all_species_ogs, total_ogs)):
        pct = 100 * n / total if total > 0 else 0
        ax3.annotate(f'{pct:.1f}%\nof OGs', (bar.get_x() + bar.get_width()/2, bar.get_height()),
                    ha='center', va='bottom', fontsize=9)
    
    ax3.set_ylabel('Number of Core Orthogroups')
    ax3.set_title('Core Genome Size (OGs present in ALL genomes)')
    ax3.set_xticklabels(species_names, rotation=15)
    ax3.grid(axis='y', linestyle=':', alpha=0.5)
    
    # 4. Summary table
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    table_data = []
    for i, sp in enumerate(species_names):
        table_data.append([
            sp,
            f"{n_genomes_list[i]}",
            f"{total_ogs[i]:,}",
            f"{single_copy[i]:,}",
            f"{all_species_ogs[i]:,}",
            f"{pct_in_og[i]:.1f}%"
        ])
    
    table = ax4.table(
        cellText=table_data,
        colLabels=['Species', 'Genomes', 'Total OGs', 'Single-copy', 'Core OGs', '% in OGs'],
        cellLoc='center',
        loc='center',
        colColours=['lightsteelblue']*6
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    ax4.set_title('OrthoFinder Summary Statistics', fontsize=12, fontweight='bold', pad=20)
    
    plt.suptitle('OrthoFinder Comparative Genomics Summary', fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    


def analyze_pangenome_classification(orthofinder_data, species_list):
    """
    Classify genes into Core/Accessory/Rare for each species using 
    the determine_core_and_rare_thresholds function.
    """
    results = {}
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, species in enumerate(species_list):
        data = orthofinder_data.get(species, {})
        pav = data.get("pav")
        
        if pav is None or pav.empty:
            print(f"[warn] No PAV matrix for {species}")
            continue
        
        # Determine thresholds using the heuristic function
        core_rate, core_n, rare_rate, rare_n = determine_core_and_rare_thresholds(pav)
        
        # Count genes in each category
        gene_freq = pav.sum(axis=1)
        n_genomes = pav.shape[1]
        
        core_cnt = int((gene_freq >= core_n).sum())
        rare_cnt = int((gene_freq <= rare_n).sum())
        accessory_cnt = int(((gene_freq < core_n) & (gene_freq > rare_n)).sum())
        total = len(gene_freq)
        
        results[species] = {
            "n_genomes": n_genomes,
            "core_n": core_n,
            "rare_n": rare_n,
            "core_rate": core_rate,
            "rare_rate": rare_rate,
            "core_count": core_cnt,
            "accessory_count": accessory_cnt,
            "rare_count": rare_cnt,
            "total_ogs": total
        }
        
        print(f"\n{species.upper()} (n={n_genomes} genomes):")
        print(f"  Core threshold: >= {core_n} genomes ({core_rate:.1%})")
        print(f"  Rare threshold: <= {rare_n} genomes ({rare_rate:.1%})")
        print(f"  Core: {core_cnt} OGs ({100*core_cnt/total:.1f}%)")
        print(f"  Accessory: {accessory_cnt} OGs ({100*accessory_cnt/total:.1f}%)")
        print(f"  Rare: {rare_cnt} OGs ({100*rare_cnt/total:.1f}%)")
        
        # Plot using the plot_gene_freq_hist function
        ax = axes[idx]
        plot_gene_freq_hist(pav, core_n, rare_n, ax=ax, annotate=True)
        ax.set_title(f"A. {species} (n={n_genomes})", fontsize=11, fontweight='bold')
    
    plt.suptitle('Gene Frequency Distribution and Pangenome Classification', 
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    
    return results, fig


def heaps_law(pav, n_boot=200, max_n=None):
    """
    Compute pangenome rarefaction curves and Heap's law fit.
    Returns bootstrap mean curves plus 2.5th/97.5th percentile confidence intervals.
    """
    import random
    species = list(pav.columns)
    if max_n is None:
        max_n = len(species)
    Ns = np.arange(1, max_n + 1)
    pangenome_size_mean = np.zeros_like(Ns, dtype=float)
    new_genes_mean = np.zeros_like(Ns, dtype=float)
    core_size_mean = np.zeros_like(Ns, dtype=float)

    # Store all bootstrap replicates for CI calculation
    all_pan = np.zeros((n_boot, len(Ns)))
    all_core = np.zeros((n_boot, len(Ns)))

    for b in range(n_boot):
        order = species[:]
        random.shuffle(order)
        seen = set()
        core = None
        acc_cumulative = []
        acc_new = []
        acc_core = []

        for i, sp in enumerate(order[:max_n], 1):
            present = set(pav.index[pav[sp] == 1])
            new = present - seen
            seen |= present

            # Track core (genes present in ALL genomes so far)
            if core is None:
                core = present.copy()
            else:
                core &= present

            acc_cumulative.append(len(seen))
            acc_new.append(len(new))
            acc_core.append(len(core))

        pangenome_size_mean += np.array(acc_cumulative)
        new_genes_mean += np.array(acc_new)
        core_size_mean += np.array(acc_core)
        all_pan[b] = acc_cumulative
        all_core[b] = acc_core

    pangenome_size_mean /= max(1, n_boot)
    new_genes_mean /= max(1, n_boot)
    core_size_mean /= max(1, n_boot)

    # Bootstrap confidence intervals (2.5th – 97.5th percentile)
    pan_ci_lo = np.percentile(all_pan, 2.5, axis=0)
    pan_ci_hi = np.percentile(all_pan, 97.5, axis=0)
    core_ci_lo = np.percentile(all_core, 2.5, axis=0)
    core_ci_hi = np.percentile(all_core, 97.5, axis=0)

    # Fit Heap's law: log(P) = log(k) + gamma * log(N)
    x = np.log(Ns)
    y = np.log(pangenome_size_mean)
    A = np.vstack([np.ones_like(x), x]).T
    klog, gamma = np.linalg.lstsq(A, y, rcond=None)[0]
    k = np.exp(klog)

    return Ns, pangenome_size_mean, new_genes_mean, core_size_mean, k, gamma, pan_ci_lo, pan_ci_hi, core_ci_lo, core_ci_hi


def analyze_heaps_law(orthofinder_data, species_list, n_boot=200,
                      results_base=None, force=False):
    """
    Analyze pangenome openness using rarefaction curves.

    If ``results_base`` is provided, per-species results are cached as
    ``{results_base}/{sp}/{sp}_heaps.npz``. Subsequent calls load from
    cache instead of re-running the bootstrap (use ``force=True`` to rerun).
    """
    import pickle
    results = {}

    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    colors = sns.color_palette("Set2", len(species_list))

    for idx, species in enumerate(species_list):
        cache_path = None
        if results_base is not None:
            cache_path = os.path.join(results_base, species, f'{species}_heaps.pkl')

        # --- Try to load from cache ---
        if not force and cache_path and os.path.exists(cache_path):
            with open(cache_path, 'rb') as f:
                cached = pickle.load(f)
            results[species] = cached
            Ns = np.arange(1, cached['n_genomes'] + 1)
            pangenome_size = cached['pangenome_size']
            new_genes = cached['new_genes']
            core_size = cached['core_size']
            k = cached['k']; gamma = cached['gamma']
            pan_ci_lo = cached['pan_ci_lo']; pan_ci_hi = cached['pan_ci_hi']
            core_ci_lo = cached['core_ci_lo']; core_ci_hi = cached['core_ci_hi']
            status = cached['status']
            n_genomes = cached['n_genomes']
            first_genome_pct = cached['first_genome_pct']
            status_color = {'OPEN': 'green', 'MODERATELY OPEN': 'orange',
                            'CLOSED': 'coral'}.get(status, 'gray')
            print(f"  {species}: loaded cached Heaps' law results (γ = {gamma:.3f}, {status})")
        else:
            data = orthofinder_data.get(species, {})
            pav = data.get("pav")
            if pav is None or pav.empty:
                print(f"[warn] No PAV matrix for {species}")
                continue

            print(f"\nAnalyzing pangenome for {species}...")
            Ns, pangenome_size, new_genes, core_size, k, gamma, pan_ci_lo, pan_ci_hi, core_ci_lo, core_ci_hi = heaps_law(pav, n_boot=n_boot)
            n_genomes = len(Ns)

            # Calculate metrics
            first_genome_pct = 100 * pangenome_size[0] / pangenome_size[-1]

            # Determine status
            if gamma > 0.8:
                status = "OPEN"; status_color = "green"
            elif gamma > 0.6:
                status = "MODERATELY OPEN"; status_color = "orange"
            else:
                status = "CLOSED"; status_color = "coral"

            results[species] = {
                "k": k, "gamma": gamma, "status": status,
                "n_genomes": n_genomes,
                "pangenome_size": pangenome_size,
                "new_genes": new_genes,
                "core_size": core_size,
                "first_genome_pct": first_genome_pct,
                "pan_ci_lo": pan_ci_lo, "pan_ci_hi": pan_ci_hi,
                "core_ci_lo": core_ci_lo, "core_ci_hi": core_ci_hi,
            }
            print(f"  {species}: γ = {gamma:.3f}, 1st genome = {first_genome_pct:.1f}% of pangenome -> {status}")

            # Save to cache
            if cache_path is not None:
                os.makedirs(os.path.dirname(cache_path), exist_ok=True)
                with open(cache_path, 'wb') as f:
                    pickle.dump(results[species], f)
                print(f"  Saved Heaps' law cache to {cache_path}")
        
        # TOP ROW: Pangenome & Core rarefaction curves
        ax_top = axes[0, idx]
        
        # Normalize to percentage of final pangenome
        pan_pct = 100 * pangenome_size / pangenome_size[-1]
        core_pct = 100 * core_size / pangenome_size[-1]
        
        ax_top.fill_between(Ns, core_pct, pan_pct, alpha=0.3, color='steelblue', label='Accessory')
        ax_top.fill_between(Ns, 0, core_pct, alpha=0.3, color='forestgreen', label='Core')
        ax_top.plot(Ns, pan_pct, 'o-', color='steelblue', markersize=2, linewidth=2, label='Pangenome')
        ax_top.plot(Ns, core_pct, 's-', color='forestgreen', markersize=2, linewidth=2, label='Core genome')
        
        ax_top.set_xlabel('Number of genomes')
        ax_top.set_ylabel('% of total pangenome')
        ax_top.set_title(f"A. {species} (n={n_genomes})", fontsize=11, fontweight='bold')
        ax_top.set_ylim(0, 105)
        ax_top.legend(loc='right', fontsize=8)
        ax_top.grid(True, linestyle=':', alpha=0.5)
        
        # Annotation
        textstr = f"γ = {gamma:.3f}\n1st genome: {first_genome_pct:.0f}%\n{status}"
        props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor=status_color, linewidth=2)
        ax_top.text(0.98, 0.02, textstr, transform=ax_top.transAxes, fontsize=9,
                   verticalalignment='bottom', horizontalalignment='right', bbox=props)
        
        # BOTTOM ROW: New genes per genome (log scale for y-axis)
        ax_bot = axes[1, idx]
        ax_bot.bar(Ns, new_genes, color=colors[idx], alpha=0.7, width=0.8)
        ax_bot.set_yscale('log')
        ax_bot.set_xlabel('Genome number')
        ax_bot.set_ylabel('New orthogroups (log scale)')
        ax_bot.set_title(f'Gene discovery rate', fontsize=10)
        ax_bot.grid(True, linestyle=':', alpha=0.5, axis='y')
        
        # Mark where 95% of pangenome is reached
        for i, pct in enumerate(pan_pct):
            if pct >= 95:
                ax_bot.axvline(x=i+1, color='red', linestyle='--', alpha=0.7)
                ax_bot.text(i+1, ax_bot.get_ylim()[1]*0.5, f'95% at\nn={i+1}', 
                           fontsize=8, ha='left', color='red')
                break
    
    plt.suptitle("Pangenome Rarefaction Analysis",
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.show()

    return results, fig


def plot_pangenome_summary(pangenome_results, heaps_results, species_list):
    """
    Create a combined summary of pangenome analysis across all species.
    """
    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.3)
    
    species_names = [f"A. {sp}" for sp in species_list]
    colors = sns.color_palette("Set2", len(species_list))
    
    # 1. Stacked bar chart: Core/Accessory/Rare proportions
    ax1 = fig.add_subplot(gs[0])
    
    core_pcts = []
    accessory_pcts = []
    rare_pcts = []
    
    for sp in species_list:
        r = pangenome_results.get(sp, {})
        total = r.get("total_ogs", 1)
        core_pcts.append(100 * r.get("core_count", 0) / total)
        accessory_pcts.append(100 * r.get("accessory_count", 0) / total)
        rare_pcts.append(100 * r.get("rare_count", 0) / total)
    
    x = np.arange(len(species_names))
    width = 0.6
    
    ax1.bar(x, core_pcts, width, label='Core', color='forestgreen')
    ax1.bar(x, accessory_pcts, width, bottom=core_pcts, label='Accessory', color='steelblue')
    ax1.bar(x, rare_pcts, width, bottom=[c+a for c,a in zip(core_pcts, accessory_pcts)], label='Rare', color='coral')
    
    ax1.set_ylabel('Percentage of Orthogroups')
    ax1.set_title('Pangenome Composition', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(species_names, rotation=15)
    ax1.legend(loc='upper right')
    ax1.set_ylim(0, 105)
    
    for i in range(len(species_names)):
        if core_pcts[i] > 5:
            ax1.text(i, core_pcts[i]/2, f'{core_pcts[i]:.0f}%', ha='center', va='center', 
                    fontsize=9, color='white', fontweight='bold')
        if accessory_pcts[i] > 5:
            ax1.text(i, core_pcts[i] + accessory_pcts[i]/2, f'{accessory_pcts[i]:.0f}%', 
                    ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        if rare_pcts[i] > 5:
            ax1.text(i, core_pcts[i] + accessory_pcts[i] + rare_pcts[i]/2, f'{rare_pcts[i]:.0f}%', 
                    ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    
    # 2. First genome contribution
    ax2 = fig.add_subplot(gs[1])
    
    first_pcts = [heaps_results.get(sp, {}).get("first_genome_pct", 0) for sp in species_list]
    bar_colors = ['coral' if p > 85 else 'orange' if p > 70 else 'forestgreen' for p in first_pcts]
    
    bars = ax2.bar(species_names, first_pcts, color=bar_colors, edgecolor='black', linewidth=1)
    ax2.axhline(y=85, color='coral', linestyle='--', linewidth=1, alpha=0.7, label='Very closed (>85%)')
    ax2.axhline(y=70, color='orange', linestyle=':', linewidth=1, alpha=0.7, label='Moderately closed')
    
    ax2.set_ylabel("% pangenome from 1st genome")
    ax2.set_title("Pangenome Saturation", fontsize=12, fontweight='bold')
    ax2.set_xticklabels(species_names, rotation=15)
    ax2.set_ylim(0, 100)
    ax2.legend(loc='lower right', fontsize=8)
    
    for bar, p in zip(bars, first_pcts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f'{p:.0f}%', 
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # 3. Summary table
    ax3 = fig.add_subplot(gs[2])
    ax3.axis('off')
    
    table_data = []
    for sp in species_list:
        pr = pangenome_results.get(sp, {})
        hr = heaps_results.get(sp, {})
        table_data.append([
            f"A. {sp}",
            f"{pr.get('n_genomes', '-')}",
            f"{pr.get('total_ogs', '-'):,}",
            f"{100*pr.get('core_count',0)/max(1,pr.get('total_ogs',1)):.0f}%",
            f"{hr.get('first_genome_pct', 0):.0f}%",
            hr.get('status', '-')
        ])
    
    table = ax3.table(
        cellText=table_data,
        colLabels=['Species', 'Genomes', 'Total OGs', 'Core %', '1st Gen %', 'Status'],
        cellLoc='center',
        loc='center',
        colColours=['lightsteelblue']*6
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.1, 1.5)
    ax3.set_title('Pangenome Summary', fontsize=12, fontweight='bold', pad=20)
    
    plt.suptitle('Comparative Pangenome Analysis Across Aspergillus Species', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    


def load_species_annotations(species, base_path, pangenome_results, orthofinder_data):
    """Load and classify annotations for a species."""
    species_dir = os.path.join(base_path, species)
    cache_path = os.path.join(species_dir, f"{species}_ortho_annot_long.tsv")
    
    print(f"\n{'='*60}")
    print(f"Loading annotations for A. {species}...")
    print(f"{'='*60}")
    
    # Build or load the annotation table
    try:
        long_annot = build_ortho_annot_table(species_dir, species, cache_path)
    except Exception as e:
        print(f"  Error loading annotations: {e}")
        return None
    
    # Get PAV matrix and classify
    pav_data = orthofinder_data.get(species, {})
    pav = pav_data.get("pav")
    
    if pav is not None:
        pr = pangenome_results.get(species, {})
        core_n = pr.get("core_n", int(0.95 * pav.shape[1]))
        rare_n = pr.get("rare_n", max(1, int(0.05 * pav.shape[1])))
        
        # Add pangenome class
        long_annot = classify_pangenome(long_annot, pav, core_n, rare_n)
        print(f"  Pangenome classification applied (Core>={core_n}, Rare<={rare_n})")
    
    # Add species column
    long_annot["Species"] = species
    
    # Stats
    n_proteins = len(long_annot)
    n_assemblies = long_annot["Assembly Accession"].nunique()
    n_ogs = long_annot["Orthogroup"].nunique()
    
    print(f"  Proteins: {n_proteins:,}")
    print(f"  Assemblies: {n_assemblies}")
    print(f"  Orthogroups: {n_ogs:,}")
    
    if "Pangenome_Class" in long_annot.columns:
        class_counts = long_annot["Pangenome_Class"].value_counts()
        print(f"  Pangenome classes: {class_counts.to_dict()}")
    
    return long_annot


def analyze_derived_categories(og_consensus_tables, species_list):
    """
    Analyze enrichment of derived categories (proteases, transporters, secretome)
    in Core vs Accessory genome.
    """
    from scipy.stats import fisher_exact
    
    print("\n" + "="*80)
    print("DERIVED CATEGORY ENRICHMENT BY PANGENOME CLASS")
    print("="*80)
    
    categories = [
        ("is_protease", "Proteases"),
        ("is_transporter", "Transporters"),
        ("is_secreted", "Secretome (SignalP + GO)"),
        ("is_secreted_signalp", "Secretome (SignalP only)"),
        ("has_secretion_GO", "Secretome (GO terms only)"),
    ]
    
    for species in species_list:
        if species not in og_consensus_tables:
            continue
        
        og = og_consensus_tables[species]
        print(f"\nA. {species}:")
        
        for col, name in categories:
            if col not in og.columns:
                continue
            
            # Count per class
            results = {}
            for cls in ["Core", "Accessory"]:
                subset = og[og["Pangenome_Class"] == cls]
                n_total = len(subset)
                n_positive = subset[col].sum()
                pct = 100 * n_positive / n_total if n_total > 0 else 0
                results[cls] = {"n": n_positive, "total": n_total, "pct": pct}
            
            # Fisher test: Accessory vs Core
            core_pos = results["Core"]["n"]
            core_neg = results["Core"]["total"] - core_pos
            acc_pos = results["Accessory"]["n"]
            acc_neg = results["Accessory"]["total"] - acc_pos
            
            if core_pos + acc_pos > 0:
                try:
                    odds_ratio, p_value = fisher_exact([[acc_pos, acc_neg], [core_pos, core_neg]])
                except:
                    odds_ratio, p_value = 1.0, 1.0
                
                sig_marker = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                direction = "up" if odds_ratio > 1 else "down" if odds_ratio < 1 else "="
                
                print(f"  {name}:")
                print(f"    Core: {results['Core']['n']:,}/{results['Core']['total']:,} ({results['Core']['pct']:.1f}%)")
                print(f"    Accessory: {results['Accessory']['n']:,}/{results['Accessory']['total']:,} ({results['Accessory']['pct']:.1f}%)")
                print(f"    Accessory vs Core: OR={odds_ratio:.2f}, p={p_value:.2e} {direction} {sig_marker}")


def debug_signalp_loading(species_annotations):
    """Debug why SignalP data might not be loading correctly."""
    print("\n" + "="*80)
    print("SIGNALP DATA DIAGNOSTIC")
    print("="*80)
    
    for species, annot_df in species_annotations.items():
        print(f"\nA. {species}:")
        
        # Check if signalp column exists
        signalp_cols = [c for c in annot_df.columns if 'signalp' in c.lower()]
        print(f"  SignalP columns found: {signalp_cols}")
        
        if 'signalp_Prediction' in annot_df.columns:
            col = annot_df['signalp_Prediction']
            n_non_empty = (col.fillna("").astype(str) != "").sum()
            unique_vals = col.fillna("").astype(str).unique()[:10]
            print(f"  Non-empty values: {n_non_empty:,} / {len(annot_df):,}")
            print(f"  Sample values: {list(unique_vals)}")
            
            # Check for SP predictions
            has_sp = col.fillna("").str.upper().str.contains("SP")
            print(f"  Proteins with 'SP' in prediction: {has_sp.sum():,}")
        else:
            print("  signalp_Prediction column NOT FOUND")
            
            # Try to find any signalp-related data
            for col_name in signalp_cols:
                col = annot_df[col_name]
                n_non_empty = (col.fillna("").astype(str) != "").sum()
                print(f"    {col_name}: {n_non_empty:,} non-empty values")


def analyze_secretome_from_consensus(og_consensus_tables, species_list):
    """
    Detailed secretome analysis using pre-computed orthogroup consensus.
    """
    from scipy.stats import fisher_exact
    
    print("\n" + "="*80)
    print("SECRETOME ANALYSIS BY PANGENOME CLASS")
    print("="*80)
    
    all_results = []
    
    for species in species_list:
        if species not in og_consensus_tables:
            continue
        
        og = og_consensus_tables[species]
        print(f"\n{'='*60}")
        print(f"A. {species}")
        print(f"{'='*60}")
        
        # Overall secretome statistics
        n_total = len(og)
        n_signalp = og["is_secreted_signalp"].sum()
        n_go_sec = og["has_secretion_GO"].sum()
        n_combined = og["is_secreted"].sum()
        
        print(f"\nSecretome Detection Summary:")
        print(f"  Total orthogroups: {n_total:,}")
        print(f"  SignalP-predicted secreted: {n_signalp:,} ({100*n_signalp/n_total:.1f}%)")
        print(f"  Secretion GO terms: {n_go_sec:,} ({100*n_go_sec/n_total:.1f}%)")
        print(f"  Combined secretome: {n_combined:,} ({100*n_combined/n_total:.1f}%)")
        
        # By pangenome class
        print(f"\nSecretome by Pangenome Class:")
        for cls in ["Core", "Accessory"]:
            subset = og[og["Pangenome_Class"] == cls]
            n_cls = len(subset)
            n_sec = subset["is_secreted"].sum()
            pct = 100 * n_sec / n_cls if n_cls > 0 else 0
            print(f"  {cls}: {n_sec:,} / {n_cls:,} ({pct:.1f}%)")
            
            all_results.append({
                "species": species,
                "pangenome_class": cls,
                "total_ogs": n_cls,
                "secreted_ogs": n_sec,
                "pct_secreted": pct
            })
        
        # Fisher's exact tests
        print(f"\nStatistical Tests:")
        core = og[og["Pangenome_Class"] == "Core"]
        acc = og[og["Pangenome_Class"] == "Accessory"]
        
        # Accessory vs Core
        if len(acc) > 0 and len(core) > 0:
            a = acc["is_secreted"].sum()
            b = len(acc) - a
            c = core["is_secreted"].sum()
            d = len(core) - c
            
            try:
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
                sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                direction = "ENRICHED" if odds_ratio > 1 else "DEPLETED"
                print(f"  Accessory vs Core: OR={odds_ratio:.2f}, p={p_value:.2e} {sig}")
                if p_value < 0.05:
                    print(f"    -> Secreted proteins {direction} in Accessory genome")
            except:
                pass
    
    return pd.DataFrame(all_results)


def plot_secretome_by_pangenome(og_consensus_tables, species_list):
    """
    Visualize secretome distribution across species and pangenome classes.
    """
    from scipy.stats import fisher_exact
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # 1. Secretome % by pangenome class
    ax1 = axes[0]
    classes = ["Core", "Accessory"]
    x = np.arange(len(species_list))
    width = 0.35
    colors_class = {"Core": "#2ca02c", "Accessory": "#ff7f0e"}
    
    for i, cls in enumerate(classes):
        pcts = []
        for species in species_list:
            if species in og_consensus_tables:
                og = og_consensus_tables[species]
                subset = og[og["Pangenome_Class"] == cls]
                pct = 100 * subset["is_secreted"].sum() / len(subset) if len(subset) > 0 else 0
                pcts.append(pct)
            else:
                pcts.append(0)
        
        ax1.bar(x + i*width, pcts, width, label=cls, color=colors_class[cls])
    
    ax1.set_xlabel("Species")
    ax1.set_ylabel("% Orthogroups Secreted")
    ax1.set_title("Secretome by Pangenome Class", fontsize=12, fontweight="bold")
    ax1.set_xticks(x + width/2)
    ax1.set_xticklabels([f"A. {sp}" for sp in species_list], rotation=15)
    ax1.legend()
    ax1.grid(axis="y", alpha=0.3)
    
    # 2. SignalP vs GO-based detection
    ax2 = axes[1]
    detection_data = {"SignalP": [], "GO terms": [], "Combined": []}
    
    for species in species_list:
        if species in og_consensus_tables:
            og = og_consensus_tables[species]
            n_total = len(og)
            detection_data["SignalP"].append(100 * og["is_secreted_signalp"].sum() / n_total)
            detection_data["GO terms"].append(100 * og["has_secretion_GO"].sum() / n_total)
            detection_data["Combined"].append(100 * og["is_secreted"].sum() / n_total)
        else:
            for key in detection_data:
                detection_data[key].append(0)
    
    x2 = np.arange(len(species_list))
    width2 = 0.25
    colors2 = ["#1f77b4", "#9467bd", "#17becf"]
    
    for i, (method, pcts) in enumerate(detection_data.items()):
        ax2.bar(x2 + i*width2, pcts, width2, label=method, color=colors2[i])
    
    ax2.set_xlabel("Species")
    ax2.set_ylabel("% Orthogroups")
    ax2.set_title("Secretome Detection Methods", fontsize=12, fontweight="bold")
    ax2.set_xticks(x2 + width2)
    ax2.set_xticklabels([f"A. {sp}" for sp in species_list], rotation=15)
    ax2.legend()
    ax2.grid(axis="y", alpha=0.3)
    
    # 3. Accessory enrichment across species
    ax3 = axes[2]
    acc_enrichment = []
    species_names = []
    p_values = []
    
    for species in species_list:
        if species in og_consensus_tables:
            og = og_consensus_tables[species]
            core = og[og["Pangenome_Class"] == "Core"]
            acc = og[og["Pangenome_Class"] == "Accessory"]
            
            if len(acc) > 0 and len(core) > 0:
                a = acc["is_secreted"].sum()
                b = len(acc) - a
                c = core["is_secreted"].sum()
                d = len(core) - c
                
                try:
                    odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
                    acc_enrichment.append(np.log2(odds_ratio) if odds_ratio > 0 else 0)
                    p_values.append(p_value)
                    species_names.append(f"A. {species}")
                except:
                    pass
    
    if acc_enrichment:
        colors3 = ["#2ca02c" if x > 0 else "#d62728" for x in acc_enrichment]
        bars = ax3.bar(species_names, acc_enrichment, color=colors3, alpha=0.7)
        
        # Add significance markers
        for i, (bar, pval) in enumerate(zip(bars, p_values)):
            if pval < 0.001:
                marker = "***"
            elif pval < 0.01:
                marker = "**"
            elif pval < 0.05:
                marker = "*"
            else:
                marker = ""
            
            y_pos = bar.get_height() + 0.05 if bar.get_height() > 0 else bar.get_height() - 0.15
            ax3.text(bar.get_x() + bar.get_width()/2, y_pos, marker, 
                    ha="center", va="bottom", fontsize=12, fontweight="bold")
        
        ax3.axhline(0, color="black", linewidth=0.8)
        ax3.set_ylabel("log2(OR) Accessory vs Core")
        ax3.set_title("Secretome Enrichment in Accessory", fontsize=12, fontweight="bold")
        ax3.set_xticklabels(species_names, rotation=15)
        ax3.grid(axis="y", alpha=0.3)
    
    plt.suptitle("Secretome Analysis Across Aspergillus Species", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    


def analyze_secreted_functions_consensus(og_consensus_tables, species_list):
    """
    What functional categories are enriched in secreted orthogroups?
    """
    from scipy.stats import fisher_exact
    
    print("\n" + "="*80)
    print("FUNCTIONAL ENRICHMENT IN SECRETED ORTHOGROUPS")
    print("="*80)
    
    for species in species_list:
        if species not in og_consensus_tables:
            continue
        
        og = og_consensus_tables[species]
        secreted = og[og["is_secreted"] == True]
        non_secreted = og[og["is_secreted"] == False]
        
        if len(secreted) < 10:
            continue
        
        print(f"\nA. {species} (n={len(secreted)} secreted orthogroups):")
        
        # Check CAZyme enrichment in secreted
        if "CAZy" in og.columns:
            sec_cazy = (secreted["CAZy"].fillna("").astype(str) != "").sum()
            nonsec_cazy = (non_secreted["CAZy"].fillna("").astype(str) != "").sum()
            
            pct_sec = 100 * sec_cazy / len(secreted)
            pct_nonsec = 100 * nonsec_cazy / len(non_secreted) if len(non_secreted) > 0 else 0
            
            # Fisher's test
            try:
                a, b = sec_cazy, len(secreted) - sec_cazy
                c, d = nonsec_cazy, len(non_secreted) - nonsec_cazy
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
                sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                print(f"  CAZymes: {pct_sec:.1f}% secreted vs {pct_nonsec:.1f}% non-secreted (OR={odds_ratio:.2f} {sig})")
            except:
                print(f"  CAZymes: {pct_sec:.1f}% secreted vs {pct_nonsec:.1f}% non-secreted")
        
        # Check protease enrichment in secreted
        if "is_protease" in og.columns:
            sec_prot = secreted["is_protease"].sum()
            nonsec_prot = non_secreted["is_protease"].sum()
            
            pct_sec = 100 * sec_prot / len(secreted)
            pct_nonsec = 100 * nonsec_prot / len(non_secreted) if len(non_secreted) > 0 else 0
            
            try:
                a, b = sec_prot, len(secreted) - sec_prot
                c, d = nonsec_prot, len(non_secreted) - nonsec_prot
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
                sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                print(f"  Proteases: {pct_sec:.1f}% secreted vs {pct_nonsec:.1f}% non-secreted (OR={odds_ratio:.2f} {sig})")
            except:
                print(f"  Proteases: {pct_sec:.1f}% secreted vs {pct_nonsec:.1f}% non-secreted")


def analyze_bgc_for_species(species, species_root_map, bigscape_base, similarity_cutoff="c0.3", ani_excluded=None):
    """
    Analyze BGC (Biosynthetic Gene Cluster) distribution for a species.
    Parses antiSMASH regions and links to BiG-SCAPE GCFs.
    """
    print(f"\n{'='*60}")
    print(f"BGC Analysis for A. {species}")
    print(f"{'='*60}")
    
    species_root = species_root_map.get(species)
    if not species_root or not species_root.exists():
        print(f"  Species root not found: {species_root}")
        return None, None
    
    antismash_dir = species_root / "antismash_output"
    bigscape_output = species_root / "bigscape_output"
    
    # Check for antiSMASH output
    if not antismash_dir.exists():
        print(f"  No antiSMASH output found at {antismash_dir}")
        return None, None
    
    # Build set of ANI-excluded accession base patterns
    import re as _re
    _exc_bases = set()
    if ani_excluded:
        for exc_acc in ani_excluded:
            _exc_bases.add(_re.match(r'(GC[AF]_\d+)', exc_acc).group(1))

    # Find assemblies from antismash output
    print("  Scanning antiSMASH output for assemblies...")
    assemblies = []

    # Try antismash_gbk_unique folder first
    unique_dir = antismash_dir / "antismash_gbk_unique"
    if unique_dir.exists():
        for f in os.listdir(unique_dir):
            if f.endswith(".gbk") and ".region" in f:
                # Extract assembly accession from filename
                acc = f.split(".")[0]
                if acc.startswith("GC") and acc not in assemblies:
                    if not any(exc in acc for exc in _exc_bases):
                        assemblies.append(acc)

    # Also check for assembly subdirectories
    if not assemblies:
        for item in os.listdir(antismash_dir):
            item_path = antismash_dir / item
            if os.path.isdir(item_path) and item.startswith("GC"):
                if not any(exc in item for exc in _exc_bases):
                    assemblies.append(item)
    
    if not assemblies:
        print(f"  No assemblies found in {antismash_dir}")
        return None, None
    
    print(f"  Found {len(assemblies)} assemblies")
    
    # Create a minimal samples dataframe for region parsing
    samples_df = pd.DataFrame({
        "Species": species,
        "Assembly Accession": assemblies,
        "UsageClass": "unknown"  # Placeholder, not used for analysis
    })
    
    # Build regions from antiSMASH output
    regions = build_regions_from_df(
        samples_df, 
        {species: species_root},
        antismash_subdir="antismash_output"
    )
    
    if regions.empty:
        print(f"  No antiSMASH regions parsed")
        return None, None
    
    print(f"  Parsed {len(regions)} BGC regions from {regions['strain'].nunique()} genomes")
    
    # Add antiSMASH labels
    regions = add_antismash_direct_labels(regions)
    
    # Summarize BGC types
    bgc_counts = regions["antismash_primary_collapsed"].value_counts()
    print(f"\n  BGC Type Distribution:")
    for bgc_type, count in bgc_counts.head(10).items():
        print(f"    {bgc_type}: {count}")
    
    # Calculate BGCs per genome
    bgcs_per_genome = regions.groupby("strain").size()
    print(f"\n  BGCs per genome: mean={bgcs_per_genome.mean():.1f}, min={bgcs_per_genome.min()}, max={bgcs_per_genome.max()}")
    
    # Load BiG-SCAPE results if available
    bigscape_df = None
    if bigscape_output.exists():
        print(f"\n  Loading BiG-SCAPE results...")
        bigscape_df = load_bigscape_full(bigscape_base, [species], similarity_cutoff)
        
        if bigscape_df is not None and len(bigscape_df) > 0:
            print(f"    BiG-SCAPE records: {len(bigscape_df)}")
            n_gcfs = bigscape_df["Family_ID"].nunique()
            print(f"    Gene Cluster Families (GCFs): {n_gcfs}")
            
            # Join regions with GCFs
            regions = join_regions_with_bigscape_full(regions, bigscape_df)
            n_with_gcf = regions["GCF"].notna().sum()
            print(f"    Regions with GCF assignment: {n_with_gcf} / {len(regions)}")
        else:
            print(f"    No BiG-SCAPE clustering found for cutoff {similarity_cutoff}")
    else:
        print(f"  No BiG-SCAPE output at {bigscape_output}")
    
    return regions, bigscape_df


def plot_bgc_summary(all_regions, species_list):
    """
    Summarize BGC distribution across species.
    """
    if not all_regions:
        print("No regions data available")
        return None

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    species_colors = sns.color_palette("Set2", len(species_list))
    species_labels = [f"A. {sp}" for sp in species_list]

    # --- Collect BGC type data for both panels 1 & 2 ---
    all_types = set()
    type_data = {}
    for sp in species_list:
        if sp in all_regions:
            counts = all_regions[sp]["antismash_primary_collapsed"].value_counts()
            type_data[sp] = counts
            all_types.update(counts.index)
    type_totals = {t: sum(type_data.get(sp, pd.Series()).get(t, 0)
                          for sp in species_list) for t in all_types}
    top_types = sorted(type_totals.keys(), key=lambda x: type_totals[x], reverse=True)
    type_colors = dict(zip(top_types, vivid_colors(len(top_types))))

    # ---- 1. BGC composition – proportional stacked bar ----
    ax1 = axes[0, 0]
    for i, sp in enumerate(species_list):
        if sp not in type_data:
            continue
        total = type_data[sp].sum()
        bottom = 0.0
        for bgc_type in top_types:
            val = type_data[sp].get(bgc_type, 0)
            pct = val / total if total else 0
            ax1.bar(i, pct, bottom=bottom, color=type_colors[bgc_type], width=0.6,
                    edgecolor='white', linewidth=0.3)
            if pct >= 0.05:
                ax1.text(i, bottom + pct / 2, f"{pct:.0%}", ha="center", va="center",
                         fontsize=7, color="white", fontweight="bold")
            bottom += pct
        ax1.text(i, 1.01, f"n={total}", ha="center", va="bottom", fontsize=8)
    ax1.set_xticks(range(len(species_list)))
    ax1.set_xticklabels(species_labels, rotation=15, fontstyle='italic')
    ax1.set_ylabel("Proportion of BGC Regions")
    ax1.set_ylim(0, 1.10)
    ax1.set_title("BGC Type Composition by Species", fontsize=12, fontweight="bold")
    # legend
    from matplotlib.patches import Patch
    legend_patches = [Patch(facecolor=type_colors[t], label=t) for t in top_types if type_totals[t] > 0]
    ax1.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)

    # ---- 2. GCF composition – proportional stacked bar ----
    ax2 = axes[0, 1]
    gcf_type_data = {}
    gcf_all_types = set()
    for sp in species_list:
        if sp in all_regions and "GCF" in all_regions[sp].columns and "antismash_primary_collapsed" in all_regions[sp].columns:
            rdf = all_regions[sp].dropna(subset=["GCF"])
            # Get one class per GCF (majority vote)
            gcf_class = rdf.groupby("GCF")["antismash_primary_collapsed"].agg(lambda x: x.value_counts().index[0])
            counts = gcf_class.value_counts()
            gcf_type_data[sp] = counts
            gcf_all_types.update(counts.index)
    gcf_type_totals = {t: sum(gcf_type_data.get(sp, pd.Series()).get(t, 0)
                              for sp in species_list) for t in gcf_all_types}
    gcf_top_types = sorted(gcf_type_totals.keys(), key=lambda x: gcf_type_totals[x], reverse=True)
    gcf_type_colors = dict(zip(gcf_top_types, vivid_colors(len(gcf_top_types))))

    for i, sp in enumerate(species_list):
        if sp not in gcf_type_data:
            continue
        total = gcf_type_data[sp].sum()
        bottom = 0.0
        for gcf_type in gcf_top_types:
            val = gcf_type_data[sp].get(gcf_type, 0)
            pct = val / total if total else 0
            ax2.bar(i, pct, bottom=bottom, color=gcf_type_colors[gcf_type], width=0.6,
                    edgecolor='white', linewidth=0.3)
            if pct >= 0.05:
                ax2.text(i, bottom + pct / 2, f"{pct:.0%}", ha="center", va="center",
                         fontsize=7, color="white", fontweight="bold")
            bottom += pct
        ax2.text(i, 1.01, f"n={total}", ha="center", va="bottom", fontsize=8)
    ax2.set_xticks(range(len(species_list)))
    ax2.set_xticklabels(species_labels, rotation=15, fontstyle='italic')
    ax2.set_ylabel("Proportion of GCFs")
    ax2.set_ylim(0, 1.10)
    ax2.set_title("GCF Type Composition by Species", fontsize=12, fontweight="bold")
    gcf_legend_patches = [Patch(facecolor=gcf_type_colors[t], label=t) for t in gcf_top_types if gcf_type_totals[t] > 0]
    ax2.legend(handles=gcf_legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)

    # ---- 3. BGC per genome distribution (boxplot with outlier annotations) ----
    ax3 = axes[1, 0]
    bgc_per_genome_data = []
    bgc_per_genome_series = {}
    for sp in species_list:
        if sp in all_regions:
            per_genome = all_regions[sp].groupby("strain").size()
            bgc_per_genome_data.append(per_genome.values if len(per_genome) > 0 else np.array([0]))
            bgc_per_genome_series[sp] = per_genome
        else:
            bgc_per_genome_data.append(np.array([0]))
            bgc_per_genome_series[sp] = pd.Series(dtype=int)

    bp = ax3.boxplot(bgc_per_genome_data, labels=species_labels, patch_artist=True,
                     flierprops=dict(marker='o', markersize=5, alpha=0.7))
    for patch, color in zip(bp['boxes'], species_colors):
        patch.set_facecolor(color)

    # Annotate outliers with strain accession
    for idx, sp in enumerate(species_list):
        per_genome = bgc_per_genome_series[sp]
        if len(per_genome) == 0:
            continue
        q1 = np.percentile(per_genome.values, 25)
        q3 = np.percentile(per_genome.values, 75)
        iqr = q3 - q1
        lo_fence = q1 - 1.5 * iqr
        hi_fence = q3 + 1.5 * iqr
        outliers = per_genome[(per_genome < lo_fence) | (per_genome > hi_fence)]
        for strain, val in outliers.items():
            # Shorten accession for label
            short = strain.split("_")[-1] if "_" in str(strain) else str(strain)
            ax3.annotate(short, (idx + 1, val), fontsize=5.5, color='red',
                         ha='left', va='center', xytext=(5, 0),
                         textcoords='offset points')

    ax3.set_ylabel("BGCs per Genome")
    ax3.set_title("BGC Count Distribution per Genome", fontsize=12, fontweight="bold")
    ax3.set_xticklabels(species_labels, rotation=15, fontstyle='italic')
    ax3.grid(axis="y", alpha=0.3)

    # ---- 4. Total BGC and GCF counts side by side ----
    ax4 = axes[1, 1]
    bgc_totals = [len(all_regions[sp]) if sp in all_regions else 0 for sp in species_list]
    gcf_totals = []
    for sp in species_list:
        if sp in all_regions and "GCF" in all_regions[sp].columns:
            gcf_totals.append(all_regions[sp]["GCF"].nunique())
        else:
            gcf_totals.append(0)
    x = np.arange(len(species_list))
    w = 0.35
    bars1 = ax4.bar(x - w/2, bgc_totals, w, label="BGC Regions", color='steelblue', edgecolor='white')
    bars2 = ax4.bar(x + w/2, gcf_totals, w, label="GCFs", color='coral', edgecolor='white')
    for bar, val in zip(bars1, bgc_totals):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, str(val),
                 ha="center", va="bottom", fontsize=9, fontweight="bold")
    for bar, val in zip(bars2, gcf_totals):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, str(val),
                 ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax4.set_xticks(x)
    ax4.set_xticklabels(species_labels, rotation=15, fontstyle='italic')
    ax4.set_ylabel("Count")
    ax4.set_title("Total BGC Regions and GCFs by Species", fontsize=12, fontweight="bold")
    ax4.legend()

    plt.suptitle("Secondary Metabolite Gene Cluster Analysis", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()



def summarize_gcf_diversity(all_regions, all_bigscape, species_list):
    """
    Summarize GCF diversity and uniqueness across species.
    """
    print("\n" + "="*80)
    print("GCF DIVERSITY SUMMARY")
    print("="*80)
    
    for species in species_list:
        if species not in all_regions:
            continue
            
        regions = all_regions[species]
        if "GCF" not in regions.columns or regions["GCF"].isna().all():
            continue
        
        print(f"\nA. {species}:")
        
        # GCF statistics
        gcf_counts = regions["GCF"].value_counts()
        n_gcfs = len(gcf_counts)
        n_singletons = (gcf_counts == 1).sum()
        largest_gcf = gcf_counts.max()
        
        print(f"  Total GCFs: {n_gcfs}")
        print(f"  Singleton GCFs (1 member): {n_singletons} ({100*n_singletons/n_gcfs:.1f}%)")
        print(f"  Largest GCF size: {largest_gcf} members")
        
        # BGC type breakdown in GCFs
        gcf_types = regions.groupby("GCF")["antismash_primary_collapsed"].agg(
            lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else "Unknown"
        )
        type_counts = gcf_types.value_counts()
        print(f"  GCFs by BGC type:")
        for bgc_type, count in type_counts.head(5).items():
            print(f"    {bgc_type}: {count} GCFs")


# =============================================================================
# SNP PCA ANALYSIS (from Gubbins output)
# =============================================================================

def load_snp_alignment(fasta_path: str) -> Dict[str, str]:
    """
    Load SNP alignment from Gubbins filtered polymorphic sites FASTA.

    Parameters
    ----------
    fasta_path : str
        Path to Gubbins .filtered_polymorphic_sites.fasta file

    Returns
    -------
    dict
        Dictionary mapping sample_id -> SNP sequence string
    """
    sequences = {}
    if not os.path.exists(fasta_path):
        print(f"[WARN] SNP alignment file not found: {fasta_path}")
        return sequences

    with open(fasta_path, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                # Extract sample ID from header
                header = line[1:]
                # Remove common suffixes
                sample_id = header.replace('_genomic_renamed.fa.ref', '').replace('_genomic_renamed.fa', '')
                sample_id = sample_id.replace('_genomic.fa', '').replace('.fa', '')
                current_id = sample_id
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    print(f"Loaded SNP alignment: {len(sequences)} samples, {len(sequences[list(sequences.keys())[0]]) if sequences else 0} SNP sites")
    return sequences


def snp_alignment_to_numeric(sequences: Dict[str, str]) -> pd.DataFrame:
    """
    Convert SNP alignment to numeric matrix for PCA.

    Encoding: A=0, C=1, G=2, T=3, N/-=np.nan

    Parameters
    ----------
    sequences : dict
        Dictionary mapping sample_id -> SNP sequence string

    Returns
    -------
    pd.DataFrame
        Numeric matrix (samples x SNP positions)
    """
    if not sequences:
        return pd.DataFrame()

    encoding = {'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3}

    sample_ids = list(sequences.keys())
    seq_len = len(sequences[sample_ids[0]])

    matrix = np.zeros((len(sample_ids), seq_len), dtype=np.float32)
    matrix[:] = np.nan

    for i, sample_id in enumerate(sample_ids):
        seq = sequences[sample_id]
        for j, base in enumerate(seq):
            if base in encoding:
                matrix[i, j] = encoding[base]

    df = pd.DataFrame(matrix, index=sample_ids)

    # Remove positions with all NaN or no variation
    valid_cols = []
    for col in range(df.shape[1]):
        col_data = df.iloc[:, col].dropna()
        if len(col_data) > 0 and col_data.nunique() > 1:
            valid_cols.append(col)

    df_filtered = df.iloc[:, valid_cols]
    print(f"SNP matrix: {df_filtered.shape[0]} samples x {df_filtered.shape[1]} informative SNPs")

    return df_filtered


def compute_snp_pca(snp_matrix: pd.DataFrame, n_components: int = 10) -> pd.DataFrame:
    """
    Compute PCA on SNP matrix.

    Parameters
    ----------
    snp_matrix : pd.DataFrame
        Numeric SNP matrix (samples x SNPs)
    n_components : int
        Number of PC components to retain (default 10)

    Returns
    -------
    pd.DataFrame
        PCA scores with columns: sample_id, PC1, PC2, ..., PCn
    """
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer

    if snp_matrix.empty:
        return pd.DataFrame()

    # Impute missing values with column mean
    imputer = SimpleImputer(strategy='mean')
    X = imputer.fit_transform(snp_matrix.values)

    # Standardize
    X_centered = X - X.mean(axis=0)
    X_std = X_centered / (X.std(axis=0) + 1e-10)

    # PCA
    n_comp = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_comp)
    pcs = pca.fit_transform(X_std)

    # Create DataFrame
    pc_cols = [f'PC{i+1}' for i in range(n_comp)]
    pcs_df = pd.DataFrame(pcs, columns=pc_cols, index=snp_matrix.index)
    pcs_df = pcs_df.reset_index().rename(columns={'index': 'sample_id'})

    # Store explained variance ratio
    pcs_df.attrs['explained_variance_ratio'] = pca.explained_variance_ratio_

    print(f"PCA computed: {n_comp} components")
    print(f"Variance explained: PC1={pca.explained_variance_ratio_[0]*100:.1f}%, PC2={pca.explained_variance_ratio_[1]*100:.1f}%")

    return pcs_df


def plot_snp_pca(pcs_df: pd.DataFrame, metadata_df: pd.DataFrame = None,
                 color_column: str = 'isolation_source', species: str = None) -> plt.Figure:
    """
    Plot PCA scatter colored by metadata label.

    Parameters
    ----------
    pcs_df : pd.DataFrame
        PCA results from compute_snp_pca()
    metadata_df : pd.DataFrame
        Metadata with sample_id and color_column
    color_column : str
        Column to use for coloring points
    species : str
        Species name for title

    Returns
    -------
    matplotlib.Figure
    """
    if pcs_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No PCA data available", ha='center', va='center', fontsize=12)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Merge with metadata if provided
    if metadata_df is not None and color_column in metadata_df.columns:
        # Try to match sample IDs
        pcs_plot = pcs_df.copy()
        pcs_plot['label'] = 'Unknown'

        for idx, row in pcs_plot.iterrows():
            sample_id = row['sample_id']
            # Try different matching strategies
            match = metadata_df[metadata_df['Assembly'].str.contains(sample_id.split('_')[0], case=False, na=False)]
            if len(match) > 0:
                pcs_plot.loc[idx, 'label'] = match.iloc[0][color_column]
    else:
        pcs_plot = pcs_df.copy()
        pcs_plot['label'] = 'Unknown'

    # Get unique labels and colors
    labels = pcs_plot['label'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(labels)))
    color_map = dict(zip(labels, colors))

    # Get variance explained
    var_exp = pcs_df.attrs.get('explained_variance_ratio', [0.1, 0.1])

    # Plot 1: PC1 vs PC2
    ax1 = axes[0]
    for label in labels:
        mask = pcs_plot['label'] == label
        ax1.scatter(pcs_plot.loc[mask, 'PC1'], pcs_plot.loc[mask, 'PC2'],
                   c=[color_map[label]], label=label, alpha=0.7, s=50, edgecolors='white', linewidth=0.5)

    ax1.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}% variance)')
    ax1.set_ylabel(f'PC2 ({var_exp[1]*100:.1f}% variance)')
    ax1.set_title(f'SNP PCA - {"A. " + species if species else "All Samples"}')
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax1.grid(alpha=0.3)

    # Plot 2: PC1 vs PC3 (if available)
    ax2 = axes[1]
    if 'PC3' in pcs_plot.columns:
        for label in labels:
            mask = pcs_plot['label'] == label
            ax2.scatter(pcs_plot.loc[mask, 'PC1'], pcs_plot.loc[mask, 'PC3'],
                       c=[color_map[label]], label=label, alpha=0.7, s=50, edgecolors='white', linewidth=0.5)
        ax2.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}% variance)')
        ax2.set_ylabel(f'PC3 ({var_exp[2]*100:.1f}% variance)')
        ax2.set_title(f'SNP PCA (PC1 vs PC3)')
    else:
        # Scree plot instead
        ax2.bar(range(1, len(var_exp)+1), var_exp * 100, color='steelblue')
        ax2.set_xlabel('Principal Component')
        ax2.set_ylabel('Variance Explained (%)')
        ax2.set_title('Scree Plot')
    ax2.grid(alpha=0.3)

    plt.tight_layout()


def run_snp_pca_analysis(gubbins_dir: str, metadata_df: pd.DataFrame = None,
                         species: str = None, n_components: int = 10,
                         exclude_patterns: set = None) -> tuple:
    """
    Complete SNP PCA analysis pipeline.

    Parameters
    ----------
    gubbins_dir : str
        Path to Gubbins output directory
    metadata_df : pd.DataFrame
        Metadata DataFrame with isolation_source column
    species : str
        Species name
    n_components : int
        Number of PCs to compute
    exclude_patterns : set
        Set of accession base patterns (e.g., {'GCA_023653635'}) to exclude

    Returns
    -------
    tuple
        (pcs_df, snp_matrix, fig)
    """
    # Find the filtered polymorphic sites file
    fasta_path = os.path.join(gubbins_dir, 'gubbins.filtered_polymorphic_sites.fasta')

    if not os.path.exists(fasta_path):
        print(f"[WARN] Gubbins SNP file not found: {fasta_path}")
        return pd.DataFrame(), pd.DataFrame(), None

    # Load and process
    sequences = load_snp_alignment(fasta_path)

    # Remove excluded samples before computing PCA
    if exclude_patterns:
        before = len(sequences)
        sequences = {k: v for k, v in sequences.items()
                     if not any(pat in k for pat in exclude_patterns)}
        removed = before - len(sequences)
        if removed > 0:
            print(f"  ANI excluded: {removed} sample(s) removed before PCA")

    snp_matrix = snp_alignment_to_numeric(sequences)
    pcs_df = compute_snp_pca(snp_matrix, n_components=n_components)

    # Plot
    fig = plot_snp_pca(pcs_df, metadata_df, species=species)

    return pcs_df, snp_matrix, fig


# =============================================================================
# SNP KINSHIP / GENETIC RELATEDNESS MATRIX (GRM)
# =============================================================================

def compute_grm(snp_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Genetic Relatedness Matrix (GRM) from SNP data.

    Uses identity-by-state (IBS) similarity.

    Parameters
    ----------
    snp_matrix : pd.DataFrame
        Numeric SNP matrix (samples x SNPs)

    Returns
    -------
    pd.DataFrame
        NxN kinship matrix
    """
    from sklearn.impute import SimpleImputer

    if snp_matrix.empty:
        return pd.DataFrame()

    # Impute missing values
    imputer = SimpleImputer(strategy='most_frequent')
    X = imputer.fit_transform(snp_matrix.values)

    # Standardize (center and scale)
    X_centered = X - X.mean(axis=0)
    X_std = X_centered / (X.std(axis=0) + 1e-10)

    # Compute GRM: (1/p) * X * X'
    # This is the realized relationship matrix
    n_snps = X_std.shape[1]
    grm = np.dot(X_std, X_std.T) / n_snps

    # Convert to DataFrame
    grm_df = pd.DataFrame(grm, index=snp_matrix.index, columns=snp_matrix.index)

    print(f"GRM computed: {grm_df.shape[0]} x {grm_df.shape[1]} matrix")

    return grm_df


def plot_grm(grm_df: pd.DataFrame, species: str = None) -> plt.Figure:
    """
    Visualize the Genetic Relatedness Matrix as a heatmap.

    Parameters
    ----------
    grm_df : pd.DataFrame
        GRM from compute_grm()
    species : str
        Species name for title

    Returns
    -------
    matplotlib.Figure
    """
    if grm_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No GRM data available", ha='center', va='center', fontsize=12)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Heatmap
    ax1 = axes[0]

    # Cluster for better visualization (if scipy available)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        # Convert similarity to distance
        dist = 1 - grm_df.values
        np.fill_diagonal(dist, 0)
        dist = np.clip(dist, 0, None)  # Ensure non-negative

        # Hierarchical clustering
        condensed = squareform(dist)
        Z = linkage(condensed, method='average')
        order = leaves_list(Z)

        grm_ordered = grm_df.iloc[order, order]
    except Exception:
        grm_ordered = grm_df

    im = ax1.imshow(grm_ordered.values, cmap='RdBu_r', aspect='auto', vmin=-0.5, vmax=1)
    plt.colorbar(im, ax=ax1, label='Relatedness')
    ax1.set_title(f'Genetic Relatedness Matrix - {"A. " + species if species else "All Samples"}', fontweight='bold')
    ax1.set_xlabel('Samples')
    ax1.set_ylabel('Samples')

    # Distribution of off-diagonal values
    ax2 = axes[1]
    off_diag = grm_df.values[np.triu_indices_from(grm_df.values, k=1)]
    ax2.hist(off_diag, bins=50, color='steelblue', edgecolor='white', alpha=0.7)
    ax2.axvline(x=0, color='red', linestyle='--', label='Expected unrelated')
    ax2.axvline(x=np.mean(off_diag), color='orange', linestyle='-', label=f'Mean: {np.mean(off_diag):.3f}')
    ax2.set_xlabel('Relatedness')
    ax2.set_ylabel('Count')
    ax2.set_title('Distribution of Pairwise Relatedness')
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()


def run_grm_analysis(snp_matrix: pd.DataFrame, species: str = None,
                     output_path: str = None) -> tuple:
    """
    Complete GRM analysis pipeline.

    Parameters
    ----------
    snp_matrix : pd.DataFrame
        SNP matrix from snp_alignment_to_numeric()
    species : str
        Species name
    output_path : str
        Path to save kinship.tsv (optional)

    Returns
    -------
    tuple
        (grm_df, fig)
    """
    grm_df = compute_grm(snp_matrix)

    if output_path and not grm_df.empty:
        grm_df.to_csv(output_path, sep='\t')
        print(f"GRM saved to: {output_path}")

    fig = plot_grm(grm_df, species=species)

    return grm_df, fig


# =============================================================================
# BGC/GCF PRESENCE-ABSENCE MATRIX
# =============================================================================

def load_bigscape_clustering(bigscape_dir: str) -> pd.DataFrame:
    """
    Load BiG-SCAPE GCF clustering results.

    Parameters
    ----------
    bigscape_dir : str
        Path to BiG-SCAPE output directory

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: Record, GBK, strain, BGC_type, GCF
    """
    all_clusters = []

    # Find clustering files
    output_dirs = glob.glob(os.path.join(bigscape_dir, 'output_files', '*_c0.*'))

    for output_dir in output_dirs:
        if not os.path.isdir(output_dir):
            continue

        # Find all clustering TSV files
        for clustering_file in glob.glob(os.path.join(output_dir, '*', '*_clustering_*.tsv')):
            try:
                df = pd.read_csv(clustering_file, sep='\t')
                if 'Family' in df.columns:
                    # Extract BGC type from path
                    bgc_type = os.path.basename(os.path.dirname(clustering_file))
                    df['BGC_type'] = bgc_type
                    all_clusters.append(df)
            except Exception as e:
                continue

    if not all_clusters:
        print(f"[WARN] No clustering files found in {bigscape_dir}")
        return pd.DataFrame()

    combined = pd.concat(all_clusters, ignore_index=True)

    # Extract strain from GBK name
    combined['strain'] = combined['GBK'].apply(
        lambda x: '_'.join(x.split('_')[:2]) if pd.notna(x) else None
    )
    combined['GCF'] = combined['Family']

    # Filter out MIBiG reference BGCs (they start with "BGC" not "GCA"/"GCF")
    n_before = len(combined)
    is_actual_strain = combined['strain'].str.startswith(('GCA', 'GCF'), na=False)
    n_mibig = (~is_actual_strain).sum()
    combined = combined[is_actual_strain].copy()

    if n_mibig > 0:
        print(f"Filtered out {n_mibig} MIBiG reference BGCs")

    print(f"Loaded BiG-SCAPE clustering: {len(combined)} BGCs, {combined['GCF'].nunique()} GCFs, {combined['strain'].nunique()} strains")

    return combined


def build_bgc_pav_matrix(bigscape_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build BGC presence-absence matrix (strains x BGC regions).

    Parameters
    ----------
    bigscape_df : pd.DataFrame
        Output from load_bigscape_clustering()

    Returns
    -------
    pd.DataFrame
        Binary PAV matrix (strains x BGC regions)
    """
    if bigscape_df.empty:
        return pd.DataFrame()

    # Create pivot: strains x BGC Records
    pav = bigscape_df.pivot_table(
        index='strain',
        columns='Record',
        values='GCF',
        aggfunc='count',
        fill_value=0
    )

    # Convert to binary
    pav = (pav > 0).astype(int)

    print(f"BGC PAV matrix: {pav.shape[0]} strains x {pav.shape[1]} BGC regions")

    return pav


def build_gcf_pav_matrix(bigscape_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build GCF presence-absence matrix (strains x GCFs).

    Parameters
    ----------
    bigscape_df : pd.DataFrame
        Output from load_bigscape_clustering()

    Returns
    -------
    pd.DataFrame
        Binary PAV matrix (strains x GCFs)
    """
    if bigscape_df.empty:
        return pd.DataFrame()

    # Create pivot: strains x GCFs
    pav = bigscape_df.pivot_table(
        index='strain',
        columns='GCF',
        values='Record',
        aggfunc='count',
        fill_value=0
    )

    # Convert to binary
    pav = (pav > 0).astype(int)

    print(f"GCF PAV matrix: {pav.shape[0]} strains x {pav.shape[1]} GCFs")

    return pav


def build_gcf_cnv_matrix(bigscape_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build GCF copy number variation matrix (strains x GCFs).

    Parameters
    ----------
    bigscape_df : pd.DataFrame
        Output from load_bigscape_clustering()

    Returns
    -------
    pd.DataFrame
        CNV matrix (strains x GCFs) - counts of BGCs per GCF per strain
    """
    if bigscape_df.empty:
        return pd.DataFrame()

    # Create pivot: strains x GCFs with counts
    cnv = bigscape_df.pivot_table(
        index='strain',
        columns='GCF',
        values='Record',
        aggfunc='count',
        fill_value=0
    )

    print(f"GCF CNV matrix: {cnv.shape[0]} strains x {cnv.shape[1]} GCFs")

    return cnv


def plot_bgc_gcf_matrices(bgc_pav: pd.DataFrame, gcf_pav: pd.DataFrame,
                          gcf_cnv: pd.DataFrame, species: str = None,
                          core_n: int = None, rare_n: int = None) -> plt.Figure:
    """
    Visualize BGC and GCF matrices.

    Parameters
    ----------
    bgc_pav : pd.DataFrame
        BGC presence-absence matrix
    gcf_pav : pd.DataFrame
        GCF presence-absence matrix
    gcf_cnv : pd.DataFrame
        GCF copy number variation matrix
    species : str
        Species name for title

    Returns
    -------
    matplotlib.Figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    title_suffix = f" - A. {species}" if species else ""

    # 1. GCF PAV summary
    ax1 = axes[0, 0]
    if not gcf_pav.empty:
        gcf_freq = gcf_pav.sum(axis=0).sort_values(ascending=False)
        gcf_freq = gcf_freq[gcf_freq > 0]  # drop GCFs absent from all strains
        n_strains = gcf_pav.shape[0]

        # Classify: Core (>95%), Accessory (5-95%), Rare (<5%)
        core = (gcf_freq >= 0.95 * n_strains).sum()
        accessory = ((gcf_freq >= 0.05 * n_strains) & (gcf_freq < 0.95 * n_strains)).sum()
        rare = (gcf_freq < 0.05 * n_strains).sum()

        ax1.bar(['Core\n(≥95%)', 'Accessory\n(5-95%)', 'Rare\n(<5%)'],
                [core, accessory, rare],
                color=['#2ecc71', '#f39c12', '#e74c3c'])
        ax1.set_ylabel('Number of GCFs')
        ax1.set_title(f'GCF Pangenome Classification{title_suffix}', fontweight='bold')
        for i, v in enumerate([core, accessory, rare]):
            ax1.text(i, v + 0.5, str(v), ha='center', fontweight='bold')
    else:
        ax1.text(0.5, 0.5, "No GCF data", ha='center', va='center')

    # 2. GCF frequency distribution
    ax2 = axes[0, 1]
    if not gcf_pav.empty:
        gcf_freq = gcf_pav.sum(axis=0)
        ax2.hist(gcf_freq, bins=30, color='steelblue', edgecolor='white', alpha=0.7)
        ax2.axvline(x=gcf_freq.median(), color='red', linestyle='--', label=f'Median: {gcf_freq.median():.0f}')
        ax2.set_xlabel('Number of Strains with GCF')
        ax2.set_ylabel('Number of GCFs')
        ax2.set_title('GCF Frequency Distribution')
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, "No GCF data", ha='center', va='center')

    # 3. GCF CNV distribution
    ax3 = axes[1, 0]
    if not gcf_cnv.empty:
        # Get max copy number per strain
        max_copies = gcf_cnv.max(axis=1)
        ax3.hist(max_copies, bins=range(1, int(max_copies.max())+2), color='purple',
                 edgecolor='white', alpha=0.7, align='left')
        ax3.set_xlabel('Max GCF Copies per Strain')
        ax3.set_ylabel('Number of Strains')
        ax3.set_title('GCF Copy Number Distribution')
        ax3.axvline(x=max_copies.mean(), color='red', linestyle='--', label=f'Mean: {max_copies.mean():.1f}')
        ax3.legend()
    else:
        ax3.text(0.5, 0.5, "No CNV data", ha='center', va='center')

    # 4. BGCs per strain
    ax4 = axes[1, 1]
    if not bgc_pav.empty:
        bgcs_per_strain = bgc_pav.sum(axis=1).sort_values()
        ax4.barh(range(len(bgcs_per_strain)), bgcs_per_strain.values, color='teal', alpha=0.7)
        ax4.set_xlabel('Number of BGC Regions')
        ax4.set_ylabel('Strain (ranked)')
        ax4.set_title(f'BGCs per Strain (n={len(bgcs_per_strain)})')
        ax4.axvline(x=bgcs_per_strain.mean(), color='red', linestyle='--', label=f'Mean: {bgcs_per_strain.mean():.1f}')
        ax4.legend()
    else:
        ax4.text(0.5, 0.5, "No BGC data", ha='center', va='center')

    plt.suptitle(f'BGC/GCF Matrix Analysis{title_suffix}', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    return fig


def plot_gcf_pangenome_combined(bgc_gcf_results, species_list):
    """
    Combined species-comparing plot for GCF-based core/accessory/rare pangenome.

    Parameters
    ----------
    bgc_gcf_results : dict
        Keyed by species short name, values are dicts with 'gcf_pav' DataFrames.
    species_list : list of str
        Species short names in display order.
    """
    species_in = [sp for sp in species_list if sp in bgc_gcf_results
                  and not bgc_gcf_results[sp]['gcf_pav'].empty]
    if not species_in:
        print("No GCF PAV data available for combined plot")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sp_labels = [f"A. {sp}" for sp in species_in]
    cat_colors = {'Core': '#2ecc71', 'Accessory': '#f39c12', 'Rare': '#e74c3c'}
    sp_colors = sns.color_palette("Set2", len(species_in))

    # Collect pangenome stats
    stats = {}
    for sp in species_in:
        gcf_pav = bgc_gcf_results[sp]['gcf_pav']
        n_strains = gcf_pav.shape[0]
        freq = gcf_pav.sum(axis=0)
        freq = freq[freq > 0]  # drop GCFs absent from all strains
        core = int((freq >= 0.95 * n_strains).sum())
        acc = int(((freq >= 0.05 * n_strains) & (freq < 0.95 * n_strains)).sum())
        rare = int((freq < 0.05 * n_strains).sum())
        stats[sp] = {'Core': core, 'Accessory': acc, 'Rare': rare,
                     'total': core + acc + rare, 'n_strains': n_strains, 'freq': freq}

    # ---- Panel 1: Stacked bar (absolute counts) ----
    ax1 = axes[0]
    x = np.arange(len(species_in))
    w = 0.5
    bottoms = np.zeros(len(species_in))
    for cat in ['Core', 'Accessory', 'Rare']:
        vals = [stats[sp][cat] for sp in species_in]
        ax1.bar(x, vals, w, bottom=bottoms, label=cat, color=cat_colors[cat], edgecolor='white')
        for i, v in enumerate(vals):
            if v > 0:
                ax1.text(i, bottoms[i] + v / 2, str(v), ha='center', va='center',
                         fontsize=9, fontweight='bold', color='white')
        bottoms += vals
    for i, sp in enumerate(species_in):
        ax1.text(i, bottoms[i] + 1, f"n={stats[sp]['total']}", ha='center', va='bottom', fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(sp_labels, rotation=15, fontstyle='italic')
    ax1.set_ylabel('Number of GCFs')
    ax1.set_title('GCF Pangenome Classification', fontsize=12, fontweight='bold')
    ax1.legend()

    # ---- Panel 2: Proportional stacked bar ----
    ax2 = axes[1]
    bottoms = np.zeros(len(species_in))
    for cat in ['Core', 'Accessory', 'Rare']:
        vals = [stats[sp][cat] / stats[sp]['total'] if stats[sp]['total'] else 0 for sp in species_in]
        ax2.bar(x, vals, w, bottom=bottoms, color=cat_colors[cat], edgecolor='white')
        for i, v in enumerate(vals):
            if v >= 0.04:
                ax2.text(i, bottoms[i] + v / 2, f"{v:.0%}", ha='center', va='center',
                         fontsize=9, fontweight='bold', color='white')
        bottoms += vals
    ax2.set_xticks(x)
    ax2.set_xticklabels(sp_labels, rotation=15, fontstyle='italic')
    ax2.set_ylabel('Proportion of GCFs')
    ax2.set_ylim(0, 1.05)
    ax2.set_title('GCF Pangenome Proportion', fontsize=12, fontweight='bold')

    # ---- Panel 3: GCF frequency distribution (overlaid histograms) ----
    ax3 = axes[2]
    for i, sp in enumerate(species_in):
        freq = stats[sp]['freq']
        n_strains = stats[sp]['n_strains']
        freq_pct = freq / n_strains * 100
        ax3.hist(freq_pct, bins=20, alpha=0.5, label=f"A. {sp} (n={n_strains})",
                 color=sp_colors[i], edgecolor='white')
    ax3.axvline(x=95, color='#2ecc71', linestyle='--', linewidth=1, alpha=0.7)
    ax3.axvline(x=5, color='#e74c3c', linestyle='--', linewidth=1, alpha=0.7)
    ax3.text(96, ax3.get_ylim()[1] * 0.9, 'Core', fontsize=8, color='#2ecc71')
    ax3.text(6, ax3.get_ylim()[1] * 0.9, 'Rare', fontsize=8, color='#e74c3c')
    ax3.set_xlabel('GCF Prevalence (% of strains)')
    ax3.set_ylabel('Number of GCFs')
    ax3.set_title('GCF Frequency Distribution', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=8)

    plt.suptitle('GCF-Based Core/Accessory/Rare Pangenome Across Species',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()


def run_bgc_gcf_matrix_analysis(bigscape_dir: str, species: str = None,
                                 output_dir: str = None) -> dict:
    """
    Complete BGC/GCF matrix analysis pipeline.

    Parameters
    ----------
    bigscape_dir : str
        Path to BiG-SCAPE output directory
    species : str
        Species name
    output_dir : str
        Directory to save output matrices (optional)

    Returns
    -------
    dict
        Dictionary with keys: bigscape_df, bgc_pav, gcf_pav, gcf_cnv, fig
    """
    # Load clustering data
    bigscape_df = load_bigscape_clustering(bigscape_dir)

    if bigscape_df.empty:
        return {
            'bigscape_df': bigscape_df,
            'bgc_pav': pd.DataFrame(),
            'gcf_pav': pd.DataFrame(),
            'gcf_cnv': pd.DataFrame(),
            'fig': None
        }

    # Build matrices
    bgc_pav = build_bgc_pav_matrix(bigscape_df)
    gcf_pav = build_gcf_pav_matrix(bigscape_df)
    gcf_cnv = build_gcf_cnv_matrix(bigscape_df)

    # Save if output_dir specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        bgc_pav.to_csv(os.path.join(output_dir, f'{species}_bgc_pav.tsv'), sep='\t')
        gcf_pav.to_csv(os.path.join(output_dir, f'{species}_gcf_pav.tsv'), sep='\t')
        gcf_cnv.to_csv(os.path.join(output_dir, f'{species}_gcf_cnv.tsv'), sep='\t')
        print(f"Matrices saved to: {output_dir}")

    # Plot
    fig = plot_bgc_gcf_matrices(bgc_pav, gcf_pav, gcf_cnv, species=species)

    return {
        'bigscape_df': bigscape_df,
        'bgc_pav': bgc_pav,
        'gcf_pav': gcf_pav,
        'gcf_cnv': gcf_cnv,
        'fig': fig
    }


# =============================================================================
# NB1 PIPELINE FUNCTIONS (extracted from notebook cells)
# =============================================================================

def load_all_orthofinder(species_list, species_root, results_base, force=False):
    """Load OrthoFinder results and build PAV/CNV matrices for all species.

    On subsequent runs, loads cached PAV/CNV TSVs from ``results_base`` if they
    exist (use ``force=True`` to rebuild from OrthoFinder output).

    Returns
    -------
    ortho_data : dict
        ``{species: orthogroups_long_df}`` (None if loaded from cache).
    pav_data : dict
        ``{species: pav_matrix}``.
    cnv_data : dict
        ``{species: cnv_matrix}``.
    """
    ortho_data = {}
    pav_data = {}
    cnv_data = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_dir = os.path.join(results_base, sp)
        pav_path = os.path.join(out_dir, f'{sp}_pav.tsv')
        cnv_path = os.path.join(out_dir, f'{sp}_cnv.tsv')

        # --- Try to load from cache first ---
        if not force and os.path.exists(pav_path) and os.path.exists(cnv_path):
            pav = pd.read_csv(pav_path, sep='\t', index_col=0)
            cnv = pd.read_csv(cnv_path, sep='\t', index_col=0)
            # Apply runtime ANI exclusion (genomes are columns)
            pav = filter_ani_excluded(pav, sp, axis='cols')
            cnv = filter_ani_excluded(cnv, sp, axis='cols')
            pav_data[sp] = pav
            cnv_data[sp] = cnv
            ortho_data[sp] = None
            print(f'  Loaded cached PAV ({pav.shape}) and CNV ({cnv.shape}) from {out_dir}')
            continue

        # --- Rebuild from OrthoFinder output ---
        sp_dir = os.path.join(species_root, sp)
        og_paths = glob.glob(os.path.join(sp_dir, 'orthofinder_output', '*',
                                          'Orthogroups', 'Orthogroups.tsv'))
        if not og_paths:
            print(f'  WARNING: No OrthoFinder results for {sp}')
            continue
        og_tsv = sorted(og_paths)[-1]

        ortho_long = load_orthogroups_long(og_tsv, species_name=sp)
        ortho_data[sp] = ortho_long
        print(f'  Orthogroups long table: {ortho_long.shape}')

        pav = build_pav(ortho_long)
        cnv = ortho_long.groupby(['Orthogroup', 'Assembly Accession']).size().unstack(
            fill_value=0)

        # Save full unfiltered matrices to disk (preserves audit trail)
        os.makedirs(out_dir, exist_ok=True)
        pav.to_csv(pav_path, sep='\t')
        cnv.to_csv(cnv_path, sep='\t')
        print(f'  Saved PAV and CNV to {out_dir}')
        print(f'  PAV matrix: {pav.shape}')
        print(f'  CNV matrix: {cnv.shape}')

        # Apply runtime ANI exclusion to the in-memory copies returned to the caller
        pav = filter_ani_excluded(pav, sp, axis='cols')
        cnv = filter_ani_excluded(cnv, sp, axis='cols')
        if sp in ANI_EXCLUDED and ANI_EXCLUDED[sp]:
            print(f'  Returned (ANI-filtered) shapes: PAV {pav.shape}, CNV {cnv.shape}')
        pav_data[sp] = pav
        cnv_data[sp] = cnv

    return ortho_data, pav_data, cnv_data


def classify_all_pangenomes(species_list, ortho_data, pav_data, results_base=None,
                            force=False):
    """Classify orthogroups into Core/Accessory/Rare for all species.

    Cached to ``{results_base}/{sp}/{sp}_pangenome_class.tsv``.

    Returns
    -------
    pangenome_class : dict
        ``{species: classified_df}`` with ``pangenome_class`` column.
    """
    pangenome_class = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        cache_path = None
        if results_base is not None:
            cache_path = os.path.join(results_base, sp, f'{sp}_pangenome_class.tsv')

        # --- Load from cache ---
        if not force and cache_path and os.path.exists(cache_path):
            classified = pd.read_csv(cache_path, sep='\t')
            pangenome_class[sp] = classified
            counts = classified['Pangenome_Class'].value_counts()
            print(f'  Loaded cached classification ({len(classified):,} OGs)')
            print(f'  Class counts:\n{counts.to_string()}')
        else:
            # --- Compute ---
            pav = pav_data[sp]
            core_rate, core_n, rare_rate, rare_n = determine_core_and_rare_thresholds(pav)
            print(f'  Core threshold (>={core_n} genomes, {core_rate:.0%}), Rare threshold (<={rare_n} genomes, {rare_rate:.0%})')

            ortho_df = ortho_data.get(sp)
            if ortho_df is None:
                ortho_df = pd.DataFrame({'Orthogroup': pav.index})
            classified = classify_pangenome(ortho_df, pav, core_n, rare_n)
            pangenome_class[sp] = classified

            counts = classified['Pangenome_Class'].value_counts()
            print(f'  Class counts:\n{counts.to_string()}')

            if cache_path:
                os.makedirs(os.path.dirname(cache_path), exist_ok=True)
                classified.to_csv(cache_path, sep='\t', index=False)
                print(f'  Saved classification to {cache_path}')

        # --- Always show plot (recompute thresholds from PAV if needed) ---
        pav = pav_data[sp]
        core_rate, core_n, rare_rate, rare_n = determine_core_and_rare_thresholds(pav)
        fig, ax = plt.subplots(figsize=(8, 4))
        plot_gene_freq_hist(pav, core_n, rare_n, ax=ax)
        ax.set_title(f'A. {sp} — Gene Frequency Distribution')
        fig.tight_layout()
        if cache_path:
            fig.savefig(os.path.join(results_base, sp, f'{sp}_gene_freq_hist.png'),
                        dpi=120, bbox_inches='tight')
        plt.show()

    return pangenome_class


def pangenome_class_summary(pangenome_class, species_list):
    """Build a summary DataFrame comparing pangenome classes across species.

    Returns
    -------
    summary_df : pd.DataFrame
    """
    summary_rows = []
    for sp in species_list:
        counts = pangenome_class[sp]['Pangenome_Class'].value_counts()
        total = counts.sum()
        row = {'Species': sp}
        for cls in ['Core', 'Accessory', 'Rare']:
            row[cls] = counts.get(cls, 0)
            row[f'{cls} (%)'] = round(100 * counts.get(cls, 0) / total, 1)
        summary_rows.append(row)
    return pd.DataFrame(summary_rows)


def load_all_annotations_and_consensus(species_list, species_root, ortho_data,
                                        pangenome_class, pav_data, results_base,
                                        force=False):
    """Load functional annotations and build OG consensus tables for all species.

    Cached to ``{results_base}/{sp}/{sp}_og_consensus.tsv``. The underlying
    ``{species_root}/{sp}/{sp}_ortho_annot_long.tsv`` is also cached by
    ``build_ortho_annot_table`` (expensive step).

    Returns
    -------
    og_consensus : dict
        ``{species: og_consensus_df}``.
    """
    og_consensus = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_path = os.path.join(results_base, sp, f'{sp}_og_consensus.tsv')

        # --- Load from cache ---
        if not force and os.path.exists(out_path):
            consensus = pd.read_csv(out_path, sep='\t')
            og_consensus[sp] = consensus
            print(f'  Loaded cached OG consensus table: {consensus.shape}')
            continue

        # --- Compute ---
        sp_dir = os.path.join(species_root, sp)
        cache_path = os.path.join(sp_dir, f'{sp}_ortho_annot_long.tsv')
        long_annot = build_ortho_annot_table(sp_dir, sp, cache_path)
        print(f'  Annotated long table: {long_annot.shape}')

        pav = pav_data[sp]
        _, core_n, _, rare_n = determine_core_and_rare_thresholds(pav)

        consensus = build_og_consensus_with_pangenome(long_annot, pav, core_n, rare_n)
        og_consensus[sp] = consensus
        print(f'  OG consensus table: {consensus.shape}')

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        consensus.to_csv(out_path, sep='\t', index=False)
        print(f'  Saved -> {out_path}')

    return og_consensus


def run_all_enrichment(species_list, og_consensus, results_base, force=False):
    """Run functional enrichment analysis for all species.

    Cached per-species to ``{results_base}/{sp}/{sp}_enrichment_*.tsv``.
    Pass ``force=True`` to rerun.

    Returns
    -------
    enrichment_results : dict
        ``{species: dict_of_enrichment_dfs}``.
    """
    enrichment_results = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_dir = os.path.join(results_base, sp)
        os.makedirs(out_dir, exist_ok=True)

        # --- Try to load from cache ---
        cached_files = glob.glob(os.path.join(out_dir, f'{sp}_enrichment_*.tsv'))
        if not force and cached_files:
            enrich_dict = {}
            for fp in cached_files:
                layer = os.path.basename(fp).replace(f'{sp}_enrichment_', '').replace('.tsv', '')
                layer = layer.replace('_', '/')  # reverse the safe_name transform
                enrich_dict[layer] = pd.read_csv(fp, sep='\t')
            enrichment_results[sp] = enrich_dict
            print(f'  Loaded {len(enrich_dict)} cached enrichment layers from {out_dir}')
            if enrich_dict:
                combined = pd.concat(enrich_dict.values(), ignore_index=True)
                plot_enrichment_heatmap_single_species(combined, title=sp)
                plt.tight_layout()
                plt.show()
            continue

        # --- Compute from scratch ---
        consensus = og_consensus[sp]
        enrich_dict = run_pangenome_enrichment_og(consensus)
        enrichment_results[sp] = enrich_dict
        print(f'  Enrichment: {len(enrich_dict)} annotation layers tested')

        if enrich_dict:
            combined = pd.concat(enrich_dict.values(), ignore_index=True)
            plot_enrichment_heatmap_single_species(combined, title=sp)
            plt.tight_layout()
            plt.show()

        summary = summarize_enrichment_single_species(enrich_dict, sp)
        if not summary.empty:
            print(summary.to_string(index=False))

        for col_name, enr_df in enrich_dict.items():
            safe_name = col_name.replace('/', '_')
            enr_df.to_csv(os.path.join(out_dir, f'{sp}_enrichment_{safe_name}.tsv'),
                          sep='\t', index=False)
        print(f'  Saved enrichment results to {out_dir}')

    return enrichment_results


def run_all_snp_pca_and_grm(species_list, species_root, results_base, force=False):
    """Compute SNP PCA and GRM/kinship for all species.

    Cached to ``{results_base}/{sp}/{sp}_snp_pcs.tsv`` and ``{sp}_kinship.tsv``.

    Returns
    -------
    snp_results : dict
        ``{species: {'snp_pcs': df, 'grm': df}}``.
    """
    snp_results = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_dir = os.path.join(results_base, sp)
        os.makedirs(out_dir, exist_ok=True)
        pcs_path = os.path.join(out_dir, f'{sp}_snp_pcs.tsv')
        grm_path = os.path.join(out_dir, f'{sp}_kinship.tsv')

        # --- Load from cache ---
        if not force and os.path.exists(pcs_path) and os.path.exists(grm_path):
            pcs_df = pd.read_csv(pcs_path, sep='\t', index_col=0)
            grm_df = pd.read_csv(grm_path, sep='\t', index_col=0)
            # Runtime ANI exclusion: PCs by row (sample_id), GRM by both axes
            pcs_df = filter_ani_excluded(pcs_df, sp, axis='auto')
            grm_df = filter_ani_excluded(grm_df, sp, axis='both')
            snp_results[sp] = {'snp_pcs': pcs_df, 'grm': grm_df}
            print(f'  Loaded cached SNP PCs ({pcs_df.shape}) and GRM ({grm_df.shape})')
        else:
            # --- Compute ---
            sp_dir = os.path.join(species_root, sp)
            gubbins_dir = os.path.join(sp_dir, 'gubbins_output')

            pcs_df, snp_matrix, _ = run_snp_pca_analysis(
                gubbins_dir, species=sp,
                exclude_patterns=ANI_EXCLUDED.get(sp, set()),
            )
            if not pcs_df.empty:
                pcs_df.to_csv(pcs_path, sep='\t')
                print(f'  SNP PCs saved: {pcs_df.shape}')
            else:
                print(f'  No SNP data available for {sp}')

            if not snp_matrix.empty:
                grm_df, _ = run_grm_analysis(snp_matrix, species=sp, output_path=grm_path)
            else:
                grm_df = pd.DataFrame()

            snp_results[sp] = {'snp_pcs': pcs_df, 'grm': grm_df}

        # --- Always plot PCA + GRM heatmap (uses original plot_grm style) ---
        pcs_df = snp_results[sp]['snp_pcs']
        grm_df = snp_results[sp]['grm']

        if not pcs_df.empty and 'PC1' in pcs_df.columns and 'PC2' in pcs_df.columns:
            fig_pca, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(pcs_df['PC1'], pcs_df['PC2'], alpha=0.7, s=30, color='steelblue')
            ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
            ax.set_title(f'A. {sp} — SNP PCA (n={len(pcs_df)})')
            ax.grid(True, linestyle=':', alpha=0.5)
            fig_pca.tight_layout()
            fig_pca.savefig(os.path.join(out_dir, f'{sp}_snp_pca.png'),
                            dpi=120, bbox_inches='tight')
            plt.show()

        if not grm_df.empty:
            plot_grm(grm_df, species=sp)
            # plot_grm creates its own figure with clustering + RdBu_r colormap
            plt.savefig(os.path.join(out_dir, f'{sp}_kinship_grm.png'),
                        dpi=120, bbox_inches='tight')
            plt.show()

    return snp_results


def run_all_bgc_gcf(species_list, species_root, results_base, force=False):
    """Load BiG-SCAPE clustering and build BGC/GCF matrices for all species.

    Cached to ``{results_base}/{sp}/{sp}_{bgc_pav,gcf_pav,gcf_cnv}.tsv``.

    Returns
    -------
    bgc_results : dict
        ``{species: {'bgc_pav': df, 'gcf_pav': df, 'gcf_cnv': df}}``.
    """
    bgc_results = {}

    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_dir = os.path.join(results_base, sp)
        os.makedirs(out_dir, exist_ok=True)
        bgc_path = os.path.join(out_dir, f'{sp}_bgc_pav.tsv')
        gcf_pav_path = os.path.join(out_dir, f'{sp}_gcf_pav.tsv')
        gcf_cnv_path = os.path.join(out_dir, f'{sp}_gcf_cnv.tsv')

        # --- Load from cache ---
        if not force and all(os.path.exists(p) for p in [bgc_path, gcf_pav_path, gcf_cnv_path]):
            bgc_pav = pd.read_csv(bgc_path, sep='\t', index_col=0)
            gcf_pav = pd.read_csv(gcf_pav_path, sep='\t', index_col=0)
            gcf_cnv = pd.read_csv(gcf_cnv_path, sep='\t', index_col=0)
            # Runtime ANI exclusion — these matrices are genomes-as-rows,
            # but BGC IDs (in columns) include the accession too, so filter both.
            bgc_pav = filter_ani_excluded(bgc_pav, sp, axis='auto')
            gcf_pav = filter_ani_excluded(gcf_pav, sp, axis='auto')
            gcf_cnv = filter_ani_excluded(gcf_cnv, sp, axis='auto')
            bgc_results[sp] = {'bgc_pav': bgc_pav, 'gcf_pav': gcf_pav, 'gcf_cnv': gcf_cnv}
            print(f'  Loaded cached BGC matrices: BGC_PAV {bgc_pav.shape}, GCF_PAV {gcf_pav.shape}, GCF_CNV {gcf_cnv.shape}')
        else:
            # --- Compute ---
            sp_dir = os.path.join(species_root, sp)
            bigscape_dir = os.path.join(sp_dir, 'bigscape_output')
            bigscape = load_bigscape_clustering(bigscape_dir)
            print(f'  BiG-SCAPE clustering loaded: {len(bigscape)} entries')

            bgc_pav = build_bgc_pav_matrix(bigscape)
            bgc_pav.to_csv(bgc_path, sep='\t')
            print(f'  BGC PAV: {bgc_pav.shape}')

            gcf_pav = build_gcf_pav_matrix(bigscape)
            gcf_pav.to_csv(gcf_pav_path, sep='\t')
            print(f'  GCF PAV: {gcf_pav.shape}')

            gcf_cnv = build_gcf_cnv_matrix(bigscape)
            gcf_cnv.to_csv(gcf_cnv_path, sep='\t')
            print(f'  GCF CNV: {gcf_cnv.shape}')

            bgc_results[sp] = {'bgc_pav': bgc_pav, 'gcf_pav': gcf_pav, 'gcf_cnv': gcf_cnv}

        # --- Always plot BGC / GCF heatmaps ---
        bgc_pav = bgc_results[sp]['bgc_pav']
        gcf_pav = bgc_results[sp]['gcf_pav']
        gcf_cnv = bgc_results[sp]['gcf_cnv']
        if not bgc_pav.empty or not gcf_pav.empty:
            fig = plot_bgc_gcf_matrices(bgc_pav, gcf_pav, gcf_cnv, species=sp)
            if fig is not None:
                fig.savefig(os.path.join(out_dir, f'{sp}_bgc_gcf.png'),
                            dpi=120, bbox_inches='tight')
            plt.show()

    return bgc_results
