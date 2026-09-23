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
from collections import defaultdict

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


def build_ortho_annot_table(base_dir: str, species_name: str, out_path: str) -> pd.DataFrame:
    """Build the *long* per-protein OrthoFinder + annotation table."""
    if os.path.isfile(out_path):
        print(f"[i] Found existing {out_path}, loading long per-protein table.")
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

    # Symmetric color scale -- use the data's actual range so extremes don't get
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


_ANTISMASH_NOISE = {"hypothetical protein", "hypothetical", "protein", "domain", "like protein", "unknown"}
_TRNA_RE = re.compile(r"\btrna-[a-z]+\b", re.I)


# =============================================================================
# ADDITIONAL NOTEBOOK-SPECIFIC FUNCTIONS
# =============================================================================

    
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

    # Bootstrap confidence intervals (2.5th - 97.5th percentile)
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
                            force=False, plot=False):
    """Classify orthogroups into Core/Accessory/Rare for all species.

    Cached to ``{results_base}/{sp}/{sp}_pangenome_class.tsv``.

    ``plot=True`` draws one gene-frequency histogram per species with the core
    and rare thresholds marked, and writes ``{sp}_gene_freq_hist.png``.

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

        # --- Optional gene-frequency histogram (recompute thresholds from PAV) ---
        if plot:
            pav = pav_data[sp]
            core_rate, core_n, rare_rate, rare_n = determine_core_and_rare_thresholds(pav)
            fig, ax = plt.subplots(figsize=(8, 4))
            plot_gene_freq_hist(pav, core_n, rare_n, ax=ax)
            ax.set_title(f'A. {sp} - Gene Frequency Distribution')
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


def run_all_enrichment(species_list, og_consensus, results_base, force=False,
                       plot=False):
    """Run functional enrichment analysis for all species.

    Cached per-species to ``{results_base}/{sp}/{sp}_enrichment_*.tsv``.
    Pass ``force=True`` to rerun.

    ``plot=True`` renders the per-species enrichment heatmap inline.

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
            if plot and enrich_dict:
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

        if plot and enrich_dict:
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


def run_all_snp_pca_and_grm(species_list, species_root, results_base, force=False,
                            plot=False):
    """Compute SNP PCA and GRM/kinship for all species.

    Cached to ``{results_base}/{sp}/{sp}_snp_pcs.tsv`` and ``{sp}_kinship.tsv``.
    The kinship matrix is the random-effect covariance in the NB2 pan-GWAS.

    ``plot=True`` writes the per-species ``{sp}_snp_pca.png`` and
    ``{sp}_kinship_grm.png``.

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

        # --- Optional per-species PCA / GRM figures ---
        if not plot:
            continue

        pcs_df = snp_results[sp]['snp_pcs']
        grm_df = snp_results[sp]['grm']

        if not pcs_df.empty and 'PC1' in pcs_df.columns and 'PC2' in pcs_df.columns:
            fig_pca, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(pcs_df['PC1'], pcs_df['PC2'], alpha=0.7, s=30, color='steelblue')
            ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
            ax.set_title(f'A. {sp}, SNP PCA (n={len(pcs_df)})')
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


def run_all_bgc_gcf(species_list, species_root, results_base, force=False,
                    plot=False):
    """Load BiG-SCAPE clustering and build BGC/GCF matrices for all species.

    Cached to ``{results_base}/{sp}/{sp}_{bgc_pav,gcf_pav,gcf_cnv}.tsv``. The
    three matrices are feature layers in the NB2 pan-GWAS.

    ``plot=True`` writes ``{sp}_bgc_gcf.png``.

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
            # Runtime ANI exclusion -- these matrices are genomes-as-rows,
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
            if plot:
                fig = plot_bgc_gcf_matrices(bgc_pav, gcf_pav, gcf_cnv, species=sp)
                if fig is not None:
                    fig.savefig(os.path.join(out_dir, f'{sp}_bgc_gcf.png'),
                                dpi=120, bbox_inches='tight')
                plt.show()

    return bgc_results


# =============================================================================
# GCF ANNOTATION (CLASS + MIBiG)
# =============================================================================

_KCB_REGION_RE = re.compile(r'(scaffold_\d+)_c(\d+)$')


def parse_knownclusterblast(filepath):
    """Parse one antiSMASH knownclusterblast .txt into scored MIBiG hits.

    Reads the ``Details`` blocks rather than the ``Significant hits`` list, so
    each hit carries the evidence antiSMASH recorded for it.

    Returns
    -------
    list of dict
        ``{'mibig_id', 'compound', 'mibig_type', 'n_genes', 'score'}``.
    """
    try:
        content = Path(filepath).read_text()
    except (OSError, UnicodeDecodeError):
        return []

    hits = []
    for block in content.split('>>')[1:]:
        m_id = re.search(r'^\s*\d+\.\s*(BGC\d+)', block, re.MULTILINE)
        if not m_id:
            continue
        m_src = re.search(r'^Source:\s*(.*)$', block, re.MULTILINE)
        m_type = re.search(r'^Type:\s*(.*)$', block, re.MULTILINE)
        m_n = re.search(r'^Number of proteins with BLAST hits to this cluster:\s*(\d+)',
                        block, re.MULTILINE)
        m_sc = re.search(r'^Cumulative BLAST score:\s*(\d+)', block, re.MULTILINE)
        hits.append({
            'mibig_id': m_id.group(1),
            'compound': m_src.group(1).strip() if m_src else '',
            'mibig_type': m_type.group(1).strip() if m_type else '',
            'n_genes': int(m_n.group(1)) if m_n else 0,
            'score': int(m_sc.group(1)) if m_sc else 0,
        })
    return hits


def build_gcf_annotation(species, species_root, max_hits=3, min_support=0.5):
    """Annotate every GCF of one species with its BGC class and MIBiG evidence.

    Two independent kinds of MIBiG evidence are reported, because they answer
    different questions and disagree often:

    ``mibig_in_family``
        MIBiG reference BGCs that BiG-SCAPE placed *inside* the family at the
        c0.3 cutoff. This is the strict reading of "has a characterized
        reference compound".
    ``kcb_*``
        antiSMASH knownclusterblast homology hits, ranked by how many member
        regions recover the hit and, within that, by cumulative BLAST score.
        ``kcb_support`` is the fraction of member regions carrying the top hit,
        so a hit seen in 1 of 5 regions is not mistaken for a family-wide one.

    ANI-excluded genomes are dropped so the counts match the NB1 GCF matrices.

    Returns
    -------
    pandas.DataFrame
        One row per GCF, or an empty frame when BiG-SCAPE output is missing.
    """
    sp_dir = Path(species_root) / species
    bigscape_dirs = sorted((sp_dir / 'bigscape_output' / 'output_files').glob('*_c0*'))
    if not bigscape_dirs:
        print(f'  No BiG-SCAPE output for {species}')
        return pd.DataFrame()
    run_dir = bigscape_dirs[-1]

    ann_path = run_dir / 'record_annotations.tsv'
    if not ann_path.exists():
        print(f'  No record_annotations.tsv in {run_dir.name}')
        return pd.DataFrame()
    ann = pd.read_csv(ann_path, sep='\t')
    rec_class = dict(zip(ann['Record'], ann['Class']))
    rec_cat = dict(zip(ann['Record'], ann['Category']))

    clus_files = sorted(run_dir.glob('*/*_clustering_*.tsv'))
    if not clus_files:
        print(f'  No clustering tables in {run_dir.name}')
        return pd.DataFrame()
    clus = pd.concat(
        [pd.read_csv(f, sep='\t').assign(bin=f.parent.name) for f in clus_files],
        ignore_index=True)

    is_ref = clus['Record'].astype(str).str.startswith('BGC')
    refs = clus[is_ref]
    members = clus[~is_ref].copy()
    members = members[~members['Record'].astype(str).map(
        lambda r: is_ani_excluded(r, species))]
    if members.empty:
        print(f'  No strain records left for {species} after ANI exclusion')
        return pd.DataFrame()

    members['strain'] = members['GBK'].astype(str).str.split('_scaffold').str[0]
    members['Class'] = members['Record'].map(rec_class)
    members['Category'] = members['Record'].map(rec_cat)

    # MIBiG references that clustered into each family
    fam_refs = {}
    for fam, g in refs.groupby('Family'):
        ids = sorted({str(r).split('.gbk')[0] for r in g['Record']})
        fam_refs[fam] = ids

    # knownclusterblast hits, keyed by BiG-SCAPE Record
    kcb_by_record = defaultdict(list)
    antismash_dir = sp_dir / 'antismash_output'
    record_index = defaultdict(list)
    for rec in members['Record']:
        record_index[str(rec).split('_scaffold')[0]].append(str(rec))
    if antismash_dir.is_dir():
        for strain_dir in antismash_dir.iterdir():
            kcb_dir = strain_dir / 'knownclusterblast'
            if not kcb_dir.is_dir():
                continue
            strain_records = record_index.get(strain_dir.name, [])
            if not strain_records:
                continue
            for txt in kcb_dir.glob('*.txt'):
                m = _KCB_REGION_RE.match(txt.stem)
                if not m:
                    continue
                scaffold, region_num = m.group(1), int(m.group(2))
                # antiSMASH zero-pads region numbers inconsistently across versions
                patterns = [f'{strain_dir.name}_{scaffold}.region{region_num:03d}',
                            f'{strain_dir.name}_{scaffold}.region{region_num:02d}',
                            f'{strain_dir.name}_{scaffold}.region{region_num}']
                matched = [r for r in strain_records if any(p in r for p in patterns)]
                if not matched:
                    continue
                hits = parse_knownclusterblast(txt)
                for rec in matched:
                    kcb_by_record[rec].extend(hits)

    rows = []
    for fam, g in members.groupby('Family'):
        n_regions = len(g)
        classes = g['Class'].dropna()
        cats = g['Category'].dropna()

        # Rank knownclusterblast hits by how many member regions recover them
        per_hit_regions = Counter()
        per_hit_best = {}
        for rec in g['Record']:
            seen = set()
            for h in kcb_by_record.get(str(rec), []):
                key = h['mibig_id']
                if key not in seen:
                    per_hit_regions[key] += 1
                    seen.add(key)
                prev = per_hit_best.get(key)
                if prev is None or h['score'] > prev['score']:
                    per_hit_best[key] = h
        ranked = sorted(per_hit_regions,
                        key=lambda k: (per_hit_regions[k], per_hit_best[k]['score']),
                        reverse=True)
        top = ranked[:max_hits]
        best = per_hit_best[ranked[0]] if ranked else None
        support = per_hit_regions[ranked[0]] / n_regions if ranked else 0.0

        rows.append({
            'GCF': fam,
            'n_regions': n_regions,
            'n_strains': g['strain'].nunique(),
            'bigscape_bin': '; '.join(sorted(set(g['bin']))),
            'antismash_class': classes.value_counts().index[0] if len(classes) else '',
            'class_purity': (classes.value_counts().iloc[0] / len(classes)) if len(classes) else float('nan'),
            'all_classes': '; '.join(sorted(set(classes))),
            'bigscape_category': cats.value_counts().index[0] if len(cats) else '',
            'mibig_in_family': '; '.join(fam_refs.get(fam, [])),
            'n_mibig_in_family': len(fam_refs.get(fam, [])),
            'kcb_top_id': best['mibig_id'] if best else '',
            'kcb_top_compound': best['compound'] if best else '',
            'kcb_top_genes': best['n_genes'] if best else 0,
            'kcb_top_score': best['score'] if best else 0,
            'kcb_support': round(support, 3),
            'kcb_confident': bool(best is not None and support >= min_support),
            'kcb_top_hits': '; '.join(
                f'{k}({per_hit_regions[k]}/{n_regions})' for k in top),
        })

    out = pd.DataFrame(rows).sort_values('GCF').reset_index(drop=True)
    print(f'  {species}: {len(out)} GCFs from {run_dir.name}, '
          f'{int((out["n_mibig_in_family"] > 0).sum())} with a MIBiG reference in-family, '
          f'{int(out["kcb_confident"].sum())} with a confident knownclusterblast hit')
    return out


def run_all_gcf_annotation(species_list, species_root, results_base, force=False,
                           max_hits=3, min_support=0.5):
    """Build and cache the GCF annotation table for every species.

    Cached to ``{results_base}/{sp}/{sp}_gcf_annotation.tsv``. The ``GCF``
    column joins onto the GCF feature columns of ``{sp}_gcf_pav.tsv`` and onto
    the ``feature`` column of the NB2 GCF association tables.

    Returns
    -------
    dict
        ``{species: DataFrame}``.
    """
    annotations = {}
    for sp in species_list:
        print(f'\n=== {sp} ===')
        out_dir = os.path.join(results_base, sp)
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, f'{sp}_gcf_annotation.tsv')

        if not force and os.path.exists(path):
            df = pd.read_csv(path, sep='\t')
            print(f'  Loaded cached GCF annotation: {df.shape}')
        else:
            df = build_gcf_annotation(sp, species_root, max_hits=max_hits,
                                      min_support=min_support)
            if df.empty:
                annotations[sp] = df
                continue
            df.to_csv(path, sep='\t', index=False)
            print(f'  Wrote {path}')
        annotations[sp] = df

    return annotations


def annotate_gcf_associations(assoc_df, annotation_df,
                              cols=('antismash_class', 'bigscape_category',
                                    'n_regions', 'n_strains', 'mibig_in_family',
                                    'kcb_top_id', 'kcb_top_compound',
                                    'kcb_support', 'kcb_confident')):
    """Left-join GCF annotations onto an association table on ``feature``.

    Returns the association table unchanged when either input is empty or the
    ``feature`` column is missing, so it is safe to map over every contrast.
    """
    if assoc_df is None or annotation_df is None:
        return assoc_df
    if assoc_df.empty or annotation_df.empty or 'feature' not in assoc_df.columns:
        return assoc_df
    keep = ['GCF'] + [c for c in cols if c in annotation_df.columns]
    return assoc_df.merge(annotation_df[keep], how='left',
                          left_on='feature', right_on='GCF').drop(columns='GCF')


def summarise_gcf_annotation(gcf_annotation):
    """Per-species MIBiG evidence counts and BiG-SCAPE category composition.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        Evidence counts per species, and GCF counts per category.
    """
    evidence = pd.DataFrame({
        sp: {'GCFs': len(a),
             'regions': int(a['n_regions'].sum()),
             'MIBiG ref in family': int((a['n_mibig_in_family'] > 0).sum()),
             'confident KCB hit': int(a['kcb_confident'].sum()),
             'no MIBiG evidence': int((~a['kcb_confident'] & (a['n_mibig_in_family'] == 0)).sum())}
        for sp, a in gcf_annotation.items() if a is not None and not a.empty}).T
    composition = pd.DataFrame(
        {sp: a['bigscape_category'].value_counts()
         for sp, a in gcf_annotation.items() if a is not None and not a.empty}
    ).fillna(0).astype(int)
    return evidence, composition


def collect_significant_gcf_associations(species_list, nb2_results):
    """Gather FDR-significant rows from every annotated GCF association table.

    Reads ``{nb2_results}/{sp}/pangwas_results/gcf_*_assoc_*_annotated.tsv``.

    Returns
    -------
    pandas.DataFrame
        One row per significant species x contrast x layer x GCF.
    """
    rows = []
    for sp in species_list:
        sp_dir = Path(nb2_results) / sp / 'pangwas_results'
        if not sp_dir.is_dir():
            continue
        for f in sorted(sp_dir.glob('gcf_*_assoc_*_annotated.tsv')):
            layer, contrast = f.stem.replace('_annotated', '').split('_assoc_')
            df = pd.read_csv(f, sep='\t')
            if 'significant' not in df.columns:
                continue
            for _, r in df[df['significant']].iterrows():
                rows.append({
                    'species': sp, 'contrast': contrast, 'layer': layer,
                    'GCF': r['feature'],
                    'class': r.get('antismash_class'),
                    'direction': 'gained' if r['beta'] > 0 else 'lost',
                    'beta': round(r['beta'], 3),
                    'qvalue': r['qvalue'],
                    'n_strains': r.get('n_strains'),
                    'mibig_in_family': r.get('mibig_in_family'),
                    'kcb_compound': r.get('kcb_top_compound'),
                    'kcb_support': r.get('kcb_support'),
                })
    if not rows:
        return pd.DataFrame()
    return (pd.DataFrame(rows)
            .sort_values(['species', 'contrast', 'layer', 'qvalue'])
            .reset_index(drop=True))


# =============================================================================
# COMBINED KINSHIP / GRM FIGURE
# =============================================================================

KINSHIP_PHENO_COLORS = {
    'Human-pathogenic':  '#e4211c',
    'Environmental':     '#2ca02c',
    'Industrial-trait':  '#1f9ede',
    'Plant-pathogenic':  '#ff7f0e',
    'Animal-pathogenic': '#8b5cb8',
    'Unknown':           '#bdbdbd',
}


def plot_kinship_grid(species_list, results_base, pheno_csv, out_path=None,
                      pheno_colors=None, figsize=(15, 14), dpi=150,
                      fontsize=13, clip_pct=(2, 98), method='average',
                      cmap_name='RdBu_r', show_left_strip=True,
                      drop_unlabelled=True):
    """Clustered kinship (GRM) heatmap per species, annotated by phenotype.

    Reads ``{results_base}/{sp}/{sp}_kinship.tsv`` for each species, orders the
    strains by hierarchical clustering of the relatedness matrix, and draws the
    phenotype of each strain as a colour strip along the matrix edges. The file
    is the one the NB2 pan-GWAS reads as random-effect covariance.

    With ``drop_unlabelled=True``, ANI-excluded genomes and strains without a
    phenotype label are removed, leaving the strains used in the pan-GWAS.
    With ``drop_unlabelled=False`` every QC-passed genome is shown.

    Parameters
    ----------
    species_list : list of str
    results_base : str
        NB1 results directory holding ``{sp}/{sp}_kinship.tsv``.
    pheno_csv : str
        CSV with ``Assembly Accession`` and ``Phenotype`` columns (NB0 output).
    out_path : str, optional
        Path of the PNG to write.
    clip_pct : (low, high)
        Percentiles of the off-diagonal values used as colour limits. The
        diagonal (self-relatedness) is masked and shown in grey.
    method : str
        scipy linkage method.

    Returns
    -------
    fig, info : matplotlib Figure and ``{species: {'order', 'labels', 'n'}}``
    """
    from matplotlib.colors import TwoSlopeNorm, to_rgb
    from matplotlib.patches import Patch
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform

    if pheno_colors is None:
        pheno_colors = KINSHIP_PHENO_COLORS

    pheno = pd.read_csv(pheno_csv)
    pheno_map = dict(zip(pheno['Assembly Accession'].astype(str), pheno['Phenotype']))

    fig = plt.figure(figsize=figsize)
    outer = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.24)
    info = {}
    seen_labels = set()

    for idx, sp in enumerate(species_list):
        path = os.path.join(results_base, sp, f'{sp}_kinship.tsv')
        if not os.path.exists(path):
            print(f'  {sp}: no kinship matrix at {path} -- skipping')
            continue

        K = pd.read_csv(path, sep='\t', index_col=0).astype(float)

        # Mirror the pan-GWAS load path (funpan_gwas.run_all_gwas)
        K = filter_ani_excluded(K, sp, axis='both')

        accs = []
        for i in K.index:
            m = ACC_RE.search(str(i))
            accs.append(m.group(1) if m else str(i))

        if drop_unlabelled:
            keep = [j for j, a in enumerate(accs) if a in pheno_map]
            n_drop = len(accs) - len(keep)
            if n_drop:
                print(f'  {sp}: dropped {n_drop} strain(s) with no phenotype label')
            K = K.iloc[keep, keep]
            accs = [accs[j] for j in keep]

        labels = [pheno_map.get(a, 'Unknown') for a in accs]

        M = K.values.copy()
        M = (M + M.T) / 2.0                       # linkage needs exact symmetry

        dist = M.max() - M
        np.fill_diagonal(dist, 0.0)
        dist[dist < 0] = 0.0
        order = leaves_list(linkage(squareform(dist, checks=False), method=method))

        Mo = M[np.ix_(order, order)]
        lab_o = [labels[i] for i in order]
        acc_o = [accs[i] for i in order]
        seen_labels.update(lab_o)
        # Per-strain diagnostics: diagonal and mean off-diagonal relatedness
        diag = np.diag(Mo)
        offd = Mo.copy(); np.fill_diagonal(offd, np.nan)
        info[sp] = {
            'order': acc_o, 'labels': lab_o, 'n': len(Mo),
            'diag': pd.Series(diag, index=acc_o),
            'mean_offdiag': pd.Series(np.nanmean(offd, axis=1), index=acc_o),
            'n_high_diag': int((diag > 3 * np.median(diag)).sum()),
        }

        off = Mo[~np.eye(len(Mo), dtype=bool)]
        lo, hi = np.percentile(off, clip_pct)
        norm = TwoSlopeNorm(vmin=min(lo, -0.05), vcenter=0.0, vmax=max(hi, 0.05))

        Mplot = Mo.copy()
        np.fill_diagonal(Mplot, np.nan)

        n = len(lab_o)
        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad('#f0f0f0')

        # aspect='equal' shrinks the heatmap axes inside its gridspec cell, so the
        # strips and colourbar attach to the heatmap axes and share its coordinates
        ax = fig.add_subplot(outer[idx])
        im = ax.imshow(Mplot, cmap=cmap, norm=norm, aspect='equal',
                       interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel('Strains (hierarchically clustered)', fontsize=fontsize - 3)

        div = make_axes_locatable(ax)
        ax_top = div.append_axes('top', size='3.2%', pad=0.05, sharex=ax)
        ax_left = (div.append_axes('left', size='3.2%', pad=0.05, sharey=ax)
                   if show_left_strip else None)
        cax = div.append_axes('right', size='3.5%', pad=0.10)

        rgb = np.array([[to_rgb(pheno_colors.get(l, '#bdbdbd')) for l in lab_o]])
        ax_top.imshow(rgb, aspect='auto', interpolation='nearest',
                      extent=(-0.5, n - 0.5, 0.0, 1.0))
        ax_top.set_xticks([]); ax_top.set_yticks([])
        for side in ax_top.spines.values():
            side.set_visible(False)
        ax_top.set_title(f'A. {sp}  (n = {n})', fontsize=fontsize + 1,
                         fontweight='bold', pad=8)

        if ax_left is not None:
            ax_left.imshow(rgb.transpose(1, 0, 2), aspect='auto',
                           interpolation='nearest',
                           extent=(0.0, 1.0, n - 0.5, -0.5))
            ax_left.set_xticks([]); ax_left.set_yticks([])
            for side in ax_left.spines.values():
                side.set_visible(False)

        cb = fig.colorbar(im, cax=cax)
        cb.set_label('Relatedness', fontsize=fontsize - 3)
        cb.ax.tick_params(labelsize=fontsize - 4)

    handles = [Patch(facecolor=c, label=l) for l, c in pheno_colors.items()
               if l in seen_labels]
    fig.legend(handles=handles, loc='lower center', ncol=min(6, len(handles)),
               frameon=False, fontsize=fontsize - 1, bbox_to_anchor=(0.5, 0.04))
    fig.suptitle('SNP-derived genetic relatedness (kinship) per species',
                 fontsize=fontsize + 4, fontweight='bold', y=0.955)

    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')

    return fig, info


# =============================================================================
# NB1 FIGURES AND FIRST-RUN SUPPORT
# =============================================================================

def pangenome_composition_counts(pangenome_class, pav_data, species_list,
                                 verbose=True):
    """Core / accessory / rare orthogroup counts per species.

    ``pangenome_class[sp]`` is the long per-protein table, so it is deduplicated
    by Orthogroup. Orthogroups present in zero analysed genomes are dropped:
    load_all_orthofinder removes ANI-excluded genomes from the PAV columns but
    leaves their orthogroup rows in place.

    Returns
    -------
    pd.DataFrame
        Columns ``species``, ``n_genomes``, ``core``, ``accessory``, ``rare``,
        ``total``.
    """
    rows = []
    for sp in species_list:
        per_og = pangenome_class[sp].drop_duplicates('Orthogroup')
        present = pav_data[sp].sum(axis=1)
        keep = set(present[present > 0].index.astype(str))
        n_drop = (~per_og['Orthogroup'].astype(str).isin(keep)).sum()
        if n_drop and verbose:
            print(f'  {sp}: dropped {n_drop} orthogroup(s) absent from all '
                  f'{pav_data[sp].shape[1]} analysed genomes')
        per_og = per_og[per_og['Orthogroup'].astype(str).isin(keep)]
        counts = per_og['Pangenome_Class'].value_counts()
        core = int(counts.get('Core', 0))
        acc = int(counts.get('Accessory', 0))
        rare = int(counts.get('Rare', 0))
        rows.append({'species': sp, 'n_genomes': pav_data[sp].shape[1],
                     'core': core, 'accessory': acc, 'rare': rare,
                     'total': core + acc + rare})
    return pd.DataFrame(rows)


def plot_pangenome_composition(pangenome_class, pav_data, species_list,
                               out_path=None, fontsize=20, figsize=(14, 10),
                               dpi=800, colors=("#0B8B08C5", "#ef990f", "#ec5e4f"),
                               xlim_factor=1.20, label_pad=0.07):
    """Stacked bar of core / accessory / rare orthogroup counts per species.

    ``pangenome_class[sp]`` is the long per-protein table, so it is deduplicated
    by Orthogroup to give per-OG counts. Orthogroups present in zero analysed
    genomes are dropped first: load_all_orthofinder removes ANI-excluded genomes
    from the PAV columns but leaves their orthogroup rows in place, and those
    rows would otherwise be counted as Rare.

    Returns
    -------
    fig, counts : Figure and a DataFrame of the per-species counts plotted.
    """
    core_c, acc_c, rare_c = colors
    tab = pangenome_composition_counts(pangenome_class, pav_data, species_list)
    core_n = tab['core'].tolist()
    acc_n = tab['accessory'].tolist()
    rare_n = tab['rare'].tolist()
    n_genomes_list = tab['n_genomes'].tolist()

    totals = [c + a + r for c, a, r in zip(core_n, acc_n, rare_n)]
    core_pcts = [100 * c / t for c, t in zip(core_n, totals)]
    acc_pcts = [100 * a / t for a, t in zip(acc_n, totals)]
    rare_pcts = [100 * r / t for r, t in zip(rare_n, totals)]

    labels = [f"A. {sp}\n(n={n})\ntotal = {t:,} OGs"
              for sp, n, t in zip(species_list, n_genomes_list, totals)]
    y = np.arange(len(species_list))
    FS = fontsize

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(y, core_n, color=core_c, edgecolor='white', label='Core')
    ax.barh(y, acc_n, left=core_n, color=acc_c, edgecolor='white', label='Accessory')
    ax.barh(y, rare_n, left=[c + a for c, a in zip(core_n, acc_n)],
            color=rare_c, edgecolor='white', label='Rare')

    xmax = max(totals)
    for i, (c, a, r, t, cp, ap, rp) in enumerate(zip(core_n, acc_n, rare_n, totals,
                                                     core_pcts, acc_pcts, rare_pcts)):
        ax.text(c / 2, i, f'{c:,}\n({cp:.1f}%)', ha='center', va='center',
                fontsize=FS, color='white', fontweight='bold')
        ax.text(c + a / 2, i, f'{a:,}\n({ap:.1f}%)', ha='center', va='center',
                fontsize=FS, color='white', fontweight='bold')
        # Rare segment is too narrow for inside text, so annotate past the bar end
        ax.text(t + xmax * label_pad, i, f'{r:,}\n({rp:.1f}%)',
                ha='center', va='center', fontsize=FS, color=rare_c, fontweight='bold')

    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=FS)
    ax.tick_params(axis='x', labelsize=FS)
    ax.set_xlabel('Number of orthogroups', fontsize=FS)
    ax.set_xlim(0, xmax * xlim_factor); ax.invert_yaxis()
    ax.set_title('Pangenome composition across Aspergillus species',
                 fontweight='bold', fontsize=FS + 2)
    ax.legend(loc='upper right', fontsize=FS, frameon=True, framealpha=0.95)
    ax.grid(True, axis='x', linestyle=':', alpha=0.5)
    plt.tight_layout()

    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')

    counts = pd.DataFrame({'species': species_list, 'n_genomes': n_genomes_list,
                           'core': core_n, 'accessory': acc_n, 'rare': rare_n,
                           'total': totals})
    return fig, counts


def plot_heaps_law_curves(species_list, results_base, out_path=None,
                          fontsize=20, figsize=(12, 10), dpi=800,
                          species_colors=None):
    """Normalised pangenome accumulation and core decay curves per species.

    Reads the cached bootstrap results from
    ``{results_base}/{sp}/{sp}_heaps.pkl``, so analyze_heaps_law must have run.
    Curves are expressed as a percentage of each species' final pangenome size
    so that species of different sizes share one axis.
    """
    import pickle
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    if species_colors is None:
        species_colors = {sp: SPECIES_COLORS.get(sp, '#333333') for sp in species_list}

    FS = fontsize
    fig, ax = plt.subplots(figsize=figsize)
    species_handles = []

    for sp in species_list:
        pkl = os.path.join(results_base, sp, f'{sp}_heaps.pkl')
        if not os.path.exists(pkl):
            print(f'  {sp}: no cached Heap\'s law result at {pkl} -- skipping')
            continue
        with open(pkl, 'rb') as fh:
            h = pickle.load(fh)

        n = h['n_genomes']
        Ns_pct = 100 * np.arange(1, n + 1) / n
        c = species_colors[sp]
        gamma = h['gamma']

        pan = np.asarray(h['pangenome_size']); core = np.asarray(h['core_size'])
        pan_lo, pan_hi = np.asarray(h['pan_ci_lo']), np.asarray(h['pan_ci_hi'])
        core_lo, core_hi = np.asarray(h['core_ci_lo']), np.asarray(h['core_ci_hi'])
        pan_final = pan[-1]

        ax.fill_between(Ns_pct, 100 * pan_lo / pan_final, 100 * pan_hi / pan_final,
                        color=c, alpha=0.15, linewidth=0)
        ax.fill_between(Ns_pct, 100 * core_lo / pan_final, 100 * core_hi / pan_final,
                        color=c, alpha=0.15, linewidth=0)
        ax.plot(Ns_pct, 100 * pan / pan_final, '-', color=c, linewidth=2.5)
        ax.plot(Ns_pct, 100 * core / pan_final, '--', color=c, linewidth=2.5)

        species_handles.append(Line2D([0], [0], color=c, lw=2.5,
                                      label=f"A. {sp}  (n={n}, $\\gamma$={gamma:.3f})"))

    ax.set_xlabel('% of genomes sampled', fontsize=FS)
    ax.set_ylabel('% of final pangenome size', fontsize=FS)
    ax.set_title("Pangenome accumulation & core decay\n(Heap's law, normalized)",
                 fontweight='bold', fontsize=FS + 2)
    ax.set_xlim(0, 100); ax.set_ylim(0, 105)
    ax.tick_params(axis='both', labelsize=FS)
    ax.grid(True, linestyle=':', alpha=0.5)

    style_handles = [
        Line2D([0], [0], color='gray', lw=2.5, linestyle='-', label='Pangenome'),
        Line2D([0], [0], color='gray', lw=2.5, linestyle='--', label='Core genome'),
        Patch(facecolor='gray', alpha=0.25, edgecolor='none', label='95% bootstrap CI'),
    ]
    leg1 = ax.legend(handles=species_handles, loc='lower left',
                     bbox_to_anchor=(0.02, 0.02), fontsize=FS - 1,
                     title='Species', title_fontsize=FS,
                     frameon=True, framealpha=0.95)
    ax.add_artist(leg1)
    ax.legend(handles=style_handles, loc='lower right',
              bbox_to_anchor=(0.98, 0.02), fontsize=FS - 1,
              title='Curve / band', title_fontsize=FS,
              frameon=True, framealpha=0.95)
    plt.tight_layout()

    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')

    return fig


def nb1_preflight(species_list, species_root, results_base, nb0_results,
                  verbose=True):
    """Check that NB1's inputs exist before anything is computed.

    Reports, per species, whether the OrthoFinder, Gubbins and BiG-SCAPE outputs
    the notebook needs are present, and whether NB0's phenotype table is there.
    Creates the per-species output directories. Nothing is computed or written
    beyond those directories.

    Returns
    -------
    status : pd.DataFrame
        One row per species with a boolean per required input.
    """
    os.makedirs(results_base, exist_ok=True)

    needed = {
        'orthofinder': 'orthofinder_output',
        'gubbins':     'gubbins_output',
        'bigscape':    'bigscape_output',
        'proteins':    'filtered_protein',
        'genomes':     'filtered_genome',
    }

    rows = []
    for sp in species_list:
        os.makedirs(os.path.join(results_base, sp), exist_ok=True)
        sp_dir = os.path.join(species_root, sp)
        row = {'species': sp, 'species_dir': os.path.isdir(sp_dir)}
        for key, sub in needed.items():
            row[key] = os.path.isdir(os.path.join(sp_dir, sub))
        rows.append(row)

    status = pd.DataFrame(rows).set_index('species')
    pheno = os.path.join(nb0_results, 'phenotype_classified_for_gwas.csv')
    pheno_ok = os.path.exists(pheno)

    if verbose:
        print('NB1 preflight')
        print('=' * 60)
        print(f'  species root : {species_root}')
        print(f'  results base : {results_base}')
        print(f'  NB0 phenotype table: {"found" if pheno_ok else "MISSING -- run NB0 first"}')
        print(f'  -> {pheno}')
        print()
        print(status.replace({True: 'ok', False: 'MISSING'}).to_string())
        missing = status.replace({True: None}).notna().sum().sum()
        print()
        if missing or not pheno_ok:
            print(f'{missing} missing input(s). Sections depending on them will fail;')
            print('see the per-section notes for which inputs each one needs.')
        else:
            print('All inputs present.')
    return status


def print_nb1_summary(pangenome_class, pav_data, species_list, results_base):
    """Print the final NB1 pangenome counts and the files written, with paths.

    Returns
    -------
    counts : pd.DataFrame
        Core / accessory / rare / total orthogroups and Heap's law gamma per species.
    files : pd.DataFrame
        Every expected output with its full path and whether it exists.
    """
    import glob
    import pickle

    counts = pangenome_composition_counts(pangenome_class, pav_data, species_list,
                                          verbose=False)
    gammas = []
    for sp in species_list:
        pkl = os.path.join(results_base, sp, f'{sp}_heaps.pkl')
        gamma = None
        if os.path.exists(pkl):
            with open(pkl, 'rb') as fh:
                gamma = pickle.load(fh).get('gamma')
        gammas.append(round(gamma, 4) if gamma is not None else None)
    counts['heaps_gamma'] = gammas
    counts = counts.set_index('species')

    top_level = ['pangenome_composition_combined.png',
                 'heaps_law_combined.png',
                 'kinship_grm_combined.png']
    per_species = ['{sp}_pav.tsv', '{sp}_cnv.tsv', '{sp}_pangenome_class.tsv',
                   '{sp}_gene_freq_hist.png', '{sp}_heaps.pkl',
                   '{sp}_og_consensus.tsv', '{sp}_enrichment_*.tsv',
                   '{sp}_snp_pcs.tsv', '{sp}_kinship.tsv',
                   '{sp}_bgc_pav.tsv', '{sp}_gcf_pav.tsv', '{sp}_gcf_cnv.tsv']

    rows = []
    for name in top_level:
        path = os.path.join(results_base, name)
        rows.append({'species': 'all', 'file': name, 'path': path,
                     'exists': os.path.exists(path)})
    for sp in species_list:
        for pattern in per_species:
            name = pattern.format(sp=sp)
            path = os.path.join(results_base, sp, name)
            if '*' in name:
                n = len(glob.glob(path))
                rows.append({'species': sp, 'file': f'{name} ({n} files)',
                             'path': path, 'exists': n > 0})
            else:
                rows.append({'species': sp, 'file': name, 'path': path,
                             'exists': os.path.exists(path)})
    files = pd.DataFrame(rows)

    print(f'\n{"=" * 70}')
    print('NB1 Pangenome Architecture complete.')
    print(f'All outputs saved to: {results_base}')
    for _, r in files.iterrows():
        mark = ' ' if r['exists'] else '!'
        label = r['path'] if '*' not in r['file'] else f"{r['path']}  ({r['file'].split('(')[1]}"
        print(f"  {mark} {label}")
    n_missing = int((~files['exists']).sum())
    if n_missing:
        print(f'\n  ! = missing ({n_missing} file(s))')
    print('=' * 70)
    print('\nOrthogroup counts per species:')

    return counts, files
