#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_convergence.py
=====================
Cross-species convergence analysis, meta-analysis, annotation transfer,
and combined pangenome matrix building.

Extracted from Theme2_Analysis_Functions.py.
"""

from funpan_utils import *
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter
from scipy import stats
from scipy.stats import fisher_exact, chi2, norm
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# =============================================================================

def load_combined_orthogroups_genecount(genecount_path: str) -> pd.DataFrame:
    """
    Load the OrthoFinder Orthogroups.GeneCount.tsv for the combined pangenome.

    The column headers have format: species__GCA_xxx.proteins
    This function parses out species and accession from column names.

    Parameters
    ----------
    genecount_path : str
        Path to Orthogroups.GeneCount.tsv

    Returns
    -------
    pd.DataFrame
        Gene count matrix with columns renamed to 'species__accession' format
    """
    df = pd.read_csv(genecount_path, sep='\t')

    # Rename columns to extract species and accession
    rename_map = {'Orthogroup': 'Orthogroup'}
    for col in df.columns:
        if col in ['Orthogroup', 'Total']:
            continue

        # Parse species__GCA_xxx.xxx_assembly.proteins format
        parts = col.split('__')
        if len(parts) >= 2:
            species = parts[0]
            acc_match = ACC_RE.search(parts[1])
            if acc_match:
                accession = acc_match.group(1)
                rename_map[col] = f"{species}__{accession}"
            else:
                rename_map[col] = col
        else:
            rename_map[col] = col

    df = df.rename(columns=rename_map)

    # Drop Total column if present
    if 'Total' in df.columns:
        df = df.drop(columns=['Total'])

    return df


def load_combined_orthogroups_long(orthogroups_tsv: str) -> pd.DataFrame:
    """
    Load and reshape OrthoFinder Orthogroups.tsv to long format for combined pangenome.

    Parameters
    ----------
    orthogroups_tsv : str
        Path to Orthogroups.tsv

    Returns
    -------
    pd.DataFrame
        Long format with columns: Orthogroup, Species, Assembly_Accession, Protein_ID
    """
    og = pd.read_csv(orthogroups_tsv, sep='\t')

    # Build mapping from original column names to (species, accession)
    col_info = {}
    for col in og.columns:
        if col == 'Orthogroup':
            continue

        parts = col.split('__')
        if len(parts) >= 2:
            species = parts[0]
            # Normalise species name: strip intermediate-analysis suffixes
            for suffix in ('_ani_filtered', '_filtered', '_old'):
                if species.endswith(suffix):
                    species = species[:-len(suffix)]
                    break
            acc_match = ACC_RE.search(parts[1])
            if acc_match:
                accession = acc_match.group(1)
                col_info[col] = (species, accession)

    # Reshape to long format
    records = []
    for _, row in og.iterrows():
        orthogroup = row['Orthogroup']
        for col, (species, accession) in col_info.items():
            proteins = row.get(col, '')
            if pd.isna(proteins) or proteins == '':
                continue

            # Split protein IDs
            for protein_id in str(proteins).split(', '):
                protein_id = protein_id.strip()
                if protein_id:
                    records.append({
                        'Orthogroup': orthogroup,
                        'Species': species,
                        'Assembly_Accession': accession,
                        'Protein_ID': protein_id
                    })

    return pd.DataFrame(records)


def build_genome_species_mapping(genecount_df: pd.DataFrame) -> Dict[str, str]:
    """
    Build mapping from genome accession to species from gene count dataframe.

    Normalises species labels stripping common suffixes added during
    intermediate analysis steps (e.g. ``niger_ani_filtered`` -> ``niger``).

    Parameters
    ----------
    genecount_df : pd.DataFrame
        Gene count matrix with columns like 'species__accession'

    Returns
    -------
    Dict[str, str]
        Mapping from accession to canonical species name
    """
    mapping = {}
    for col in genecount_df.columns:
        if col == 'Orthogroup':
            continue
        if '__' in col:
            parts = col.split('__')
            species = parts[0]
            # Normalise: strip the '_ani_filtered' / '_filtered' / '_old'
            # suffixes that were used for intermediate analyses
            for suffix in ('_ani_filtered', '_filtered', '_old'):
                if species.endswith(suffix):
                    species = species[:-len(suffix)]
                    break
            accession = parts[1] if len(parts) > 1 else col
            mapping[accession] = species

    return mapping


def build_pav_from_genecount(genecount_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert gene count matrix to presence/absence (PAV) matrix.

    Parameters
    ----------
    genecount_df : pd.DataFrame
        Gene count matrix with Orthogroup as first column

    Returns
    -------
    pd.DataFrame
        PAV matrix (0/1) indexed by Orthogroup
    """
    pav = genecount_df.set_index('Orthogroup')
    pav = (pav > 0).astype(int)
    return pav


def split_pav_by_species(pav: pd.DataFrame, genome_mapping: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """
    Split combined PAV matrix into per-species PAV matrices.

    Parameters
    ----------
    pav : pd.DataFrame
        Combined PAV matrix
    genome_mapping : Dict[str, str]
        Mapping from genome column to species

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary of species -> PAV matrix
    """
    species_pav = {}

    for species in SPECIES_LIST:
        # Find columns belonging to this species
        species_cols = [col for col in pav.columns
                       if col.startswith(f"{species}__")]

        if species_cols:
            sp_pav = pav[species_cols].copy()
            # Rename columns to just accession
            sp_pav.columns = [col.replace(f"{species}__", "") for col in sp_pav.columns]
            species_pav[species] = sp_pav

    return species_pav


def load_species_metadata(species_dirs: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """
    Load metadata for each species.

    Parameters
    ----------
    species_dirs : Dict[str, str]
        Mapping from species name to metadata CSV path

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary of species -> metadata DataFrame
    """
    metadata = {}
    for species, path in species_dirs.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            metadata[species] = df
            print(f"Loaded {len(df)} samples for {species}")
        else:
            print(f"[WARN] Metadata not found: {path}")

    return metadata


def load_classified_metadata(isolation_csv: str) -> pd.DataFrame:
    """
    Load pre-classified isolation source metadata.

    Parameters
    ----------
    isolation_csv : str
        Path to isolation_categorized_5_class.csv or similar

    Returns
    -------
    pd.DataFrame
        Metadata with IsolationClass column
    """
    df = pd.read_csv(isolation_csv)

    # Standardize column names if needed
    col_map = {}
    for col in df.columns:
        col_lower = col.lower().replace(' ', '_').replace('-', '_')
        # Match exact 'isolationclass' or 'isolation_class' but NOT 'isolationsubclass'
        if col_lower in ('isolationclass', 'isolation_class'):
            col_map[col] = 'IsolationClass'
        elif col_lower in ('assembly_accession', 'assemblyaccession'):
            col_map[col] = 'Assembly_Accession'

    if col_map:
        df = df.rename(columns=col_map)

    # Standardize Assembly_Accession column name
    if 'Assembly Accession' in df.columns and 'Assembly_Accession' not in df.columns:
        df = df.rename(columns={'Assembly Accession': 'Assembly_Accession'})

    return df

# =============================================================================
# CROSS-SPECIES CONVERGENT HIT DETECTION (Section 6.1)
# =============================================================================

def identify_convergent_hits(
    species_results: Dict[str, pd.DataFrame],
    fdr_threshold: float = 0.1,
    min_species_significant: int = 2,
    min_species_consistent: int = 3,
    moderate_p_threshold: float = 0.1
) -> pd.DataFrame:
    """
    Identify cross-species convergent hits based on:
    - Significant (FDR < threshold) in >= min_species_significant species with consistent direction
    - OR consistent direction in >= min_species_consistent species with moderate evidence

    Parameters
    ----------
    species_results : Dict[str, pd.DataFrame]
        Per-species association results
    fdr_threshold : float
        FDR threshold for significance
    min_species_significant : int
        Minimum species with significant results (default 2)
    min_species_consistent : int
        Minimum species with consistent direction for moderate evidence (default 3)
    moderate_p_threshold : float
        P-value threshold for moderate evidence

    Returns
    -------
    pd.DataFrame
        Convergent hits with summary statistics
    """
    # Collect all orthogroups
    all_orthogroups = set()
    for res in species_results.values():
        all_orthogroups.update(res['Orthogroup'].values)

    records = []

    for og in all_orthogroups:
        species_data = []

        for species, res in species_results.items():
            og_data = res[res['Orthogroup'] == og]
            if len(og_data) == 0:
                continue

            row = og_data.iloc[0]
            species_data.append({
                'species': species,
                'p_value': row['p_value'],
                'q_value': row['q_value'],
                'odds_ratio': row['odds_ratio'],
                'direction': row['direction']
            })

        if len(species_data) < 2:
            continue

        # Count significant and direction-consistent species
        enriched_sig = sum(1 for d in species_data
                         if d['q_value'] < fdr_threshold and d['direction'] == 'enriched_in_case')
        depleted_sig = sum(1 for d in species_data
                          if d['q_value'] < fdr_threshold and d['direction'] == 'depleted_in_case')

        # Count direction-consistent with moderate evidence
        enriched_mod = sum(1 for d in species_data
                          if d['p_value'] < moderate_p_threshold and d['direction'] == 'enriched_in_case')
        depleted_mod = sum(1 for d in species_data
                          if d['p_value'] < moderate_p_threshold and d['direction'] == 'depleted_in_case')

        # Determine convergence
        is_convergent = False
        convergence_type = 'none'
        convergent_direction = 'none'

        # Criterion 1: Significant in >= min_species with consistent direction
        if enriched_sig >= min_species_significant:
            is_convergent = True
            convergence_type = 'multi_significant'
            convergent_direction = 'enriched_in_case'
        elif depleted_sig >= min_species_significant:
            is_convergent = True
            convergence_type = 'multi_significant'
            convergent_direction = 'depleted_in_case'

        # Criterion 2: Consistent direction with moderate evidence
        elif enriched_mod >= min_species_consistent:
            is_convergent = True
            convergence_type = 'consistent_moderate'
            convergent_direction = 'enriched_in_case'
        elif depleted_mod >= min_species_consistent:
            is_convergent = True
            convergence_type = 'consistent_moderate'
            convergent_direction = 'depleted_in_case'

        if is_convergent:
            # Calculate summary statistics
            p_values = [d['p_value'] for d in species_data if not np.isnan(d['p_value'])]
            q_values = [d['q_value'] for d in species_data if not np.isnan(d['q_value'])]
            odds_ratios = [d['odds_ratio'] for d in species_data
                          if not np.isnan(d['odds_ratio']) and np.isfinite(d['odds_ratio'])]

            records.append({
                'Orthogroup': og,
                'n_species_tested': len(species_data),
                'n_species_significant': enriched_sig + depleted_sig,
                'n_enriched_sig': enriched_sig,
                'n_depleted_sig': depleted_sig,
                'n_enriched_mod': enriched_mod,
                'n_depleted_mod': depleted_mod,
                'convergence_type': convergence_type,
                'convergent_direction': convergent_direction,
                'min_p_value': min(p_values) if p_values else np.nan,
                'median_p_value': np.median(p_values) if p_values else np.nan,
                'mean_log_odds': np.mean(np.log(odds_ratios)) if odds_ratios else np.nan,
                'species_significant': ';'.join([d['species'] for d in species_data
                                                 if d['q_value'] < fdr_threshold]),
                'species_tested': ';'.join([d['species'] for d in species_data])
            })

    if not records:
        return pd.DataFrame()

    convergent_df = pd.DataFrame(records)
    convergent_df = convergent_df.sort_values('min_p_value')

    return convergent_df


# =============================================================================
# META-ANALYSIS (Section 6.2)
# =============================================================================

def fisher_combine_pvalues(p_values: List[float]) -> Tuple[float, float]:
    """
    Combine p-values using Fisher's method.

    chi2 = -2 * sum(log(p_i))

    Parameters
    ----------
    p_values : List[float]
        List of p-values

    Returns
    -------
    Tuple[float, float]
        (chi2_statistic, combined_p_value)
    """
    p_values = [p for p in p_values if not np.isnan(p) and p > 0]
    if len(p_values) == 0:
        return np.nan, np.nan

    # Clip very small p-values to avoid -inf in log
    p_values = [max(p, 1e-300) for p in p_values]

    chi2_stat = -2 * sum(np.log(p) for p in p_values)
    df = 2 * len(p_values)

    combined_p = 1 - chi2.cdf(chi2_stat, df)

    return chi2_stat, combined_p


def stouffer_combine_pvalues(p_values: List[float], weights: Optional[List[float]] = None) -> Tuple[float, float]:
    """
    Combine p-values using Stouffer's Z-score method.

    Z = sum(w_i * Z_i) / sqrt(sum(w_i^2))

    Parameters
    ----------
    p_values : List[float]
        List of p-values
    weights : List[float], optional
        Weights for each p-value (default: equal weights)

    Returns
    -------
    Tuple[float, float]
        (z_statistic, combined_p_value)
    """
    p_values = [p for p in p_values if not np.isnan(p) and 0 < p < 1]
    if len(p_values) == 0:
        return np.nan, np.nan

    if weights is None:
        weights = [1.0] * len(p_values)
    else:
        weights = [w for w, p in zip(weights, p_values) if not np.isnan(p) and 0 < p < 1]

    # Convert p-values to Z-scores
    z_scores = [norm.ppf(1 - p) for p in p_values]

    # Weighted combination
    z_combined = sum(w * z for w, z in zip(weights, z_scores)) / np.sqrt(sum(w**2 for w in weights))

    # Two-tailed p-value
    combined_p = 2 * (1 - norm.cdf(abs(z_combined)))

    return z_combined, combined_p


def run_meta_analysis(
    species_results: Dict[str, pd.DataFrame],
    method: str = 'fisher',
    weights: Optional[Dict[str, float]] = None
) -> pd.DataFrame:
    """
    Run meta-analysis combining per-species p-values for each orthogroup.

    Parameters
    ----------
    species_results : Dict[str, pd.DataFrame]
        Per-species association results
    method : str
        'fisher' or 'stouffer'
    weights : Dict[str, float], optional
        Per-species weights for Stouffer method

    Returns
    -------
    pd.DataFrame
        Meta-analysis results with combined p-values
    """
    # Collect all orthogroups
    all_orthogroups = set()
    for res in species_results.values():
        all_orthogroups.update(res['Orthogroup'].values)

    records = []

    for og in all_orthogroups:
        species_data = {}

        for species, res in species_results.items():
            og_data = res[res['Orthogroup'] == og]
            if len(og_data) == 0:
                continue

            row = og_data.iloc[0]
            species_data[species] = {
                'p_value': row['p_value'],
                'odds_ratio': row['odds_ratio'],
                'direction': row['direction']
            }

        if len(species_data) < 2:
            continue

        # Get p-values and directions
        p_values = [d['p_value'] for d in species_data.values()]
        directions = [d['direction'] for d in species_data.values()]
        odds_ratios = [d['odds_ratio'] for d in species_data.values()]

        # Check direction consistency
        n_enriched = sum(1 for d in directions if d == 'enriched_in_case')
        n_depleted = sum(1 for d in directions if d == 'depleted_in_case')
        direction_consistent = (n_enriched == len(directions) or n_depleted == len(directions))

        # Check for heterogeneity (one species driving signal)
        valid_p = [p for p in p_values if not np.isnan(p)]
        if len(valid_p) >= 2:
            min_p = min(valid_p)
            second_min_p = sorted(valid_p)[1]
            heterogeneity_flag = (min_p < 0.01 and second_min_p > 0.1)
        else:
            heterogeneity_flag = False

        # Combine p-values
        if method == 'fisher':
            stat, meta_p = fisher_combine_pvalues(p_values)
        else:
            w = [weights.get(s, 1.0) for s in species_data.keys()] if weights else None
            stat, meta_p = stouffer_combine_pvalues(p_values, w)

        # Calculate mean effect size
        valid_or = [or_ for or_ in odds_ratios if not np.isnan(or_) and np.isfinite(or_)]
        mean_log_or = np.mean([np.log(or_) for or_ in valid_or]) if valid_or else np.nan

        records.append({
            'Orthogroup': og,
            'n_species': len(species_data),
            'meta_p_value': meta_p,
            'meta_statistic': stat,
            'n_enriched': n_enriched,
            'n_depleted': n_depleted,
            'direction_consistent': direction_consistent,
            'heterogeneity_flag': heterogeneity_flag,
            'mean_log_odds_ratio': mean_log_or,
            'consensus_direction': 'enriched' if n_enriched > n_depleted else ('depleted' if n_depleted > n_enriched else 'mixed'),
            'species_tested': ';'.join(species_data.keys()),
            'per_species_p': ';'.join([f"{s}:{d['p_value']:.2e}" for s, d in species_data.items()])
        })

    if not records:
        return pd.DataFrame()

    meta_df = pd.DataFrame(records)

    # Calculate FDR
    valid_mask = ~meta_df['meta_p_value'].isna()
    meta_df['meta_q_value'] = 1.0
    if valid_mask.sum() > 0:
        _, q_values, _, _ = multipletests(
            meta_df.loc[valid_mask, 'meta_p_value'].values,
            method='fdr_bh'
        )
        meta_df.loc[valid_mask, 'meta_q_value'] = q_values

    meta_df = meta_df.sort_values('meta_p_value')

    return meta_df


# =============================================================================
# LEAVE-ONE-SPECIES-OUT ROBUSTNESS (Section 6.3)
# =============================================================================

def run_leave_one_out_analysis(
    species_results: Dict[str, pd.DataFrame],
    convergent_hits: pd.DataFrame,
    meta_method: str = 'fisher'
) -> pd.DataFrame:
    """
    Run leave-one-species-out analysis for robustness testing.

    For each convergent hit, recalculate meta p-value leaving out each species.
    Hits that remain significant are the most robust.

    Parameters
    ----------
    species_results : Dict[str, pd.DataFrame]
        Per-species association results
    convergent_hits : pd.DataFrame
        Identified convergent hits
    meta_method : str
        Meta-analysis method ('fisher' or 'stouffer')

    Returns
    -------
    pd.DataFrame
        Robustness analysis results
    """
    if convergent_hits.empty:
        return pd.DataFrame()

    species_list = list(species_results.keys())
    records = []

    for _, hit in convergent_hits.iterrows():
        og = hit['Orthogroup']

        # Get original meta p-value
        original_meta = run_meta_analysis({k: v for k, v in species_results.items()}, method=meta_method)
        og_meta = original_meta[original_meta['Orthogroup'] == og]

        if len(og_meta) == 0:
            continue

        original_p = og_meta.iloc[0]['meta_p_value']
        original_q = og_meta.iloc[0]['meta_q_value']

        # Leave-one-out analysis
        loo_results = {}
        for species in species_list:
            if species not in species_results:
                continue

            # Create subset without this species
            subset_results = {k: v for k, v in species_results.items() if k != species}

            if len(subset_results) < 2:
                continue

            # Run meta-analysis
            loo_meta = run_meta_analysis(subset_results, method=meta_method)
            loo_og = loo_meta[loo_meta['Orthogroup'] == og]

            if len(loo_og) > 0:
                loo_results[species] = {
                    'p_value': loo_og.iloc[0]['meta_p_value'],
                    'q_value': loo_og.iloc[0]['meta_q_value']
                }

        # Determine robustness
        n_still_significant = sum(1 for r in loo_results.values() if r['q_value'] < 0.1)
        max_loo_p = max([r['p_value'] for r in loo_results.values()]) if loo_results else np.nan

        # Find which species, when removed, causes biggest change
        species_impact = {}
        for species, r in loo_results.items():
            if not np.isnan(r['p_value']) and not np.isnan(original_p):
                impact = np.log10(r['p_value']) - np.log10(original_p)
                species_impact[species] = impact

        driving_species = max(species_impact, key=lambda k: abs(species_impact.get(k, 0))) if species_impact else 'none'

        records.append({
            'Orthogroup': og,
            'original_meta_p': original_p,
            'original_meta_q': original_q,
            'n_loo_tests': len(loo_results),
            'n_still_significant': n_still_significant,
            'max_loo_p': max_loo_p,
            'robust': n_still_significant == len(loo_results),
            'driving_species': driving_species,
            'convergence_type': hit.get('convergence_type', 'unknown'),
            'convergent_direction': hit.get('convergent_direction', 'unknown')
        })

    robustness_df = pd.DataFrame(records)

    if not robustness_df.empty:
        robustness_df = robustness_df.sort_values('original_meta_p')

    return robustness_df


def get_tier_classification(robustness_df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify hits into tiers based on robustness.

    Tier 1: Robust (significant after removing any species)
    Tier 2: Partially robust (significant after removing most species)
    Tier 3: Species-dependent (loses significance when key species removed)

    Parameters
    ----------
    robustness_df : pd.DataFrame
        Robustness analysis results

    Returns
    -------
    pd.DataFrame
        Results with tier classification
    """
    if robustness_df.empty:
        return robustness_df

    df = robustness_df.copy()

    def classify_tier(row):
        if row['robust']:
            return 'Tier1_Robust'
        elif row['n_still_significant'] >= row['n_loo_tests'] - 1:
            return 'Tier2_PartiallyRobust'
        else:
            return 'Tier3_SpeciesDependent'

    df['tier'] = df.apply(classify_tier, axis=1)

    return df


# =============================================================================
# ANNOTATION AND INTERPRETATION
# =============================================================================

def annotate_convergent_hits(
    convergent_hits: pd.DataFrame,
    og_annotations: pd.DataFrame,
    annotation_cols: List[str] = None
) -> pd.DataFrame:
    """
    Add functional annotations to convergent hits.

    Parameters
    ----------
    convergent_hits : pd.DataFrame
        Convergent hit results
    og_annotations : pd.DataFrame
        Orthogroup annotations (e.g., from build_og_consensus_table)
    annotation_cols : List[str], optional
        Columns to include from annotations

    Returns
    -------
    pd.DataFrame
        Annotated convergent hits
    """
    if convergent_hits.empty:
        return convergent_hits

    if annotation_cols is None:
        annotation_cols = ['Description', 'GOs', 'PFAMs', 'CAZy', 'COG_category',
                          'KEGG_ko', 'KEGG_Pathway']

    # Filter to available columns
    available_cols = [c for c in annotation_cols if c in og_annotations.columns]

    if not available_cols:
        return convergent_hits

    merge_cols = ['Orthogroup'] + available_cols
    annot_subset = og_annotations[merge_cols].drop_duplicates(subset=['Orthogroup'])

    annotated = convergent_hits.merge(annot_subset, on='Orthogroup', how='left')

    return annotated


def summarize_convergent_functions(
    convergent_hits: pd.DataFrame,
    annotation_col: str = 'COG_category',
    top_n: int = 10
) -> pd.DataFrame:
    """
    Summarize functional categories of convergent hits.

    Parameters
    ----------
    convergent_hits : pd.DataFrame
        Annotated convergent hits
    annotation_col : str
        Annotation column to summarize
    top_n : int
        Number of top categories to return

    Returns
    -------
    pd.DataFrame
        Summary of functional categories
    """
    if annotation_col not in convergent_hits.columns:
        print(f"Column {annotation_col} not found")
        return pd.DataFrame()

    # Count categories
    all_terms = []
    for val in convergent_hits[annotation_col].dropna():
        terms = str(val).split(';')
        all_terms.extend([t.strip() for t in terms if t.strip()])

    term_counts = Counter(all_terms)

    summary = pd.DataFrame([
        {'Term': term, 'Count': count, 'Percentage': count / len(convergent_hits) * 100}
        for term, count in term_counts.most_common(top_n)
    ])

    return summary


def generate_cross_species_report(
    convergent_hits: pd.DataFrame,
    meta_results: pd.DataFrame,
    robustness_results: pd.DataFrame,
    contrast_name: str
) -> str:
    """
    Generate a text report summarizing cross-species analysis results.

    Parameters
    ----------
    convergent_hits : pd.DataFrame
        Convergent hit results
    meta_results : pd.DataFrame
        Meta-analysis results
    robustness_results : pd.DataFrame
        Robustness analysis results
    contrast_name : str
        Name of the contrast

    Returns
    -------
    str
        Formatted report text
    """
    report = []
    report.append("=" * 80)
    report.append(f"CROSS-SPECIES ANALYSIS REPORT: {contrast_name}")
    report.append("=" * 80)

    # Convergent hits summary
    report.append("\n1. CONVERGENT HITS (Section 6.1)")
    report.append("-" * 40)
    if not convergent_hits.empty:
        n_total = len(convergent_hits)
        report.append(f"Total convergent hits: {n_total}")

        # Optional richer breakdown when extended columns are available
        if 'convergence_type' in convergent_hits.columns:
            n_multi_sig = (convergent_hits['convergence_type'] == 'multi_significant').sum()
            n_consistent = (convergent_hits['convergence_type'] == 'consistent_moderate').sum()
            report.append(f"  - Multi-species significant: {n_multi_sig}")
            report.append(f"  - Consistent with moderate evidence: {n_consistent}")
        elif 'category' in convergent_hits.columns:
            for cat, n in convergent_hits['category'].value_counts().items():
                report.append(f"  - {cat}: {n}")

        if 'convergent_direction' in convergent_hits.columns:
            n_enriched = (convergent_hits['convergent_direction'] == 'enriched_in_case').sum()
            n_depleted = (convergent_hits['convergent_direction'] == 'depleted_in_case').sum()
            report.append(f"  - Enriched in case group: {n_enriched}")
            report.append(f"  - Depleted in case group: {n_depleted}")
    else:
        report.append("No convergent hits identified")

    # Meta-analysis summary
    report.append("\n2. META-ANALYSIS RESULTS (Section 6.2)")
    report.append("-" * 40)
    if not meta_results.empty:
        n_tested = len(meta_results)
        report.append(f"Orthogroups tested: {n_tested}")
        if 'meta_q_value' in meta_results.columns:
            n_sig_01 = (meta_results['meta_q_value'] < 0.1).sum()
            n_sig_05 = (meta_results['meta_q_value'] < 0.05).sum()
            report.append(f"Significant at FDR < 0.1: {n_sig_01}")
            report.append(f"Significant at FDR < 0.05: {n_sig_05}")
        if 'direction_consistent' in meta_results.columns:
            report.append(f"Direction-consistent across species: {meta_results['direction_consistent'].sum()}")
        if 'heterogeneity_flag' in meta_results.columns:
            report.append(f"Heterogeneous (single species driving): {meta_results['heterogeneity_flag'].sum()}")
    else:
        report.append("No meta-analysis results available")

    # Robustness summary
    report.append("\n3. ROBUSTNESS ANALYSIS (Section 6.3)")
    report.append("-" * 40)
    if not robustness_results.empty:
        if 'tier' in robustness_results.columns:
            tier_counts = robustness_results['tier'].value_counts()
            for tier, count in tier_counts.items():
                report.append(f"{tier}: {count}")

        n_robust = robustness_results['robust'].sum() if 'robust' in robustness_results.columns else 0
        report.append(f"\nFully robust hits (persist after removing any species): {n_robust}")

        if 'driving_species' in robustness_results.columns:
            driver_counts = robustness_results['driving_species'].value_counts()
            report.append("\nSignal-driving species distribution:")
            for species, count in driver_counts.items():
                report.append(f"  {species}: {count}")
    else:
        report.append("No robustness results available")

    # Key findings
    report.append("\n4. KEY FINDINGS")
    report.append("-" * 40)

    if not robustness_results.empty and 'tier' in robustness_results.columns:
        tier1 = robustness_results[robustness_results['tier'] == 'Tier1_Robust']
        if len(tier1) > 0:
            report.append("\nTOP-TIER ROBUST HITS (strongest cross-species evidence):")
            for _, hit in tier1.head(10).iterrows():
                report.append(f"  {hit['Orthogroup']}: direction={hit.get('convergent_direction', 'unknown')}")

    report.append("\n" + "=" * 80)

    return "\n".join(report)


# =============================================================================
# MAIN WORKFLOW FUNCTIONS
# =============================================================================

def run_full_cross_species_analysis(
    genecount_path: str,
    metadata_path: str,
    output_dir: str,
    contrasts: Optional[List[Dict]] = None,
    fdr_threshold: float = 0.1
) -> Dict[str, Any]:
    """
    Run complete cross-species analysis pipeline.

    Parameters
    ----------
    genecount_path : str
        Path to combined Orthogroups.GeneCount.tsv
    metadata_path : str
        Path to classified metadata CSV
    output_dir : str
        Output directory for results
    contrasts : List[Dict], optional
        Contrast definitions (uses defaults if None)
    fdr_threshold : float
        FDR threshold for significance

    Returns
    -------
    Dict[str, Any]
        Dictionary containing all results
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("CROSS-SPECIES PANGENOME ANALYSIS")
    print("=" * 60)

    # Load data
    print("\n1. Loading combined orthogroup data...")
    genecount = load_combined_orthogroups_genecount(genecount_path)
    pav = build_pav_from_genecount(genecount)
    genome_mapping = build_genome_species_mapping(genecount)
    print(f"   Loaded {len(pav)} orthogroups across {len(pav.columns)} genomes")

    # Split by species
    species_pav = split_pav_by_species(pav, genome_mapping)
    print(f"   Species: {', '.join(species_pav.keys())}")

    # Load metadata
    print("\n2. Loading metadata...")
    metadata = load_classified_metadata(metadata_path)
    print(f"   Loaded {len(metadata)} samples")

    # Create phenotypes
    print("\n3. Creating phenotype data...")
    phenotype = create_phenotype_from_isolation_class(pav, metadata)

    # Split phenotypes by species
    species_phenotypes = {}
    for species in species_pav.keys():
        sp_pheno = phenotype[phenotype['species'] == species]
        if len(sp_pheno) > 0:
            species_phenotypes[species] = sp_pheno

    # Get contrasts
    if contrasts is None:
        contrasts = define_isolation_contrasts()

    # Results storage
    all_results = {}

    for contrast in contrasts:
        contrast_name = contrast['name']
        print(f"\n{'='*60}")
        print(f"CONTRAST: {contrast_name}")
        print(f"{'='*60}")

        # 1. Per-species association tests
        species_results = run_all_species_association(
            species_pav, species_phenotypes, contrast
        )

        if len(species_results) < 2:
            print("  Insufficient species with results, skipping cross-species analysis")
            continue

        # 2. Identify convergent hits (Section 6.1)
        print("\n  Identifying convergent hits...")
        convergent = identify_convergent_hits(
            species_results,
            fdr_threshold=fdr_threshold
        )
        print(f"    Found {len(convergent)} convergent hits")

        # 3. Meta-analysis (Section 6.2)
        print("\n  Running meta-analysis...")
        meta = run_meta_analysis(species_results, method='fisher')
        n_sig = (meta['meta_q_value'] < fdr_threshold).sum()
        print(f"    {n_sig} significant at FDR < {fdr_threshold}")

        # 4. Robustness analysis (Section 6.3)
        print("\n  Running leave-one-out robustness analysis...")
        if not convergent.empty:
            robustness = run_leave_one_out_analysis(species_results, convergent)
            robustness = get_tier_classification(robustness)
            n_robust = robustness['robust'].sum() if 'robust' in robustness.columns else 0
            print(f"    {n_robust} robust hits identified")
        else:
            robustness = pd.DataFrame()

        # Store results
        all_results[contrast_name] = {
            'species_results': species_results,
            'convergent_hits': convergent,
            'meta_analysis': meta,
            'robustness': robustness
        }

        # Generate report
        report = generate_cross_species_report(
            convergent, meta, robustness, contrast_name
        )
        print(report)

        # Save results
        contrast_dir = os.path.join(output_dir, contrast_name)
        os.makedirs(contrast_dir, exist_ok=True)

        if not convergent.empty:
            convergent.to_csv(os.path.join(contrast_dir, 'convergent_hits.csv'), index=False)
        meta.to_csv(os.path.join(contrast_dir, 'meta_analysis.csv'), index=False)
        if not robustness.empty:
            robustness.to_csv(os.path.join(contrast_dir, 'robustness.csv'), index=False)

        for species, res in species_results.items():
            res.to_csv(os.path.join(contrast_dir, f'{species}_association.csv'), index=False)

        with open(os.path.join(contrast_dir, 'report.txt'), 'w') as f:
            f.write(report)

    return all_results


# =============================================================================
# ANNOTATION MAPPING FROM INDIVIDUAL SPECIES TO COMBINED PANGENOME
# =============================================================================
#
# Key insight: OrthoFinder assigns NEW orthogroup IDs each run, so we cannot
# directly map orthogroup IDs between individual species runs and the combined run.
# Instead, we use PROTEIN IDs as the linking key - these are consistent across runs.
#
# Approach:
# 1. Load per-protein annotations from each species (EggNOG, dbCAN, InterProScan, SignalP)
# 2. Use the combined OrthoFinder long format to know which proteins are in which orthogroup
# 3. Aggregate protein-level annotations up to the combined orthogroup level
# =============================================================================

def load_species_protein_annotations(
    species_base_dirs: Dict[str, str],
    annotation_sources: List[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    Load per-protein annotation tables from each species.

    These are the raw annotation outputs (EggNOG, dbCAN, InterProScan, SignalP)
    which have Protein_ID as the key column.

    Parameters
    ----------
    species_base_dirs : Dict[str, str]
        Mapping from species name to base directory path
    annotation_sources : List[str], optional
        Which annotation sources to load ['eggnog', 'dbcan', 'interproscan', 'signalp']

    Returns
    -------
    Dict[str, pd.DataFrame]
        Per-species protein annotation tables (merged from all sources)
    """
    if annotation_sources is None:
        annotation_sources = ['eggnog', 'dbcan', 'interproscan', 'signalp']

    species_annotations = {}

    for species, base_dir in species_base_dirs.items():
        print(f"\nLoading annotations for {species}...")

        protein_annots = []

        # EggNOG annotations
        if 'eggnog' in annotation_sources:
            eggnog_dir = os.path.join(base_dir, 'eggnog_output')
            if os.path.isdir(eggnog_dir):
                eggnog_files = [f for f in os.listdir(eggnog_dir)
                               if f.endswith('.emapper.annotations')]
                for ef in eggnog_files[:50]:  # Limit for memory
                    try:
                        eggnog_path = os.path.join(eggnog_dir, ef)
                        df = _read_eggnog_simple(eggnog_path)
                        if not df.empty:
                            protein_annots.append(df)
                    except Exception as e:
                        pass
                print(f"  Loaded {len(eggnog_files)} EggNOG files")

        # dbCAN annotations
        if 'dbcan' in annotation_sources:
            dbcan_dir = os.path.join(base_dir, 'dbcan_output')
            if os.path.isdir(dbcan_dir):
                dbcan_count = 0
                for subdir in os.listdir(dbcan_dir)[:50]:
                    overview_path = os.path.join(dbcan_dir, subdir, 'overview.tsv')
                    if os.path.exists(overview_path):
                        try:
                            df = _read_dbcan_simple(overview_path)
                            if not df.empty:
                                protein_annots.append(df)
                                dbcan_count += 1
                        except Exception:
                            pass
                print(f"  Loaded {dbcan_count} dbCAN files")

        # InterProScan annotations
        if 'interproscan' in annotation_sources:
            interpro_dir = os.path.join(base_dir, 'interproscan_output')
            if os.path.isdir(interpro_dir):
                interpro_count = 0
                for subdir in os.listdir(interpro_dir)[:50]:
                    subdir_path = os.path.join(interpro_dir, subdir)
                    if os.path.isdir(subdir_path):
                        tsv_files = [f for f in os.listdir(subdir_path) if f.endswith('.tsv')]
                        for tsv in tsv_files[:1]:
                            try:
                                df = _read_interproscan_simple(os.path.join(subdir_path, tsv))
                                if not df.empty:
                                    protein_annots.append(df)
                                    interpro_count += 1
                            except Exception:
                                pass
                print(f"  Loaded {interpro_count} InterProScan files")

        # SignalP annotations
        if 'signalp' in annotation_sources:
            signalp_dir = os.path.join(base_dir, 'signalp_output')
            if os.path.isdir(signalp_dir):
                signalp_count = 0
                for subdir in os.listdir(signalp_dir)[:50]:
                    pred_path = os.path.join(signalp_dir, subdir, 'prediction_results.txt')
                    if os.path.exists(pred_path):
                        try:
                            df = _read_signalp_simple(pred_path)
                            if not df.empty:
                                protein_annots.append(df)
                                signalp_count += 1
                        except Exception:
                            pass
                print(f"  Loaded {signalp_count} SignalP files")

        # Merge all protein annotations for this species
        if protein_annots:
            merged = pd.concat(protein_annots, ignore_index=True)
            # Group by Protein_ID and take first non-empty value for each column
            merged = merged.groupby('Protein_ID', as_index=False).first()
            merged['Species'] = species
            species_annotations[species] = merged
            print(f"  Total: {len(merged)} proteins with annotations")
        else:
            print(f"  [WARN] No annotations found for {species}")

    return species_annotations


def _read_eggnog_simple(path: str) -> pd.DataFrame:
    """Simple EggNOG reader - extracts key columns."""
    # Find header line
    hdr_idx = None
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if line.startswith('#query') or line.startswith('query'):
                hdr_idx = i
                break

    if hdr_idx is None:
        return pd.DataFrame()

    df = pd.read_csv(path, sep='\t', header=0, skiprows=hdr_idx, dtype=str, na_filter=False)

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = c.lstrip('#').strip().lower()
        if cl in ('query', 'query_name'):
            col_map[c] = 'Protein_ID'
        elif cl == 'description':
            col_map[c] = 'Description_EggNog'
        elif cl in ('cog_category', 'cog', 'cogs'):
            col_map[c] = 'COG_category'
        elif cl in ('go', 'gos', 'go_terms'):
            col_map[c] = 'GOs_EggNog'
        elif cl in ('kegg_ko', 'ko'):
            col_map[c] = 'KEGG_ko'
        elif cl in ('kegg_pathway', 'pathway'):
            col_map[c] = 'KEGG_Pathway'
        elif cl in ('pfams', 'pfam'):
            col_map[c] = 'PFAMs_EggNog'
        elif cl in ('ec', 'ec_number'):
            col_map[c] = 'EC_EggNog'

    df = df.rename(columns=col_map)

    keep_cols = ['Protein_ID', 'Description_EggNog', 'COG_category', 'GOs_EggNog',
                 'KEGG_ko', 'KEGG_Pathway', 'PFAMs_EggNog', 'EC_EggNog']
    available = [c for c in keep_cols if c in df.columns]

    if 'Protein_ID' not in available:
        return pd.DataFrame()

    return df[available]


def _read_dbcan_simple(path: str) -> pd.DataFrame:
    """Simple dbCAN reader - extracts CAZy families."""
    df = pd.read_csv(path, sep='\t', dtype=str, na_filter=False)

    # Normalize column names
    col_map = {'Gene ID': 'Protein_ID', 'Recommend Results': 'CAZy_dbCAN'}
    df = df.rename(columns=col_map)

    if 'Protein_ID' not in df.columns:
        return pd.DataFrame()

    keep_cols = ['Protein_ID', 'CAZy_dbCAN']
    available = [c for c in keep_cols if c in df.columns]

    return df[available]


def _read_interproscan_simple(path: str) -> pd.DataFrame:
    """Simple InterProScan reader - extracts IPR and GO terms."""
    try:
        df = pd.read_csv(
            path, sep='\t', header=None, dtype=str, na_filter=False,
            names=['Protein_ID', 'MD5', 'Length', 'Analysis', 'Signature_Acc',
                   'Signature_Desc', 'Start', 'Stop', 'Score', 'Status', 'Date',
                   'IPR_Acc', 'IPR_Desc', 'GO', 'Pathways']
        )
    except Exception:
        return pd.DataFrame()

    # Aggregate per protein
    grouped = df.groupby('Protein_ID', as_index=False).agg({
        'IPR_Acc': lambda x: ';'.join([v for v in x if v and v != '-']),
        'GO': lambda x: ';'.join([v for v in x if v and v != '-']),
        'IPR_Desc': lambda x: ';'.join(set([v for v in x if v and v != '-']))[:200]
    })

    grouped = grouped.rename(columns={
        'IPR_Acc': 'interpro_IPR',
        'GO': 'GOs_Interpro',
        'IPR_Desc': 'Description_Interpro'
    })

    return grouped


def _read_signalp_simple(path: str) -> pd.DataFrame:
    """Simple SignalP reader."""
    rows = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                prot_id = parts[0].split()[0]
                prediction = parts[1]
                sp_prob = parts[3] if len(parts) > 3 else ''
                rows.append({
                    'Protein_ID': prot_id,
                    'signalp_Prediction': prediction,
                    'signalp_SP_Probability': sp_prob
                })

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def aggregate_protein_annotations_to_orthogroups(
    combined_og_long: pd.DataFrame,
    species_protein_annotations: Dict[str, pd.DataFrame],
    annotation_cols: List[str] = None
) -> pd.DataFrame:
    """
    Aggregate protein-level annotations up to orthogroup level.

    Uses the combined OrthoFinder long format to know which proteins belong
    to which orthogroup, then aggregates annotations using consensus.

    Parameters
    ----------
    combined_og_long : pd.DataFrame
        Combined pangenome long format with Orthogroup, Species, Protein_ID
    species_protein_annotations : Dict[str, pd.DataFrame]
        Per-species protein annotation tables
    annotation_cols : List[str], optional
        Annotation columns to aggregate

    Returns
    -------
    pd.DataFrame
        Orthogroup-level annotation table
    """
    if annotation_cols is None:
        annotation_cols = [
            'Description_EggNog', 'Description_Interpro',
            'COG_category', 'GOs_EggNog', 'GOs_Interpro',
            'KEGG_ko', 'KEGG_Pathway', 'PFAMs_EggNog',
            'EC_EggNog', 'CAZy_dbCAN', 'interpro_IPR',
            'signalp_Prediction', 'signalp_SP_Probability'
        ]

    print("\nAggregating protein annotations to orthogroups...")

    # Merge all species protein annotations
    all_protein_annots = []
    for species, annot_df in species_protein_annotations.items():
        if not annot_df.empty:
            all_protein_annots.append(annot_df)

    if not all_protein_annots:
        print("  [WARN] No protein annotations to aggregate")
        return pd.DataFrame({'Orthogroup': combined_og_long['Orthogroup'].unique()})

    protein_annots = pd.concat(all_protein_annots, ignore_index=True)

    # Remove duplicates (same protein from multiple sources)
    protein_annots = protein_annots.groupby('Protein_ID', as_index=False).first()

    print(f"  Total proteins with annotations: {len(protein_annots)}")

    # Join with combined OG long format
    og_with_annots = combined_og_long.merge(
        protein_annots,
        on='Protein_ID',
        how='left'
    )

    print(f"  Proteins matched to orthogroups: {og_with_annots['Description_EggNog'].notna().sum() if 'Description_EggNog' in og_with_annots.columns else 0}")

    # Aggregate by orthogroup
    def consensus_agg(series):
        """Get most common non-empty value."""
        vals = series.dropna().astype(str)
        vals = vals[vals != '']
        vals = vals[vals != '-']
        vals = vals[vals != 'nan']
        if len(vals) == 0:
            return ''
        counts = Counter(vals)
        return counts.most_common(1)[0][0]

    # Filter to available annotation columns
    available_cols = [c for c in annotation_cols if c in og_with_annots.columns]

    # Build aggregation dict
    agg_dict = {col: consensus_agg for col in available_cols}
    agg_dict['Protein_ID'] = 'count'  # Count proteins per OG
    if 'Species_x' in og_with_annots.columns:
        agg_dict['Species_x'] = lambda x: ';'.join(sorted(set(x.dropna())))
    elif 'Species' in og_with_annots.columns:
        agg_dict['Species'] = lambda x: ';'.join(sorted(set(x.dropna())))

    og_annotations = og_with_annots.groupby('Orthogroup', as_index=False).agg(agg_dict)

    # Rename columns
    og_annotations = og_annotations.rename(columns={
        'Protein_ID': 'n_proteins',
        'Species_x': 'species_contributing',
        'Species': 'species_contributing'
    })

    # Create combined columns
    if 'Description_EggNog' in og_annotations.columns:
        og_annotations['Description'] = og_annotations['Description_EggNog']
    if 'Description_Interpro' in og_annotations.columns and 'Description' in og_annotations.columns:
        # Fill missing with Interpro description
        mask = og_annotations['Description'] == ''
        og_annotations.loc[mask, 'Description'] = og_annotations.loc[mask, 'Description_Interpro']

    if 'GOs_EggNog' in og_annotations.columns and 'GOs_Interpro' in og_annotations.columns:
        og_annotations['GOs'] = og_annotations.apply(
            lambda r: ';'.join(set(filter(None, [str(r.get('GOs_EggNog', '')), str(r.get('GOs_Interpro', ''))]))),
            axis=1
        )
    elif 'GOs_EggNog' in og_annotations.columns:
        og_annotations['GOs'] = og_annotations['GOs_EggNog']

    if 'PFAMs_EggNog' in og_annotations.columns:
        og_annotations['PFAMs'] = og_annotations['PFAMs_EggNog']

    if 'CAZy_dbCAN' in og_annotations.columns:
        og_annotations['CAZy'] = og_annotations['CAZy_dbCAN']

    if 'EC_EggNog' in og_annotations.columns:
        og_annotations['EC'] = og_annotations['EC_EggNog']

    # Report coverage
    print("\nAnnotation coverage:")
    for col in ['Description', 'COG_category', 'GOs', 'PFAMs', 'CAZy', 'KEGG_ko',
                'interpro_IPR', 'signalp_Prediction']:
        if col in og_annotations.columns:
            n_annot = (og_annotations[col].fillna('').astype(str) != '').sum()
            pct = 100 * n_annot / len(og_annotations)
            print(f"  {col}: {n_annot}/{len(og_annotations)} ({pct:.1f}%)")

    return og_annotations


def transfer_annotations_to_combined_via_proteins(
    combined_og_long: pd.DataFrame,
    species_base_dirs: Dict[str, str],
    output_path: str = None
) -> pd.DataFrame:
    """
    Main function to transfer annotations from individual species to combined pangenome.

    Uses protein IDs as the linking key (NOT orthogroup IDs).

    Parameters
    ----------
    combined_og_long : pd.DataFrame
        Combined pangenome in long format
    species_base_dirs : Dict[str, str]
        Species directory mapping
    output_path : str, optional
        Path to save output

    Returns
    -------
    pd.DataFrame
        Combined orthogroup annotation table
    """
    print("=" * 60)
    print("TRANSFERRING ANNOTATIONS VIA PROTEIN IDS")
    print("=" * 60)

    # Step 1: Load protein-level annotations from each species
    print("\nStep 1: Loading protein-level annotations from each species...")
    species_protein_annotations = load_species_protein_annotations(species_base_dirs)

    if not species_protein_annotations:
        print("[WARN] No protein annotations loaded")
        return pd.DataFrame({'Orthogroup': combined_og_long['Orthogroup'].unique()})

    # Step 2: Aggregate to orthogroup level
    print("\nStep 2: Aggregating to orthogroup level...")
    og_annotations = aggregate_protein_annotations_to_orthogroups(
        combined_og_long,
        species_protein_annotations
    )

    # Save if requested
    if output_path:
        og_annotations.to_csv(output_path, sep='\t', index=False)
        print(f"\nSaved to: {output_path}")

    return og_annotations


def load_species_bgc_gcf_data(species_base_dirs: Dict[str, str]) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Load BGC (antiSMASH) and GCF (BiG-SCAPE) data from each species.

    Parses BiG-SCAPE output structure:
    - record_annotations.tsv for BGC types
    - *_clustering_*.tsv files for GCF assignments

    Parameters
    ----------
    species_base_dirs : Dict[str, str]
        Mapping from species name to base directory path

    Returns
    -------
    Dict[str, Dict[str, pd.DataFrame]]
        Nested dict: species -> {'bgc': bgc_df, 'gcf': gcf_df}
    """
    species_bgc_data = {}

    for species, base_dir in species_base_dirs.items():
        bgc_data = {}

        # Find BiG-SCAPE output directory
        bigscape_dir = os.path.join(base_dir, 'bigscape_output', 'output_files')

        if os.path.isdir(bigscape_dir):
            # Find the clustering output subdirectory (e.g., 2025-12-18_14-33-21_c0.3)
            subdirs = [d for d in os.listdir(bigscape_dir)
                      if os.path.isdir(os.path.join(bigscape_dir, d)) and 'c0.' in d]

            if subdirs:
                # Use most recent
                clustering_dir = os.path.join(bigscape_dir, sorted(subdirs)[-1])

                # Load record_annotations.tsv for BGC data
                annot_path = os.path.join(clustering_dir, 'record_annotations.tsv')
                if os.path.exists(annot_path):
                    try:
                        bgc_df = pd.read_csv(annot_path, sep='\t', dtype=str).fillna('')
                        # Extract assembly accession from Record column
                        bgc_df['Assembly_Accession'] = bgc_df['Record'].apply(
                            lambda x: ACC_RE.search(x).group(1) if ACC_RE.search(x) else ''
                        )
                        bgc_data['bgc'] = bgc_df
                        print(f"Loaded BGC data for {species}: {len(bgc_df)} regions from BiG-SCAPE")
                    except Exception as e:
                        print(f"[WARN] Failed to load BGC data for {species}: {e}")

                # Load GCF clustering files
                gcf_dfs = []
                for subdir in os.listdir(clustering_dir):
                    subdir_path = os.path.join(clustering_dir, subdir)
                    if os.path.isdir(subdir_path):
                        # Look for clustering file
                        for f in os.listdir(subdir_path):
                            if f.endswith('_clustering_c0.3.tsv'):
                                try:
                                    gcf_df = pd.read_csv(os.path.join(subdir_path, f), sep='\t', dtype=str).fillna('')
                                    gcf_df['BGC_Category'] = subdir
                                    gcf_df['Assembly_Accession'] = gcf_df['Record'].apply(
                                        lambda x: ACC_RE.search(x).group(1) if ACC_RE.search(x) else ''
                                    )
                                    gcf_dfs.append(gcf_df)
                                except Exception:
                                    pass

                if gcf_dfs:
                    combined_gcf = pd.concat(gcf_dfs, ignore_index=True)
                    bgc_data['gcf'] = combined_gcf
                    print(f"Loaded GCF data for {species}: {len(combined_gcf)} mappings")

        # Fallback: check for pre-existing summary files
        if 'bgc' not in bgc_data:
            for p in [os.path.join(base_dir, 'bgc_annotations.tsv'),
                     os.path.join(base_dir, 'antismash_summary.tsv')]:
                if os.path.exists(p):
                    try:
                        bgc_data['bgc'] = pd.read_csv(p, sep='\t', dtype=str).fillna('')
                        print(f"Loaded BGC data for {species}: {len(bgc_data['bgc'])} regions")
                        break
                    except Exception as e:
                        print(f"[WARN] Failed to load BGC data for {species}: {e}")

        if bgc_data:
            species_bgc_data[species] = bgc_data

    return species_bgc_data


def map_bgc_to_orthogroups(
    bgc_df: pd.DataFrame,
    og_long_df: pd.DataFrame,
    gene_col: str = 'genes',
    bgc_id_col: str = 'bgc_id'
) -> pd.DataFrame:
    """
    Map BGCs to orthogroups based on gene content.

    Parameters
    ----------
    bgc_df : pd.DataFrame
        BGC annotation table with gene lists
    og_long_df : pd.DataFrame
        Orthogroup long format with Protein_ID
    gene_col : str
        Column containing gene/protein IDs in BGC df
    bgc_id_col : str
        Column containing BGC identifiers

    Returns
    -------
    pd.DataFrame
        BGC to orthogroup mapping
    """
    # Create protein to OG lookup
    protein_to_og = dict(zip(og_long_df['Protein_ID'], og_long_df['Orthogroup']))

    records = []

    for _, row in bgc_df.iterrows():
        bgc_id = row.get(bgc_id_col, '')
        genes_str = row.get(gene_col, '')

        if not genes_str:
            continue

        # Parse gene list
        genes = [g.strip() for g in str(genes_str).split(';') if g.strip()]

        # Map to orthogroups
        mapped_ogs = set()
        for gene in genes:
            og = protein_to_og.get(gene)
            if og:
                mapped_ogs.add(og)

        if mapped_ogs:
            records.append({
                'bgc_id': bgc_id,
                'n_genes': len(genes),
                'n_ogs_mapped': len(mapped_ogs),
                'orthogroups': ';'.join(sorted(mapped_ogs)),
                'product': row.get('product', row.get('raw_product', '')),
                'bgc_type': row.get('bgc_type', row.get('type', ''))
            })

    if not records:
        return pd.DataFrame()

    return pd.DataFrame(records)


def create_og_to_bgc_mapping(
    bgc_og_mapping: pd.DataFrame
) -> pd.DataFrame:
    """
    Invert BGC-to-OG mapping to get OG-to-BGC mapping.

    Parameters
    ----------
    bgc_og_mapping : pd.DataFrame
        BGC to orthogroup mapping

    Returns
    -------
    pd.DataFrame
        Orthogroup to BGC mapping
    """
    records = []

    for _, row in bgc_og_mapping.iterrows():
        bgc_id = row['bgc_id']
        ogs = str(row['orthogroups']).split(';')
        product = row.get('product', '')
        bgc_type = row.get('bgc_type', '')

        for og in ogs:
            if og:
                records.append({
                    'Orthogroup': og,
                    'bgc_id': bgc_id,
                    'bgc_product': product,
                    'bgc_type': bgc_type
                })

    if not records:
        return pd.DataFrame()

    # Aggregate multiple BGCs per orthogroup
    df = pd.DataFrame(records)

    agg_df = df.groupby('Orthogroup').agg({
        'bgc_id': lambda x: ';'.join(sorted(set(x))),
        'bgc_product': lambda x: ';'.join(sorted(set(x))),
        'bgc_type': lambda x: ';'.join(sorted(set(x)))
    }).reset_index()

    agg_df['n_bgcs'] = agg_df['bgc_id'].str.split(';').str.len()

    return agg_df


def build_combined_annotation_table(
    combined_og_genecount: pd.DataFrame,
    combined_og_long: pd.DataFrame,
    species_base_dirs: Dict[str, str],
    output_path: str = None
) -> pd.DataFrame:
    """
    Build comprehensive annotation table for combined pangenome.

    Uses PROTEIN IDs as the linking key (not orthogroup IDs, which differ
    between individual species runs and the combined run).

    Integrates:
    - Functional annotations (EggNOG, dbCAN, InterProScan, SignalP) via protein IDs
    - BGC/GCF mappings from antiSMASH/BiG-SCAPE
    - Pangenome statistics (presence counts, species coverage)

    Parameters
    ----------
    combined_og_genecount : pd.DataFrame
        Combined pangenome gene count matrix
    combined_og_long : pd.DataFrame
        Combined pangenome long format (Orthogroup, Species, Protein_ID)
    species_base_dirs : Dict[str, str]
        Species directory mapping
    output_path : str, optional
        Path to save output

    Returns
    -------
    pd.DataFrame
        Comprehensive annotation table
    """
    print("\n" + "=" * 60)
    print("BUILDING COMBINED ANNOTATION TABLE")
    print("(Using protein IDs as linking key)")
    print("=" * 60)

    # Step 1-2: Transfer annotations via protein IDs
    print("\n1-2. Transferring annotations via protein IDs...")
    combined_annot = transfer_annotations_to_combined_via_proteins(
        combined_og_long=combined_og_long,
        species_base_dirs=species_base_dirs,
        output_path=None  # Don't save yet - we'll add more columns
    )

    if combined_annot.empty:
        # Fallback: create minimal table
        combined_ogs = combined_og_genecount['Orthogroup'].tolist() if 'Orthogroup' in combined_og_genecount.columns else combined_og_genecount.index.tolist()
        combined_annot = pd.DataFrame({'Orthogroup': combined_ogs})

    # Step 3: Load and map BGC/GCF data
    print("\n3. Loading BGC/GCF data...")
    bgc_data = load_species_bgc_gcf_data(species_base_dirs)

    if bgc_data:
        print("\n4. Mapping BGCs to orthogroups...")
        all_bgc_mappings = []

        for species, data in bgc_data.items():
            if 'bgc' in data:
                # Filter combined_og_long to this species
                species_long = combined_og_long[combined_og_long['Species'] == species]

                bgc_og_map = map_bgc_to_orthogroups(
                    data['bgc'], species_long
                )

                if not bgc_og_map.empty:
                    bgc_og_map['species'] = species
                    all_bgc_mappings.append(bgc_og_map)

        if all_bgc_mappings:
            combined_bgc = pd.concat(all_bgc_mappings, ignore_index=True)
            og_bgc = create_og_to_bgc_mapping(combined_bgc)

            # Merge with annotation table
            combined_annot = combined_annot.merge(
                og_bgc[['Orthogroup', 'n_bgcs', 'bgc_type', 'bgc_product']],
                on='Orthogroup', how='left'
            )
            combined_annot['n_bgcs'] = combined_annot['n_bgcs'].fillna(0).astype(int)

            print(f"  {(combined_annot['n_bgcs'] > 0).sum()} orthogroups mapped to BGCs")

    # Step 4: Add pangenome statistics
    print("\n5. Adding pangenome statistics...")

    # Calculate presence across genomes
    pav_cols = [c for c in combined_og_genecount.columns if c != 'Orthogroup' and '__' in c]
    if pav_cols:
        og_col = 'Orthogroup'
        if og_col in combined_og_genecount.columns:
            presence_counts = (combined_og_genecount[pav_cols] > 0).sum(axis=1)
            og_presence = pd.DataFrame({
                'Orthogroup': combined_og_genecount[og_col],
                'n_genomes_present': presence_counts.values,
                'presence_fraction': (presence_counts / len(pav_cols)).values
            })

            combined_annot = combined_annot.merge(og_presence, on='Orthogroup', how='left')

    # Count species representation from long format
    species_counts = combined_og_long.groupby('Orthogroup')['Species'].nunique().reset_index()
    species_counts.columns = ['Orthogroup', 'n_species_present']

    if 'n_species_present' not in combined_annot.columns:
        combined_annot = combined_annot.merge(species_counts, on='Orthogroup', how='left')

    # Save if path provided
    if output_path:
        combined_annot.to_csv(output_path, sep='\t', index=False)
        print(f"\nSaved combined annotation table to: {output_path}")

    print(f"\nFinal annotation table: {len(combined_annot)} orthogroups")

    return combined_annot


def get_annotation_summary(combined_annot: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for combined annotation table.

    Parameters
    ----------
    combined_annot : pd.DataFrame
        Combined annotation table

    Returns
    -------
    pd.DataFrame
        Summary statistics
    """
    summary_rows = []

    # Overall counts
    n_total = len(combined_annot)

    # Annotation coverage per column
    annotation_cols = [
        'Description', 'GOs', 'PFAMs', 'CAZy', 'EC',
        'COG_category', 'KEGG_ko', 'KEGG_Pathway',
        'interpro_IPR', 'dbcan_Substrate'
    ]

    for col in annotation_cols:
        if col in combined_annot.columns:
            n_annotated = (combined_annot[col].fillna('').astype(str) != '').sum()
            pct = 100 * n_annotated / n_total
            summary_rows.append({
                'Category': col,
                'N_Annotated': n_annotated,
                'Percentage': f"{pct:.1f}%"
            })

    # BGC coverage
    if 'n_bgcs' in combined_annot.columns:
        n_bgc = (combined_annot['n_bgcs'] > 0).sum()
        pct = 100 * n_bgc / n_total
        summary_rows.append({
            'Category': 'BGC_associated',
            'N_Annotated': n_bgc,
            'Percentage': f"{pct:.1f}%"
        })

    # Species coverage
    if 'n_species_present' in combined_annot.columns:
        for n in range(1, 5):
            n_multi = (combined_annot['n_species_present'] >= n).sum()
            summary_rows.append({
                'Category': f'Present_in_{n}+_species',
                'N_Annotated': n_multi,
                'Percentage': f"{100*n_multi/n_total:.1f}%"
            })

    return pd.DataFrame(summary_rows)

# =============================================================================
# KINSHIP / GENOMIC RELATIONSHIP MATRIX (GRM)
# =============================================================================

def calculate_kinship_matrix(
    pav: pd.DataFrame,
    method: str = 'realized',
    center: bool = True
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Calculate kinship/genomic relationship matrix from PAV data.

    Parameters
    ----------
    pav : pd.DataFrame
        Presence/absence matrix (binary 0/1) with samples as columns, OGs as rows
    method : str
        Kinship method: 'realized' (VanRaden Method 1), 'ibs' (identity by state)
    center : bool
        Whether to center the marker matrix (subtract mean frequency)

    Returns
    -------
    Tuple[pd.DataFrame, np.ndarray]
        Kinship matrix as DataFrame and raw numpy array
    """
    # Get sample columns (exclude Orthogroup column)
    sample_cols = [c for c in pav.columns if c != 'Orthogroup']

    if not sample_cols:
        raise ValueError("No sample columns found in PAV matrix")

    # Create binary marker matrix (samples x markers)
    # Transpose so samples are rows, markers are columns
    X = pav[sample_cols].T.values.astype(float)
    n_samples, n_markers = X.shape

    print(f"Calculating kinship matrix: {n_samples} samples x {n_markers} markers")

    if method == 'realized':
        # VanRaden Method 1 (centered, scaled)
        if center:
            # Calculate allele frequencies
            p = X.mean(axis=0)
            # Center the matrix
            Z = X - p
            # Scaling factor
            scale = 2 * np.sum(p * (1 - p))
            if scale > 0:
                K = Z @ Z.T / scale
            else:
                K = Z @ Z.T / n_markers
        else:
            K = X @ X.T / n_markers

    elif method == 'ibs':
        # Identity by state - proportion of shared markers
        K = np.zeros((n_samples, n_samples))
        for i in range(n_samples):
            for j in range(i, n_samples):
                # Count matching alleles
                matches = np.sum(X[i] == X[j])
                K[i, j] = matches / n_markers
                K[j, i] = K[i, j]

    else:
        raise ValueError(f"Unknown method: {method}. Use 'realized' or 'ibs'")

    # Create DataFrame with sample names
    kinship_df = pd.DataFrame(K, index=sample_cols, columns=sample_cols)

    return kinship_df, K


def plot_kinship_heatmap(
    kinship_df: pd.DataFrame,
    species_mapping: Dict[str, str] = None,
    isolation_mapping: Dict[str, str] = None,
    title: str = "Kinship Matrix",
    figsize: Tuple[int, int] = (12, 10),
    output_path: str = None
):
    """
    Plot kinship matrix as annotated heatmap.

    Parameters
    ----------
    kinship_df : pd.DataFrame
        Kinship matrix
    species_mapping : Dict[str, str], optional
        Sample -> species mapping for color annotation
    isolation_mapping : Dict[str, str], optional
        Sample -> isolation class mapping for color annotation
    title : str
        Plot title
    figsize : Tuple[int, int]
        Figure size
    output_path : str, optional
        Path to save figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Cluster the kinship matrix
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform

    # Convert kinship to distance
    K = kinship_df.values
    # Ensure symmetry and no negative values for distance
    K_sym = (K + K.T) / 2
    K_min = K_sym.min()
    if K_min < 0:
        K_sym = K_sym - K_min
    K_max = K_sym.max()
    if K_max > 0:
        D = 1 - (K_sym / K_max)  # Convert similarity to distance
    else:
        D = np.zeros_like(K_sym)
    np.fill_diagonal(D, 0)

    # Cluster
    try:
        condensed = squareform(D, checks=False)
        Z = linkage(condensed, method='average')
        order = leaves_list(Z)
    except Exception:
        order = list(range(len(kinship_df)))

    # Reorder
    ordered_samples = kinship_df.index[order]
    K_ordered = kinship_df.loc[ordered_samples, ordered_samples]

    # Plot heatmap
    sns.heatmap(
        K_ordered,
        ax=ax,
        cmap='RdBu_r',
        center=0,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'Kinship coefficient'}
    )

    ax.set_title(title)

    # Add species color bar if mapping provided
    if species_mapping:
        # Create color annotations
        species_list = [species_mapping.get(s, 'Unknown') for s in ordered_samples]
        unique_species = sorted(set(species_list))
        species_palette = dict(zip(unique_species, sns.color_palette('Set1', len(unique_species))))

        # Add as colored bar on top
        for i, species in enumerate(species_list):
            ax.add_patch(plt.Rectangle(
                (i, len(ordered_samples) + 0.5), 1, 2,
                color=species_palette.get(species, 'gray'),
                transform=ax.get_xaxis_transform(),
                clip_on=False
            ))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved kinship heatmap to: {output_path}")



# =============================================================================
# PCA ANALYSIS
# =============================================================================

def run_pca_analysis(
    pav: pd.DataFrame,
    n_components: int = 10,
    scale: bool = True
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Run PCA on presence/absence variation matrix.

    Parameters
    ----------
    pav : pd.DataFrame
        PAV matrix (orthogroups x samples)
    n_components : int
        Number of principal components to compute
    scale : bool
        Whether to scale variables to unit variance

    Returns
    -------
    Tuple[pd.DataFrame, np.ndarray, np.ndarray]
        - PC scores DataFrame (samples x PCs)
        - Explained variance ratios
        - Loadings (PCs x features)
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Get sample columns
    sample_cols = [c for c in pav.columns if c != 'Orthogroup']

    if not sample_cols:
        raise ValueError("No sample columns found in PAV matrix")

    # Transpose: samples as rows, features (OGs) as columns
    X = pav[sample_cols].T.values.astype(float)

    # Handle missing values
    X = np.nan_to_num(X, nan=0)

    # Scale if requested
    if scale:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Fit PCA
    n_components = min(n_components, min(X_scaled.shape) - 1)
    pca = PCA(n_components=n_components)
    pc_scores = pca.fit_transform(X_scaled)

    # Create scores DataFrame
    pc_columns = [f'PC{i+1}' for i in range(n_components)]
    scores_df = pd.DataFrame(pc_scores, index=sample_cols, columns=pc_columns)

    explained_var = pca.explained_variance_ratio_
    loadings = pca.components_

    print(f"PCA completed: {n_components} components")
    print(f"Variance explained: PC1={explained_var[0]:.2%}, PC2={explained_var[1]:.2%}, "
          f"PC3={explained_var[2]:.2%}" if len(explained_var) >= 3 else "")
    print(f"Cumulative: {sum(explained_var):.2%}")

    return scores_df, explained_var, loadings


def plot_pca(
    scores_df: pd.DataFrame,
    explained_var: np.ndarray,
    color_mapping: Dict[str, str] = None,
    color_label: str = "Group",
    shape_mapping: Dict[str, str] = None,
    shape_label: str = "Species",
    pc_x: int = 1,
    pc_y: int = 2,
    figsize: Tuple[int, int] = (10, 8),
    title: str = "PCA of PAV",
    output_path: str = None
):
    """
    Plot PCA scatter plot with optional grouping.

    Parameters
    ----------
    scores_df : pd.DataFrame
        PC scores DataFrame from run_pca_analysis
    explained_var : np.ndarray
        Explained variance ratios
    color_mapping : Dict[str, str], optional
        Sample -> group mapping for coloring
    color_label : str
        Legend label for color groups
    shape_mapping : Dict[str, str], optional
        Sample -> group mapping for marker shapes
    shape_label : str
        Legend label for shape groups
    pc_x, pc_y : int
        Which PCs to plot (1-indexed)
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    output_path : str, optional
        Path to save figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    x_col = f'PC{pc_x}'
    y_col = f'PC{pc_y}'

    if color_mapping is not None:
        # Get unique groups and assign colors
        groups = [color_mapping.get(s, 'Unknown') for s in scores_df.index]
        unique_groups = sorted(set(groups))
        palette = dict(zip(unique_groups, sns.color_palette('Set2', len(unique_groups))))

        # Plot by group
        for group in unique_groups:
            mask = [g == group for g in groups]
            ax.scatter(
                scores_df.loc[mask, x_col],
                scores_df.loc[mask, y_col],
                c=[palette[group]],
                label=group,
                alpha=0.7,
                s=60,
                edgecolors='white',
                linewidths=0.5
            )
        ax.legend(title=color_label, bbox_to_anchor=(1.02, 1), loc='upper left')
    else:
        ax.scatter(
            scores_df[x_col],
            scores_df[y_col],
            alpha=0.7,
            s=60,
            edgecolors='white',
            linewidths=0.5
        )

    ax.set_xlabel(f'{x_col} ({explained_var[pc_x-1]:.1%} variance)')
    ax.set_ylabel(f'{y_col} ({explained_var[pc_y-1]:.1%} variance)')
    ax.set_title(title)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved PCA plot to: {output_path}")



def plot_pca_grid(
    scores_df: pd.DataFrame,
    explained_var: np.ndarray,
    species_mapping: Dict[str, str] = None,
    isolation_mapping: Dict[str, str] = None,
    n_pcs: int = 4,
    figsize: Tuple[int, int] = (14, 12),
    output_path: str = None
):
    """
    Plot grid of PCA plots colored by species and isolation class.

    Parameters
    ----------
    scores_df : pd.DataFrame
        PC scores DataFrame
    explained_var : np.ndarray
        Explained variance ratios
    species_mapping : Dict[str, str], optional
        Sample -> species mapping
    isolation_mapping : Dict[str, str], optional
        Sample -> isolation class mapping
    n_pcs : int
        Number of PCs to show (plots PC1 vs PC2, PC3, PC4, etc.)
    figsize : Tuple[int, int]
        Figure size
    output_path : str, optional
        Path to save figure
    """
    n_rows = 2  # Species row and isolation class row
    n_cols = n_pcs - 1  # PC1 vs PC2, PC1 vs PC3, etc.

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)

    # Species palette
    if species_mapping:
        species_list = [species_mapping.get(s, 'Unknown') for s in scores_df.index]
        unique_species = sorted(set(species_list))
        species_palette = dict(zip(unique_species, sns.color_palette('Set1', len(unique_species))))

    # Isolation palette
    if isolation_mapping:
        iso_list = [isolation_mapping.get(s, 'Unknown') for s in scores_df.index]
        unique_iso = sorted(set(iso_list))
        iso_palette = dict(zip(unique_iso,
                              [ISOLATION_CLASS_COLORS.get(c, '#bdc3c7') for c in unique_iso]))

    for col_idx in range(n_cols):
        pc_x = 1
        pc_y = col_idx + 2

        x_col = f'PC{pc_x}'
        y_col = f'PC{pc_y}'

        # Row 0: Species coloring
        ax = axes[0, col_idx]
        if species_mapping:
            for species in unique_species:
                mask = [s == species for s in species_list]
                ax.scatter(
                    scores_df.loc[mask, x_col],
                    scores_df.loc[mask, y_col],
                    c=[species_palette[species]],
                    label=species if col_idx == 0 else None,
                    alpha=0.7, s=40, edgecolors='white', linewidths=0.3
                )
        else:
            ax.scatter(scores_df[x_col], scores_df[y_col], alpha=0.7, s=40)

        ax.set_xlabel(f'{x_col} ({explained_var[pc_x-1]:.1%})')
        ax.set_ylabel(f'{y_col} ({explained_var[pc_y-1]:.1%})')
        if col_idx == 0:
            ax.set_title('By Species')
            if species_mapping:
                ax.legend(fontsize=8, loc='upper right')

        # Row 1: Isolation class coloring
        ax = axes[1, col_idx]
        if isolation_mapping:
            for iso in unique_iso:
                mask = [i == iso for i in iso_list]
                ax.scatter(
                    scores_df.loc[mask, x_col],
                    scores_df.loc[mask, y_col],
                    c=[iso_palette[iso]],
                    label=iso if col_idx == 0 else None,
                    alpha=0.7, s=40, edgecolors='white', linewidths=0.3
                )
        else:
            ax.scatter(scores_df[x_col], scores_df[y_col], alpha=0.7, s=40)

        ax.set_xlabel(f'{x_col} ({explained_var[pc_x-1]:.1%})')
        ax.set_ylabel(f'{y_col} ({explained_var[pc_y-1]:.1%})')
        if col_idx == 0:
            ax.set_title('By Isolation Class')
            if isolation_mapping:
                ax.legend(fontsize=7, loc='upper right')

    plt.suptitle('PCA of Combined Pangenome', fontsize=14, y=1.02)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved PCA grid to: {output_path}")



# =============================================================================
# PAV AND CNV MATRICES
# =============================================================================

def build_pav_from_genecount(genecount_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build presence/absence variation matrix from gene count data.

    Parameters
    ----------
    genecount_df : pd.DataFrame
        Gene count matrix (orthogroups x samples)

    Returns
    -------
    pd.DataFrame
        Binary PAV matrix (1=present, 0=absent)
    """
    sample_cols = [c for c in genecount_df.columns if c != 'Orthogroup']

    pav = genecount_df.copy()
    pav[sample_cols] = (pav[sample_cols] > 0).astype(int)

    return pav


def build_cnv_from_genecount(genecount_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build copy number variation matrix from gene count data.

    This is essentially the raw gene count matrix, but can include
    normalization or categorization.

    Parameters
    ----------
    genecount_df : pd.DataFrame
        Gene count matrix (orthogroups x samples)

    Returns
    -------
    pd.DataFrame
        CNV matrix (raw counts)
    """
    return genecount_df.copy()


def categorize_cnv(cnv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Categorize CNV into discrete categories: absent, single, multi-copy.

    Parameters
    ----------
    cnv_df : pd.DataFrame
        CNV matrix (raw counts)

    Returns
    -------
    pd.DataFrame
        Categorized CNV matrix (0=absent, 1=single, 2=multi-copy)
    """
    sample_cols = [c for c in cnv_df.columns if c != 'Orthogroup']

    cnv_cat = cnv_df.copy()

    for col in sample_cols:
        cnv_cat[col] = cnv_df[col].apply(
            lambda x: 0 if x == 0 else (1 if x == 1 else 2)
        )

    return cnv_cat


def get_pav_cnv_summary(
    pav: pd.DataFrame,
    cnv: pd.DataFrame,
    species_mapping: Dict[str, str] = None
) -> pd.DataFrame:
    """
    Generate summary statistics for PAV and CNV.

    Parameters
    ----------
    pav : pd.DataFrame
        PAV matrix
    cnv : pd.DataFrame
        CNV matrix
    species_mapping : Dict[str, str], optional
        Sample -> species mapping

    Returns
    -------
    pd.DataFrame
        Summary statistics
    """
    sample_cols = [c for c in pav.columns if c != 'Orthogroup']

    summary_rows = []

    # Overall statistics
    n_og = len(pav)
    n_samples = len(sample_cols)

    # Core (present in all) vs accessory vs singleton
    presence_per_og = (pav[sample_cols] > 0).sum(axis=1)
    n_core = (presence_per_og == n_samples).sum()
    n_singleton = (presence_per_og == 1).sum()
    n_accessory = n_og - n_core - n_singleton

    summary_rows.append({
        'Category': 'Total orthogroups',
        'Count': n_og,
        'Percentage': '100%'
    })
    summary_rows.append({
        'Category': 'Core (all samples)',
        'Count': n_core,
        'Percentage': f'{100*n_core/n_og:.1f}%'
    })
    summary_rows.append({
        'Category': 'Accessory',
        'Count': n_accessory,
        'Percentage': f'{100*n_accessory/n_og:.1f}%'
    })
    summary_rows.append({
        'Category': 'Singletons',
        'Count': n_singleton,
        'Percentage': f'{100*n_singleton/n_og:.1f}%'
    })

    # Multi-copy statistics from CNV
    cnv_vals = cnv[sample_cols].values
    n_multicopy = (cnv_vals > 1).any(axis=1).sum()

    summary_rows.append({
        'Category': 'With multi-copy in any sample',
        'Count': n_multicopy,
        'Percentage': f'{100*n_multicopy/n_og:.1f}%'
    })

    # Per-species statistics if mapping provided
    if species_mapping:
        species_to_samples = {}
        for sample, species in species_mapping.items():
            if species not in species_to_samples:
                species_to_samples[species] = []
            if sample in sample_cols:
                species_to_samples[species].append(sample)

        for species, samples in sorted(species_to_samples.items()):
            if samples:
                n_present = (pav[samples].sum(axis=1) > 0).sum()
                summary_rows.append({
                    'Category': f'{species}: present OGs',
                    'Count': n_present,
                    'Percentage': f'{100*n_present/n_og:.1f}%'
                })

    return pd.DataFrame(summary_rows)


# =============================================================================
# BGC AND GCF PAV/CNV MATRICES
# =============================================================================

def build_bgc_pav_matrix(
    antismash_summary_df: pd.DataFrame,
    sample_col: str = 'Assembly_Accession',
    bgc_type_col: str = 'product_class'
) -> pd.DataFrame:
    """
    Build BGC presence/absence matrix per sample.

    Parameters
    ----------
    antismash_summary_df : pd.DataFrame
        antiSMASH summary table with BGC annotations per genome
    sample_col : str
        Column containing sample/genome identifiers
    bgc_type_col : str
        Column containing BGC type/product class

    Returns
    -------
    pd.DataFrame
        BGC PAV matrix (BGC types x samples)
    """
    if sample_col not in antismash_summary_df.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in DataFrame")

    # Count BGCs per type per sample
    if bgc_type_col in antismash_summary_df.columns:
        bgc_counts = antismash_summary_df.groupby(
            [sample_col, bgc_type_col]
        ).size().unstack(fill_value=0)
    else:
        # If no type column, just count regions
        bgc_counts = antismash_summary_df.groupby(sample_col).size()
        bgc_counts = pd.DataFrame({'BGC_total': bgc_counts})

    # Transpose to have BGC types as rows, samples as columns
    pav = (bgc_counts > 0).astype(int).T
    pav.insert(0, 'BGC_Type', pav.index)
    pav = pav.reset_index(drop=True)

    return pav


def build_bgc_cnv_matrix(
    antismash_summary_df: pd.DataFrame,
    sample_col: str = 'Assembly_Accession',
    bgc_type_col: str = 'product_class'
) -> pd.DataFrame:
    """
    Build BGC copy number matrix per sample.

    Parameters
    ----------
    antismash_summary_df : pd.DataFrame
        antiSMASH summary table
    sample_col : str
        Column containing sample/genome identifiers
    bgc_type_col : str
        Column containing BGC type/product class

    Returns
    -------
    pd.DataFrame
        BGC CNV matrix (BGC types x samples, values = count)
    """
    if sample_col not in antismash_summary_df.columns:
        raise ValueError(f"Sample column '{sample_col}' not found")

    if bgc_type_col in antismash_summary_df.columns:
        bgc_counts = antismash_summary_df.groupby(
            [sample_col, bgc_type_col]
        ).size().unstack(fill_value=0)
    else:
        bgc_counts = antismash_summary_df.groupby(sample_col).size()
        bgc_counts = pd.DataFrame({'BGC_total': bgc_counts})

    # Transpose
    cnv = bgc_counts.T.copy()
    cnv.insert(0, 'BGC_Type', cnv.index)
    cnv = cnv.reset_index(drop=True)

    return cnv


def build_gcf_pav_matrix(
    gcf_membership_df: pd.DataFrame,
    sample_col: str = 'Assembly_Accession',
    gcf_col: str = 'GCF_ID'
) -> pd.DataFrame:
    """
    Build GCF (Gene Cluster Family) presence/absence matrix.

    Parameters
    ----------
    gcf_membership_df : pd.DataFrame
        GCF membership table mapping BGCs to GCFs
    sample_col : str
        Column containing sample/genome identifiers
    gcf_col : str
        Column containing GCF identifiers

    Returns
    -------
    pd.DataFrame
        GCF PAV matrix (GCFs x samples)
    """
    if sample_col not in gcf_membership_df.columns:
        raise ValueError(f"Sample column '{sample_col}' not found")
    if gcf_col not in gcf_membership_df.columns:
        raise ValueError(f"GCF column '{gcf_col}' not found")

    # Count GCF occurrences per sample
    gcf_counts = gcf_membership_df.groupby(
        [gcf_col, sample_col]
    ).size().unstack(fill_value=0)

    # Convert to binary
    pav = (gcf_counts > 0).astype(int)
    pav = pav.reset_index()
    pav = pav.rename(columns={gcf_col: 'GCF_ID'})

    return pav


def build_gcf_cnv_matrix(
    gcf_membership_df: pd.DataFrame,
    sample_col: str = 'Assembly_Accession',
    gcf_col: str = 'GCF_ID'
) -> pd.DataFrame:
    """
    Build GCF copy number matrix.

    Parameters
    ----------
    gcf_membership_df : pd.DataFrame
        GCF membership table
    sample_col : str
        Column containing sample identifiers
    gcf_col : str
        Column containing GCF identifiers

    Returns
    -------
    pd.DataFrame
        GCF CNV matrix (GCFs x samples, values = counts)
    """
    if sample_col not in gcf_membership_df.columns:
        raise ValueError(f"Sample column '{sample_col}' not found")
    if gcf_col not in gcf_membership_df.columns:
        raise ValueError(f"GCF column '{gcf_col}' not found")

    gcf_counts = gcf_membership_df.groupby(
        [gcf_col, sample_col]
    ).size().unstack(fill_value=0)

    cnv = gcf_counts.reset_index()
    cnv = cnv.rename(columns={gcf_col: 'GCF_ID'})

    return cnv


def get_bgc_gcf_summary(
    bgc_pav: pd.DataFrame = None,
    bgc_cnv: pd.DataFrame = None,
    gcf_pav: pd.DataFrame = None,
    gcf_cnv: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Generate summary statistics for BGC and GCF matrices.

    Parameters
    ----------
    bgc_pav : pd.DataFrame, optional
        BGC PAV matrix
    bgc_cnv : pd.DataFrame, optional
        BGC CNV matrix
    gcf_pav : pd.DataFrame, optional
        GCF PAV matrix
    gcf_cnv : pd.DataFrame, optional
        GCF CNV matrix

    Returns
    -------
    pd.DataFrame
        Summary statistics
    """
    summary_rows = []

    if bgc_pav is not None:
        id_col = 'BGC_Type' if 'BGC_Type' in bgc_pav.columns else bgc_pav.columns[0]
        sample_cols = [c for c in bgc_pav.columns if c != id_col]

        n_types = len(bgc_pav)
        n_samples = len(sample_cols)

        # Core BGC types (present in all samples)
        presence = bgc_pav[sample_cols].sum(axis=1)
        n_core = (presence == n_samples).sum()
        n_accessory = (presence > 0).sum() - n_core

        summary_rows.extend([
            {'Category': 'BGC: Total types', 'Count': n_types, 'Percentage': '-'},
            {'Category': 'BGC: Core types', 'Count': n_core, 'Percentage': f'{100*n_core/max(1,n_types):.1f}%'},
            {'Category': 'BGC: Accessory types', 'Count': n_accessory, 'Percentage': f'{100*n_accessory/max(1,n_types):.1f}%'}
        ])

    if gcf_pav is not None:
        id_col = 'GCF_ID' if 'GCF_ID' in gcf_pav.columns else gcf_pav.columns[0]
        sample_cols = [c for c in gcf_pav.columns if c != id_col]

        n_gcf = len(gcf_pav)
        n_samples = len(sample_cols)

        presence = gcf_pav[sample_cols].sum(axis=1)
        n_core = (presence == n_samples).sum()
        n_accessory = (presence > 0).sum() - n_core

        summary_rows.extend([
            {'Category': 'GCF: Total families', 'Count': n_gcf, 'Percentage': '-'},
            {'Category': 'GCF: Core families', 'Count': n_core, 'Percentage': f'{100*n_core/max(1,n_gcf):.1f}%'},
            {'Category': 'GCF: Accessory families', 'Count': n_accessory, 'Percentage': f'{100*n_accessory/max(1,n_gcf):.1f}%'}
        ])

    if bgc_cnv is not None:
        id_col = 'BGC_Type' if 'BGC_Type' in bgc_cnv.columns else bgc_cnv.columns[0]
        sample_cols = [c for c in bgc_cnv.columns if c != id_col]

        avg_per_sample = bgc_cnv[sample_cols].sum().mean()
        max_per_sample = bgc_cnv[sample_cols].sum().max()

        summary_rows.extend([
            {'Category': 'BGC: Avg per genome', 'Count': f'{avg_per_sample:.1f}', 'Percentage': '-'},
            {'Category': 'BGC: Max per genome', 'Count': max_per_sample, 'Percentage': '-'}
        ])

    if gcf_cnv is not None:
        id_col = 'GCF_ID' if 'GCF_ID' in gcf_cnv.columns else gcf_cnv.columns[0]
        sample_cols = [c for c in gcf_cnv.columns if c != id_col]

        avg_per_sample = gcf_cnv[sample_cols].sum().mean()
        max_per_sample = gcf_cnv[sample_cols].sum().max()

        summary_rows.extend([
            {'Category': 'GCF: Avg per genome', 'Count': f'{avg_per_sample:.1f}', 'Percentage': '-'},
            {'Category': 'GCF: Max per genome', 'Count': max_per_sample, 'Percentage': '-'}
        ])

    return pd.DataFrame(summary_rows)


# =============================================================================
# COMBINED POPULATION STRUCTURE ANALYSIS
# =============================================================================

def run_population_structure_analysis(
    genecount_df: pd.DataFrame,
    species_mapping: Dict[str, str] = None,
    isolation_mapping: Dict[str, str] = None,
    output_dir: str = None
) -> Dict[str, Any]:
    """
    Run complete population structure analysis: PAV, CNV, kinship, PCA.

    Parameters
    ----------
    genecount_df : pd.DataFrame
        Gene count matrix from OrthoFinder
    species_mapping : Dict[str, str], optional
        Sample -> species mapping
    isolation_mapping : Dict[str, str], optional
        Sample -> isolation class mapping
    output_dir : str, optional
        Directory to save outputs

    Returns
    -------
    Dict[str, Any]
        Dictionary containing all results:
        - pav: PAV DataFrame
        - cnv: CNV DataFrame
        - kinship: Kinship DataFrame
        - pca_scores: PCA scores DataFrame
        - pca_variance: Explained variance array
        - summary: Summary DataFrame
    """
    print("=" * 60)
    print("POPULATION STRUCTURE ANALYSIS")
    print("=" * 60)

    results = {}

    # Build PAV matrix
    print("\n1. Building PAV matrix...")
    pav = build_pav_from_genecount(genecount_df)
    results['pav'] = pav
    print(f"   PAV matrix: {len(pav)} orthogroups x {len(pav.columns)-1} samples")

    # Build CNV matrix
    print("\n2. Building CNV matrix...")
    cnv = build_cnv_from_genecount(genecount_df)
    results['cnv'] = cnv

    # Calculate kinship
    print("\n3. Calculating kinship matrix...")
    try:
        kinship_df, K = calculate_kinship_matrix(pav, method='realized')
        results['kinship'] = kinship_df
    except Exception as e:
        print(f"   [WARN] Kinship calculation failed: {e}")
        results['kinship'] = None

    # Run PCA
    print("\n4. Running PCA...")
    try:
        scores_df, explained_var, loadings = run_pca_analysis(pav, n_components=10)
        results['pca_scores'] = scores_df
        results['pca_variance'] = explained_var
        results['pca_loadings'] = loadings
    except Exception as e:
        print(f"   [WARN] PCA failed: {e}")
        results['pca_scores'] = None
        results['pca_variance'] = None

    # Generate summary
    print("\n5. Generating summary...")
    summary = get_pav_cnv_summary(pav, cnv, species_mapping)
    results['summary'] = summary
    print(summary.to_string(index=False))

    # Save outputs if directory provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        pav.to_csv(os.path.join(output_dir, 'pav_matrix.tsv'), sep='\t', index=False)
        cnv.to_csv(os.path.join(output_dir, 'cnv_matrix.tsv'), sep='\t', index=False)

        if results['kinship'] is not None:
            results['kinship'].to_csv(os.path.join(output_dir, 'kinship_matrix.tsv'), sep='\t')

        if results['pca_scores'] is not None:
            results['pca_scores'].to_csv(os.path.join(output_dir, 'pca_scores.tsv'), sep='\t')

            # Save PCA plots
            if species_mapping or isolation_mapping:
                plot_pca_grid(
                    results['pca_scores'],
                    results['pca_variance'],
                    species_mapping=species_mapping,
                    isolation_mapping=isolation_mapping,
                    output_path=os.path.join(output_dir, 'pca_grid.png')
                )

        if results['kinship'] is not None:
            plot_kinship_heatmap(
                results['kinship'],
                species_mapping=species_mapping,
                output_path=os.path.join(output_dir, 'kinship_heatmap.png')
            )

        summary.to_csv(os.path.join(output_dir, 'population_structure_summary.tsv'), sep='\t', index=False)

        print(f"\nOutputs saved to: {output_dir}")

    return results


# if __name__ == "__main__":
    # Example usage
#     print("Convergence functions loaded successfully!")
#     print("\nKey functions:")
#     print("  - load_combined_orthogroups_genecount(): Load OrthoFinder gene counts")
#     print("  - identify_convergent_hits(): Find cross-species convergent markers")
#     print("  - run_meta_analysis(): Combine p-values across species")
#     print("  - run_leave_one_out_analysis(): Test robustness")
#     print("  - build_combined_annotation_table(): Create unified annotation table")
#     print("  - run_full_cross_species_analysis(): Complete pipeline")
#     print("\nNew population structure functions:")
#     print("  - calculate_kinship_matrix(): Compute GRM from PAV")
#     print("  - run_pca_analysis(): PCA on PAV data")
#     print("  - build_pav_from_genecount(): Build PAV matrix")
#     print("  - build_cnv_from_genecount(): Build CNV matrix")
#     print("  - build_bgc_pav_matrix(): BGC presence/absence")
#     print("  - build_gcf_pav_matrix(): GCF presence/absence")
#     print("  - run_population_structure_analysis(): Complete analysis")


# --- Functions extracted from NB3_Convergence notebook ---

# --- Transfer pan-GWAS significant hits to combined OG space ---
# For each species, map the significant species-OGs to combined OGs.

def transfer_gwas_hits_to_combined(species_gwas, species_og_to_combined, layer_pattern='pav_assoc'):
    """Map per-species GWAS significant OGs into combined OG space.
    
    Returns dict: {species: DataFrame with columns [species_og, combined_og, pvalue, qvalue, beta/effect, ...]}
    """
    mapped_results = {}

    for sp in SPECIES_LIST:
        if sp not in species_gwas or sp not in species_og_to_combined:
            continue

        og_map = species_og_to_combined[sp]
        sp_data = species_gwas[sp]

        # Find relevant result files matching the layer pattern
        matching_keys = [k for k in sp_data.keys() if layer_pattern in k]

        for key in matching_keys:
            res_df = sp_data[key].copy()
            # Use the 'feature' column (actual orthogroup ID) if present,
            # else fall back to the index. Cast to str so OG IDs like
            # 'OG0000123' match the combined-OG keys.
            if 'feature' in res_df.columns:
                res_df['species_og'] = res_df['feature'].astype(str)
            else:
                res_df['species_og'] = res_df.index.astype(str)

            # Map to combined OGs
            records = []
            for _, row in res_df.iterrows():
                sp_og = row['species_og']
                combined_ogs = og_map.get(sp_og, set())
                for cog in combined_ogs:
                    rec = row.to_dict()
                    rec['combined_og'] = cog
                    records.append(rec)

            if records:
                mapped_df = pd.DataFrame(records)
                result_key = f'{sp}__{key}'
                mapped_results[result_key] = mapped_df
                n_sig = int(mapped_df['significant'].sum()) if 'significant' in mapped_df.columns else 0
                print(f'  {result_key}: {len(mapped_df)} mapped rows, {n_sig} significant')

    return mapped_results

# Transfer PAV results
# print('Transferring PAV GWAS hits to combined OG space:')
# mapped_pav = transfer_gwas_hits_to_combined(species_gwas, species_og_to_combined, 'pav_assoc')

# Transfer CNV results
# print('\nTransferring CNV GWAS hits to combined OG space:')
# mapped_cnv = transfer_gwas_hits_to_combined(species_gwas, species_og_to_combined, 'cnv_assoc')

# Transfer functional results
# print('\nTransferring functional GWAS hits to combined OG space:')
# mapped_func = transfer_gwas_hits_to_combined(species_gwas, species_og_to_combined, 'functional_')

# --- Helper: find convergent hits between two species in mapped results ---

def find_convergent_hits(sp1, sp2, mapped_results, contrast_pattern, fdr_col='qvalue', fdr_thresh=0.1):
    """Find orthogroups significant in BOTH species with consistent direction.
    
    Parameters
    ----------
    sp1, sp2 : str
        Species names
    mapped_results : dict
        Mapped GWAS results from transfer_gwas_hits_to_combined()
    contrast_pattern : str
        Pattern to match the contrast (e.g. 'human_pathogenic_vs_rest')
    
    Returns
    -------
    pd.DataFrame
        Convergence table with columns for both species' statistics
    """
    # Find keys for each species matching the contrast
    sp1_keys = [k for k in mapped_results if k.startswith(f'{sp1}__') and contrast_pattern in k]
    sp2_keys = [k for k in mapped_results if k.startswith(f'{sp2}__') and contrast_pattern in k]

    if not sp1_keys or not sp2_keys:
        print(f'  No matching results for {sp1} or {sp2} with pattern "{contrast_pattern}"')
        return pd.DataFrame()

    # Combine results per species (may be multiple layers)
    sp1_df = pd.concat([mapped_results[k] for k in sp1_keys], ignore_index=True)
    sp2_df = pd.concat([mapped_results[k] for k in sp2_keys], ignore_index=True)

    # Get significant combined OGs for each species. The fallback only
    # accepts a column that exists in BOTH dataframes — otherwise sp2 would
    # silently return empty when the column found in sp1 is absent from sp2.
    if fdr_col not in sp1_df.columns or fdr_col not in sp2_df.columns:
        for alt in ['qvalue', 'q_value', 'fdr', 'FDR']:
            if alt in sp1_df.columns and alt in sp2_df.columns:
                fdr_col = alt
                break

    sp1_sig = sp1_df[sp1_df[fdr_col] < fdr_thresh] if fdr_col in sp1_df.columns else pd.DataFrame()
    sp2_sig = sp2_df[sp2_df[fdr_col] < fdr_thresh] if fdr_col in sp2_df.columns else pd.DataFrame()

    sp1_sig_ogs = set(sp1_sig['combined_og'].unique()) if not sp1_sig.empty else set()
    sp2_sig_ogs = set(sp2_sig['combined_og'].unique()) if not sp2_sig.empty else set()

    print(f'  {sp1}: {len(sp1_sig_ogs)} significant combined OGs')
    print(f'  {sp2}: {len(sp2_sig_ogs)} significant combined OGs')

    shared_sig = sp1_sig_ogs & sp2_sig_ogs
    print(f'  Shared significant OGs: {len(shared_sig)}')

    if not shared_sig:
        return pd.DataFrame()

    # Check direction consistency
    # Determine effect direction from beta or odds_ratio column
    beta_col = None
    for candidate in ['beta', 'effect', 'log_odds', 'odds_ratio']:
        if candidate in sp1_df.columns and candidate in sp2_df.columns:
            beta_col = candidate
            break

    convergent_records = []
    for og in shared_sig:
        r1 = sp1_sig[sp1_sig['combined_og'] == og].iloc[0]
        r2 = sp2_sig[sp2_sig['combined_og'] == og].iloc[0]

        if beta_col:
            b1 = float(r1[beta_col])
            b2 = float(r2[beta_col])
            same_dir = (b1 > 0) == (b2 > 0)
            category = 'CONVERGENT' if same_dir else 'DISCORDANT'
        else:
            b1, b2 = np.nan, np.nan
            category = 'UNKNOWN_DIRECTION'

        convergent_records.append({
            'combined_og': og,
            f'{sp1}_species_og': r1.get('species_og', ''),
            f'{sp2}_species_og': r2.get('species_og', ''),
            f'{sp1}_pvalue': r1.get('pvalue', np.nan),
            f'{sp2}_pvalue': r2.get('pvalue', np.nan),
            f'{sp1}_qvalue': r1.get(fdr_col, np.nan),
            f'{sp2}_qvalue': r2.get(fdr_col, np.nan),
            f'{sp1}_beta': b1,
            f'{sp2}_beta': b2,
            'category': category,
        })

    conv_df = pd.DataFrame(convergent_records)
    n_conv = (conv_df['category'] == 'CONVERGENT').sum()
    n_disc = (conv_df['category'] == 'DISCORDANT').sum()
    print(f'  CONVERGENT (same direction): {n_conv}')
    print(f'  DISCORDANT (opposite direction): {n_disc}')

    return conv_df


# --- Pathogenic convergence: fumigatus vs flavus ---
# print('=' * 80)
# print('PATHOGENIC CONVERGENCE: A. fumigatus vs A. flavus')
# print('  Contrast: human_pathogenic_vs_rest')
# print('=' * 80)

# PAV convergence
# print('\n--- PAV OG Convergence ---')
# path_pav_conv = find_convergent_hits('fumigatus', 'flavus', mapped_pav, 'human_pathogenic')

# CNV convergence
# print('\n--- CNV OG Convergence ---')
# path_cnv_conv = find_convergent_hits('fumigatus', 'flavus', mapped_cnv, 'human_pathogenic')

# --- Functional-level convergence for pathogenic species ---
# Compare functional GWAS results (COG, Pfam, CAZy, KEGG, GO, InterPro) across fumigatus and flavus.

def cross_species_functional_comparison(sp1, sp2, species_gwas, contrast_pattern,
                                         fdr_thresh=0.1):
    """Compare functional GWAS results between two species.
    
    Finds functional terms tested in both species, classifies shared significant ones
    as CONVERGENT (same direction) or DISCORDANT (opposite direction).
    """
    sp1_results = species_gwas.get(sp1, {})
    sp2_results = species_gwas.get(sp2, {})

    # Find functional result keys matching the contrast
    sp1_func_keys = [k for k in sp1_results if 'functional_' in k and contrast_pattern in k]
    sp2_func_keys = [k for k in sp2_results if 'functional_' in k and contrast_pattern in k]

    # Group by functional layer
    def extract_layer(key):
        # e.g. 'functional_cog_assoc_human_pathogenic_vs_rest' -> 'cog'
        parts = key.split('functional_')[1].split('_assoc')[0] if 'functional_' in key else key
        # Handle PAV variants
        parts = parts.replace('_pav', '').replace('_cnv', '')
        return parts

    all_records = []

    for sp1_key in sp1_func_keys:
        layer = extract_layer(sp1_key)
        # Find matching sp2 key for same layer and contrast
        sp2_matches = [k for k in sp2_func_keys if layer in k]
        if not sp2_matches:
            continue

        sp2_key = sp2_matches[0]
        r1 = sp1_results[sp1_key]
        r2 = sp2_results[sp2_key]

        # Identify the feature/term column (index or a named column)
        r1_features = set(r1.index)
        r2_features = set(r2.index)
        shared_features = r1_features & r2_features

        # Determine FDR and beta columns
        fdr_col = 'qvalue' if 'qvalue' in r1.columns else ('fdr' if 'fdr' in r1.columns else None)
        beta_col = 'beta' if 'beta' in r1.columns else ('effect' if 'effect' in r1.columns else None)

        if fdr_col is None:
            continue

        for feat in shared_features:
            row1 = r1.loc[feat]
            row2 = r2.loc[feat]

            q1 = float(row1[fdr_col]) if not pd.isna(row1[fdr_col]) else 1.0
            q2 = float(row2[fdr_col]) if not pd.isna(row2[fdr_col]) else 1.0

            b1 = float(row1[beta_col]) if beta_col and beta_col in r1.columns else 0.0
            b2 = float(row2[beta_col]) if beta_col and beta_col in r2.columns else 0.0

            sig1 = q1 < fdr_thresh
            sig2 = q2 < fdr_thresh

            if sig1 and sig2:
                category = 'CONVERGENT' if (b1 > 0) == (b2 > 0) else 'DISCORDANT'
            elif sig1:
                category = f'{sp1}_ONLY'
            elif sig2:
                category = f'{sp2}_ONLY'
            else:
                category = 'NOT_SIGNIFICANT'

            all_records.append({
                'feature': feat,
                'layer': layer,
                f'{sp1}_qvalue': q1,
                f'{sp2}_qvalue': q2,
                f'{sp1}_beta': b1,
                f'{sp2}_beta': b2,
                'category': category,
            })

    comp_df = pd.DataFrame(all_records)
    return comp_df


# print('=' * 60)
# print('PATHOGENIC FUNCTIONAL CONVERGENCE: fumigatus vs flavus')
# print('=' * 60)

# pathogenic_func_conv = cross_species_functional_comparison(
#     'fumigatus', 'flavus', species_gwas, 'human_pathogenic')

# if not pathogenic_func_conv.empty:
#     sig_df = pathogenic_func_conv[pathogenic_func_conv['category'].isin(['CONVERGENT', 'DISCORDANT'])]
#     n_conv = (sig_df['category'] == 'CONVERGENT').sum()
#     n_disc = (sig_df['category'] == 'DISCORDANT').sum()
#     n_fum = pathogenic_func_conv['category'].str.contains('fumigatus').sum()
#     n_fla = pathogenic_func_conv['category'].str.contains('flavus').sum()
# 
#     print(f'\n  Shared significant terms: {len(sig_df)}')
#     print(f'    CONVERGENT (same direction): {n_conv}')
#     print(f'    DISCORDANT (opposite direction): {n_disc}')
#     print(f'  fumigatus-only: {n_fum}')
#     print(f'  flavus-only: {n_fla}')
# 
#     if n_conv > 0:
#         print('\n  Convergent terms by layer:')
#         for layer, grp in sig_df[sig_df['category'] == 'CONVERGENT'].groupby('layer'):
#             print(f'    {layer}: {len(grp)} terms')
# 
    # Save
#     pathogenic_func_conv.to_csv(NB3_RESULTS / 'convergence' / 'pathogenic_functional_convergence.csv',
#                                  index=False)
#     print(f'\n  Saved: pathogenic_functional_convergence.csv')
# else:
#     print('  No functional GWAS results to compare.')

# --- Meta-analysis and leave-one-out for pathogenic convergence ---
# Use funpan_convergence module functions for formal meta-analysis.
# We need to reformulate per-species results with standardised columns.

def prepare_species_results_for_meta(species_gwas, species_og_to_combined,
                                      species_list, contrast_pattern, layer_pattern='pav_assoc'):
    """Prepare per-species results in the format expected by fpc.identify_convergent_hits().
    
    Each species result DataFrame needs: Orthogroup, p_value, q_value, odds_ratio, direction
    """
    species_results = {}

    for sp in species_list:
        if sp not in species_gwas or sp not in species_og_to_combined:
            continue

        sp_data = species_gwas[sp]
        og_map = species_og_to_combined[sp]

        # Find matching keys
        matching_keys = [k for k in sp_data if layer_pattern in k and contrast_pattern in k]
        if not matching_keys:
            continue

        all_rows = []
        for key in matching_keys:
            res = sp_data[key]
            for sp_og in res.index:
                combined_ogs = og_map.get(sp_og, set())
                if not combined_ogs:
                    continue

                row = res.loc[sp_og]
                p_val = float(row.get('pvalue', np.nan))
                q_val = float(row.get('qvalue', row.get('q_value', row.get('fdr', np.nan))))
                beta = float(row.get('beta', row.get('effect', row.get('log_odds', 0))))
                or_val = float(row.get('odds_ratio', np.exp(beta) if not np.isnan(beta) else np.nan))

                direction = 'enriched_in_case' if beta > 0 else 'depleted_in_case'

                for cog in combined_ogs:
                    all_rows.append({
                        'Orthogroup': cog,
                        'p_value': p_val,
                        'q_value': q_val,
                        'odds_ratio': or_val,
                        'direction': direction,
                    })

        if all_rows:
            df = pd.DataFrame(all_rows)
            # Keep best p-value per orthogroup
            df = df.sort_values('p_value').drop_duplicates(subset='Orthogroup', keep='first')
            species_results[sp] = df
            n_sig = (df['q_value'] < FDR_THRESHOLD).sum()
            print(f'  {sp}: {len(df)} orthogroups, {n_sig} significant')

    return species_results


# print('Preparing pathogenic species results for meta-analysis...')
# pathogenic_species = ['fumigatus', 'flavus']
# path_meta_input = prepare_species_results_for_meta(
#     species_gwas, species_og_to_combined,
#     pathogenic_species, 'human_pathogenic', 'pav_assoc')

# Identify convergent hits
# print('\nIdentifying convergent hits (PAV level)...')
# path_convergent = fpc.identify_convergent_hits(path_meta_input, fdr_threshold=FDR_THRESHOLD)
# print(f'  Convergent hits found: {len(path_convergent)}')

# if not path_convergent.empty:
#     display(path_convergent.head(10))

# Meta-analysis (Fisher method)
# print('\nRunning Fisher meta-analysis...')
# path_meta = fpc.run_meta_analysis(path_meta_input, method='fisher')
# n_meta_sig = (path_meta['meta_q_value'] < FDR_THRESHOLD).sum() if not path_meta.empty else 0
# print(f'  Meta-analysis significant (FDR < {FDR_THRESHOLD}): {n_meta_sig}')

# Leave-one-out robustness
# print('\nRunning leave-one-out analysis...')
# if not path_convergent.empty:
#     path_robust = fpc.run_leave_one_out_analysis(path_meta_input, path_convergent)
#     path_robust = fpc.get_tier_classification(path_robust)
#     if not path_robust.empty:
#         print(path_robust['tier'].value_counts().to_string())
# else:
#     path_robust = pd.DataFrame()
#     print('  No convergent hits to test for robustness.')

# --- Multi-layer convergence summary ---

# convergence_summary = []

# Helper to count convergent hits in a convergence comparison DataFrame
def count_convergent(conv_df, layer_name, comparison_name):
    if conv_df is None or conv_df.empty:
        return {'layer': layer_name, 'comparison': comparison_name,
                'n_tested': 0, 'n_convergent': 0, 'n_discordant': 0,
                'n_sp1_only': 0, 'n_sp2_only': 0}
    
    n_tested = len(conv_df[conv_df['category'].isin(['CONVERGENT', 'DISCORDANT'])])
    n_conv = (conv_df['category'] == 'CONVERGENT').sum()
    n_disc = (conv_df['category'] == 'DISCORDANT').sum()
    n_sp1 = conv_df['category'].str.endswith('_ONLY').sum() // 2  # rough count
    n_sp2 = n_sp1  # symmetric
    
    return {'layer': layer_name, 'comparison': comparison_name,
            'n_tested': n_tested, 'n_convergent': n_conv, 'n_discordant': n_disc,
            'n_sp1_only': int(conv_df['category'].str.contains('fumigatus|niger').sum()),
            'n_sp2_only': int(conv_df['category'].str.contains('flavus|oryzae').sum())}


# --- PATHOGENIC ---
# OG PAV
# convergence_summary.append(count_convergent(path_pav_conv, 'OG_PAV', 'pathogenic'))
# OG CNV
# convergence_summary.append(count_convergent(path_cnv_conv, 'OG_CNV', 'pathogenic'))
# Functional terms
# if not pathogenic_func_conv.empty:
#     for layer in pathogenic_func_conv['layer'].unique():
#         layer_df = pathogenic_func_conv[pathogenic_func_conv['layer'] == layer]
#         convergence_summary.append(count_convergent(layer_df, f'Functional_{layer}', 'pathogenic'))

# --- INDUSTRIAL ---
# convergence_summary.append(count_convergent(ind_pav_conv, 'OG_PAV', 'industrial'))
# convergence_summary.append(count_convergent(ind_cnv_conv, 'OG_CNV', 'industrial'))
# if not industrial_func_conv.empty:
#     for layer in industrial_func_conv['layer'].unique():
#         layer_df = industrial_func_conv[industrial_func_conv['layer'] == layer]
#         convergence_summary.append(count_convergent(layer_df, f'Functional_{layer}', 'industrial'))

# --- BGC/GCF layer ---
# Test GCF convergence across species
# for comparison_name, sp1, sp2, contrast in [
#     ('pathogenic', 'fumigatus', 'flavus', 'human_pathogenic'),
#     ('industrial', 'niger', 'oryzae', 'industrial')]:
# 
#     sp1_gcf_keys = [k for k in species_gwas.get(sp1, {}) if 'gcf' in k and contrast in k]
#     sp2_gcf_keys = [k for k in species_gwas.get(sp2, {}) if 'gcf' in k and contrast in k]
# 
#     if sp1_gcf_keys and sp2_gcf_keys:
        # Direct comparison of GCF features (shared family IDs)
#         r1 = species_gwas[sp1][sp1_gcf_keys[0]]
#         r2 = species_gwas[sp2][sp2_gcf_keys[0]]
# 
#         shared = set(r1.index) & set(r2.index)
#         fdr_col = 'qvalue' if 'qvalue' in r1.columns else 'fdr'
#         beta_col = 'beta' if 'beta' in r1.columns else 'effect'
# 
#         n_conv_gcf, n_disc_gcf = 0, 0
#         for feat in shared:
#             q1 = float(r1.loc[feat].get(fdr_col, 1))
#             q2 = float(r2.loc[feat].get(fdr_col, 1))
#             if q1 < FDR_THRESHOLD and q2 < FDR_THRESHOLD:
#                 if beta_col in r1.columns and beta_col in r2.columns:
#                     b1 = float(r1.loc[feat][beta_col])
#                     b2 = float(r2.loc[feat][beta_col])
#                     if (b1 > 0) == (b2 > 0):
#                         n_conv_gcf += 1
#                     else:
#                         n_disc_gcf += 1
# 
#         convergence_summary.append({
#             'layer': 'BGC_GCF', 'comparison': comparison_name,
#             'n_tested': len(shared), 'n_convergent': n_conv_gcf,
#             'n_discordant': n_disc_gcf, 'n_sp1_only': 0, 'n_sp2_only': 0
#         })
#     else:
#         convergence_summary.append({
#             'layer': 'BGC_GCF', 'comparison': comparison_name,
#             'n_tested': 0, 'n_convergent': 0, 'n_discordant': 0,
#             'n_sp1_only': 0, 'n_sp2_only': 0
#         })


# Build summary table
# convergence_summary_df = pd.DataFrame(convergence_summary)
# convergence_summary_df = convergence_summary_df.sort_values(['comparison', 'layer'])

# print('=' * 80)
# print('CONVERGENCE SUMMARY BY ANALYTICAL LAYER')
# print('=' * 80)
# display(convergence_summary_df)

# Save
# convergence_summary_df.to_csv(NB3_RESULTS / 'convergence_summary_by_layer.csv', index=False)
# print(f'\nSaved: convergence_summary_by_layer.csv')

# Print key finding
# total_conv = convergence_summary_df['n_convergent'].sum()
# total_tested = convergence_summary_df['n_tested'].sum()
# print(f'\nTotal convergent hits across all layers: {total_conv}/{total_tested}')

# --- Permutation null model for convergence ---

def run_convergence_null_model(sp1_name, sp2_name,
                                sp1_sig_ogs, sp2_sig_ogs,
                                sp1_background, sp2_background,
                                n_permutations=10000, seed=42):
    """Permutation null model for convergence testing.
    
    Repeatedly draw random sets of the same size as the observed significant sets
    from each species' background, count the overlap.
    
    Returns
    -------
    observed : int
        Observed overlap count
    null_distribution : np.ndarray
        Null distribution of overlap counts
    p_value : float
        Empirical p-value (fraction of null >= observed)
    """
    rng = np.random.RandomState(seed)

    sp1_sig_list = list(sp1_sig_ogs)
    sp2_sig_list = list(sp2_sig_ogs)
    n1 = len(sp1_sig_list)
    n2 = len(sp2_sig_list)

    # Observed overlap
    observed = len(sp1_sig_ogs & sp2_sig_ogs)

    # Background (all OGs testable in both)
    bg1 = list(sp1_background)
    bg2 = list(sp2_background)

    null_overlaps = np.zeros(n_permutations, dtype=int)

    for i in range(n_permutations):
        perm1 = set(rng.choice(bg1, size=min(n1, len(bg1)), replace=False))
        perm2 = set(rng.choice(bg2, size=min(n2, len(bg2)), replace=False))
        null_overlaps[i] = len(perm1 & perm2)

    # Empirical p-value
    p_value = (np.sum(null_overlaps >= observed) + 1) / (n_permutations + 1)

    return observed, null_overlaps, p_value


# --- Run permutation test for pathogenic convergence ---
# print('=' * 60)
# print('PERMUTATION NULL MODEL: Pathogenic Convergence')
# print('=' * 60)

# Gather significant OGs and backgrounds in combined OG space
# path_sp1_keys = [k for k in mapped_pav if k.startswith('fumigatus__') and 'human_pathogenic' in k]
# path_sp2_keys = [k for k in mapped_pav if k.startswith('flavus__') and 'human_pathogenic' in k]

# fdr_col_detect = 'qvalue'

# if path_sp1_keys and path_sp2_keys:
#     sp1_mapped = pd.concat([mapped_pav[k] for k in path_sp1_keys], ignore_index=True)
#     sp2_mapped = pd.concat([mapped_pav[k] for k in path_sp2_keys], ignore_index=True)
# 
    # Detect FDR column
#     for alt in ['qvalue', 'q_value', 'fdr']:
#         if alt in sp1_mapped.columns:
#             fdr_col_detect = alt
#             break
# 
#     fum_sig_ogs = set(sp1_mapped[sp1_mapped[fdr_col_detect] < FDR_THRESHOLD]['combined_og'].unique())
#     fla_sig_ogs = set(sp2_mapped[sp2_mapped[fdr_col_detect] < FDR_THRESHOLD]['combined_og'].unique())
#     fum_bg = set(sp1_mapped['combined_og'].unique())
#     fla_bg = set(sp2_mapped['combined_og'].unique())
# 
#     print(f'  fumigatus: {len(fum_sig_ogs)} sig OGs from {len(fum_bg)} background')
#     print(f'  flavus: {len(fla_sig_ogs)} sig OGs from {len(fla_bg)} background')
# 
#     obs_path, null_path, pval_path = run_convergence_null_model(
#         'fumigatus', 'flavus',
#         fum_sig_ogs, fla_sig_ogs,
#         fum_bg, fla_bg,
#         n_permutations=10000
#     )
#     print(f'\n  Observed overlap: {obs_path}')
#     print(f'  Null distribution: mean={null_path.mean():.2f}, median={np.median(null_path):.1f}, '
#           f'max={null_path.max()}')
#     print(f'  Empirical p-value: {pval_path:.4f}')
# else:
#     obs_path, null_path, pval_path = 0, np.zeros(10000), 1.0
#     print('  No mapped PAV results for pathogenic convergence test.')

# --- Run permutation test for industrial convergence ---
# print(f'\n{"=" * 60}')
# print('PERMUTATION NULL MODEL: Industrial Convergence')
# print('=' * 60)

# ind_sp1_keys = [k for k in mapped_pav if k.startswith('niger__') and 'industrial' in k]
# ind_sp2_keys = [k for k in mapped_pav if k.startswith('oryzae__') and 'industrial' in k]

# if ind_sp1_keys and ind_sp2_keys:
#     sp1_mapped = pd.concat([mapped_pav[k] for k in ind_sp1_keys], ignore_index=True)
#     sp2_mapped = pd.concat([mapped_pav[k] for k in ind_sp2_keys], ignore_index=True)
# 
#     for alt in ['qvalue', 'q_value', 'fdr']:
#         if alt in sp1_mapped.columns:
#             fdr_col_detect = alt
#             break
# 
#     nig_sig_ogs = set(sp1_mapped[sp1_mapped[fdr_col_detect] < FDR_THRESHOLD]['combined_og'].unique())
#     ory_sig_ogs = set(sp2_mapped[sp2_mapped[fdr_col_detect] < FDR_THRESHOLD]['combined_og'].unique())
#     nig_bg = set(sp1_mapped['combined_og'].unique())
#     ory_bg = set(sp2_mapped['combined_og'].unique())
# 
#     print(f'  niger: {len(nig_sig_ogs)} sig OGs from {len(nig_bg)} background')
#     print(f'  oryzae: {len(ory_sig_ogs)} sig OGs from {len(ory_bg)} background')
# 
#     obs_ind, null_ind, pval_ind = run_convergence_null_model(
#         'niger', 'oryzae',
#         nig_sig_ogs, ory_sig_ogs,
#         nig_bg, ory_bg,
#         n_permutations=10000
#     )
#     print(f'\n  Observed overlap: {obs_ind}')
#     print(f'  Null distribution: mean={null_ind.mean():.2f}, median={np.median(null_ind):.1f}, '
#           f'max={null_ind.max()}')
#     print(f'  Empirical p-value: {pval_ind:.4f}')
# else:
#     obs_ind, null_ind, pval_ind = 0, np.zeros(10000), 1.0
#     print('  No mapped PAV results for industrial convergence test (expected: A. oryzae has 0 sig hits).')



# --- From NB3 cell 25 ---
# --- CAZy family cross-species comparison ---

# CAZY_RE = re.compile(r'\b(?:GH|GT|PL|CE|CBM|AA)\d+(?:_\d+)?\b')

def build_family_profile(og_consensus_dict, family_col, family_regex=None):
    """Build a per-species family count profile from OG consensus annotations."""
    species_profiles = {}

    for sp, df in og_consensus_dict.items():
        # Try alternate column names
        col_actual = None
        for alt in [family_col, f'{family_col}_dbCAN', f'{family_col}_EggNog',
                    'CAZy_dbCAN', 'CAZy', 'PFAMs_EggNog', 'PFAMs']:
            if alt in df.columns:
                col_actual = alt
                break
        if col_actual is None:
            continue

        counts = Counter()
        for val in df[col_actual].dropna():
            val_str = str(val).strip()
            if val_str in ('', '-', 'nan'):
                continue
            if family_regex:
                families = family_regex.findall(val_str)
            else:
                families = [t.strip() for t in re.split(r'[;,|]', val_str) if t.strip()]
            for fam in families:
                counts[fam] += 1
        species_profiles[sp] = counts

    all_families = sorted(set().union(*[set(c.keys()) for c in species_profiles.values()]))
    present_species = [sp for sp in SPECIES_LIST if sp in species_profiles]
    rows = [{fam: species_profiles[sp].get(fam, 0) for fam in all_families} for sp in present_species]
    return pd.DataFrame(rows, index=present_species)


# cazy_profile = build_family_profile(og_consensus, 'CAZy', CAZY_RE)
# print(f'CAZy family profile: {cazy_profile.shape}')

# if not cazy_profile.empty:
#     top_n = min(30, cazy_profile.shape[1])
#     top_families = cazy_profile.sum().nlargest(top_n).index
#     plot_df = cazy_profile[top_families]
# 
#     fig, ax = plt.subplots(figsize=(14, 4))
#     sns.heatmap(plot_df, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
#                 linewidths=0.5, cbar_kws={'label': 'OG count'})
#     ax.set_title('CAZy Family Distribution Across Species (Top 30)')
#     ax.set_ylabel('Species')
#     plt.tight_layout()
#     fig.savefig(NB3_RESULTS / 'cazy_cross_species_comparison.png', dpi=150, bbox_inches='tight')
# 
#     cazy_profile.to_csv(NB3_RESULTS / 'cazy_family_comparison.csv')
#     print('Saved: cazy_cross_species_comparison.png, cazy_family_comparison.csv')
# else:
#     print('No CAZy annotation data available.')



# --- From NB3 cell 26 ---
# --- Protease family cross-species comparison ---
# Use PFAM families known to be proteases (from funpan_utils.PROTEASE_PFAMS)

# PROTEASE_PFAM_IDS = set(fpu.PROTEASE_PFAMS.keys())

def build_protease_profile(og_consensus_dict):
    """Count protease PFAM families per species."""
    species_profiles = {}

    for sp, df in og_consensus_dict.items():
        pfam_col = None
        for alt in ['PFAMs', 'PFAMs_EggNog', 'PFAMs_Interpro']:
            if alt in df.columns:
                pfam_col = alt
                break
        if pfam_col is None:
            continue

        counts = Counter()
        for val in df[pfam_col].dropna():
            val_str = str(val).strip()
            if val_str in ('', '-', 'nan'):
                continue
            pfams = [t.strip() for t in re.split(r'[;,|]', val_str) if t.strip()]
            for pf in pfams:
                if pf in PROTEASE_PFAM_IDS:
                    counts[pf] += 1
        species_profiles[sp] = counts

    all_pfams = sorted(set().union(*[set(c.keys()) for c in species_profiles.values()]))
    present_species = [sp for sp in SPECIES_LIST if sp in species_profiles]
    rows = [{pf: species_profiles[sp].get(pf, 0) for pf in all_pfams} for sp in present_species]
    profile_df = pd.DataFrame(rows, index=present_species)

    # Add human-readable descriptions as column suffixes
    renamed = {}
    for col in profile_df.columns:
        desc = fpu.PROTEASE_PFAMS.get(col, '')
        renamed[col] = f'{col} ({desc})' if desc else col
    profile_df = profile_df.rename(columns=renamed)

    return profile_df


# protease_profile = build_protease_profile(og_consensus)
# print(f'Protease family profile: {protease_profile.shape}')

# if not protease_profile.empty:
#     fig, ax = plt.subplots(figsize=(14, 4))
#     sns.heatmap(protease_profile, annot=True, fmt='d', cmap='YlGnBu', ax=ax,
#                 linewidths=0.5, cbar_kws={'label': 'OG count'})
#     ax.set_title('Protease PFAM Family Distribution Across Species')
#     ax.set_ylabel('Species')
#     plt.xticks(rotation=45, ha='right', fontsize=8)
#     plt.tight_layout()
#     fig.savefig(NB3_RESULTS / 'protease_cross_species_comparison.png', dpi=150, bbox_inches='tight')
# 
#     protease_profile.to_csv(NB3_RESULTS / 'protease_family_comparison.csv')
#     print('Saved: protease_cross_species_comparison.png, protease_family_comparison.csv')
# else:
#     print('No protease PFAM data available.')



# --- From NB3 cell 27 ---
# --- Secretome cross-species comparison ---
# Secreted proteins: identified by SignalP prediction (SP/Sec) or secretion-related GO terms

def build_secretome_profile(og_consensus_dict):
    """Count secreted protein OGs per species using SignalP and GO annotations."""
    records = []

    for sp, df in og_consensus_dict.items():
        n_total = len(df)

        # SignalP-based secretome
        n_signalp = 0
        for col in ['signalp_Prediction', 'signalp_prediction']:
            if col in df.columns:
                n_signalp = df[col].fillna('').astype(str).str.contains('SP|Sec', case=False).sum()
                break

        # GO-based secretome (extracellular/secretion GO terms)
        n_go_secreted = 0
        secretion_go_ids = set(fpu.SECRETION_GO_TERMS.keys())
        for col in ['GOs', 'GOs_EggNog', 'GOs_Interpro']:
            if col in df.columns:
                for val in df[col].dropna():
                    gos = set(re.findall(r'GO:\d+', str(val)))
                    if gos & secretion_go_ids:
                        n_go_secreted += 1
                break

        # Combined (union of both criteria)
        n_either = 0
        for _, row in df.iterrows():
            is_sp = False
            for col in ['signalp_Prediction', 'signalp_prediction']:
                if col in df.columns:
                    val = str(row.get(col, '')).strip()
                    if 'SP' in val.upper() or 'Sec' in val:
                        is_sp = True
                    break

            is_go = False
            for col in ['GOs', 'GOs_EggNog', 'GOs_Interpro']:
                if col in df.columns:
                    gos = set(re.findall(r'GO:\d+', str(row.get(col, ''))))
                    if gos & secretion_go_ids:
                        is_go = True
                    break

            if is_sp or is_go:
                n_either += 1

        records.append({
            'Species': sp,
            'Total_OGs': n_total,
            'SignalP_secreted': n_signalp,
            'GO_secreted': n_go_secreted,
            'Either_secreted': n_either,
            'Secretome_fraction': round(n_either / n_total * 100, 1) if n_total > 0 else 0,
        })

    return pd.DataFrame(records)


# secretome_df = build_secretome_profile(og_consensus)
# print(f'Secretome comparison:')

# if not secretome_df.empty:
#     display(secretome_df)
# 
    # Bar plot
#     fig, axes = plt.subplots(1, 2, figsize=(12, 4))
# 
#     ax = axes[0]
#     x = range(len(secretome_df))
#     colors = [fpu.SPECIES_COLORS.get(sp, '#888888') for sp in secretome_df['Species']]
#     ax.bar(x, secretome_df['Either_secreted'], color=colors, edgecolor='white')
#     ax.set_xticks(x)
#     ax.set_xticklabels([f'A. {sp}' for sp in secretome_df['Species']], rotation=30, ha='right')
#     ax.set_ylabel('Number of secreted OGs')
#     ax.set_title('Secretome Size (OGs with signal peptide or secretion GO)')
# 
#     ax = axes[1]
#     ax.bar(x, secretome_df['Secretome_fraction'], color=colors, edgecolor='white')
#     ax.set_xticks(x)
#     ax.set_xticklabels([f'A. {sp}' for sp in secretome_df['Species']], rotation=30, ha='right')
#     ax.set_ylabel('Secretome fraction (%)')
#     ax.set_title('Secretome as Fraction of Total Pangenome')
# 
#     plt.tight_layout()
#     fig.savefig(NB3_RESULTS / 'secretome_cross_species_comparison.png', dpi=150, bbox_inches='tight')
# 
#     secretome_df.to_csv(NB3_RESULTS / 'secretome_comparison.csv', index=False)
#     print('Saved: secretome_cross_species_comparison.png, secretome_comparison.csv')
# else:
#     print('No secretome data available.')


# =============================================================================
# Functions extracted from NB3_Convergence notebook
# =============================================================================


def load_gwas_results_from_nb2(species_list, nb2_results_dir):
    """Load per-species pan-GWAS result tables from NB2_Results.

    Parameters
    ----------
    species_list : list of str
        Species keys, e.g. ['fumigatus', 'flavus', 'niger', 'oryzae'].
    nb2_results_dir : pathlib.Path or str
        Root directory of NB2 results (contains per-species sub-dirs).

    Returns
    -------
    dict
        ``{species: {result_stem: pd.DataFrame}}``
    """
    from pathlib import Path
    nb2_results_dir = Path(nb2_results_dir)
    species_gwas = {}

    for sp in species_list:
        sp_results_dir = nb2_results_dir / sp / 'pangwas_results'
        if not sp_results_dir.exists():
            print(f'[WARN] No pan-GWAS results directory for {sp}: {sp_results_dir}')
            continue

        sp_results = {}
        result_files = sorted(sp_results_dir.glob('*.tsv'))

        for fpath in result_files:
            fname = fpath.stem
            try:
                df = pd.read_csv(fpath, sep='\t', index_col=0)
                sp_results[fname] = df
            except Exception as e:
                print(f'  [WARN] Failed to load {fname}: {e}')

        species_gwas[sp] = sp_results
        print(f'{sp}: loaded {len(sp_results)} result files')
        for k in sorted(sp_results.keys()):
            df = sp_results[k]
            n_sig = int(df['significant'].sum()) if 'significant' in df.columns else 0
            print(f'  {k}: {len(df)} features, {n_sig} significant')

    return species_gwas


def build_species_to_combined_og_map(species_list, species_root, combined_og_long):
    """Build a mapping from per-species OG IDs to combined-pangenome OG IDs.

    Uses shared protein IDs as the linking key between species-level and
    combined OrthoFinder runs.

    Parameters
    ----------
    species_list : list of str
        Species keys.
    species_root : pathlib.Path or str
        Root directory containing per-species OrthoFinder outputs
        (e.g. ``/datadrive/Species/Aspergillus``).
    combined_og_long : pd.DataFrame
        Long-format combined OrthoFinder table with columns
        ``Protein_ID`` and ``Orthogroup``.

    Returns
    -------
    tuple
        ``(species_og_to_combined, SPECIES_OG_PATHS)`` where
        ``species_og_to_combined = {species: {species_og: set(combined_ogs)}}``
    """
    from pathlib import Path
    from collections import defaultdict
    species_root = Path(species_root)

    SPECIES_OG_PATHS = {
        sp: species_root / sp / 'orthofinder_output' / 'Orthogroups' / 'Orthogroups.tsv'
        for sp in species_list
    }

    for sp in species_list:
        if not SPECIES_OG_PATHS[sp].exists():
            alt_paths = list((species_root / sp).glob(
                'orthofinder_output/Results_*/Orthogroups/Orthogroups.tsv'))
            if alt_paths:
                SPECIES_OG_PATHS[sp] = sorted(alt_paths)[-1]
                print(f'  {sp}: using {SPECIES_OG_PATHS[sp]}')

    protein_to_combined_og = dict(zip(
        combined_og_long['Protein_ID'],
        combined_og_long['Orthogroup']
    ))
    print(f'Protein-to-combined-OG lookup: {len(protein_to_combined_og)} proteins')

    species_og_to_combined = {}

    for sp in species_list:
        og_path = SPECIES_OG_PATHS.get(sp)
        if og_path is None or not og_path.exists():
            print(f'[WARN] No species OG file for {sp}')
            continue

        print(f'\n{sp}: loading species OrthoFinder from {og_path}')
        sp_og = pd.read_csv(og_path, sep='\t')

        mapping = defaultdict(set)
        n_mapped = 0

        for _, row in sp_og.iterrows():
            sp_og_id = row['Orthogroup']
            for col in sp_og.columns:
                if col == 'Orthogroup':
                    continue
                proteins = row.get(col, '')
                if pd.isna(proteins) or proteins == '':
                    continue
                for prot in str(proteins).split(', '):
                    prot = prot.strip()
                    if prot and prot in protein_to_combined_og:
                        mapping[sp_og_id].add(protein_to_combined_og[prot])
                        n_mapped += 1

        species_og_to_combined[sp] = dict(mapping)
        print(f'  {len(mapping)} species OGs mapped to combined OGs ({n_mapped} protein links)')

    return species_og_to_combined, SPECIES_OG_PATHS


def test_industrial_convergence(species_gwas, mapped_pav, mapped_cnv,
                                species_og_to_combined, fdr_threshold=0.1,
                                nb3_results=None):
    """Run industrial convergence tests (A. niger vs A. oryzae).

    Parameters
    ----------
    species_gwas : dict
        Per-species GWAS results from :func:`load_gwas_results_from_nb2`.
    mapped_pav, mapped_cnv : dict
        Mapped GWAS results from :func:`transfer_gwas_hits_to_combined`.
    species_og_to_combined : dict
        From :func:`build_species_to_combined_og_map`.
    fdr_threshold : float
    nb3_results : pathlib.Path or None
        Output directory.

    Returns
    -------
    dict
        Keys: ind_pav_conv, ind_cnv_conv, industrial_func_conv,
        ind_convergent, ind_report.
    """
    from pathlib import Path
    import funpan_convergence as fpc

    print('=' * 80)
    print('INDUSTRIAL CONVERGENCE: A. niger vs A. oryzae')
    print('  Contrast: industrial_trait_vs_rest')
    print('=' * 80)

    # PAV convergence
    print('\n--- PAV OG Convergence ---')
    ind_pav_conv = find_convergent_hits('niger', 'oryzae', mapped_pav, 'industrial')

    # CNV convergence
    print('\n--- CNV OG Convergence ---')
    ind_cnv_conv = find_convergent_hits('niger', 'oryzae', mapped_cnv, 'industrial')

    # Functional convergence
    print('\n--- Functional Convergence ---')
    industrial_func_conv = cross_species_functional_comparison(
        'niger', 'oryzae', species_gwas, 'industrial')

    if not industrial_func_conv.empty:
        sig_df = industrial_func_conv[
            industrial_func_conv['category'].isin(['CONVERGENT', 'DISCORDANT'])]
        n_conv = (sig_df['category'] == 'CONVERGENT').sum()
        n_disc = (sig_df['category'] == 'DISCORDANT').sum()
        print(f'\n  Shared significant terms: {len(sig_df)}')
        print(f'    CONVERGENT: {n_conv}')
        print(f'    DISCORDANT: {n_disc}')

        if nb3_results is not None:
            industrial_func_conv.to_csv(
                Path(nb3_results) / 'convergence' / 'industrial_functional_convergence.csv',
                index=False)
            print('  Saved: industrial_functional_convergence.csv')
    else:
        print('  No functional GWAS results to compare '
              '(expected: A. oryzae has 0 significant hits).')

    # Meta-analysis for industrial
    print('\nPreparing industrial species results for meta-analysis...')
    industrial_species = ['niger', 'oryzae']
    ind_meta_input = prepare_species_results_for_meta(
        species_gwas, species_og_to_combined,
        industrial_species, 'industrial', 'pav_assoc')

    print('\nIdentifying convergent hits (PAV level)...')
    ind_convergent = fpc.identify_convergent_hits(ind_meta_input, fdr_threshold=fdr_threshold)
    print(f'  Convergent hits found: {len(ind_convergent)}')

    ind_meta = (fpc.run_meta_analysis(ind_meta_input, method='fisher')
                if len(ind_meta_input) >= 2 else pd.DataFrame())
    ind_report = fpc.generate_cross_species_report(
        ind_convergent if not ind_convergent.empty else pd.DataFrame(),
        ind_meta,
        pd.DataFrame(),
        'Industrial (niger vs oryzae)'
    )
    print(ind_report)

    return {
        'ind_pav_conv': ind_pav_conv,
        'ind_cnv_conv': ind_cnv_conv,
        'industrial_func_conv': industrial_func_conv,
        'ind_convergent': ind_convergent,
        'ind_report': ind_report,
    }


def build_gcf_convergence_table(species_gwas, fdr_threshold=0.1, nb3_results=None):
    """Build a cross-species GCF convergence table.

    Parameters
    ----------
    species_gwas : dict
        Per-species GWAS results.
    fdr_threshold : float
    nb3_results : pathlib.Path or None
        Output directory.

    Returns
    -------
    pd.DataFrame
        GCF convergence comparison table.
    """
    from pathlib import Path

    gcf_records = []
    for comparison_name, sp1, sp2, contrast in [
        ('pathogenic', 'fumigatus', 'flavus', 'human_pathogenic'),
        ('industrial', 'niger', 'oryzae', 'industrial'),
    ]:
        sp1_gcf_keys = [k for k in species_gwas.get(sp1, {}) if 'gcf' in k and contrast in k]
        sp2_gcf_keys = [k for k in species_gwas.get(sp2, {}) if 'gcf' in k and contrast in k]

        if not sp1_gcf_keys or not sp2_gcf_keys:
            continue

        r1 = species_gwas[sp1][sp1_gcf_keys[0]]
        r2 = species_gwas[sp2][sp2_gcf_keys[0]]
        shared = set(r1.index) & set(r2.index)
        fdr_col = 'qvalue' if 'qvalue' in r1.columns else 'fdr'
        beta_col = 'beta' if 'beta' in r1.columns else 'effect'

        for feat in shared:
            q1 = float(r1.loc[feat].get(fdr_col, 1))
            q2 = float(r2.loc[feat].get(fdr_col, 1))
            b1 = (float(r1.loc[feat].get(beta_col, 0))
                  if beta_col in r1.columns else 0)
            b2 = (float(r2.loc[feat].get(beta_col, 0))
                  if beta_col in r2.columns else 0)
            sig1 = q1 < fdr_threshold
            sig2 = q2 < fdr_threshold

            if sig1 and sig2:
                cat = 'CONVERGENT' if (b1 > 0) == (b2 > 0) else 'DISCORDANT'
            elif sig1:
                cat = f'{sp1}_ONLY'
            elif sig2:
                cat = f'{sp2}_ONLY'
            else:
                cat = 'NOT_SIGNIFICANT'

            gcf_records.append({
                'gcf_family': feat, 'comparison': comparison_name,
                f'{sp1}_qvalue': q1, f'{sp2}_qvalue': q2,
                f'{sp1}_beta': b1, f'{sp2}_beta': b2,
                'category': cat,
            })

    if gcf_records:
        gcf_conv_df = pd.DataFrame(gcf_records)
        if nb3_results is not None:
            gcf_conv_df.to_csv(
                Path(nb3_results) / 'cross_species_gcf_convergence.csv', index=False)
        print(f'Saved cross_species_gcf_convergence.csv: '
              f'{len(gcf_conv_df)} GCF families compared')
        return gcf_conv_df
    else:
        print('No GCF data available for cross-species comparison.')
        return pd.DataFrame()


def plot_convergence_null_model(null_path, obs_path, pval_path,
                                null_ind, obs_ind, pval_ind,
                                nb3_results=None):
    """Plot two-panel null-model figure for convergence permutation test.

    Parameters
    ----------
    null_path, null_ind : array-like
        Null distributions (pathogenic, industrial).
    obs_path, obs_ind : int
        Observed overlaps.
    pval_path, pval_ind : float
        P-values.
    nb3_results : pathlib.Path or None
        Output directory.
    """
    import matplotlib.pyplot as plt
    from pathlib import Path

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(null_path,
            bins=range(0, max(int(null_path.max()) + 3, 5)),
            color='steelblue', alpha=0.7, edgecolor='white',
            density=True, label='Null distribution')
    ax.axvline(obs_path, color='red', linewidth=2, linestyle='--',
               label=f'Observed = {obs_path}')
    ax.set_xlabel('Number of overlapping significant OGs')
    ax.set_ylabel('Density')
    ax.set_title(f'Pathogenic Convergence Null Model\n'
                 f'(fumigatus vs flavus, p={pval_path:.4f})')
    ax.legend()

    ax = axes[1]
    ax.hist(null_ind,
            bins=range(0, max(int(null_ind.max()) + 3, 5)),
            color='darkorange', alpha=0.7, edgecolor='white',
            density=True, label='Null distribution')
    ax.axvline(obs_ind, color='red', linewidth=2, linestyle='--',
               label=f'Observed = {obs_ind}')
    ax.set_xlabel('Number of overlapping significant OGs')
    ax.set_ylabel('Density')
    ax.set_title(f'Industrial Convergence Null Model\n'
                 f'(niger vs oryzae, p={pval_ind:.4f})')
    ax.legend()

    plt.tight_layout()
    if nb3_results is not None:
        fig.savefig(Path(nb3_results) / 'convergence_null_model.png',
                    dpi=150, bbox_inches='tight')
    print('Saved: convergence_null_model.png')


def load_og_consensus_annotations(species_list, nb1_results):
    """Load OG consensus annotation tables for each species from NB1 results.

    Parameters
    ----------
    species_list : list of str
    nb1_results : pathlib.Path or str

    Returns
    -------
    dict
        ``{species: pd.DataFrame}``
    """
    from pathlib import Path
    nb1_results = Path(nb1_results)
    og_consensus = {}

    for sp in species_list:
        for pattern in [f'{sp}_og_consensus.tsv', 'og_consensus_annotations.tsv',
                        f'{sp}_og_consensus_annotations.tsv']:
            path = nb1_results / sp / pattern
            if path.exists():
                og_consensus[sp] = pd.read_csv(path, sep='\t')
                print(f'{sp}: loaded OG consensus from {path.name} '
                      f'({og_consensus[sp].shape})')
                break
        else:
            print(f'{sp}: no OG consensus file found in {nb1_results / sp}')

    if og_consensus:
        sample_sp = list(og_consensus.keys())[0]
        print(f'\nAvailable columns ({sample_sp}): '
              f'{list(og_consensus[sample_sp].columns[:20])}')

    return og_consensus

