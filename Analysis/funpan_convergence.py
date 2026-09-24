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


# =============================================================================
# ANNOTATION AND INTERPRETATION
# =============================================================================


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


# =============================================================================
# KINSHIP / GENOMIC RELATIONSHIP MATRIX (GRM)
# =============================================================================


# =============================================================================
# PCA ANALYSIS
# =============================================================================


# =============================================================================
# PAV AND CNV MATRICES
# =============================================================================


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


# =============================================================================
# COMBINED POPULATION STRUCTURE ANALYSIS
# =============================================================================


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
    # accepts a column that exists in BOTH dataframes, otherwise sp2 would
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
        (e.g. ``<repo>/Species/Aspergillus``).
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

    # Funannotate assigns generic per-genome locus tags (FUN_000001-T1), so a
    # protein ID alone is NOT unique across genomes -- keying on it collapses
    # ~2.3M protein rows onto ~13k keys and makes the mapping arbitrary. Key on
    # (assembly accession, protein ID) instead.
    protein_to_combined_og = dict(zip(
        combined_og_long['Assembly_Accession'].astype(str)
        + '|' + combined_og_long['Protein_ID'].astype(str),
        combined_og_long['Orthogroup']
    ))
    if len(protein_to_combined_og) < 0.99 * len(combined_og_long):
        raise ValueError(
            f'Protein key is not unique: {len(combined_og_long)} rows collapsed to '
            f'{len(protein_to_combined_og)} keys. Check Assembly_Accession parsing.')
    print(f'Protein-to-combined-OG lookup: {len(protein_to_combined_og)} proteins')

    species_og_to_combined = {}

    for sp in species_list:
        og_path = SPECIES_OG_PATHS.get(sp)
        if og_path is None or not og_path.exists():
            print(f'[WARN] No species OG file for {sp}')
            continue

        print(f'\n{sp}: loading species OrthoFinder from {og_path}')
        sp_og = pd.read_csv(og_path, sep='\t')

        # Species-level column headers carry the same assembly accession as the
        # combined run, so resolve each column to its accession once up front.
        col_acc = {}
        for col in sp_og.columns:
            if col == 'Orthogroup':
                continue
            m = ACC_RE.search(str(col))
            if m:
                col_acc[col] = m.group(1)
        if not col_acc:
            print(f'[WARN] {sp}: no assembly accessions parsed from column headers')

        mapping = defaultdict(set)
        n_mapped = 0

        for _, row in sp_og.iterrows():
            sp_og_id = row['Orthogroup']
            for col, acc in col_acc.items():
                proteins = row.get(col, '')
                if pd.isna(proteins) or proteins == '':
                    continue
                for prot in str(proteins).split(', '):
                    prot = prot.strip()
                    if not prot:
                        continue
                    cog = protein_to_combined_og.get(f'{acc}|{prot}')
                    if cog is not None:
                        mapping[sp_og_id].add(cog)
                        n_mapped += 1

        species_og_to_combined[sp] = dict(mapping)
        _fan = np.array([len(v) for v in mapping.values()]) if mapping else np.array([0])
        print(f'  {len(mapping)} species OGs mapped to combined OGs ({n_mapped} protein links)')
        print(f'    combined OGs per species OG: median {np.median(_fan):.0f}, '
              f'mean {_fan.mean():.2f}, max {_fan.max()}')

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


def gcf_orthogroup_overlap(species_pair, contrast, species_root, nb1_results,
                          nb2_results, combined_og_tsv, nb3_results=None):
    """Compare two species' significant GCFs through genus-level orthogroups.

    BiG-SCAPE families are defined per species and share no identifier across
    runs, so families are resolved to the orthogroups of their constituent
    genes in the combined multi-species OrthoFinder run and compared there.

    Returns
    -------
    (dict, pandas.DataFrame)
        ``{species: {GCF: set(orthogroups)}}`` for the significant families,
        and one row per cross-species family pair with the shared orthogroup
        count.
    """
    import glob
    from pathlib import Path

    comb = pd.read_csv(combined_og_tsv, sep='\t', dtype=str)
    acc_re = re.compile(r'(GC[AF]_\d+\.\d+)')
    gene2og = {}
    for col in comb.columns[1:]:
        m = acc_re.search(col)
        if not m:
            continue
        acc = m.group(1)
        cells = comb[col].dropna()
        for og, cell in zip(comb.loc[cells.index, 'Orthogroup'], cells):
            for g in cell.split(', '):
                b = re.match(r'(FUN_\d+)', g.strip())
                if b:
                    gene2og[(acc, b.group(1))] = og

    fam_ogs, sig = {}, {}
    for sp in species_pair:
        runs = sorted((Path(species_root) / sp / 'bigscape_output' /
                       'output_files').glob('*_c0*'))
        clus = pd.concat([pd.read_csv(f, sep='\t')
                          for f in glob.glob(str(runs[-1] / '*' / '*_clustering_*.tsv'))])
        clus = clus[~clus.Record.astype(str).str.startswith('BGC')].copy()
        ex = clus.Record.str.extract(r'(GC[AF]_\d+\.\d+).*?(scaffold_\d+)\.region(\d+)')
        clus['acc'], clus['scaf'], clus['reg'] = ex[0], ex[1], ex[2]
        clus = clus.dropna(subset=['acc', 'scaf', 'reg'])

        gmap = pd.read_csv(Path(nb1_results) / sp / 'bgc_og_mapping.tsv', sep='\t',
                           usecols=['gene_id', 'accession', 'bgc_id'])
        parts = gmap.bgc_id.str.split('__', expand=True)
        gmap['scaf'] = parts[1]
        gmap['reg'] = parts[2].str.extract(r'region(\d+)')
        gmap['base'] = gmap.gene_id.str.extract(r'(FUN_\d+)')
        gmap = gmap.dropna(subset=['base', 'reg'])
        by_region = gmap.groupby(['accession', 'scaf', 'reg'])['base'].apply(list)

        assoc = pd.read_csv(Path(nb2_results) / sp / 'pangwas_results' /
                            f'gcf_pav_assoc_{contrast}.tsv', sep='\t')
        sig[sp] = assoc.loc[assoc['significant'], 'feature'].tolist()

        fam_ogs[sp] = {}
        for fam in sig[sp]:
            ogs = set()
            g = clus[clus.Family == fam]
            for acc, scaf, reg in zip(g.acc, g.scaf, g.reg):
                for b in by_region.get((acc, scaf, str(int(reg))), []):
                    og = gene2og.get((acc, b))
                    if og:
                        ogs.add(og)
            fam_ogs[sp][fam] = ogs

    sp1, sp2 = species_pair
    rows = []
    for f1, o1 in fam_ogs[sp1].items():
        for f2, o2 in fam_ogs[sp2].items():
            inter = o1 & o2
            rows.append({f'{sp1}_GCF': f1, f'{sp2}_GCF': f2,
                         f'{sp1}_n_orthogroups': len(o1),
                         f'{sp2}_n_orthogroups': len(o2),
                         'shared_orthogroups': len(inter),
                         'shared_ids': '; '.join(sorted(inter))})
    pairs = pd.DataFrame(rows)

    if nb3_results is not None and not pairs.empty:
        pairs.to_csv(Path(nb3_results) / 'gcf_orthogroup_overlap.csv', index=False)

    u1 = set().union(*fam_ogs[sp1].values()) if fam_ogs[sp1] else set()
    u2 = set().union(*fam_ogs[sp2].values()) if fam_ogs[sp2] else set()
    print(f'{sp1}: {len(sig[sp1])} significant GCFs -> {len(u1)} genus-level orthogroups')
    print(f'{sp2}: {len(sig[sp2])} significant GCFs -> {len(u2)} genus-level orthogroups')
    print(f'shared orthogroups between the two sets: {len(u1 & u2)}')
    return fam_ogs, pairs


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


# =============================================================================
# NB3 ANALYSIS AND FIGURE FUNCTIONS
# =============================================================================

def _show(obj):
    """Print a DataFrame or object; stands in for the notebook's display()."""
    try:
        print(obj.to_string())
    except AttributeError:
        print(obj)


def method3_directional_convergence(species_og_to_combined, nb2_results, nb3_results, pair_m3=['fumigatus', 'flavus'], contrast_m3='human_pathogenic_vs_rest', pangenome_size_m3={'fumigatus': 10881, 'flavus': 14450}):
    """Directional convergence across orthogroups tested in both species.

    Projects each species' pan-GWAS results into genus-level orthogroup
    coordinates, keeps the orthogroups tested in both, and classifies each by
    whether the sign of beta agrees. Writes the projection accounting,
    concordance and merge-rule sensitivity tables under ``nb3_results``.

    Returns a dict of the three DataFrames.
    """
    PAIR_M3 = pair_m3
    CONTRAST_M3 = contrast_m3
    PANGENOME_SIZE_M3 = pangenome_size_m3
    NB2_RESULTS = Path(nb2_results)
    NB3_RESULTS = Path(nb3_results)

    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from scipy import stats

    if 'species_og_to_combined' not in dir():
        raise NameError(
            "species_og_to_combined is not defined - run the setup cells and "
            "build_species_to_combined_og_map(...) above before this cell.")

    (NB3_RESULTS / 'convergence').mkdir(parents=True, exist_ok=True)

    _proj, _conc, _sens = [], [], []

    for _layer in ['pav', 'cnv']:
        _per_sp = {}
        for _sp in PAIR_M3:
            _d = pd.read_csv(
                NB2_RESULTS / _sp / 'pangwas_results' / f'{_layer}_assoc_{CONTRAST_M3}.tsv',
                sep='\t')
            _d['feature'] = _d['feature'].astype(str)
            _m = species_og_to_combined[_sp]
            _L = pd.DataFrame(
                [(str(r.feature), c, r.beta, r.pvalue, r.qvalue)
                 for r in _d.itertuples() for c in _m.get(str(r.feature), set())],
                columns=['species_og', 'genus_og', 'beta', 'p', 'q'])

            _fan = _L.groupby('species_og')['genus_og'].nunique()
            _rev = _L.groupby('genus_og')['species_og'].nunique()
            _split_extra = int((_fan - 1).clip(lower=0).sum())
            _merge_lost = int((_rev - 1).clip(lower=0).sum())
            assert len(_d) + _split_extra == len(_L)
            assert len(_L) - _merge_lost == _L['genus_og'].nunique()

            _proj.append({
                'layer': _layer.upper(), 'species': _sp,
                'og_in_pangenome': PANGENOME_SIZE_M3[_sp],
                'og_tested': len(_d),
                'og_sig_q05': int((_d.qvalue < 0.05).sum()),
                'species_og_split': int((_fan > 1).sum()),
                'links_added_by_splits': _split_extra,
                'total_links': len(_L),
                'genus_og_receiving_multiple': int((_rev > 1).sum()),
                'links_lost_to_merges': _merge_lost,
                'distinct_genus_og': _L['genus_og'].nunique(),
            })
            _per_sp[_sp] = _L

        _A, _B = _per_sp['fumigatus'], _per_sp['flavus']
        _shared = sorted(set(_A.genus_og) & set(_B.genus_og))

        def _collapse(L, rule, shared=None):
            shared = _shared if shared is None else shared
            L = L[L.genus_og.isin(shared)]
            if rule == 'most_significant':
                return L.sort_values('p').drop_duplicates('genus_og').set_index('genus_og')
            if rule == 'median_beta':
                return L.groupby('genus_og').agg(beta=('beta', 'median'),
                                                 p=('p', 'min'), q=('q', 'min'))
            if rule == 'drop_ambiguous':
                _k = L.groupby('genus_og')['species_og'].nunique()
                return L[L.genus_og.isin(_k[_k == 1].index)].set_index('genus_og')
            raise ValueError(rule)

        for _rule in ['most_significant', 'median_beta', 'drop_ambiguous']:
            _a, _b = _collapse(_A, _rule), _collapse(_B, _rule)
            _i = _a.index.intersection(_b.index)
            _j = pd.DataFrame({'b1': _a.loc[_i, 'beta'], 'b2': _b.loc[_i, 'beta'],
                               'q1': _a.loc[_i, 'q'], 'q2': _b.loc[_i, 'q']})
            _j = _j[(_j.b1 != 0) & (_j.b2 != 0)]
            _same = (_j.b1 > 0) == (_j.b2 > 0)
            _n, _k = len(_j), int(_same.sum())
            _bt = stats.binomtest(_k, _n, 0.5)
            _rho, _prho = stats.spearmanr(_j.b1, _j.b2)
            _r, _pr = stats.pearsonr(_j.b1, _j.b2)
            _rec = {'layer': _layer.upper(), 'merge_rule': _rule, 'n_shared_tested': _n,
                    'n_concordant': _k, 'n_discordant': _n - _k,
                    'pct_concordant': round(100 * _k / _n, 1),
                    'sign_test_p': _bt.pvalue,
                    'ci_low': _bt.proportion_ci().low, 'ci_high': _bt.proportion_ci().high,
                    'spearman_rho': _rho, 'spearman_p': _prho,
                    'pearson_r': _r, 'pearson_p': _pr}
            _sens.append(_rec)

            if _rule == 'most_significant':
                _rec = dict(_rec)
                _rec.update({
                    'n_concordant_up': int(((_j.b1 > 0) & (_j.b2 > 0)).sum()),
                    'n_concordant_down': int(((_j.b1 < 0) & (_j.b2 < 0)).sum()),
                    'n_sig_fumigatus_only': int(((_j.q1 < .05) & (_j.q2 >= .05)).sum()),
                    'n_sig_flavus_only': int(((_j.q2 < .05) & (_j.q1 >= .05)).sum()),
                    'n_sig_both': int(((_j.q1 < .05) & (_j.q2 < .05)).sum()),
                    'n_sig_neither': int(((_j.q1 >= .05) & (_j.q2 >= .05)).sum())})
                _conc.append(_rec)

                _out = _j.reset_index().rename(columns={
                    'genus_og': 'combined_og', 'b1': 'fumigatus_beta', 'b2': 'flavus_beta',
                    'q1': 'fumigatus_qvalue', 'q2': 'flavus_qvalue'})
                _out['category'] = np.where(
                    (_out.fumigatus_beta > 0) == (_out.flavus_beta > 0),
                    'CONVERGENT', 'DISCORDANT')
                _out.to_csv(NB3_RESULTS / 'convergence' /
                            f'pathogenic_method3_{_layer}_convergence.csv', index=False)

        _sa = set(_A[_A.q < 0.05].genus_og); _sb = set(_B[_B.q < 0.05].genus_og)
        print(f'{_layer.upper()}: sig genus OGs {len(_sa)} vs {len(_sb)}, '
              f'jointly significant {len(_sa & _sb)}; significant but untested in the '
              f'partner species: fumigatus {len(_sa - set(_B.genus_og))}, '
              f'flavus {len(_sb - set(_A.genus_og))}')

    method3_projection = pd.DataFrame(_proj)
    method3_concordance = pd.DataFrame(_conc)
    method3_sensitivity = pd.DataFrame(_sens)
    method3_projection.to_csv(NB3_RESULTS / 'method3_projection_accounting.csv', index=False)
    method3_concordance.to_csv(NB3_RESULTS / 'method3_directional_concordance.csv', index=False)
    method3_sensitivity.to_csv(NB3_RESULTS / 'method3_merge_rule_sensitivity.csv', index=False)

    print('\n--- projection accounting ---'); _show(method3_projection)
    print('--- directional concordance ---'); _show(method3_concordance)
    print('--- merge-rule sensitivity ---'); _show(method3_sensitivity)

    return {'projection': method3_projection, 'concordance': method3_concordance,
            'sensitivity': method3_sensitivity}


def run_permutation_null(mapped_pav, mapped_cnv, nb3_results, n_perms=1000,
                         fdr_threshold=0.1):
    """Permutation null for the overlap of per-species significant orthogroups.

    Draws random orthogroup sets of the observed sizes and records how often
    they overlap at least as much as the real sets, per layer and species pair.

    Returns a dict with the summary table and, per pair and layer, the null
    array, the observed overlap and the p-value.
    """
    N_PERMS = n_perms
    FDR_THRESHOLD = fdr_threshold
    NB3_RESULTS = Path(nb3_results)

    import numpy as _np


    def _per_species_sig_sets(mapped_dict, contrast_pattern, fdr_thresh, fdr_col='qvalue'):
        out = {}
        for key, df in mapped_dict.items():
            if contrast_pattern not in key:
                continue
            sp = key.split('__')[0]
            if df.empty or fdr_col not in df.columns:
                continue
            sig = set(df.loc[df[fdr_col] < fdr_thresh, 'combined_og'].dropna().astype(str))
            bg = set(df['combined_og'].dropna().astype(str))
            out.setdefault(sp, {'sig': set(), 'bg': set()})
            out[sp]['sig'].update(sig)
            out[sp]['bg'].update(bg)
        return out

    def _null_overlap(mapped_dict, sp1, sp2, contrast_pattern, fdr_thresh,
                      n_perms=1000, seed=42):
        rng_local = _np.random.default_rng(seed)
        sets = _per_species_sig_sets(mapped_dict, contrast_pattern, fdr_thresh)
        if sp1 not in sets or sp2 not in sets:
            return 0, _np.zeros(n_perms), 0, 0
        n1, n2 = len(sets[sp1]['sig']), len(sets[sp2]['sig'])
        bg1 = _np.array(sorted(sets[sp1]['bg']))
        bg2 = _np.array(sorted(sets[sp2]['bg']))
        observed = len(sets[sp1]['sig'] & sets[sp2]['sig'])
        nulls = _np.zeros(n_perms)
        if min(n1, n2, len(bg1), len(bg2)) == 0:
            return observed, nulls, n1, n2
        for i in range(n_perms):
            s1 = set(rng_local.choice(bg1, size=min(n1, len(bg1)), replace=False))
            s2 = set(rng_local.choice(bg2, size=min(n2, len(bg2)), replace=False))
            nulls[i] = len(s1 & s2)
        return observed, nulls, n1, n2

    obs_path_pav, null_path_pav, n_fum_pav, n_fla_pav = _null_overlap(
        mapped_pav, 'fumigatus', 'flavus', 'human_pathogenic_vs_rest',
        fdr_thresh=FDR_THRESHOLD, n_perms=N_PERMS, seed=42)
    obs_path_cnv, null_path_cnv, n_fum_cnv, n_fla_cnv = _null_overlap(
        mapped_cnv, 'fumigatus', 'flavus', 'human_pathogenic_vs_rest',
        fdr_thresh=FDR_THRESHOLD, n_perms=N_PERMS, seed=43)
    obs_ind_pav, null_ind_pav, n_nig_pav, n_ory_pav = _null_overlap(
        mapped_pav, 'niger', 'oryzae', 'industrial_trait_vs_rest',
        fdr_thresh=FDR_THRESHOLD, n_perms=N_PERMS, seed=44)
    obs_ind_cnv, null_ind_cnv, n_nig_cnv, n_ory_cnv = _null_overlap(
        mapped_cnv, 'niger', 'oryzae', 'industrial_trait_vs_rest',
        fdr_thresh=FDR_THRESHOLD, n_perms=N_PERMS, seed=45)

    p_path_pav = (null_path_pav >= obs_path_pav).sum() / max(1, len(null_path_pav))
    p_path_cnv = (null_path_cnv >= obs_path_cnv).sum() / max(1, len(null_path_cnv))
    p_ind_pav  = (null_ind_pav  >= obs_ind_pav ).sum() / max(1, len(null_ind_pav))
    p_ind_cnv  = (null_ind_cnv  >= obs_ind_cnv ).sum() / max(1, len(null_ind_cnv))

    null_summary = pd.DataFrame([
        {'contrast':'pathogenic (fumigatus vs flavus)', 'layer':'PAV',
         'n_sp1_sig': n_fum_pav, 'n_sp2_sig': n_fla_pav,
         'observed': obs_path_pav, 'null_mean': float(null_path_pav.mean()),
         'null_max': int(null_path_pav.max()), 'null_q95': float(_np.quantile(null_path_pav, 0.95)),
         'p_value': p_path_pav, 'n_perms': N_PERMS},
        {'contrast':'pathogenic (fumigatus vs flavus)', 'layer':'CNV',
         'n_sp1_sig': n_fum_cnv, 'n_sp2_sig': n_fla_cnv,
         'observed': obs_path_cnv, 'null_mean': float(null_path_cnv.mean()),
         'null_max': int(null_path_cnv.max()), 'null_q95': float(_np.quantile(null_path_cnv, 0.95)),
         'p_value': p_path_cnv, 'n_perms': N_PERMS},
        {'contrast':'industrial (niger vs oryzae)', 'layer':'PAV',
         'n_sp1_sig': n_nig_pav, 'n_sp2_sig': n_ory_pav,
         'observed': obs_ind_pav, 'null_mean': float(null_ind_pav.mean()),
         'null_max': int(null_ind_pav.max()), 'null_q95': float(_np.quantile(null_ind_pav, 0.95)),
         'p_value': p_ind_pav, 'n_perms': N_PERMS},
        {'contrast':'industrial (niger vs oryzae)', 'layer':'CNV',
         'n_sp1_sig': n_nig_cnv, 'n_sp2_sig': n_ory_cnv,
         'observed': obs_ind_cnv, 'null_mean': float(null_ind_cnv.mean()),
         'null_max': int(null_ind_cnv.max()), 'null_q95': float(_np.quantile(null_ind_cnv, 0.95)),
         'p_value': p_ind_cnv, 'n_perms': N_PERMS},
    ])
    null_summary.to_csv(os.path.join(str(NB3_RESULTS), 'permutation_null_model.csv'), index=False)
    _np.savez(os.path.join(str(NB3_RESULTS), 'permutation_null_arrays.npz'),
             null_path_pav=null_path_pav, null_path_cnv=null_path_cnv,
             null_ind_pav=null_ind_pav, null_ind_cnv=null_ind_cnv)

    # Backwards-compatible aliases for cells that already reference these names
    null_path = null_path_pav
    null_ind  = null_ind_pav
    obs_path  = obs_path_pav
    obs_ind   = obs_ind_pav
    pval_path = p_path_pav
    pval_ind  = p_ind_pav

    print(null_summary.to_string(index=False))

    return {'summary': null_summary,
            'path_pav': (null_path_pav, obs_path_pav, p_path_pav),
            'path_cnv': (null_path_cnv, obs_path_cnv, p_path_cnv),
            'ind_pav': (null_ind_pav, obs_ind_pav, p_ind_pav),
            'ind_cnv': (null_ind_cnv, obs_ind_cnv, p_ind_cnv)}


def plot_null_model(null_result, nb3_results, n_perms=1000):
    """Null distribution and observed overlap, one panel per pair and layer.
    """
    NB3_RESULTS = Path(nb3_results)
    N_PERMS = n_perms
    null_path_pav, obs_path_pav, p_path_pav = null_result['path_pav']
    null_path_cnv, obs_path_cnv, p_path_cnv = null_result['path_cnv']
    null_ind_pav, obs_ind_pav, p_ind_pav = null_result['ind_pav']
    null_ind_cnv, obs_ind_cnv, p_ind_cnv = null_result['ind_cnv']
    fig = None

    import matplotlib.pyplot as _plt
    import numpy as _np

    fig, axes = _plt.subplots(2, 2, figsize=(13, 8))
    panels = [
        (axes[0, 0], null_path_pav, obs_path_pav, p_path_pav,
            'a  Pathogenic (fumigatus vs flavus) - PAV layer'),
        (axes[0, 1], null_path_cnv, obs_path_cnv, p_path_cnv,
            'b  Pathogenic (fumigatus vs flavus) - CNV layer'),
        (axes[1, 0], null_ind_pav, obs_ind_pav, p_ind_pav,
            'c  Industrial (niger vs oryzae) - PAV layer'),
        (axes[1, 1], null_ind_cnv, obs_ind_cnv, p_ind_cnv,
            'd  Industrial (niger vs oryzae) - CNV layer'),
    ]
    for ax, nulls, obs, pval, title in panels:
        if nulls.max() <= 0 and obs == 0:
            ax.text(0.5, 0.5, 'No data\n(contrast inestimable)',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=10, color='#888')
            ax.set_xticks([]); ax.set_yticks([])
        else:
            ax.hist(nulls, bins=max(15, int(nulls.max() + 1)),
                    color='#cccccc', edgecolor='white', alpha=0.85)
            ax.axvline(obs, color='#d62728', linewidth=2.5,
                       label=f'observed = {int(obs)}')
            ax.axvline(_np.quantile(nulls, 0.95), color='#444',
                       linestyle='--', linewidth=1.0, label='null 95%')
            ax.set_xlabel('Cross-species OG overlap (random resampling)')
            ax.set_ylabel('Permutation count')
            ax.legend(frameon=False, fontsize=8)
        ax.set_title(title + f'\n(null mean={nulls.mean():.1f}, p={pval:.3f}, n_perm={N_PERMS})',
                     loc='left', fontsize=9, fontweight='bold')
    fig.suptitle('Figure 3 (revised). Permutation null model -- per-layer, per-species sig sizes',
                 fontsize=11, fontweight='bold', y=1.005)
    _plt.tight_layout()
    fig.savefig(os.path.join(str(NB3_RESULTS), 'convergence_null_model.png'),
                dpi=200, bbox_inches='tight')
    print('Saved: convergence_null_model.png')

    return fig


def compare_cazy_families(og_consensus, species_list, nb3_results, plot=False):
    """CAZy family counts per species, written to cazy_family_comparison.csv.

    ``plot=True`` also draws a heatmap of the top families.
    """
    SPECIES_LIST = species_list
    NB3_RESULTS = Path(nb3_results)
    fig = None

    import collections, matplotlib.pyplot as _plt

    cazy_rows = []
    for sp in SPECIES_LIST:
        if sp not in og_consensus or og_consensus[sp].empty:
            continue
        ogc = og_consensus[sp]
        pclass_col = 'Pangenome_Class' if 'Pangenome_Class' in ogc.columns else None
        cazy = ogc['CAZy'].dropna()
        for og_idx, terms in cazy.items():
            pclass = ogc.loc[og_idx, pclass_col] if pclass_col else 'Unknown'
            for term in str(terms).split(';'):
                term = term.strip()
                if term and term != 'nan':
                    cazy_rows.append({'species': sp, 'family': term,
                                      'pangenome_class': pclass})
    cazy_df = pd.DataFrame(cazy_rows)
    cazy_pivot = cazy_df.groupby(['family', 'species']).size().unstack(fill_value=0)
    cazy_pivot['total'] = cazy_pivot.sum(axis=1)
    cazy_pivot = cazy_pivot.sort_values('total', ascending=False)
    cazy_pivot.to_csv(os.path.join(str(NB3_RESULTS), 'cazy_family_comparison.csv'))
    print(f'CAZy: {len(cazy_pivot)} families across {len(SPECIES_LIST)} species')

    if plot:
        # Heatmap of top 30
        top_n = 30
        top = cazy_pivot.drop(columns='total').head(top_n)
        fig, ax = _plt.subplots(figsize=(8, max(8, top_n * 0.32)))
        import numpy as _np
        im = ax.imshow(top.values, aspect='auto', cmap='YlOrBr')
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=7)
        ax.set_xticks(range(len(top.columns))); ax.set_xticklabels(
            [f'A. {c}' for c in top.columns], rotation=45, ha='right')
        for i in range(top.shape[0]):
            for j in range(top.shape[1]):
                ax.text(j, i, int(top.iloc[i, j]), ha='center', va='center',
                        fontsize=6,
                        color='white' if top.iloc[i, j] > top.values.max() / 2 else 'black')
        _plt.colorbar(im, ax=ax, label='OGs in family')
        ax.set_title(f'CAZy family OG counts (top {top_n} by total)',
                     loc='left', fontweight='bold')
        _plt.tight_layout()
        fig.savefig(os.path.join(str(NB3_RESULTS), 'cazy_cross_species_comparison.png'),
                    dpi=200, bbox_inches='tight')
        print('Saved: cazy_cross_species_comparison.png')

    return cazy_pivot, fig


def compare_protease_families(og_consensus, species_list, nb3_results, plot=False):
    """Protease family counts per species, written to protease_family_comparison.csv.

    ``plot=True`` also draws a heatmap of the top families.
    """
    SPECIES_LIST = species_list
    NB3_RESULTS = Path(nb3_results)
    fig = None
    prot_pivot = pd.DataFrame()
    import matplotlib.pyplot as _plt

    prot_rows = []
    for sp in SPECIES_LIST:
        if sp not in og_consensus or og_consensus[sp].empty:
            continue
        ogc = og_consensus[sp]
        if 'is_protease' not in ogc.columns:
            continue
        proteases = ogc[ogc['is_protease'] == True]
        for og_idx, fam in proteases['protease_families'].dropna().items():
            pclass = ogc.loc[og_idx, 'Pangenome_Class'] if 'Pangenome_Class' in ogc.columns else 'Unknown'
            for f in str(fam).split(';'):
                f = f.strip()
                if f and f != 'nan':
                    prot_rows.append({'species': sp, 'family': f, 'pangenome_class': pclass})
        n_total = len(proteases)
        print(f'  {sp}: {n_total} protease OGs')
    prot_df = pd.DataFrame(prot_rows)
    if prot_df.empty:
        print('No protease annotations available.')
    else:
        prot_pivot = prot_df.groupby(['family', 'species']).size().unstack(fill_value=0)
        prot_pivot['total'] = prot_pivot.sum(axis=1)
        prot_pivot = prot_pivot.sort_values('total', ascending=False)
        prot_pivot.to_csv(os.path.join(str(NB3_RESULTS), 'protease_family_comparison.csv'))
        if plot:
            top_n = 25
            top = prot_pivot.drop(columns='total').head(top_n)
            fig, ax = _plt.subplots(figsize=(8, max(6, top_n * 0.3)))
            im = ax.imshow(top.values, aspect='auto', cmap='Blues')
            ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=7)
            ax.set_xticks(range(len(top.columns))); ax.set_xticklabels(
                [f'A. {c}' for c in top.columns], rotation=45, ha='right')
            for i in range(top.shape[0]):
                for j in range(top.shape[1]):
                    ax.text(j, i, int(top.iloc[i, j]), ha='center', va='center', fontsize=6,
                            color='white' if top.iloc[i, j] > top.values.max() / 2 else 'black')
            _plt.colorbar(im, ax=ax, label='Protease OGs in family')
            ax.set_title(f'Protease family OG counts (top {top_n} by total)',
                         loc='left', fontweight='bold')
            _plt.tight_layout()
            fig.savefig(os.path.join(str(NB3_RESULTS), 'protease_cross_species_comparison.png'),
                        dpi=200, bbox_inches='tight')
            print('Saved: protease_cross_species_comparison.png and protease_family_comparison.csv')

    return prot_pivot, fig


def compare_secretomes(og_consensus, species_list, nb3_results, plot=False):
    """Predicted secretome size and pangenome class split per species.

    Written to secretome_comparison.csv. ``plot=True`` also draws a bar chart.
    """
    SPECIES_LIST = species_list
    NB3_RESULTS = Path(nb3_results)
    fig = None
    import matplotlib.pyplot as _plt

    import numpy as _np
    secretome_rows = []
    for sp in SPECIES_LIST:
        if sp not in og_consensus or og_consensus[sp].empty:
            continue
        ogc = og_consensus[sp]
        n_total = len(ogc)
        n_secreted = int(ogc['is_secreted_signalp'].fillna(False).astype(bool).sum()) \
            if 'is_secreted_signalp' in ogc.columns else 0
        n_secreted_GO = int(ogc['has_secretion_GO'].fillna(False).astype(bool).sum()) \
            if 'has_secretion_GO' in ogc.columns else 0
        pclass_breakdown = {}
        if 'Pangenome_Class' in ogc.columns and 'is_secreted_signalp' in ogc.columns:
            sub = ogc[ogc['is_secreted_signalp'].fillna(False).astype(bool)]
            pclass_breakdown = sub['Pangenome_Class'].value_counts().to_dict()
        secretome_rows.append({
            'species': sp, 'n_total_OGs': n_total,
            'n_secreted_signalp': n_secreted,
            'pct_secreted_signalp': 100.0 * n_secreted / max(1, n_total),
            'n_secretion_GO': n_secreted_GO,
            'n_secreted_in_Core': pclass_breakdown.get('Core', 0),
            'n_secreted_in_Accessory': pclass_breakdown.get('Accessory', 0),
            'n_secreted_in_Rare': pclass_breakdown.get('Rare', 0),
        })
    secretome_df = pd.DataFrame(secretome_rows)
    if secretome_df.empty:
        print('No secretome annotations available; skipping.')
    else:
        secretome_df.to_csv(os.path.join(str(NB3_RESULTS), 'secretome_comparison.csv'),
                            index=False)
        if plot:
            fig, axes = _plt.subplots(1, 2, figsize=(11, 4.5))
            sp_lbls = [f'A. {s}' for s in secretome_df['species']]
            axes[0].bar(sp_lbls, secretome_df['n_secreted_signalp'],
                        color=['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4'])
            for i, v in enumerate(secretome_df['n_secreted_signalp']):
                axes[0].text(i, v + 5, str(int(v)), ha='center', fontweight='bold', fontsize=9)
            axes[0].set_title('a  Secretome size (SignalP-positive OGs)',
                              loc='left', fontweight='bold')
            axes[0].set_ylabel('Secreted OGs')
            bottom = _np.zeros(len(secretome_df))
            for cls, color in [('Core', '#4c72b0'),
                               ('Accessory', '#dd8452'),
                               ('Rare', '#c44e52')]:
                vals = secretome_df[f'n_secreted_in_{cls}'].values
                axes[1].bar(sp_lbls, vals, bottom=bottom, label=cls, color=color)
                bottom += vals
            axes[1].set_title('b  Secretome by pangenome class',
                              loc='left', fontweight='bold')
            axes[1].set_ylabel('Secreted OGs'); axes[1].legend(frameon=False)
            _plt.suptitle('Cross-species secretome comparison', fontweight='bold')
            _plt.tight_layout()
            fig.savefig(os.path.join(str(NB3_RESULTS),
                                      'secretome_cross_species_comparison.png'),
                        dpi=200, bbox_inches='tight')
            print('Saved: secretome_cross_species_comparison.png and secretome_comparison.csv')

    return secretome_df, fig


def run_functional_convergence(og_consensus, nb2_results, nb3_results, func_conv_n_perm=1000,
                               functional_layers=None):
    """Jaccard overlap of each species pair's enriched annotation terms, against a
    permutation null, per annotation layer and contrast.
    """
    FUNC_CONV_N_PERM = func_conv_n_perm
    FUNCTIONAL_LAYERS = functional_layers if functional_layers is not None else [
        'COG_category', 'PFAMs', 'CAZy', 'KEGG_ko', 'KEGG_Pathway',
        'KEGG_TC', 'GOs', 'interpro_IPR', 'EC', 'dbcan_Substrate']
    NB2_RESULTS = Path(nb2_results)
    NB3_RESULTS = Path(nb3_results)

    import glob, time
    from statsmodels.stats.multitest import multipletests as _mt


    def _precompute_og_to_terms(og_consensus, layers):
        """Build OG -> {term} per layer, dropping uninformative annotations
        (COG-S/R, Pfam DUFs, GO root nodes, EC -.-.-.-) so they are excluded from
        Jaccard set construction."""
        from funpan_utils import is_informative_term as _info
        out = {}
        for layer in layers:
            if layer not in og_consensus.columns:
                continue
            out[layer] = {}
            for _, row in og_consensus.iterrows():
                v = row.get(layer)
                if pd.isna(v):
                    continue
                terms = {t.strip() for t in str(v).split(';')
                         if t.strip() and t.strip() != 'nan'}
                terms = {t for t in terms if _info(layer, t)}
                if terms:
                    out[layer][row['Orthogroup']] = terms
        return out

    def _terms_for_og_set(og_set, og_to_terms_layer):
        out = set()
        for og in og_set:
            ts = og_to_terms_layer.get(og)
            if ts:
                out |= ts
        return out

    def _jaccard(a, b):
        if not a and not b:
            return 0.0
        union = a | b
        return len(a & b) / len(union) if union else 0.0

    def _get_sig_and_bg_ogs(species, contrast, layers=('pav', 'cnv')):
        sig, bg = set(), set()
        for layer in layers:
            fp = f'{NB2_RESULTS}/{species}/pangwas_results/{layer}_assoc_{contrast}.tsv'
            if not glob.glob(fp):
                continue
            df = pd.read_csv(fp, sep='\t')
            if 'feature' not in df.columns:
                continue
            bg.update(df['feature'].astype(str).tolist())
            if 'significant' in df.columns:
                sig.update(df.loc[df['significant'], 'feature'].astype(str).tolist())
        return sig, bg

    def _functional_convergence_pair(sp1, sp2, contrast, layers,
                                     n_perm=FUNC_CONV_N_PERM, seed=42):
        sig1, bg1 = _get_sig_and_bg_ogs(sp1, contrast)
        sig2, bg2 = _get_sig_and_bg_ogs(sp2, contrast)
        print(f'\n  {sp1} vs {sp2} | {contrast}: '
              f'{len(sig1)} vs {len(sig2)} sig OGs')
        if len(sig1) == 0 or len(sig2) == 0:
            print('    -> NOT TESTABLE (one species has 0 sig OGs)')
            return None
        o2t1 = _precompute_og_to_terms(og_consensus[sp1], layers)
        o2t2 = _precompute_og_to_terms(og_consensus[sp2], layers)
        bg1_arr = _np.array(sorted(bg1))
        bg2_arr = _np.array(sorted(bg2))
        rows = []
        for layer in layers:
            if layer not in o2t1 or layer not in o2t2:
                continue
            sig1_terms = _terms_for_og_set(sig1, o2t1[layer])
            sig2_terms = _terms_for_og_set(sig2, o2t2[layer])
            observed_J = _jaccard(sig1_terms, sig2_terms)
            observed_shared = sig1_terms & sig2_terms
            n1 = min(len(sig1), len(bg1_arr)); n2 = min(len(sig2), len(bg2_arr))
            rng_layer = _np.random.default_rng(seed + (hash(layer) % 10000))
            null_J = _np.zeros(n_perm)
            for k in range(n_perm):
                s1_rand = set(rng_layer.choice(bg1_arr, size=n1, replace=False))
                s2_rand = set(rng_layer.choice(bg2_arr, size=n2, replace=False))
                t1 = _terms_for_og_set(s1_rand, o2t1[layer])
                t2 = _terms_for_og_set(s2_rand, o2t2[layer])
                null_J[k] = _jaccard(t1, t2)
            p = float((null_J >= observed_J).sum()) / max(1, len(null_J))
            rows.append({
                'sp1': sp1, 'sp2': sp2, 'contrast': contrast, 'layer': layer,
                'n_sig_sp1': len(sig1), 'n_sig_sp2': len(sig2),
                'n_terms_sp1_sig': len(sig1_terms),
                'n_terms_sp2_sig': len(sig2_terms),
                'n_shared_terms': len(observed_shared),
                'observed_jaccard': observed_J,
                'null_mean_jaccard': float(null_J.mean()),
                'null_q95_jaccard': float(_np.quantile(null_J, 0.95)),
                'p_value': p,
                'shared_terms': ';'.join(sorted(observed_shared)[:30]),
            })
        df = pd.DataFrame(rows)
        return df

    # --- Run for all testable contrasts ---
    _pairs = [
        ('fumigatus','flavus','human_pathogenic_vs_rest'),
        ('fumigatus','flavus','environmental_vs_rest'),
        ('niger','oryzae','industrial_trait_vs_rest'),
        ('niger','oryzae','environmental_vs_rest'),
        ('flavus','niger','environmental_vs_rest'),
    ]
    import numpy as _np
    _t0 = time.time()
    _results = []
    for sp1, sp2, contrast in _pairs:
        res = _functional_convergence_pair(sp1, sp2, contrast, FUNCTIONAL_LAYERS)
        if res is not None and len(res):
            _results.append(res)
    _elapsed = time.time() - _t0
    print(f'\n  Completed {len(_results)} pair x contrast in {_elapsed:.0f}s')

    if _results:
        func_conv_df = pd.concat(_results, ignore_index=True)
        func_conv_df['q_value'] = 1.0
        for (sp1, sp2, c), grp in func_conv_df.groupby(['sp1', 'sp2', 'contrast']):
            func_conv_df.loc[grp.index, 'q_value'] = _mt(grp['p_value'], method='fdr_bh')[1]
        func_conv_df['significant'] = func_conv_df['q_value'] < 0.05
        func_conv_df.to_csv(os.path.join(str(NB3_RESULTS),
                                         'functional_convergence_results.csv'),
                            index=False)
        print(f'\n  Saved: NB3_Results/functional_convergence_results.csv')
        n_sig = int(func_conv_df['significant'].sum())
        print(f'\n  Significant functional convergence (q<0.05 per pair/contrast): {n_sig}')
        if n_sig:
            print(func_conv_df[func_conv_df['significant']][
                ['sp1','sp2','contrast','layer','n_shared_terms',
                 'observed_jaccard','null_mean_jaccard','p_value','q_value',
                 'shared_terms']].to_string(index=False))
        else:
            print('\n  Top 5 trends by p_value (none survive FDR):')
            print(func_conv_df.sort_values('p_value').head(5)[
                ['sp1','sp2','contrast','layer','n_shared_terms',
                 'observed_jaccard','null_mean_jaccard','p_value']].to_string(index=False))
    else:
        print('\n  All pairs have at least one species with 0 sig OGs; functional convergence not testable.')

    return func_conv_df, _results


def plot_functional_convergence_jaccard(func_conv_df, nb3_results):
    """Observed against null Jaccard per species pair and contrast.
    """
    NB3_RESULTS = Path(nb3_results)
    _results = func_conv_df is not None and len(func_conv_df) > 0
    fig = None
    import numpy as _np

    if _results:
        import matplotlib.pyplot as _plt
        pairs_with_data = func_conv_df.groupby(['sp1','sp2','contrast']).size().index.tolist()
        n_panels = len(pairs_with_data)
        if n_panels:
            ncols = 2; nrows = (n_panels + ncols - 1) // ncols
            fig, axes = _plt.subplots(nrows, ncols, figsize=(13, 3.5 * nrows),
                                      squeeze=False)
            for ax, (sp1, sp2, c) in zip(axes.flatten(), pairs_with_data):
                sub = func_conv_df[(func_conv_df['sp1']==sp1) &
                                    (func_conv_df['sp2']==sp2) &
                                    (func_conv_df['contrast']==c)].copy()
                x = _np.arange(len(sub))
                ax.bar(x - 0.18, sub['null_mean_jaccard'], width=0.36,
                       color='#cccccc', edgecolor='white', label='Null mean Jaccard')
                ax.bar(x + 0.18, sub['observed_jaccard'], width=0.36,
                       color=['#d62728' if s else '#1f77b4' for s in sub['significant']],
                       edgecolor='white', label='Observed Jaccard')
                ax.set_xticks(x)
                ax.set_xticklabels(sub['layer'], rotation=45, ha='right', fontsize=7)
                ax.set_ylabel('Jaccard similarity')
                ax.set_title(f'{sp1} vs {sp2}  |  {c}', loc='left',
                             fontsize=9, fontweight='bold')
                for xi, q in zip(x, sub['q_value']):
                    if q < 0.05:
                        ax.text(xi, sub['observed_jaccard'].iloc[xi] + 0.01,
                                f'q={q:.2f}', fontsize=7, ha='center')
                ax.legend(frameon=False, fontsize=7)
            # Hide unused panels
            for ax in axes.flatten()[n_panels:]:
                ax.axis('off')
            fig.suptitle('Functional convergence: observed vs null Jaccard per layer',
                         fontsize=11, fontweight='bold', y=1.005)
            _plt.tight_layout()
            fig.savefig(os.path.join(str(NB3_RESULTS),
                                     'functional_convergence.png'),
                        dpi=200, bbox_inches='tight')
            print('Saved: NB3_Results/functional_convergence.png')
    else:
        print('No functional convergence data to plot.')

    return fig


def plot_convergence_pies(nb3_results, fs_cs=15, dpi_cs=400, fig_w_per_pie_cs=6.5, fig_h_cs=5.6):
    """Concordant up, concordant down and discordant shares as pie charts.
    """
    FS_CS = fs_cs
    DPI_CS = dpi_cs
    FIG_W_PER_PIE_CS = fig_w_per_pie_cs
    FIG_H_CS = fig_h_cs
    NB3_RESULTS = Path(nb3_results)
    fig = None


    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import funpan_utils as fpu

    _panels = [
        ('PAV (pathogenic: fumigatus vs flavus)',
         NB3_RESULTS / 'convergence' / 'pathogenic_method3_pav_convergence.csv'),
        ('CNV (pathogenic: fumigatus vs flavus)',
         NB3_RESULTS / 'convergence' / 'pathogenic_method3_cnv_convergence.csv'),
    ]
    _panels = [(t, p) for t, p in _panels if p.exists()]
    if not _panels:
        print('No per-OG convergence CSVs in NB3_Results/convergence/')
    else:
        fig, axes = plt.subplots(1, len(_panels),
                                 figsize=(FIG_W_PER_PIE_CS * len(_panels), FIG_H_CS),
                                 squeeze=False)
        for ax, (title, path) in zip(axes.ravel(), _panels):
            df = pd.read_csv(path)
            sig = df[df['category'].isin(['CONVERGENT', 'DISCORDANT'])].copy()
            conv = sig[sig['category'] == 'CONVERGENT']
            disc = sig[sig['category'] == 'DISCORDANT']
            up_up   = ((conv['fumigatus_beta'] > 0) & (conv['flavus_beta'] > 0)).sum()
            dn_dn   = ((conv['fumigatus_beta'] < 0) & (conv['flavus_beta'] < 0)).sum()
            disc_n  = len(disc)
            n_total = up_up + dn_dn + disc_n
            sizes  = [up_up, dn_dn, disc_n]
            labels = [f'Concordant UP\n(both β > 0)\nn = {up_up:,}',
                      f'Concordant DOWN\n(both β < 0)\nn = {dn_dn:,}',
                      f'Discordant\n(opposite β)\nn = {disc_n:,}']
            colors = ['#2ca02c', '#1f4e79', '#c0392b']

            wedges, _txt, autotxts = ax.pie(
                sizes, labels=labels, colors=colors,
                autopct=lambda p: f'{p:.1f}%' if p >= 1 else '',
                startangle=90,
                wedgeprops=dict(linewidth=1.0, edgecolor='white'),
                textprops=dict(fontsize=FS_CS - 1,))# fontweight='bold'))
            for t in autotxts:
                t.set_color('white'); t.set_fontsize(FS_CS); t.set_fontweight('bold')

            pct_conv = 100 * (up_up + dn_dn) / max(n_total, 1)
            ax.set_title(
                f"{title}\n"
                f"tested in both = {n_total:,}, concordant = "
                f"{up_up + dn_dn:,} ({pct_conv:.1f}%)",
                fontsize=FS_CS, pad=4)
            ax.set_aspect('equal')

        fig.suptitle('Cross-species Directional Convergence',
                     fontsize=FS_CS + 3, fontweight='bold', y=0.99)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_png = NB3_RESULTS / 'convergence_summary_combined.png'
        fig.savefig(out_png, dpi=DPI_CS, bbox_inches='tight')
        print(f'Saved: {out_png}')

    return fig


def plot_family_comparison_combined(nb3_results, fs_fam=20, dpi_fam=400, top_n_cazy=10, top_n_prot=10, fig_w_fam=12.0, fig_h_fam=8.0, species=SPECIES_LIST):
    """Top CAZy and protease families per species, as two heatmaps.
    """
    FS_FAM = fs_fam
    DPI_FAM = dpi_fam
    TOP_N_CAZY = top_n_cazy
    TOP_N_PROT = top_n_prot
    FIG_W_FAM = fig_w_fam
    FIG_H_FAM = fig_h_fam
    SPECIES = species
    NB3_RESULTS = Path(nb3_results)
    fig = None


    import funpan_utils as fpu

    cazy_path = NB3_RESULTS / 'cazy_family_comparison.csv'
    prot_path = NB3_RESULTS / 'protease_family_comparison.csv'
    if not (cazy_path.exists() and prot_path.exists()):
        print('CAZy or protease comparison CSV missing.')
    else:
        cazy = (pd.read_csv(cazy_path)
                  .sort_values('total', ascending=False).head(TOP_N_CAZY))
        prot = (pd.read_csv(prot_path)
                  .sort_values('total', ascending=False).head(TOP_N_PROT))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W_FAM, FIG_H_FAM))

        def _heat(ax, df, title):
            M = df[SPECIES].values
            im = ax.imshow(M, aspect='auto', cmap='YlOrRd')
            ax.set_xticks(np.arange(len(SPECIES)))
            ax.set_xticklabels([f"A. {s}" for s in SPECIES],
                               fontsize=FS_FAM - 1, rotation=20, ha='right')
            ax.set_yticks(np.arange(len(df)))
            ax.set_yticklabels(df['family'], fontsize=FS_FAM - 2)
            ax.set_title(title, fontsize=FS_FAM + 1, fontweight='bold')
            # Remove any inner gridlines / minor ticks that would slash through annotations
            ax.grid(False)
            ax.tick_params(which='both', length=0)
            ax.set_xticks(np.arange(M.shape[1] + 1) - 0.5, minor=True)
            ax.set_yticks(np.arange(M.shape[0] + 1) - 0.5, minor=True)
            ax.tick_params(which='minor', length=0)
            for spine in ('top', 'right', 'left', 'bottom'):
                ax.spines[spine].set_visible(False)
            vmax = M.max() if M.size else 1
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    v = M[i, j]
                    ax.text(j, i, f'{int(v)}',
                            ha='center', va='center',
                            fontsize=FS_FAM - 2, #fontweight='bold',
                            color='white' if v > vmax * 0.55 else '#222')
            cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
            cb.ax.tick_params(labelsize=FS_FAM - 3)
            cb.set_label('Genes / species', fontsize=FS_FAM - 2)
        _heat(ax1, cazy, f'CAZy families (top {len(cazy)})')
        _heat(ax2, prot, f'Protease families (top {len(prot)})')
        plt.tight_layout()
        out_png = NB3_RESULTS / 'family_comparison_combined.png'
        fig.savefig(out_png, dpi=DPI_FAM, bbox_inches='tight')
        print(f'Saved: {out_png}')

    return fig


def plot_secretome_stacked(nb3_results, fs_sec=18, dpi_sec=400, fig_w_sec=5.5, fig_h_sec=7.5, core_c='#3b8a86', acc_c='#d6a14a', rare_c='#a85d7a', species=SPECIES_LIST):
    """Secretome size per species, split by pangenome class.
    """
    FS_SEC = fs_sec
    DPI_SEC = dpi_sec
    FIG_W_SEC = fig_w_sec
    FIG_H_SEC = fig_h_sec
    CORE_C = core_c
    ACC_C = acc_c
    RARE_C = rare_c
    SPECIES = species
    NB3_RESULTS = Path(nb3_results)
    fig = None


    import funpan_utils as fpu

    _path = NB3_RESULTS / 'secretome_comparison.csv'
    if not _path.exists():
        print(f'No secretome CSV at {_path}.')
    else:
        sec = pd.read_csv(_path).set_index('species').reindex(SPECIES).reset_index()
        core_v = sec['n_secreted_in_Core'].values.astype(int)
        acc_v  = sec['n_secreted_in_Accessory'].values.astype(int)
        rare_v = sec['n_secreted_in_Rare'].values.astype(int)
        tot    = core_v + acc_v + rare_v

        x = np.arange(len(SPECIES))
        fig, ax = plt.subplots(figsize=(FIG_W_SEC, FIG_H_SEC))
        ax.bar(x, core_v, color=CORE_C, edgecolor='white', linewidth=0.8, label='Core')
        ax.bar(x, acc_v,  bottom=core_v, color=ACC_C, edgecolor='white', linewidth=0.8,
               label='Accessory')
        ax.bar(x, rare_v, bottom=core_v + acc_v, color=RARE_C, edgecolor='white',
               linewidth=0.8, label='Rare')

        ymax = tot.max()
        for i in range(len(SPECIES)):
            t = int(tot[i])
            # Inline labels for Core / Accessory (always wide enough vertically)
            for centre, val in [(core_v[i] / 2,             core_v[i]),
                                (core_v[i] + acc_v[i] / 2,  acc_v[i])]:
                pct = 100 * val / t if t else 0
                ax.text(i, centre, f'{int(val):,}\n({pct:.1f}%)',
                        ha='center', va='center',
                        fontsize=FS_SEC - 3, color='white')
            # Rare segment is too thin to label inline -> count + percent above bar
            r = int(rare_v[i]); rpct = 100 * r / t if t else 0
            ax.text(i, t + ymax * 0.015, f'{r:,} ({rpct:.1f}%)',
                    ha='center', va='bottom',
                    fontsize=FS_SEC - 3, color=RARE_C)

        ax.set_xticks(x)
        ax.set_xticklabels([f"A. {s}\nn = {int(t):,}"
                            for s, t in zip(SPECIES, tot)],
                           fontsize=FS_SEC - 2,
                           rotation=40, ha='right', rotation_mode='anchor')
        ax.set_ylabel('Number of secreted orthogroups', fontsize=FS_SEC)
        ax.set_ylim(0, ymax * 1.13)
        ax.tick_params(axis='y', labelsize=FS_SEC - 3)
        ax.set_title('Secretome Composition',
                     fontsize=FS_SEC + 1, fontweight='bold')
        ax.legend(loc='lower right', fontsize=FS_SEC - 3, frameon=True, framealpha=0.95)
        ax.grid(True, axis='y', linestyle=':', alpha=0.5)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        plt.tight_layout()

        out_png = NB3_RESULTS / 'secretome_combined.png'
        fig.savefig(out_png, dpi=DPI_SEC, bbox_inches='tight')
        print(f'Saved: {out_png}')

    return fig


def plot_functional_null_result(nb3_results, fs_fc=15, dpi_fc=400, fig_w_fc=13.0, fig_h_fc=6.0, q_threshold_fc=0.05, sp1_fc='fumigatus', sp2_fc='flavus', contrast_fc='human_pathogenic_vs_rest'):
    """Observed against null Jaccard per annotation layer, with the q-value per layer.
    """
    FS_FC = fs_fc
    DPI_FC = dpi_fc
    FIG_W_FC = fig_w_fc
    FIG_H_FC = fig_h_fc
    Q_THRESHOLD_FC = q_threshold_fc
    SP1_FC = sp1_fc
    SP2_FC = sp2_fc
    CONTRAST_FC = contrast_fc
    NB3_RESULTS = Path(nb3_results)
    fig = None


    import funpan_utils as fpu

    # Focal contrast for this panel (used only if the results file holds several)

    _path = NB3_RESULTS / 'functional_convergence_results_FILTERED.csv'
    if not _path.exists():
        _path = NB3_RESULTS / 'functional_convergence_results.csv'   # current results
    if not _path.exists():
        print(f'No functional convergence results at {_path}.')
    else:
        fc = pd.read_csv(_path)
        if {'sp1', 'sp2', 'contrast'} <= set(fc.columns):
            fc = fc[(fc['sp1'] == SP1_FC) & (fc['sp2'] == SP2_FC) & (fc['contrast'] == CONTRAST_FC)]
            print(f'{SP1_FC} x {SP2_FC} | {CONTRAST_FC}: {len(fc)} annotation layers')
        # Sort by observed jaccard so the visually strongest layers go first
        fc = fc.sort_values('observed_jaccard', ascending=False).reset_index(drop=True)

        fig, ax = plt.subplots(figsize=(FIG_W_FC, FIG_H_FC))
        x = np.arange(len(fc))
        w = 0.36

        # Color the observed bars by whether obs > null (trending convergent) or not
        obs_colors = ['#27ae60' if o > n else '#7f8c8d'
                      for o, n in zip(fc['observed_jaccard'], fc['null_mean_jaccard'])]
        ax.bar(x - w/2, fc['observed_jaccard'], width=w, color=obs_colors,
               edgecolor='white', linewidth=0.6,
               label='Observed Jaccard\n(green = obs > null trend; grey = below)')
        ax.bar(x + w/2, fc['null_mean_jaccard'], width=w, color='#bdc3c7',
               edgecolor='white', linewidth=0.6,
               label='Null mean (1000 permutations)')

        # Per-layer q-value annotation above the higher of the two bars
        y_room = max(fc[['observed_jaccard', 'null_mean_jaccard']].max().max() * 1.45, 0.1)
        for i, row in fc.iterrows():
            tip = max(row['observed_jaccard'], row['null_mean_jaccard'])
            sig = bool(row['significant'])
            q   = row['q_value']
            col = '#c0392b' if sig else '#666'
            marker = '★' if sig else 'ns'
            ax.text(i, tip + y_room * 0.04,
                    f"q = {q:.2f}\n{marker}",
                    ha='center', va='bottom',
                    fontsize=FS_FC - 4, color=col, fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(fc['layer'], rotation=20, ha='right',
                           fontsize=FS_FC - 2)
        ax.set_ylabel('Jaccard index', fontsize=FS_FC)
        ax.set_ylim(0, y_room)
        n_sig    = int(fc['significant'].sum())
        n_total  = len(fc)
        n_trend  = int((fc['observed_jaccard'] > fc['null_mean_jaccard']).sum())
        ax.set_title('Cross-species functional convergence\n'
                     'A. fumigatus × A. flavus, human-pathogenic vs rest',
                     fontsize=FS_FC + 1, fontweight='bold')

        ax.legend(loc='upper left', fontsize=FS_FC - 3, frameon=True, framealpha=0.95)
        ax.tick_params(axis='y', labelsize=FS_FC - 2)
        ax.grid(True, axis='y', linestyle=':', alpha=0.5)
        plt.tight_layout()
        out_png = NB3_RESULTS / 'functional_convergence_combined.png'
        fig.savefig(out_png, dpi=DPI_FC, bbox_inches='tight')
        print(f'Saved: {out_png}')

    return fig


def load_combined_orthofinder(genecount_path, og_tsv_path, species_list,
                              ani_excluded=None):
    """Load the combined OrthoFinder run and restrict it to the analysed genomes.

    The combined run predates ANI and contamination QC, so excluded genomes are
    still present as columns and are dropped after loading.

    Returns
    -------
    combined_genecount, combined_og_long, genome_mapping
    """
    if ani_excluded is None:
        ani_excluded = ANI_EXCLUDED
    excl = set().union(*ani_excluded.values())

    print('Loading combined OrthoFinder gene count matrix...')
    combined_genecount = load_combined_orthogroups_genecount(str(genecount_path))
    print(f'  Gene count matrix: {combined_genecount.shape}')

    print('\nLoading combined OrthoFinder long format...')
    combined_og_long = load_combined_orthogroups_long(str(og_tsv_path))
    print(f'  Long format: {combined_og_long.shape}')
    print(f'  Unique orthogroups: {combined_og_long["Orthogroup"].nunique()}')
    print(f'  Unique species: {combined_og_long["Species"].nunique()}')

    dropped = [c for c in combined_genecount.columns if any(a in str(c) for a in excl)]
    combined_genecount = combined_genecount[
        [c for c in combined_genecount.columns if c not in dropped]]
    combined_og_long = combined_og_long[
        ~combined_og_long['Assembly_Accession'].astype(str).apply(
            lambda a: any(x in a for x in excl))]
    print(f'  Dropped {len(dropped)} excluded genome column(s): {dropped}')

    genome_mapping = build_genome_species_mapping(combined_genecount)
    print(f'\nGenome-species mapping: {len(genome_mapping)} genomes')
    for sp in species_list:
        n = sum(1 for v in genome_mapping.values() if v == sp)
        print(f'  {sp}: {n} genomes')
    return combined_genecount, combined_og_long, genome_mapping


def load_og_consensus_filtered(species_list, nb1_results, species_root,
                               ani_excluded=None):
    """Load per-species orthogroup annotations, without contaminant-private rows.

    ``{sp}_og_consensus.tsv`` comes from the full OrthoFinder run, so it carries
    orthogroups whose only member proteins sit on an excluded genome. Dropping a
    genome column cannot remove those rows, so they are filtered here.

    Returns
    -------
    dict
        ``{species: DataFrame}``.
    """
    import glob as _glob

    if ani_excluded is None:
        ani_excluded = ANI_EXCLUDED
    excl = set().union(*ani_excluded.values())

    og_consensus = load_og_consensus_annotations(species_list, nb1_results)
    for sp in species_list:
        if sp not in og_consensus or og_consensus[sp].empty:
            continue
        priv = set()
        pattern = f'{species_root}/{sp}/orthofinder_output/*/Orthogroups/'
        for f in _glob.glob(pattern + 'Orthogroups.tsv') + \
                 _glob.glob(pattern + 'Orthogroups_UnassignedGenes.tsv'):
            with open(f) as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                ex = {i for i, h in enumerate(hdr) if any(a in h for a in excl)}
                if not ex:
                    continue
                for line in fh:
                    row = line.rstrip('\n').split('\t')
                    present = [i for i in range(1, len(row)) if row[i].strip()]
                    if present and all(i in ex for i in present):
                        priv.add(row[0])
        if priv:
            df = og_consensus[sp]
            og_consensus[sp] = df[~df['Orthogroup'].isin(priv)].copy()
            print(f'{sp}: dropped {len(priv)} contaminant-private OGs from og_consensus')
    return og_consensus


def nb3_preflight(species_list, species_root, nb0_results, nb1_results,
                  nb2_results, nb3_results, combined_of_dir, verbose=True):
    """Check NB3's inputs exist and create the output directories.

    Looks for the combined OrthoFinder run, NB0's phenotype table, and per
    species the NB1 annotation table and the NB2 association results. Nothing is
    computed or written beyond directories.

    Returns
    -------
    pd.DataFrame
        One row per species with a boolean per required input.
    """
    import glob as _glob

    nb3_results = Path(nb3_results)
    (nb3_results / 'convergence').mkdir(parents=True, exist_ok=True)

    rows = []
    for sp in species_list:
        rows.append({
            'species': sp,
            'og_consensus': (Path(nb1_results) / sp / f'{sp}_og_consensus.tsv').exists(),
            'orthofinder': bool(_glob.glob(f'{species_root}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups.tsv')),
            'gwas_results': bool(_glob.glob(f'{nb2_results}/{sp}/pangwas_results/*_assoc_*.tsv')),
            'enrichment': (Path(nb2_results) / sp / f'{sp}_functional_enrichment.csv').exists(),
        })
    status = pd.DataFrame(rows).set_index('species')

    combined = Path(combined_of_dir) / 'Orthogroups' / 'Orthogroups.tsv'
    pheno = Path(nb0_results) / 'phenotype_classified_for_gwas.csv'
    if verbose:
        print('NB3 inputs')
        print('=' * 60)
        print(f'  combined OrthoFinder run: {"found" if combined.exists() else "MISSING"}')
        print(f'  -> {combined}')
        print(f'  NB0 phenotype table: {"found" if pheno.exists() else "MISSING, run NB0 first"}')
        print()
        print(status.replace({True: 'ok', False: 'MISSING'}).to_string())
        missing = int((~status).sum().sum())
        print()
        if missing or not combined.exists() or not pheno.exists():
            print(f'{missing} missing per-species input(s); run NB1 and NB2 first.')
        else:
            print('All inputs present.')
    return status


def print_nb3_summary(species_list, nb3_results):
    """Print every NB3 output with its full path, and the convergence counts.

    Returns
    -------
    counts : pd.DataFrame
        Directional concordance per layer, empty if Method 3 has not run.
    files : pd.DataFrame
        Every expected output with its full path and whether it exists.
    """
    nb3_results = Path(nb3_results)
    expected = [
        'method3_projection_accounting.csv',
        'method3_directional_concordance.csv',
        'method3_merge_rule_sensitivity.csv',
        'convergence/pathogenic_method3_pav_convergence.csv',
        'convergence/pathogenic_method3_cnv_convergence.csv',
        'cross_species_gcf_convergence.csv',
        'permutation_null_model.csv',
        'convergence_null_model.png',
        'cazy_family_comparison.csv',
        'protease_family_comparison.csv',
        'secretome_comparison.csv',
        'functional_convergence_results.csv',
        'functional_convergence.png',
        'convergence_summary_combined.png',
        'family_comparison_combined.png',
        'secretome_combined.png',
        'functional_convergence_combined.png',
    ]
    rows = [{'path': str(nb3_results / n), 'exists': (nb3_results / n).exists()}
            for n in expected]
    files = pd.DataFrame(rows)

    conc = nb3_results / 'method3_directional_concordance.csv'
    counts = pd.read_csv(conc) if conc.exists() else pd.DataFrame()
    if not counts.empty:
        keep = [c for c in ['layer', 'n_shared_tested', 'n_concordant', 'n_discordant',
                            'pct_concordant', 'sign_test_p', 'spearman_rho', 'n_sig_both']
                if c in counts.columns]
        counts = counts[keep].set_index('layer')

    print('=' * 70)
    print('NB3 Cross-Species Convergence complete.')
    print(f'All outputs saved to: {nb3_results}')
    for _, r in files.iterrows():
        print(f"  {' ' if r['exists'] else '!'} {r['path']}")
    n_missing = int((~files['exists']).sum())
    if n_missing:
        print(f'\n  ! = missing ({n_missing} path(s))')
    print('=' * 70)
    print('\nDirectional concordance per layer:')
    return counts, files
