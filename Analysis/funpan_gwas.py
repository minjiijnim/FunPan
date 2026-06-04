#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_gwas.py
==============
Pan-GWAS association testing, diagnostic plotting,
functional enrichment of significant hits, and power analysis.

Extracted from Theme2_Analysis_Functions.py and power_analysis.py.
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
# PHENOTYPE AND CONTRAST SETUP
# =============================================================================

def create_phenotype_from_isolation_class(
    pav: pd.DataFrame,
    metadata: pd.DataFrame,
    isolation_class_col: str = 'IsolationClass',
    accession_col: str = 'Assembly_Accession',
    species_col: str = 'Species'
) -> pd.DataFrame:
    """
    Create phenotype DataFrame from metadata matching PAV matrix samples.

    Parameters
    ----------
    pav : pd.DataFrame
        PAV matrix with genome columns
    metadata : pd.DataFrame
        Metadata with isolation class
    isolation_class_col : str
        Column name for isolation class
    accession_col : str
        Column name for accession
    species_col : str
        Column name for species

    Returns
    -------
    pd.DataFrame
        Phenotype DataFrame with columns: sample_id, species, isolation_class
    """
    # Handle different column name formats
    if accession_col not in metadata.columns:
        for col in metadata.columns:
            if 'assembly' in col.lower() and 'accession' in col.lower():
                accession_col = col
                break

    # Try to find the isolation class column - prefer Phenotype from NB2
    if isolation_class_col not in metadata.columns:
        # Priority order: Phenotype > IsolationClass > other variants
        for col_name in ['Phenotype', 'IsolationClass', 'isolation_class']:
            if col_name in metadata.columns:
                isolation_class_col = col_name
                break
        else:
            # Fallback: search for any column with isolation/class
            for col in metadata.columns:
                if 'isolation' in col.lower() and 'class' in col.lower():
                    isolation_class_col = col
                    break

    # Build lookup - extract just the GCA/GCF accession from metadata
    acc_to_class = {}
    acc_to_species_meta = {}

    # Valid classes from NB2 Phenotype categories AND legacy IsolationClass
    valid_phenotype_classes = {
        'Human-pathogenic', 'Animal-pathogenic', 'Plant-pathogenic',
        'Industrial-trait', 'Environmental', 'Lab'
    }
    valid_legacy_classes = {
        'Human', 'Animal', 'Plant', 'Industrial', 'Environmental', 'Lab'
    }
    valid_classes = valid_phenotype_classes | valid_legacy_classes

    for _, row in metadata.iterrows():
        acc_raw = str(row.get(accession_col, ''))
        # Extract GCA/GCF accession
        match = ACC_RE.search(acc_raw)
        if match:
            acc = match.group(1)
        else:
            acc = acc_raw

        iso_class = row.get(isolation_class_col, 'Unknown')
        # Only use valid isolation classes
        if iso_class not in valid_classes:
            iso_class = 'Unknown'

        acc_to_class[acc] = iso_class

        if species_col in metadata.columns:
            acc_to_species_meta[acc] = row.get(species_col, 'unknown')

    records = []
    for col in pav.columns:
        if col == 'Orthogroup':
            continue

        # Parse species and accession from column name
        # Format: species__GCA_xxx.x_assembly.proteins or species__GCA_xxx.x
        if '__' in col:
            parts = col.split('__')
            species = parts[0]
            rest = parts[1]
            # Extract just the GCA/GCF accession using regex
            match = ACC_RE.search(rest)
            if match:
                accession = match.group(1)
            else:
                accession = rest
        else:
            # Try to extract accession directly
            match = ACC_RE.search(col)
            if match:
                accession = match.group(1)
            else:
                accession = col
            species = acc_to_species_meta.get(accession, 'unknown')

        iso_class = acc_to_class.get(accession, 'Unknown')

        records.append({
            'sample_id': col,
            'species': species,
            'accession': accession,
            'isolation_class': iso_class
        })

    phenotype_df = pd.DataFrame(records)

    # Report matching stats
    n_matched = (phenotype_df['isolation_class'] != 'Unknown').sum()
    n_total = len(phenotype_df)
    print(f"Matched {n_matched}/{n_total} samples to isolation class")

    return phenotype_df


def define_isolation_contrasts() -> List[Dict[str, Any]]:
    """
    Define contrasts for isolation class comparisons.

    Returns
    -------
    List[Dict]
        List of contrast definitions
    """
    # Using NB2 Phenotype categories
    # Pathogenic classes: Human-pathogenic, Animal-pathogenic, Plant-pathogenic
    # Trait classes: Industrial-trait, Environmental, Lab
    contrasts = [
        {
            'name': 'Industrial_vs_Others',
            'case_class': 'Industrial-trait',
            'control_classes': ['Human-pathogenic', 'Animal-pathogenic',
                               'Plant-pathogenic', 'Environmental'],
            'description': 'Industrial-trait strains vs pathogenic/environmental'
        },
        {
            'name': 'Human_Pathogenic_vs_Others',
            'case_class': 'Human-pathogenic',
            'control_classes': ['Industrial-trait', 'Plant-pathogenic',
                               'Environmental', 'Animal-pathogenic'],
            'description': 'Human-pathogenic strains vs others'
        },
        {
            'name': 'Industrial_vs_Human_Pathogenic',
            'case_class': 'Industrial-trait',
            'control_classes': ['Human-pathogenic'],
            'description': 'Industrial-trait vs human-pathogenic (direct comparison)'
        },
        {
            'name': 'Plant_Pathogenic_vs_Others',
            'case_class': 'Plant-pathogenic',
            'control_classes': ['Human-pathogenic', 'Animal-pathogenic',
                               'Industrial-trait', 'Environmental'],
            'description': 'Plant-pathogenic strains vs others'
        },
        {
            'name': 'Environmental_vs_Adapted',
            'case_class': 'Environmental',
            'control_classes': ['Human-pathogenic', 'Animal-pathogenic',
                               'Industrial-trait', 'Plant-pathogenic'],
            'description': 'Environmental vs host/trait-adapted strains'
        },
        {
            'name': 'Pathogenic_vs_NonPathogenic',
            'case_class': ['Human-pathogenic', 'Animal-pathogenic', 'Plant-pathogenic'],
            'control_classes': ['Industrial-trait', 'Environmental', 'Lab'],
            'description': 'All pathogenic vs non-pathogenic strains'
        }
    ]

    return contrasts


# =============================================================================
# ASSOCIATION TESTING (PER-SPECIES)
# =============================================================================

def run_fisher_test_orthogroup(
    pav_row: pd.Series,
    case_samples: List[str],
    control_samples: List[str]
) -> Tuple[float, float, str]:
    """
    Run Fisher's exact test for a single orthogroup.

    Parameters
    ----------
    pav_row : pd.Series
        Presence/absence row for one orthogroup
    case_samples : List[str]
        Sample IDs in case group
    control_samples : List[str]
        Sample IDs in control group

    Returns
    -------
    Tuple[float, float, str]
        (odds_ratio, p_value, direction)
    """
    # Count presence in case and control
    case_present = sum(pav_row.get(s, 0) > 0 for s in case_samples if s in pav_row.index)
    case_absent = len([s for s in case_samples if s in pav_row.index]) - case_present

    ctrl_present = sum(pav_row.get(s, 0) > 0 for s in control_samples if s in pav_row.index)
    ctrl_absent = len([s for s in control_samples if s in pav_row.index]) - ctrl_present

    # Build contingency table
    # [[case_present, case_absent], [ctrl_present, ctrl_absent]]
    table = [[case_present, case_absent], [ctrl_present, ctrl_absent]]

    try:
        odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
    except Exception:
        return np.nan, 1.0, 'none'

    # Determine direction
    if odds_ratio > 1:
        direction = 'enriched_in_case'
    elif odds_ratio < 1:
        direction = 'depleted_in_case'
    else:
        direction = 'none'

    return odds_ratio, p_value, direction


def run_species_association_test(
    pav: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    contrast: Dict[str, Any],
    min_samples_per_group: int = 3
) -> Optional[pd.DataFrame]:
    """
    Run association test for all orthogroups within a single species.

    Parameters
    ----------
    pav : pd.DataFrame
        PAV matrix for one species (rows=orthogroups, cols=accessions)
    phenotype_df : pd.DataFrame
        Phenotype data with sample_id, isolation_class columns
    contrast : Dict
        Contrast definition with case_class and control_classes
    min_samples_per_group : int
        Minimum samples required per group

    Returns
    -------
    pd.DataFrame or None
        Results with columns: Orthogroup, odds_ratio, p_value, direction, q_value
    """
    case_class = contrast['case_class']
    control_classes = contrast['control_classes']

    # Handle case_class as list or string
    if isinstance(case_class, str):
        case_classes = [case_class]
    else:
        case_classes = list(case_class)

    # Map accession to isolation classes (use 'accession' column, not 'sample_id')
    # This handles the case where PAV columns are just accessions (after split_pav_by_species)
    if 'accession' in phenotype_df.columns:
        sample_class = dict(zip(phenotype_df['accession'], phenotype_df['isolation_class']))
    else:
        # Fallback: extract accession from sample_id
        sample_class = {}
        for _, row in phenotype_df.iterrows():
            sid = row['sample_id']
            iso = row['isolation_class']
            # Add both full sample_id and extracted accession
            sample_class[sid] = iso
            match = ACC_RE.search(sid)
            if match:
                sample_class[match.group(1)] = iso

    # Find case and control samples that exist in PAV
    case_samples = []
    control_samples = []

    for col in pav.columns:
        # Try direct lookup first
        iso = sample_class.get(col)

        # If not found, try extracting accession from column name
        if iso is None:
            match = ACC_RE.search(col)
            if match:
                iso = sample_class.get(match.group(1))

        if iso in case_classes:
            case_samples.append(col)
        elif iso in control_classes:
            control_samples.append(col)

    if len(case_samples) < min_samples_per_group or len(control_samples) < min_samples_per_group:
        print(f"  Insufficient samples: {len(case_samples)} cases, {len(control_samples)} controls")
        return None

    results = []
    for orthogroup in pav.index:
        pav_row = pav.loc[orthogroup]
        odds_ratio, p_value, direction = run_fisher_test_orthogroup(
            pav_row, case_samples, control_samples
        )

        # Calculate frequencies
        case_freq = sum(pav_row.get(s, 0) > 0 for s in case_samples) / len(case_samples) if case_samples else 0
        ctrl_freq = sum(pav_row.get(s, 0) > 0 for s in control_samples) / len(control_samples) if control_samples else 0

        results.append({
            'Orthogroup': orthogroup,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'direction': direction,
            'case_freq': case_freq,
            'control_freq': ctrl_freq,
            'n_case': len(case_samples),
            'n_control': len(control_samples)
        })

    results_df = pd.DataFrame(results)

    # Calculate FDR
    valid_mask = ~results_df['p_value'].isna()
    results_df['q_value'] = 1.0
    if valid_mask.sum() > 0:
        _, q_values, _, _ = multipletests(
            results_df.loc[valid_mask, 'p_value'].values,
            method='fdr_bh'
        )
        results_df.loc[valid_mask, 'q_value'] = q_values

    return results_df


def run_all_species_association(
    species_pav: Dict[str, pd.DataFrame],
    species_phenotypes: Dict[str, pd.DataFrame],
    contrast: Dict[str, Any]
) -> Dict[str, pd.DataFrame]:
    """
    Run association tests for all species.

    Parameters
    ----------
    species_pav : Dict[str, pd.DataFrame]
        Per-species PAV matrices
    species_phenotypes : Dict[str, pd.DataFrame]
        Per-species phenotype data
    contrast : Dict
        Contrast definition

    Returns
    -------
    Dict[str, pd.DataFrame]
        Per-species association results
    """
    results = {}

    print(f"\nRunning contrast: {contrast['name']}")
    print(f"  {contrast['description']}")

    for species in species_pav.keys():
        if species not in species_phenotypes:
            print(f"  {species}: No phenotype data")
            continue

        print(f"  Testing {species}...")
        pav = species_pav[species]
        pheno = species_phenotypes[species]

        res = run_species_association_test(pav, pheno, contrast)
        if res is not None:
            res['species'] = species
            results[species] = res

            n_sig = (res['q_value'] < 0.1).sum()
            print(f"    {len(res)} orthogroups tested, {n_sig} significant (FDR < 0.1)")

    return results


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_species_results_comparison(
    species_results: Dict[str, pd.DataFrame],
    contrast_name: str,
    fdr_threshold: float = 0.1,
    figsize: Tuple[int, int] = (14, 10)
) -> plt.Figure:
    """
    Plot comparison of association results across species.

    Parameters
    ----------
    species_results : Dict[str, pd.DataFrame]
        Per-species association results
    contrast_name : str
        Name of the contrast for plot title
    fdr_threshold : float
        FDR threshold for significance
    figsize : Tuple
        Figure size

    Returns
    -------
    plt.Figure
        The figure object
    """
    n_species = len(species_results)
    if n_species == 0:
        return None

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    # Plot 1: Number of significant hits per species
    ax = axes[0]
    species_names = list(species_results.keys())
    n_sig = [sum(res['q_value'] < fdr_threshold) for res in species_results.values()]
    colors = [SPECIES_COLORS.get(s, 'gray') for s in species_names]

    ax.bar(species_names, n_sig, color=colors, alpha=0.8)
    ax.set_ylabel('Number of significant orthogroups')
    ax.set_title(f'Significant hits per species (FDR < {fdr_threshold})')
    ax.tick_params(axis='x', rotation=45)

    # Plot 2: P-value distributions
    ax = axes[1]
    for species, res in species_results.items():
        p_vals = res['p_value'].dropna()
        ax.hist(p_vals, bins=50, alpha=0.5, label=SPECIES_DISPLAY.get(species, species),
               color=SPECIES_COLORS.get(species, 'gray'))
    ax.set_xlabel('P-value')
    ax.set_ylabel('Count')
    ax.set_title('P-value distributions')
    ax.legend()

    # Plot 3: Effect sizes (odds ratios)
    ax = axes[2]
    for species, res in species_results.items():
        sig = res[res['q_value'] < fdr_threshold]
        if len(sig) > 0:
            log_or = np.log2(sig['odds_ratio'].replace([np.inf, -np.inf], np.nan).dropna())
            ax.boxplot([log_or], positions=[list(species_results.keys()).index(species)],
                      widths=0.6, patch_artist=True,
                      boxprops=dict(facecolor=SPECIES_COLORS.get(species, 'gray'), alpha=0.7))

    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names, rotation=45)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_ylabel('log2(Odds Ratio)')
    ax.set_title('Effect sizes (significant hits)')

    # Plot 4: Overlap of significant hits
    ax = axes[3]
    sig_sets = {}
    for species, res in species_results.items():
        sig_sets[species] = set(res[res['q_value'] < fdr_threshold]['Orthogroup'])

    # Count shared hits
    if len(sig_sets) >= 2:
        shared_counts = {}
        for i, s1 in enumerate(species_names):
            for s2 in species_names[i+1:]:
                overlap = len(sig_sets.get(s1, set()) & sig_sets.get(s2, set()))
                shared_counts[f"{s1[:3]}-{s2[:3]}"] = overlap

        ax.bar(shared_counts.keys(), shared_counts.values(), color='steelblue', alpha=0.7)
        ax.set_ylabel('Shared significant orthogroups')
        ax.set_title('Overlap between species')
        ax.tick_params(axis='x', rotation=45)

    fig.suptitle(f'Cross-species comparison: {contrast_name}', fontsize=14, fontweight='bold')
    plt.tight_layout()



def plot_convergent_hits_summary(
    convergent_hits: pd.DataFrame,
    contrast_name: str,
    figsize: Tuple[int, int] = (12, 8)
) -> plt.Figure:
    """
    Plot summary of convergent hits.

    Parameters
    ----------
    convergent_hits : pd.DataFrame
        Convergent hit results
    contrast_name : str
        Name of the contrast
    figsize : Tuple
        Figure size

    Returns
    -------
    plt.Figure
        The figure object
    """
    if convergent_hits.empty:
        print("No convergent hits to plot")
        return None

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: Convergence type distribution
    ax = axes[0, 0]
    type_counts = convergent_hits['convergence_type'].value_counts()
    ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%',
          colors=['#3498db', '#e74c3c', '#27ae60'])
    ax.set_title('Convergence type distribution')

    # Plot 2: Direction distribution
    ax = axes[0, 1]
    dir_counts = convergent_hits['convergent_direction'].value_counts()
    colors_dir = {'enriched_in_case': '#27ae60', 'depleted_in_case': '#e74c3c', 'none': 'gray'}
    ax.bar(dir_counts.index, dir_counts.values,
          color=[colors_dir.get(d, 'gray') for d in dir_counts.index], alpha=0.8)
    ax.set_ylabel('Number of orthogroups')
    ax.set_title('Direction of effect')
    ax.tick_params(axis='x', rotation=45)

    # Plot 3: Number of species contributing
    ax = axes[1, 0]
    n_species = convergent_hits['n_species_tested'].value_counts().sort_index()
    ax.bar(n_species.index, n_species.values, color='steelblue', alpha=0.8)
    ax.set_xlabel('Number of species tested')
    ax.set_ylabel('Number of orthogroups')
    ax.set_title('Species coverage of convergent hits')

    # Plot 4: P-value distribution
    ax = axes[1, 1]
    ax.hist(convergent_hits['min_p_value'].dropna(), bins=30, color='purple', alpha=0.7)
    ax.set_xlabel('Minimum p-value across species')
    ax.set_ylabel('Count')
    ax.set_title('P-value distribution of convergent hits')
    ax.set_yscale('log')

    fig.suptitle(f'Convergent hits summary: {contrast_name}\n({len(convergent_hits)} total hits)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()



def plot_meta_analysis_results(
    meta_df: pd.DataFrame,
    contrast_name: str,
    fdr_threshold: float = 0.1,
    top_n: int = 30,
    figsize: Tuple[int, int] = (14, 10)
) -> plt.Figure:
    """
    Plot meta-analysis results.

    Parameters
    ----------
    meta_df : pd.DataFrame
        Meta-analysis results
    contrast_name : str
        Name of the contrast
    fdr_threshold : float
        FDR threshold for significance
    top_n : int
        Number of top hits to show
    figsize : Tuple
        Figure size

    Returns
    -------
    plt.Figure
        The figure object
    """
    if meta_df.empty:
        print("No meta-analysis results to plot")
        return None

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: Q-Q plot
    ax = axes[0, 0]
    observed = -np.log10(meta_df['meta_p_value'].dropna().sort_values())
    n = len(observed)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))

    ax.scatter(expected, observed, alpha=0.5, s=20)
    max_val = max(max(expected), max(observed))
    ax.plot([0, max_val], [0, max_val], 'r--', label='Expected')
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')
    ax.set_title('Q-Q Plot')
    ax.legend()

    # Plot 2: Volcano-style plot
    ax = axes[0, 1]
    x = meta_df['mean_log_odds_ratio'].replace([np.inf, -np.inf], np.nan)
    y = -np.log10(meta_df['meta_p_value'])

    colors = ['#e74c3c' if q < fdr_threshold else '#bdc3c7'
             for q in meta_df['meta_q_value']]

    ax.scatter(x, y, c=colors, alpha=0.6, s=30)
    ax.axhline(-np.log10(fdr_threshold), color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Mean log(Odds Ratio)')
    ax.set_ylabel('-log10(Meta p-value)')
    ax.set_title('Volcano plot')

    # Plot 3: Top hits by p-value
    ax = axes[1, 0]
    top_hits = meta_df.nsmallest(min(top_n, len(meta_df)), 'meta_p_value')

    y_pos = range(len(top_hits))
    colors_bar = ['#27ae60' if d == 'enriched' else '#e74c3c' if d == 'depleted' else 'gray'
                 for d in top_hits['consensus_direction']]

    ax.barh(y_pos, -np.log10(top_hits['meta_p_value']), color=colors_bar, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_hits['Orthogroup'], fontsize=8)
    ax.set_xlabel('-log10(Meta p-value)')
    ax.set_title(f'Top {min(top_n, len(top_hits))} hits by meta p-value')
    ax.invert_yaxis()

    # Plot 4: Direction consistency
    ax = axes[1, 1]
    consistent = meta_df['direction_consistent'].value_counts()
    ax.pie(consistent.values, labels=['Consistent', 'Mixed'] if True in consistent.index else ['Mixed', 'Consistent'],
          autopct='%1.1f%%', colors=['#27ae60', '#e74c3c'])
    ax.set_title('Direction consistency across species')

    n_sig = (meta_df['meta_q_value'] < fdr_threshold).sum()
    fig.suptitle(f'Meta-analysis results: {contrast_name}\n({n_sig} significant at FDR < {fdr_threshold})',
                fontsize=14, fontweight='bold')
    plt.tight_layout()



def plot_robustness_analysis(
    robustness_df: pd.DataFrame,
    contrast_name: str,
    figsize: Tuple[int, int] = (12, 8)
) -> plt.Figure:
    """
    Plot leave-one-out robustness analysis results.

    Parameters
    ----------
    robustness_df : pd.DataFrame
        Robustness analysis results
    contrast_name : str
        Name of the contrast
    figsize : Tuple
        Figure size

    Returns
    -------
    plt.Figure
        The figure object
    """
    if robustness_df.empty:
        print("No robustness results to plot")
        return None

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: Tier distribution
    ax = axes[0, 0]
    if 'tier' in robustness_df.columns:
        tier_counts = robustness_df['tier'].value_counts()
        colors_tier = {
            'Tier1_Robust': '#27ae60',
            'Tier2_PartiallyRobust': '#f39c12',
            'Tier3_SpeciesDependent': '#e74c3c'
        }
        ax.bar(tier_counts.index, tier_counts.values,
              color=[colors_tier.get(t, 'gray') for t in tier_counts.index], alpha=0.8)
        ax.set_ylabel('Number of orthogroups')
        ax.set_title('Robustness tier distribution')
        ax.tick_params(axis='x', rotation=45)

    # Plot 2: Driving species
    ax = axes[0, 1]
    driver_counts = robustness_df['driving_species'].value_counts()
    colors_sp = [SPECIES_COLORS.get(s, 'gray') for s in driver_counts.index]
    ax.bar(driver_counts.index, driver_counts.values, color=colors_sp, alpha=0.8)
    ax.set_ylabel('Number of orthogroups')
    ax.set_title('Signal-driving species')
    ax.tick_params(axis='x', rotation=45)

    # Plot 3: Original vs max LOO p-value
    ax = axes[1, 0]
    x = -np.log10(robustness_df['original_meta_p'].replace(0, 1e-300))
    y = -np.log10(robustness_df['max_loo_p'].replace(0, 1e-300))

    colors_robust = ['#27ae60' if r else '#e74c3c' for r in robustness_df['robust']]
    ax.scatter(x, y, c=colors_robust, alpha=0.6, s=40)

    max_val = max(max(x.dropna()), max(y.dropna()))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax.set_xlabel('-log10(Original meta p-value)')
    ax.set_ylabel('-log10(Max LOO p-value)')
    ax.set_title('Robustness: Original vs LOO')

    # Plot 4: LOO significance retention
    ax = axes[1, 1]
    retention = robustness_df['n_still_significant'] / robustness_df['n_loo_tests']
    ax.hist(retention.dropna(), bins=20, color='steelblue', alpha=0.7)
    ax.axvline(1.0, color='green', linestyle='--', label='Fully robust')
    ax.set_xlabel('Fraction of LOO tests still significant')
    ax.set_ylabel('Count')
    ax.set_title('LOO significance retention')
    ax.legend()

    n_robust = robustness_df['robust'].sum() if 'robust' in robustness_df.columns else 0
    fig.suptitle(f'Leave-one-out robustness: {contrast_name}\n({n_robust} robust hits)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()




# =============================================================================
# POWER ANALYSIS
# =============================================================================

# Statistical Power Analysis: A. oryzae (31:1) vs A. niger (13:15) pan-GWAS designs
# Addresses reviewer concern: Is the oryzae null result due to extreme case-control imbalance?
import pandas as pd
from scipy import stats
from itertools import product
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================================
# PART 1: Analytical minimum detectable effect sizes
# ============================================================================
# print("=" * 80)
# print("PART 1: ANALYTICAL POWER CALCULATIONS (Fisher's exact test framework)")
# print("=" * 80)

def power_fisher(n_case, n_control, gene_freq, odds_ratio, alpha=0.05, n_sim=50000):
    """
    Estimate power for Fisher's exact test via simulation.
    gene_freq = baseline frequency in controls.
    odds_ratio = OR for cases vs controls.
    """
    # Probability of gene presence in controls
    p0 = gene_freq
    # Probability of gene presence in cases (derived from OR)
    p1 = (odds_ratio * p0) / (1 - p0 + odds_ratio * p0)
    p1 = min(p1, 1.0)

    sig_count = 0
    for _ in range(n_sim):
        # Simulate gene presence/absence
        case_present = np.random.binomial(n_case, p1)
        control_present = np.random.binomial(n_control, p0)

        # Build 2x2 table
        table = np.array([
            [case_present, n_case - case_present],
            [control_present, n_control - control_present]
        ])

        _, pval = stats.fisher_exact(table)
        if pval < alpha:
            sig_count += 1

    return sig_count / n_sim

def find_min_detectable_or(n_case, n_control, gene_freq, target_power=0.80, alpha=0.05):
    """Binary search for minimum detectable OR at given power."""
    lo, hi = 1.0, 1000.0

    # First check if even OR=1000 gives enough power
    p = power_fisher(n_case, n_control, gene_freq, hi, alpha, n_sim=5000)
    if p < target_power:
        return float('inf')

    for _ in range(15):  # ~15 iterations of binary search
        mid = (lo + hi) / 2
        p = power_fisher(n_case, n_control, gene_freq, mid, alpha, n_sim=10000)
        if p < target_power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

# Designs
# designs = {
#     'A. niger (13 vs 15)': (13, 15),
#     'A. oryzae (31 vs 1)': (31, 1),
# }
# gene_freqs = [0.10, 0.20, 0.50]

# print("\nMinimum detectable Odds Ratio at 80% power, alpha=0.05:")
# print("-" * 70)
# print(f"{'Design':<25} {'Gene Freq':>10} {'Min OR':>15}")
# print("-" * 70)

# mde_results = []
# for design_name, (nc, nctrl) in designs.items():
#     for gf in gene_freqs:
#         min_or = find_min_detectable_or(nc, nctrl, gf)
#         or_str = f"{min_or:.1f}" if min_or < 999 else "INF (>1000)"
#         print(f"{design_name:<25} {gf:>10.0%} {or_str:>15}")
#         mde_results.append({
#             'analysis': 'minimum_detectable_OR',
#             'design': design_name,
#             'n_case': nc,
#             'n_control': nctrl,
#             'gene_freq': gf,
#             'min_detectable_OR': min_or if min_or < 999 else np.inf,
#             'target_power': 0.80,
#             'alpha': 0.05
#         })

# print("\nInterpretation:")
# print("  With only 1 control, no effect size is detectable — Fisher's exact test")
# print("  cannot distinguish signal from noise with a single observation in one group.")

# ============================================================================
# PART 2: Simulation-based power analysis
# ============================================================================
# print("\n" + "=" * 80)
# print("PART 2: SIMULATION-BASED POWER ANALYSIS (10,000 simulations each)")
# print("=" * 80)

def simulate_power(n_case, n_control, gene_freq, odds_ratio, n_sim=10000, alpha=0.05):
    """Simulate datasets and compute power as fraction with p < alpha."""
    p0 = gene_freq
    p1 = (odds_ratio * p0) / (1 - p0 + odds_ratio * p0)
    p1 = min(p1, 1.0)

    sig = 0
    for _ in range(n_sim):
        case_present = np.random.binomial(n_case, p1)
        control_present = np.random.binomial(n_control, p0)
        table = np.array([
            [case_present, n_case - case_present],
            [control_present, n_control - control_present]
        ])
        _, pval = stats.fisher_exact(table)
        if pval < alpha:
            sig += 1
    return sig / n_sim

# Three designs to compare


# =============================================================================
# SECTION 1.2: CREATE PHENOTYPE FILES (STRICT, BROAD, LENIENT)
# =============================================================================

def create_phenotype_files(samples, contrasts, output_dir):
    """
    Create binary phenotype files for each contrast.
    
    strict: high-confidence only
    broad: high + medium confidence
    lenient: all confidence levels (high + medium + low)
    """
    phenotype_files = {}
    
    for case_label, control_label, contrast_name in contrasts:
        for mode in ['strict', 'broad', 'lenient']:
            # Filter by confidence
            if mode == 'strict':
                conf_mask = samples['confidence'] == 'high'
            elif mode == 'broad':
                conf_mask = samples['confidence'].isin(['high', 'medium'])
            else:  # lenient - include all
                conf_mask = pd.Series([True] * len(samples), index=samples.index)
            
            # Filter by labels
            label_mask = samples['label_short'].isin([case_label, control_label])
            subset = samples[conf_mask & label_mask].copy()
            
            # Create binary phenotype (1 = case, 0 = control)
            subset['phenotype'] = (subset['label_short'] == case_label).astype(int)
            
            # Create phenotype DataFrame
            pheno_df = subset[['sample_id', 'phenotype']].copy()
            
            # Save
            filename = f"pheno_{contrast_name}_{mode}.tsv"
            filepath = os.path.join(output_dir, filename)
            pheno_df.to_csv(filepath, sep='\t', index=False)
            
            n_case = (pheno_df['phenotype'] == 1).sum()
            n_control = (pheno_df['phenotype'] == 0).sum()
            
            phenotype_files[f"{contrast_name}_{mode}"] = {
                'path': filepath,
                'n_case': n_case,
                'n_control': n_control,
                'n_total': len(pheno_df),
                'case_label': case_label,
                'control_label': control_label,
                'samples': pheno_df['sample_id'].tolist()
            }
            
            print(f"  {filename}: {case_label}={n_case}, {control_label}={n_control}, total={len(pheno_df)}")
    
    return phenotype_files

# Create output directory
# pheno_dir = f"{RESULTS_DIR}/phenotypes"
# os.makedirs(pheno_dir, exist_ok=True)

# Create phenotype files
# print("Creating phenotype files...")
# print("  strict = high confidence only")
# print("  broad = high + medium confidence")
# print("  lenient = all confidence levels (recommended for small sample sizes)")
# print()
# phenotype_files = create_phenotype_files(gwas_samples, contrasts, pheno_dir)
# print(f"\nCreated {len(phenotype_files)} phenotype files in {pheno_dir}")

# =============================================================================
# SECTION 2.2: SAMPLE ID NORMALIZATION (Handle format mismatches)
# =============================================================================

import re

def normalize_sample_id(sample_id):
    """
    Normalize sample ID to base accession format (e.g., GCA_000002855).
    
    Handles formats like:
    - GCA_000002855.2 (metadata format)
    - GCF_000002855.4_ASM285v2 (PCA/kinship format)
    - GCA_000230395.2_ASPNI_v3.0.proteins (OrthoFinder format)
    - GCA_000230395.2_ASPNI_v3.0.faa (FASTA format)
    """
    if not isinstance(sample_id, str):
        return str(sample_id)
    
    # Remove common suffixes
    s = sample_id.replace('.proteins', '').replace('.faa', '').replace('.fasta', '')
    
    # Extract GCA/GCF accession pattern (e.g., GCA_000002855 or GCF_000002855)
    match = re.match(r'(GC[AF]_\d+)', s)
    if match:
        return match.group(1)
    
    return s

def create_sample_id_mappings(phenotype_ids, pca_ids, feature_ids):
    """
    Create mappings between different sample ID formats.
    
    Returns:
    --------
    dict with:
        - 'pheno_to_pca': {pheno_id: pca_id}
        - 'pheno_to_feature': {pheno_id: feature_id}
        - 'common_samples': list of phenotype IDs with matches in both PCA and feature data
    """
    # Normalize all IDs
    pheno_normalized = {normalize_sample_id(s): s for s in phenotype_ids}
    pca_normalized = {normalize_sample_id(s): s for s in pca_ids}
    feature_normalized = {normalize_sample_id(s): s for s in feature_ids}
    
    # Create mappings
    pheno_to_pca = {}
    pheno_to_feature = {}
    
    for norm_id, pheno_id in pheno_normalized.items():
        if norm_id in pca_normalized:
            pheno_to_pca[pheno_id] = pca_normalized[norm_id]
        if norm_id in feature_normalized:
            pheno_to_feature[pheno_id] = feature_normalized[norm_id]
    
    # Find common samples (in phenotype, PCA, and feature matrices)
    common_samples = [s for s in phenotype_ids if s in pheno_to_pca and s in pheno_to_feature]
    
    return {
        'pheno_to_pca': pheno_to_pca,
        'pheno_to_feature': pheno_to_feature,
        'common_samples': common_samples
    }

# Test the normalization function
# test_ids = [
#     'GCA_000002855.2',
#     'GCF_000002855.4_ASM285v2',
#     'GCA_000230395.2_ASPNI_v3.0.proteins',
#     'GCA_001741905.1_ASM174190v1.proteins'
# ]
# print("Sample ID normalization test:")
# for s in test_ids:
#     print(f"  {s} -> {normalize_sample_id(s)}")

# Show PCA sample IDs for reference
# if pcs_df is not None:
#     pc_col = 'sample_id' if 'sample_id' in pcs_df.columns else pcs_df.columns[0]
#     pca_sample_ids = pcs_df[pc_col].tolist()
#     print(f"\nPCA samples: {len(pca_sample_ids)}")
#     print(f"  Example: {pca_sample_ids[0]} -> {normalize_sample_id(pca_sample_ids[0])}")

# Show phenotype sample IDs
# pheno_sample_ids = gwas_samples['sample_id'].tolist()
# print(f"\nPhenotype samples: {len(pheno_sample_ids)}")
# print(f"  Example: {pheno_sample_ids[0]} -> {normalize_sample_id(pheno_sample_ids[0])}")

# print("\nSample ID normalization functions loaded."
#       "\nFull mapping will be created after feature matrices are loaded.")

# =============================================================================
# SECTION 4: LINEAR MIXED MODEL (LMM) ASSOCIATION TESTING FRAMEWORK
# =============================================================================
# 
# Uses EMMA/EMMAX-style LMM when kinship matrix is available:
#   y = Xβ + u + ε,  where u ~ N(0, σ²_g K) and ε ~ N(0, σ²_e I)
#
# The kinship matrix K directly models genetic relatedness, making post-hoc
# genomic control (GC) correction redundant. This is the standard approach
# in microbial pan-GWAS (cf. pyseer, GEMMA, FaST-LMM).
#
# For binary phenotypes (case/control), using LMM is an approximation that
# is well-validated in GWAS literature — under the null, the score test from
# LMM is equivalent to that from a logistic mixed model.
#
# Falls back to logistic regression + PCs (or Fisher's exact) when kinship
# is not available.
# =============================================================================

from scipy import stats
from scipy.optimize import minimize_scalar
from statsmodels.stats.multitest import multipletests
import warnings
import time

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*Inverting hessian.*')
warnings.filterwarnings('ignore', message='.*Maximum Likelihood optimization failed.*')
warnings.filterwarnings('ignore', message='.*ConvergenceWarning.*')


# ---------------------------------------------------------------------------
# LMM core: eigendecomposition + REML + per-feature Wald test
# ---------------------------------------------------------------------------

def _eigendecompose_kinship(K):
    """Eigendecompose kinship matrix, clipping negative eigenvalues."""
    eigenvalues, eigenvectors = np.linalg.eigh(K)
    # Clip small/negative eigenvalues (numerical noise)
    eigenvalues = np.maximum(eigenvalues, 1e-10)
    return eigenvalues, eigenvectors


def _reml_loglik(log_delta, Uty, UtX, S):
    """
    Profile REML log-likelihood as a function of log(delta),
    where delta = sigma_e^2 / sigma_g^2.
    """
    delta = np.exp(log_delta)
    n = len(Uty)
    p = UtX.shape[1]
    
    D = S + delta  # diagonal of covariance in rotated space
    D_inv = 1.0 / D
    
    # Weighted least squares: beta = (X'V^-1 X)^-1 X'V^-1 y
    UtX_w = UtX * D_inv[:, None]
    XtVX = UtX.T @ UtX_w
    XtVy = UtX_w.T @ Uty
    
    try:
        beta = np.linalg.solve(XtVX, XtVy)
    except np.linalg.LinAlgError:
        return 1e10
    
    residuals = Uty - UtX @ beta
    RSS = np.sum(residuals**2 * D_inv)
    
    # Profile REML log-likelihood (up to constant)
    # L = -0.5 * [sum(log(D)) + log|X'V^-1 X| + (n-p)*log(RSS/(n-p)) + (n-p)]
    ll = -0.5 * np.sum(np.log(D))
    sign, logdet = np.linalg.slogdet(XtVX)
    if sign <= 0:
        return 1e10
    ll -= 0.5 * logdet
    ll -= 0.5 * (n - p) * np.log(RSS / (n - p))
    
    return -ll  # minimize negative log-likelihood


def fit_null_lmm(y, X, K, log_delta_bounds=(-10, 15), boundary_tol=0.05):
    """
    Fit null LMM (no feature effect) via REML.

    Returns
    -------
    S, U, delta_opt, Uty, UtX, h2, h2_status
        ``h2_status`` is one of:
          * ``'interior'``    -- REML found an interior optimum, h2 trustworthy
          * ``'upper_bound'`` -- log_delta hit upper wall, h2 clamped near 0
          * ``'lower_bound'`` -- log_delta hit lower wall, h2 clamped near 1
        Use ``h2_status`` (not ``h2``) when deciding whether to report h2.
    """
    S, U = _eigendecompose_kinship(K)
    Ut = U.T
    Uty = Ut @ y
    UtX = Ut @ X

    # Optimize delta = sigma_e^2 / sigma_g^2
    result = minimize_scalar(
        _reml_loglik, bounds=log_delta_bounds, method='bounded',
        args=(Uty, UtX, S)
    )
    log_delta = float(result.x)
    delta_opt = float(np.exp(log_delta))
    h2 = 1.0 / (1.0 + delta_opt)

    if abs(log_delta - log_delta_bounds[0]) < boundary_tol:
        h2_status = 'lower_bound'   # log_delta -> -inf, h2 -> 1
    elif abs(log_delta - log_delta_bounds[1]) < boundary_tol:
        h2_status = 'upper_bound'   # log_delta -> +inf, h2 -> 0
    else:
        h2_status = 'interior'

    return S, U, delta_opt, Uty, UtX, h2, h2_status


def test_feature_lmm(Uty, UtX, Utg, S, delta):
    """
    Wald test for a single feature in the rotated (eigendecomposed) space.
    
    Tests H0: beta_g = 0 in the model y = X*beta + g*beta_g + u + e.
    
    Returns: beta, se, pvalue
    """
    n = len(Uty)
    D = S + delta
    D_inv = 1.0 / D
    
    # Full model design matrix: [X, g] in rotated space
    UtXg = np.column_stack([UtX, Utg])
    p_full = UtXg.shape[1]
    
    # Weighted least squares
    UtXg_w = UtXg * D_inv[:, None]
    XgVXg = UtXg.T @ UtXg_w
    XgVy = UtXg_w.T @ Uty
    
    try:
        beta_full = np.linalg.solve(XgVXg, XgVy)
        XgVXg_inv = np.linalg.inv(XgVXg)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    
    # Residual variance
    residuals = Uty - UtXg @ beta_full
    sigma2 = np.sum(residuals**2 * D_inv) / (n - p_full)
    
    # Extract gene effect (last coefficient)
    beta_g = beta_full[-1]
    se_g = np.sqrt(max(sigma2 * XgVXg_inv[-1, -1], 0))
    
    if se_g == 0 or np.isnan(se_g):
        return beta_g, np.nan, np.nan
    
    # Wald test (chi-squared with 1 df)
    wald = (beta_g / se_g) ** 2
    pval = 1.0 - stats.chi2.cdf(wald, df=1)
    
    return beta_g, se_g, pval


# ---------------------------------------------------------------------------
# Main association function
# ---------------------------------------------------------------------------

def run_association_test(feature_matrix, phenotype_df, pcs_df=None, grm_df=None,
                         n_pcs=0, min_maf=0.05, use_kinship=True):
    """
    Run association testing for each feature against a binary phenotype.
    
    When grm_df (kinship matrix) is provided and use_kinship=True:
        Uses LMM with EMMA/EMMAX approach — eigendecomposes the kinship matrix,
        estimates variance components via REML under the null, then runs a Wald
        test for each feature. This directly models population structure from
        genetic relatedness, making post-hoc GC correction unnecessary.
    
    Fallback (no kinship):
        Uses logistic regression with PC covariates for large samples,
        or Fisher's exact test for small/imbalanced samples.
    """
    # Get sample IDs from each data source
    pheno_samples = phenotype_df['sample_id'].tolist()
    feature_samples = feature_matrix.columns.tolist()
    pc_col = 'sample_id' if 'sample_id' in pcs_df.columns else pcs_df.columns[0]
    pca_samples = pcs_df[pc_col].tolist()
    
    # Create normalized sample ID mappings
    mappings = create_sample_id_mappings(pheno_samples, pca_samples, feature_samples)
    common_samples = mappings['common_samples']
    pheno_to_pca = mappings['pheno_to_pca']
    pheno_to_feature = mappings['pheno_to_feature']
    
    # If kinship available, also intersect with kinship samples
    pheno_to_kinship = {}
    if grm_df is not None and use_kinship:
        kinship_ids = grm_df.index.tolist()
        kinship_normalized = {normalize_sample_id(s): s for s in kinship_ids}
        for pheno_id in common_samples:
            norm = normalize_sample_id(pheno_id)
            if norm in kinship_normalized:
                pheno_to_kinship[pheno_id] = kinship_normalized[norm]
        common_samples = [s for s in common_samples if s in pheno_to_kinship]
    
    n_samples = len(common_samples)
    print(f"  Sample matching:")
    print(f"    Phenotype samples: {len(pheno_samples)}")
    print(f"    PCA samples: {len(pca_samples)}")
    print(f"    Feature samples: {len(feature_samples)}")
    if grm_df is not None:
        print(f"    Kinship samples: {len(grm_df)}")
    print(f"    Common samples: {n_samples}")
    
    if n_samples == 0:
        print(f"  ERROR: No common samples found!")
        return pd.DataFrame()
    
    # Prepare phenotype vector
    pheno_dict = dict(zip(phenotype_df['sample_id'], phenotype_df['phenotype']))
    y = np.array([pheno_dict[s] for s in common_samples], dtype=float)
    n_cases = int(y.sum())
    n_controls = int(len(y) - y.sum())
    print(f"    Cases: {n_cases}, Controls: {n_controls}")
    
    # =========================================================================
    # LMM PATH (kinship available)
    # =========================================================================
    if grm_df is not None and use_kinship:
        print(f"\n  Using LMM with kinship matrix (EMMA/EMMAX approach)")
        
        # Build kinship submatrix for common samples
        kin_ids = [pheno_to_kinship[s] for s in common_samples]
        K = grm_df.loc[kin_ids, kin_ids].values.astype(float)
        
        # Build covariate matrix (intercept + PCs)
        pca_dict = {}
        for i in range(1, n_pcs + 1):
            pc_name = f'PC{i}'
            if pc_name in pcs_df.columns:
                pca_dict[pc_name] = dict(zip(pcs_df[pc_col], pcs_df[pc_name]))
        
        X = np.ones((n_samples, 1 + len(pca_dict)))
        for i, (pc_name, pc_values) in enumerate(pca_dict.items()):
            for j, pheno_sample in enumerate(common_samples):
                pca_sample = pheno_to_pca[pheno_sample]
                if pca_sample in pc_values:
                    X[j, i + 1] = pc_values[pca_sample]
        
        print(f"  Covariates: intercept + {len(pca_dict)} PCs")
        
        # Fit null model (REML)
        t0 = time.time()
        S, U, delta, Uty, UtX, h2, h2_status = fit_null_lmm(y, X, K)
        flag = '' if h2_status == 'interior' else f' [BOUNDARY: {h2_status}]'
        print(f"  Null model: delta={delta:.4f}, h²={h2:.4f}{flag} (time: {time.time()-t0:.1f}s)")
        
        Ut = U.T  # for rotating features
        
        # Test each feature
        results = []
        n_features = feature_matrix.shape[0]
        n_tested = 0
        t0 = time.time()
        
        for idx, (feature_name, row) in enumerate(feature_matrix.iterrows()):
            if idx % 2000 == 0 and idx > 0:
                elapsed = time.time() - t0
                rate = idx / elapsed
                print(f"    [{idx}/{n_features}] {rate:.0f}/s, ETA {(n_features-idx)/rate:.0f}s")
            
            # Get feature values for common samples
            x = np.array([row[pheno_to_feature[s]] for s in common_samples], dtype=float)
            
            # MAF filter (binarise for count data: any nonzero = present)
            x_binary = (x > 0).astype(float)
            maf = min(x_binary.mean(), 1 - x_binary.mean())
            if min_maf > 0 and maf < min_maf:
                continue
            
            n_tested += 1
            
            # Rotate feature into eigenspace
            Utg = Ut @ x
            
            # Wald test
            beta, se, pval = test_feature_lmm(Uty, UtX, Utg, S, delta)
            
            if np.isnan(pval):
                continue
            
            results.append({
                'feature': feature_name,
                'beta': beta,
                'se': se,
                'pvalue': pval,
                'n_samples': n_samples,
                'maf': maf,
                'n_case': n_cases,
                'n_control': n_controls,
                'method': 'lmm',
                'h2': h2,
                'h2_status': h2_status,
            })
        
        elapsed = time.time() - t0
        print(f"    Tested: {n_tested} features ({n_tested/max(elapsed,1):.0f}/s)")
    
    # =========================================================================
    # FALLBACK PATH (no kinship — logistic regression or Fisher's exact)
    # =========================================================================
    else:
        print(f"\n  No kinship matrix — falling back to logistic/Fisher")
        
        pca_dict = {}
        for i in range(1, n_pcs + 1):
            pc_name = f'PC{i}'
            if pc_name in pcs_df.columns:
                pca_dict[pc_name] = dict(zip(pcs_df[pc_col], pcs_df[pc_name]))
        
        covariates = np.ones((n_samples, 1 + len(pca_dict)))
        for i, (pc_name, pc_values) in enumerate(pca_dict.items()):
            for j, pheno_sample in enumerate(common_samples):
                pca_sample = pheno_to_pca[pheno_sample]
                if pca_sample in pc_values:
                    covariates[j, i + 1] = pc_values[pca_sample]
        
        results = []
        n_features = feature_matrix.shape[0]
        n_tested = 0
        n_fisher = 0
        n_logistic = 0
        
        for idx, (feature_name, row) in enumerate(feature_matrix.iterrows()):
            if idx % 2000 == 0 and idx > 0:
                print(f"    Processed {idx}/{n_features} features...")
            
            x = np.array([row[pheno_to_feature[s]] for s in common_samples])
            # MAF filter (binarise for count data: any nonzero = present)
            x_binary = (x > 0).astype(float)
            maf = min(x_binary.mean(), 1 - x_binary.mean())
            if min_maf > 0 and maf < min_maf:
                continue
            
            n_tested += 1
            is_binary = set(np.unique(x)).issubset({0, 1})
            use_fisher = is_binary and (n_samples < 30 or n_cases < 10 or n_controls < 10)
            
            if use_fisher:
                try:
                    table = [
                        [int(((y == 1) & (x > 0)).sum()), int(((y == 1) & (x == 0)).sum())],
                        [int(((y == 0) & (x > 0)).sum()), int(((y == 0) & (x == 0)).sum())]
                    ]
                    odds_ratio, pval = stats.fisher_exact(table)
                    beta = np.log(odds_ratio) if odds_ratio > 0 else 0
                    n_fisher += 1
                    results.append({
                        'feature': feature_name, 'beta': beta, 'se': np.nan,
                        'pvalue': pval, 'n_samples': n_samples, 'maf': maf,
                        'n_case': n_cases, 'n_control': n_controls, 'method': 'fisher',
                    })
                except:
                    continue
            else:
                try:
                    from statsmodels.api import Logit
                    X = np.column_stack([covariates, x])
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        result = Logit(y, X).fit(disp=0, method='bfgs', maxiter=100,
                                                  warn_convergence=False)
                    pval = result.pvalues[-1]
                    if np.isnan(pval):
                        continue
                    n_logistic += 1
                    results.append({
                        'feature': feature_name, 'beta': result.params[-1],
                        'se': result.bse[-1], 'pvalue': pval,
                        'n_samples': n_samples, 'maf': maf,
                        'n_case': n_cases, 'n_control': n_controls, 'method': 'logistic',
                    })
                except:
                    continue
        
        print(f"    Tested: {n_tested} features (Fisher: {n_fisher}, Logistic: {n_logistic})")
    
    # =========================================================================
    # FDR correction and return
    # =========================================================================
    # Note: FDR/q-value computation is handled centrally in
    # run_all_associations (which writes a single 'qvalue' column). We do not
    # add an 'fdr' column here to avoid the redundant fdr+qvalue pair.
    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        results_df = results_df.sort_values('pvalue')

    return results_df

# print("Association testing framework loaded (LMM with kinship / Fisher+Logistic fallback).")

# =============================================================================
# SECTION 4.2: DIAGNOSTIC PLOTTING FUNCTIONS
# =============================================================================

def calculate_genomic_inflation(pvalues):
    """Calculate genomic inflation factor (lambda)."""
    pvalues = np.array(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues > 0]
    
    chi2_observed = stats.chi2.ppf(1 - pvalues, df=1)
    lambda_gc = np.median(chi2_observed) / stats.chi2.ppf(0.5, df=1)
    return lambda_gc

def plot_qq(pvalues, title="QQ Plot", ax=None):
    """Create QQ plot of p-values."""
    pvalues = np.array(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues > 0]
    pvalues = np.sort(pvalues)
    
    n = len(pvalues)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(pvalues)
    
    lambda_gc = calculate_genomic_inflation(pvalues)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    
    # Plot diagonal line
    max_val = max(max(expected), max(observed))
    ax.plot([0, max_val], [0, max_val], 'r--', lw=1, label='Expected')
    
    # Plot observed vs expected
    ax.scatter(expected, observed, s=10, alpha=0.5, c='blue')
    
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')
    ax.set_title(f'{title}\nλ = {lambda_gc:.3f}')
    ax.legend()
    
    return lambda_gc

def plot_pvalue_histogram(pvalues, title="P-value Distribution", ax=None):
    """Create histogram of p-values."""
    pvalues = np.array(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    
    ax.hist(pvalues, bins=50, edgecolor='black', alpha=0.7)
    ax.axhline(y=len(pvalues)/50, color='r', linestyle='--', label='Uniform expectation')
    ax.set_xlabel('P-value')
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.legend()

def plot_manhattan(results_df, title="Manhattan Plot", fdr_threshold=0.05, ax=None):
    """Create Manhattan-style plot for association results."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    
    results = results_df.copy()
    results['neg_log_p'] = -np.log10(results['pvalue'])
    results['idx'] = range(len(results))

    # Color by significance — accept either 'qvalue' (current schema) or
    # legacy 'fdr' for backwards compatibility.
    fdr_col = 'qvalue' if 'qvalue' in results.columns else 'fdr'
    sig_mask = results[fdr_col] < fdr_threshold
    
    ax.scatter(results.loc[~sig_mask, 'idx'], results.loc[~sig_mask, 'neg_log_p'], 
               s=10, alpha=0.5, c='gray', label='Not significant')
    ax.scatter(results.loc[sig_mask, 'idx'], results.loc[sig_mask, 'neg_log_p'], 
               s=20, alpha=0.8, c='red', label=f'FDR < {fdr_threshold}')
    
    # Add significance threshold line
    if sig_mask.any():
        min_sig_p = results.loc[sig_mask, 'pvalue'].max()
        ax.axhline(y=-np.log10(min_sig_p), color='red', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Feature index')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(title)
    ax.legend()

def create_diagnostic_plots(pvals_or_df, title=None, contrast_name=None, output_dir=None):
    """Create diagnostic plots (QQ + p-value histogram + optional Manhattan).

    Flexible input:
      - ``pvals_or_df`` can be a pandas Series/array of p-values OR a
        DataFrame with a ``pvalue`` column (and optional ``feature`` index
        for the Manhattan plot).
      - ``title`` (or legacy ``contrast_name``) is used in the figure titles.
      - If ``output_dir`` is given, saves ``{contrast_name}_diagnostics.png``.

    Returns
    -------
    matplotlib.Figure
    """
    import numpy as _np
    if title is None:
        title = contrast_name or ''

    if isinstance(pvals_or_df, pd.DataFrame):
        results_df = pvals_or_df
        pvals = _np.asarray(results_df['pvalue'].dropna().values)
    else:
        results_df = None
        pvals = _np.asarray(pd.Series(pvals_or_df).dropna().values)

    if results_df is not None and len(results_df) > 0:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    plot_qq(pvals, title=f'QQ Plot: {title}', ax=axes[0])
    plot_pvalue_histogram(pvals, title=f'P-values: {title}', ax=axes[1])
    if results_df is not None and len(axes) > 2:
        plot_manhattan(results_df, title=f'Manhattan: {title}', ax=axes[2])

    plt.tight_layout()
    if output_dir is not None and contrast_name is not None:
        fig.savefig(f'{output_dir}/{contrast_name}_diagnostics.png',
                    dpi=150, bbox_inches='tight')
    return fig

# print("Diagnostic plotting functions loaded.")

# =============================================================================
# SECTION 6.4: CAZy FAMILY EXPANSION PAN-GWAS
# =============================================================================
# Test if phenotype groups differ in CAZyme repertoire SIZE (not individual families)
# Hypothesis: Industrial/pathogenic strains may have expanded specific enzyme classes
#
# This tests AGGREGATE counts - biologically meaningful genome-level traits:
#   - Total CAZymes per genome
#   - Count per CAZy class (GH, GT, PL, CE, AA, CBM)

import re
from scipy.stats import mannwhitneyu
from collections import defaultdict

# print("="*70)
# print("SECTION 6.4: CAZy FAMILY EXPANSION PAN-GWAS")
# print("="*70)

# Load dbCAN data and compute aggregate counts per genome
# dbcan_dir = f"{SPECIES_DIR}/dbcan_output"
# cazy_counts = {}

# if os.path.exists(dbcan_dir):
#     print(f"\nLoading CAZy annotations from dbCAN...")
#     
#     for acc_dir in os.listdir(dbcan_dir):
#         overview_path = f"{dbcan_dir}/{acc_dir}/overview.tsv"
#         if os.path.exists(overview_path):
#             try:
#                 df = pd.read_csv(overview_path, sep='\t', dtype=str, na_filter=False)
#                 acc_match = re.search(r'(GCA_\d+)', acc_dir)
#                 if acc_match:
#                     acc = acc_match.group(1)
#                     
                    # Initialize counts
#                     counts = {'total': 0, 'GH': 0, 'GT': 0, 'PL': 0, 'CE': 0, 'AA': 0, 'CBM': 0, 'other': 0}
#                     
                    # Get CAZy families from "Recommend Results" column
#                     rec_col = [c for c in df.columns if 'Recommend' in c]
#                     if rec_col:
#                         for val in df[rec_col[0]]:
#                             if val and val != '-':
                                # Count each CAZy annotation
#                                 families = re.findall(r'([A-Z]+)\d+', str(val))
#                                 for fam_class in families:
#                                     counts['total'] += 1
#                                     if fam_class in counts:
#                                         counts[fam_class] += 1
#                                     else:
#                                         counts['other'] += 1
#                     
#                     cazy_counts[acc] = counts
#             except Exception as e:
#                 pass
# 
#     print(f"  Loaded CAZy data for {len(cazy_counts)} genomes")
#     
#     if cazy_counts:
        # Build dataframe
#         cazy_df = pd.DataFrame.from_dict(cazy_counts, orient='index')
#         cazy_df.index.name = 'accession'
#         
#         print(f"\n  CAZy class summary (mean ± std per genome):")
#         for col in ['total', 'GH', 'GT', 'PL', 'CE', 'AA', 'CBM']:
#             if col in cazy_df.columns:
#                 print(f"    {col}: {cazy_df[col].mean():.1f} ± {cazy_df[col].std():.1f}")
#         
        # Normalize sample IDs for matching
#         cazy_df_norm = cazy_df.copy()
#         cazy_df_norm.index = [normalize_sample_id(str(idx)) for idx in cazy_df_norm.index]
#         
        # Run association tests
#         print("\n" + "-"*70)
#         print("CAZy Expansion Association Testing (Mann-Whitney U)")
#         print("-"*70)
#         
#         cazy_expansion_results = {}
#         
#         for contrast_key, contrast_info in phenotype_files.items():
#             n_case = contrast_info['n_case']
#             n_control = contrast_info['n_control']
#             
#             if n_case < 3 or n_control < 3:
#                 continue
#             
#             print(f"\n  {contrast_key} ({n_case} cases, {n_control} controls)")
#             
#             pheno_df = pd.read_csv(contrast_info['path'], sep='\t')
#             pheno_df['sample_id_norm'] = pheno_df['sample_id'].apply(lambda x: normalize_sample_id(str(x)))
#             
#             cases = pheno_df[pheno_df['phenotype'] == 1]['sample_id_norm'].tolist()
#             controls = pheno_df[pheno_df['phenotype'] == 0]['sample_id_norm'].tolist()
#             
#             cases_matched = [s for s in cases if s in cazy_df_norm.index]
#             controls_matched = [s for s in controls if s in cazy_df_norm.index]
#             
#             if len(cases_matched) < 3 or len(controls_matched) < 3:
#                 print(f"    [SKIP] Insufficient matched samples")
#                 continue
#             
#             results_list = []
#             
#             for feature in ['total', 'GH', 'GT', 'PL', 'CE', 'AA', 'CBM']:
#                 if feature not in cazy_df_norm.columns:
#                     continue
#                     
#                 case_vals = cazy_df_norm.loc[cases_matched, feature].values
#                 control_vals = cazy_df_norm.loc[controls_matched, feature].values
#                 
#                 try:
#                     stat, pval = mannwhitneyu(case_vals, control_vals, alternative='two-sided')
#                 except:
#                     pval = 1.0
#                 
#                 case_mean = np.mean(case_vals)
#                 control_mean = np.mean(control_vals)
#                 
#                 results_list.append({
#                     'feature': f'CAZy_{feature}',
#                     'case_mean': case_mean,
#                     'control_mean': control_mean,
#                     'diff': case_mean - control_mean,
#                     'fold_change': case_mean / control_mean if control_mean > 0 else np.nan,
#                     'pvalue': pval,
#                     'direction': 'expanded_in_case' if case_mean > control_mean else 'reduced_in_case'
#                 })
#             
#             results = pd.DataFrame(results_list)
#             
            # FDR correction
#             if len(results) > 0:
#                 from statsmodels.stats.multitest import multipletests
#                 _, fdr, _, _ = multipletests(results['pvalue'], method='fdr_bh')
#                 results['fdr'] = fdr
#             
#             output_file = f"{results_dir}/cazy_expansion_assoc_{contrast_key}.tsv"
#             results.to_csv(output_file, sep='\t', index=False)
#             cazy_expansion_results[contrast_key] = results
#             
            # Print significant results
#             for _, row in results.iterrows():
#                 sig = "*" if row['pvalue'] < 0.05 else ""
#                 sig2 = "**" if row['fdr'] < 0.1 else ""
#                 print(f"    {row['feature']}: case={row['case_mean']:.1f}, ctrl={row['control_mean']:.1f}, "
#                       f"p={row['pvalue']:.4f}{sig}{sig2}")
# else:
#     print("\n[SKIP] dbCAN output directory not found")
#     cazy_expansion_results = {}

# print("\n" + "="*70)
# print("CAZy Expansion Pan-GWAS complete!")
# =============================================================================
# AUTOMATIC INTERPRETATION: CAZy EXPANSION RESULTS
# =============================================================================

# print("\n" + "="*70)
# print("INTERPRETATION: CAZy FAMILY EXPANSION ANALYSIS")
# print("="*70)

def interpret_cazy(feature):
    interp = {
        'GH': 'Glycoside hydrolases → polysaccharide degradation',
        'AA': 'Auxiliary activities → lignin/cellulose oxidation',
        'GT': 'Glycosyl transferases → cell wall biosynthesis', 
        'CBM': 'Carbohydrate-binding modules → substrate targeting',
        'CE': 'Carbohydrate esterases → hemicellulose processing',
        'PL': 'Polysaccharide lyases → pectin degradation',
        'total': 'Overall CAZyme repertoire size'
    }
    for key, val in interp.items():
        if key in feature:
            return val
    return ''

# if cazy_expansion_results:
#     for contrast_key, results in cazy_expansion_results.items():
#         print(f"\n{'='*60}")
#         print(f"Contrast: {contrast_key}")
#         print(f"{'='*60}")
#         
        # Sort by p-value
#         results_sorted = results.sort_values('pvalue')
#         
        # Significant at p < 0.05
#         sig_05 = results_sorted[results_sorted['pvalue'] < 0.05]
        # Trending at p < 0.1
#         trending = results_sorted[(results_sorted['pvalue'] >= 0.05) & (results_sorted['pvalue'] < 0.1)]
#         
#         if len(sig_05) > 0:
#             print(f"\n  ★ SIGNIFICANT FINDINGS (p < 0.05):")
#             for _, row in sig_05.iterrows():
#                 direction = "EXPANDED" if row['case_mean'] > row['control_mean'] else "REDUCED"
#                 fold = row['case_mean'] / row['control_mean'] if row['control_mean'] > 0 else 0
#                 print(f"\n    {row['feature']}: {direction} in case phenotype")
#                 print(f"      Case mean: {row['case_mean']:.1f}")
#                 print(f"      Control mean: {row['control_mean']:.1f}")
#                 print(f"      Fold change: {fold:.2f}x")
#                 print(f"      P-value: {row['pvalue']:.4f}")
#                 bio = interpret_cazy(row['feature'])
#                 if bio:
#                     print(f"      Biology: {bio}")
#         
#         if len(trending) > 0:
#             print(f"\n  ○ TRENDING (0.05 ≤ p < 0.10):")
#             for _, row in trending.iterrows():
#                 direction = "expanded" if row['case_mean'] > row['control_mean'] else "reduced"
#                 print(f"    • {row['feature']}: {direction} (p={row['pvalue']:.4f})")
#         
#         if len(sig_05) == 0 and len(trending) == 0:
#             print(f"\n  No significant or trending associations (p < 0.1)")
#         
        # Overall pattern
#         print(f"\n  SUMMARY:")
#         n_expanded = (results['case_mean'] > results['control_mean']).sum()
#         n_reduced = (results['case_mean'] < results['control_mean']).sum()
#         print(f"    {n_expanded} CAZy features higher in cases, {n_reduced} higher in controls")
# else:
#     print("\nNo CAZy results available")


# =============================================================================
# SECTION 6.5: PROTEASE FAMILY EXPANSION PAN-GWAS
# =============================================================================
# Test if phenotype groups differ in protease repertoire SIZE
# Hypothesis: Pathogenic strains may have expanded protease arsenals for host invasion
#
# This tests AGGREGATE counts - biologically meaningful genome-level traits:
#   - Total proteases per genome
#   - Count per MEROPS class (Serine, Aspartic, Metallo, Cysteine, Threonine)

# print("="*70)
# print("SECTION 6.5: PROTEASE FAMILY EXPANSION PAN-GWAS")
# print("="*70)

# Load InterProScan data and count proteases by class
# interpro_dir = f"{SPECIES_DIR}/interproscan_output"
# protease_counts = {}

# MEROPS protease class patterns
# PROTEASE_PATTERNS = {
#     'S': [r'[Ss]erine[- ]?(?:protease|peptidase|endopeptidase)', r'[Pp]eptidase[_ ]S\d+', r'[Ss]ubtilis'],
#     'A': [r'[Aa]spartic[- ]?(?:protease|peptidase)', r'[Pp]eptidase[_ ]A\d+', r'[Pp]epsin'],
#     'M': [r'[Mm]etallo[- ]?(?:protease|peptidase)', r'[Pp]eptidase[_ ]M\d+', r'[Zz]inc[- ]?peptidase'],
#     'C': [r'[Cc]ysteine[- ]?(?:protease|peptidase)', r'[Pp]eptidase[_ ]C\d+', r'[Pp]apain'],
#     'T': [r'[Tt]hreonine[- ]?(?:protease|peptidase)', r'[Pp]eptidase[_ ]T\d+', r'[Pp]roteasome'],
# }

# if os.path.exists(interpro_dir):
#     print(f"\nLoading protease annotations from InterProScan...")
#     
#     for acc_dir in os.listdir(interpro_dir):
#         tsv_files = [f for f in os.listdir(f"{interpro_dir}/{acc_dir}") if f.endswith('.tsv')]
#         if tsv_files:
#             tsv_path = f"{interpro_dir}/{acc_dir}/{tsv_files[0]}"
#             try:
#                 df = pd.read_csv(tsv_path, sep='\t', header=None, dtype=str, na_filter=False)
#                 
#                 acc_match = re.search(r'(GCA_\d+)', acc_dir)
#                 if acc_match:
#                     acc = acc_match.group(1)
#                     
                    # Initialize counts
#                     counts = {'total': 0, 'Serine': 0, 'Aspartic': 0, 'Metallo': 0, 'Cysteine': 0, 'Threonine': 0, 'Other': 0}
#                     counted_proteins = set()  # Avoid double-counting same protein
#                     
#                     for _, row in df.iterrows():
#                         protein_id = str(row[0]) if len(row) > 0 else ""
#                         desc = str(row[5]) if len(row) > 5 else ""
#                         ipr_desc = str(row[12]) if len(row) > 12 else ""
#                         combined = f"{desc} {ipr_desc}"
#                         
                        # Check if this is a protease
#                         is_protease = False
#                         protease_class = None
#                         
#                         for pclass, patterns in PROTEASE_PATTERNS.items():
#                             for pattern in patterns:
#                                 if re.search(pattern, combined, re.IGNORECASE):
#                                     is_protease = True
#                                     protease_class = pclass
#                                     break
#                             if is_protease:
#                                 break
#                         
                        # Also check for generic peptidase mentions
#                         if not is_protease and re.search(r'[Pp]eptidase|[Pp]rotease|[Pp]roteinase', combined):
#                             is_protease = True
#                             protease_class = 'Other'
#                         
#                         if is_protease and protein_id not in counted_proteins:
#                             counted_proteins.add(protein_id)
#                             counts['total'] += 1
#                             
#                             class_names = {'S': 'Serine', 'A': 'Aspartic', 'M': 'Metallo', 
#                                           'C': 'Cysteine', 'T': 'Threonine'}
#                             if protease_class in class_names:
#                                 counts[class_names[protease_class]] += 1
#                             else:
#                                 counts['Other'] += 1
#                     
#                     if counts['total'] > 0:
#                         protease_counts[acc] = counts
#             except Exception as e:
#                 pass
# 
#     print(f"  Loaded protease data for {len(protease_counts)} genomes")
#     
#     if protease_counts:
        # Build dataframe
#         protease_df = pd.DataFrame.from_dict(protease_counts, orient='index')
#         protease_df.index.name = 'accession'
#         
#         print(f"\n  Protease class summary (mean ± std per genome):")
#         for col in ['total', 'Serine', 'Aspartic', 'Metallo', 'Cysteine', 'Threonine']:
#             if col in protease_df.columns:
#                 print(f"    {col}: {protease_df[col].mean():.1f} ± {protease_df[col].std():.1f}")
#         
        # Normalize sample IDs
#         protease_df_norm = protease_df.copy()
#         protease_df_norm.index = [normalize_sample_id(str(idx)) for idx in protease_df_norm.index]
#         
        # Run association tests
#         print("\n" + "-"*70)
#         print("Protease Expansion Association Testing (Mann-Whitney U)")
#         print("-"*70)
#         
#         protease_expansion_results = {}
#         
#         for contrast_key, contrast_info in phenotype_files.items():
#             n_case = contrast_info['n_case']
#             n_control = contrast_info['n_control']
#             
#             if n_case < 3 or n_control < 3:
#                 continue
#             
#             print(f"\n  {contrast_key} ({n_case} cases, {n_control} controls)")
#             
#             pheno_df = pd.read_csv(contrast_info['path'], sep='\t')
#             pheno_df['sample_id_norm'] = pheno_df['sample_id'].apply(lambda x: normalize_sample_id(str(x)))
#             
#             cases = pheno_df[pheno_df['phenotype'] == 1]['sample_id_norm'].tolist()
#             controls = pheno_df[pheno_df['phenotype'] == 0]['sample_id_norm'].tolist()
#             
#             cases_matched = [s for s in cases if s in protease_df_norm.index]
#             controls_matched = [s for s in controls if s in protease_df_norm.index]
#             
#             if len(cases_matched) < 3 or len(controls_matched) < 3:
#                 print(f"    [SKIP] Insufficient matched samples")
#                 continue
#             
#             results_list = []
#             
#             for feature in ['total', 'Serine', 'Aspartic', 'Metallo', 'Cysteine', 'Threonine']:
#                 if feature not in protease_df_norm.columns:
#                     continue
#                     
#                 case_vals = protease_df_norm.loc[cases_matched, feature].values
#                 control_vals = protease_df_norm.loc[controls_matched, feature].values
#                 
#                 try:
#                     stat, pval = mannwhitneyu(case_vals, control_vals, alternative='two-sided')
#                 except:
#                     pval = 1.0
#                 
#                 case_mean = np.mean(case_vals)
#                 control_mean = np.mean(control_vals)
#                 
#                 results_list.append({
#                     'feature': f'Protease_{feature}',
#                     'case_mean': case_mean,
#                     'control_mean': control_mean,
#                     'diff': case_mean - control_mean,
#                     'fold_change': case_mean / control_mean if control_mean > 0 else np.nan,
#                     'pvalue': pval,
#                     'direction': 'expanded_in_case' if case_mean > control_mean else 'reduced_in_case'
#                 })
#             
#             results = pd.DataFrame(results_list)
#             
            # FDR correction
#             if len(results) > 0:
#                 from statsmodels.stats.multitest import multipletests
#                 _, fdr, _, _ = multipletests(results['pvalue'], method='fdr_bh')
#                 results['fdr'] = fdr
#             
#             output_file = f"{results_dir}/protease_expansion_assoc_{contrast_key}.tsv"
#             results.to_csv(output_file, sep='\t', index=False)
#             protease_expansion_results[contrast_key] = results
#             
            # Print results
#             for _, row in results.iterrows():
#                 sig = "*" if row['pvalue'] < 0.05 else ""
#                 sig2 = "**" if row['fdr'] < 0.1 else ""
#                 print(f"    {row['feature']}: case={row['case_mean']:.1f}, ctrl={row['control_mean']:.1f}, "
#                       f"p={row['pvalue']:.4f}{sig}{sig2}")
# else:
#     print("\n[SKIP] InterProScan output directory not found")
#     protease_expansion_results = {}

# print("\n" + "="*70)
# print("Protease Expansion Pan-GWAS complete!")
# =============================================================================
# AUTOMATIC INTERPRETATION: PROTEASE EXPANSION RESULTS
# =============================================================================

# print("\n" + "="*70)
# print("INTERPRETATION: PROTEASE FAMILY EXPANSION ANALYSIS")
# print("="*70)

def interpret_protease(feature):
    interp = {
        'Serine': 'Serine proteases → secreted enzymes, extracellular degradation',
        'Metallo': 'Metalloproteases → collagen degradation, virulence',
        'Aspartic': 'Aspartic proteases → food processing, host proteins',
        'Cysteine': 'Cysteine proteases → intracellular turnover',
        'Threonine': 'Threonine proteases → proteasome activity',
        'total': 'Overall protease repertoire size'
    }
    for key, val in interp.items():
        if key in feature:
            return val
    return ''

# if protease_expansion_results:
#     for contrast_key, results in protease_expansion_results.items():
#         print(f"\n{'='*60}")
#         print(f"Contrast: {contrast_key}")
#         print(f"{'='*60}")
#         
#         results_sorted = results.sort_values('pvalue')
#         sig_05 = results_sorted[results_sorted['pvalue'] < 0.05]
#         trending = results_sorted[(results_sorted['pvalue'] >= 0.05) & (results_sorted['pvalue'] < 0.1)]
#         
#         if len(sig_05) > 0:
#             print(f"\n  ★ SIGNIFICANT FINDINGS (p < 0.05):")
#             for _, row in sig_05.iterrows():
#                 direction = "EXPANDED" if row['case_mean'] > row['control_mean'] else "REDUCED"
#                 fold = row['case_mean'] / row['control_mean'] if row['control_mean'] > 0 else 0
#                 print(f"\n    {row['feature']}: {direction} in case phenotype")
#                 print(f"      Case: {row['case_mean']:.1f}, Control: {row['control_mean']:.1f}")
#                 print(f"      Fold: {fold:.2f}x, P={row['pvalue']:.4f}")
#                 bio = interpret_protease(row['feature'])
#                 if bio:
#                     print(f"      → {bio}")
#         
#         if len(trending) > 0:
#             print(f"\n  ○ TRENDING (0.05 ≤ p < 0.10):")
#             for _, row in trending.iterrows():
#                 direction = "expanded" if row['case_mean'] > row['control_mean'] else "reduced"
#                 print(f"    • {row['feature']}: {direction} (p={row['pvalue']:.4f})")
#         
#         if len(sig_05) == 0 and len(trending) == 0:
#             print(f"\n  No significant associations (p < 0.1)")
#         
#         print(f"\n  BIOLOGICAL CONTEXT:")
#         print(f"    • Protease expansion is a hallmark of pathogenic fungi")
#         print(f"    • Secreted proteases enable host tissue invasion")
#         print(f"    • Industrial strains may expand proteases for protein hydrolysis")
# else:
#     print("\nNo protease results available")


# =============================================================================
# SECTION 6.7: BIOLOGICAL INTERPRETATION OF BGC/GCF RESULTS
# =============================================================================
# Annotate significant BGC/GCF hits with detailed biosynthetic pathway information

from pathlib import Path
from collections import Counter, defaultdict

def load_bgc_annotations(species_dir):
    """Load BGC annotations from BiG-SCAPE output."""
    bigscape_dir = Path(species_dir) / "bigscape_output" / "output_files"
    
    output_dirs = list(bigscape_dir.glob("*_c0*"))
    if not output_dirs:
        print("No BiG-SCAPE output found")
        return None, None
    
    output_dir = output_dirs[0]
    
    # Load record annotations (BGC details)
    annot_path = output_dir / "record_annotations.tsv"
    if annot_path.exists():
        annot_df = pd.read_csv(annot_path, sep='\t')
        n_total = len(annot_df)
        
        # Filter out MIBiG reference BGCs (their Record starts with "BGC")
        is_actual = annot_df['Record'].str.startswith(('GCA', 'GCF'), na=False)
        n_mibig = (~is_actual).sum()
        annot_df = annot_df[is_actual].copy()
        
        if n_mibig > 0:
            print(f"Filtered out {n_mibig} MIBiG reference BGCs")
        print(f"Loaded {len(annot_df)} BGC annotations from actual strains")
    else:
        print(f"Annotation file not found: {annot_path}")
        return None, None
    
    # Load all clustering files to map BGCs to GCFs
    gcf_mapping = {}
    for cluster_file in output_dir.glob("*/*_clustering_*.tsv"):
        df = pd.read_csv(cluster_file, sep='\t')
        for _, row in df.iterrows():
            record = row['Record']
            # Skip MIBiG references
            if not str(record).startswith(('GCA', 'GCF')):
                continue
            family = row['Family']
            gcf_mapping[record] = family
    
    print(f"Mapped {len(gcf_mapping)} BGCs to GCFs")
    
    return annot_df, gcf_mapping

def parse_knownclusterblast_file(filepath):
    """Parse a knownclusterblast txt file and extract MIBiG hits."""
    import re
    hits = []
    with open(filepath, 'r') as f:
        content = f.read()
    
    sig_match = re.search(r'Significant hits:\s*\n(.*?)(?:\n\n|Details:)', content, re.DOTALL)
    if not sig_match:
        return hits
    
    sig_section = sig_match.group(1)
    for line in sig_section.strip().split('\n'):
        match = re.match(r'\d+\.\s+(BGC\d+)(?:\.\d+)?\s+(.+)', line.strip())
        if match:
            bgc_id = match.group(1)
            compound = match.group(2).strip()
            hits.append((bgc_id, compound))
    return hits

def build_mibig_annotation_table(species_dir):
    """
    Build MIBiG annotation table by parsing antiSMASH knownclusterblast results.
    Links antiSMASH regions to BiG-SCAPE GCFs via compound predictions.
    
    Returns:
    - gcf_mibig: dict mapping GCF -> list of (mibig_id, compound_name)
    """
    import re
    antismash_dir = Path(species_dir) / "antismash_output"
    bigscape_dir = Path(species_dir) / "bigscape_output" / "output_files"
    
    output_dirs = list(bigscape_dir.glob("*_c0*"))
    if not output_dirs:
        print("No BiG-SCAPE output found for MIBiG annotation")
        return {}
    
    output_dir = output_dirs[0]
    
    # Build GCF mapping from BiG-SCAPE clustering
    gcf_mapping = {}
    for cluster_file in output_dir.glob("*/*_clustering_*.tsv"):
        df = pd.read_csv(cluster_file, sep='\t')
        for _, row in df.iterrows():
            if str(row['Record']).startswith(('GCA', 'GCF')):
                gcf_mapping[row['Record']] = row['Family']
    
    # Parse all knownclusterblast results
    bgc_mibig = defaultdict(list)
    
    for strain_dir in antismash_dir.iterdir():
        if not strain_dir.is_dir():
            continue
        
        strain_id = strain_dir.name
        kcb_dir = strain_dir / 'knownclusterblast'
        
        if not kcb_dir.exists():
            continue
        
        for txt_file in kcb_dir.glob('*.txt'):
            match = re.match(r'(scaffold_\d+)_c(\d+)', txt_file.stem)
            if not match:
                continue
            
            scaffold = match.group(1)
            region_num = int(match.group(2))
            
            # Match to BiG-SCAPE records
            patterns = [
                f"{strain_id}_{scaffold}.region{region_num:03d}",
                f"{strain_id}_{scaffold}.region{region_num:02d}",
                f"{strain_id}_{scaffold}.region{region_num}"
            ]
            
            matching_records = []
            for pattern in patterns:
                matching_records.extend([r for r in gcf_mapping.keys() if pattern in r])
            
            if matching_records:
                hits = parse_knownclusterblast_file(txt_file)
                for record in set(matching_records):
                    bgc_mibig[record].extend(hits)
    
    # Aggregate at GCF level
    gcf_mibig = defaultdict(list)
    for bgc_record, hits in bgc_mibig.items():
        if bgc_record in gcf_mapping:
            gcf = gcf_mapping[bgc_record]
            gcf_mibig[gcf].extend(hits)
    
    # Deduplicate - keep unique MIBiG IDs with their compound names
    for gcf in gcf_mibig:
        seen = {}
        for bgc_id, compound in gcf_mibig[gcf]:
            if bgc_id not in seen:
                seen[bgc_id] = compound
        gcf_mibig[gcf] = [(k, v) for k, v in seen.items()]
    
    print(f"Loaded MIBiG annotations: {len([g for g in gcf_mibig if gcf_mibig[g]])} GCFs with known compound hits")
    
    return dict(gcf_mibig)

def build_gcf_annotation_table(bgc_annot_df, gcf_mapping, phenotype_df=None, gcf_mibig=None):
    """
    Build detailed GCF annotation table with:
    - GCF family ID
    - Number of BGCs in family
    - Number of strains with the family
    - BGC types (Class, Category)
    - Strain list
    - Phenotype distribution (if provided)
    """
    if bgc_annot_df is None:
        return pd.DataFrame()
    
    # Build reverse mapping: GCF -> list of BGC records
    gcf_to_bgcs = defaultdict(list)
    for record, family in gcf_mapping.items():
        gcf_to_bgcs[family].append(record)
    
    # Create annotation lookup
    bgc_info = {}
    for _, row in bgc_annot_df.iterrows():
        record = row['Record']
        # Extract strain from record name (e.g., GCA_023625355.1_ASM2362535v1_scaffold...)
        strain = '_'.join(record.split('_')[:3]) if '_' in record else record.split('.')[0]
        bgc_info[record] = {
            'class': row.get('Class', ''),
            'category': row.get('Category', ''),
            'strain': strain,
        }
    
    # Build phenotype lookup if provided
    pheno_lookup = {}
    if phenotype_df is not None:
        for _, row in phenotype_df.iterrows():
            # Normalize sample ID
            sample = row['sample_id']
            norm_sample = normalize_sample_id(sample) if 'normalize_sample_id' in dir() else sample
            pheno_lookup[norm_sample] = row['phenotype']
    
    # Build GCF table
    gcf_rows = []
    for gcf, bgc_list in gcf_to_bgcs.items():
        # Get BGC info
        classes = []
        categories = []
        strains = set()
        
        for bgc in bgc_list:
            info = bgc_info.get(bgc, {})
            if info.get('class'):
                classes.append(info['class'])
            if info.get('category'):
                categories.append(info['category'])
            if info.get('strain'):
                strains.add(info['strain'])
        
        # Get most common type
        main_class = Counter(classes).most_common(1)[0][0] if classes else 'Unknown'
        main_category = Counter(categories).most_common(1)[0][0] if categories else 'Unknown'
        
        # Count phenotypes if available
        n_case = 0
        n_control = 0
        if pheno_lookup:
            for strain in strains:
                norm_strain = normalize_sample_id(strain) if 'normalize_sample_id' in dir() else strain
                if norm_strain in pheno_lookup:
                    if pheno_lookup[norm_strain] == 1:
                        n_case += 1
                    else:
                        n_control += 1
        
        # Get MIBiG annotations if available
        mibig_ids = ''
        mibig_compounds = ''
        if gcf_mibig and gcf in gcf_mibig:
            mibig_hits = gcf_mibig[gcf]
            mibig_ids = '; '.join([h[0] for h in mibig_hits[:3]])
            mibig_compounds = '; '.join([h[1][:40] for h in mibig_hits[:2]])
            if len(mibig_hits) > 3:
                mibig_ids += '...'
        
        gcf_rows.append({
            'GCF': gcf,
            'n_BGCs': len(bgc_list),
            'n_strains': len(strains),
            'main_class': main_class,
            'main_category': main_category,
            'all_classes': '; '.join(sorted(set(classes))),
            'strains': ', '.join(sorted(strains)[:5]) + ('...' if len(strains) > 5 else ''),
            'n_case': n_case,
            'n_control': n_control,
            'mibig_ids': mibig_ids,
            'mibig_compounds': mibig_compounds,
        })
    
    gcf_table = pd.DataFrame(gcf_rows)
    gcf_table = gcf_table.sort_values('n_strains', ascending=False)
    
    return gcf_table

def annotate_gcf_results_detailed(results_df, gcf_table, feature_col='feature'):
    """Add detailed GCF annotations to results."""
    if gcf_table.empty or results_df.empty:
        return results_df
    
    # Create lookup
    gcf_lookup = gcf_table.set_index('GCF').to_dict('index')
    
    results_df = results_df.copy()
    results_df['gcf_category'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('main_category', 'Unknown'))
    results_df['gcf_class'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('main_class', 'Unknown'))
    results_df['gcf_n_strains'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('n_strains', 0))
    results_df['gcf_n_bgcs'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('n_BGCs', 0))
    results_df['mibig_ids'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('mibig_ids', ''))
    results_df['mibig_compounds'] = results_df[feature_col].map(
        lambda x: gcf_lookup.get(x, {}).get('mibig_compounds', ''))
    
    return results_df

def display_significant_gcfs(results_df, gcf_table, contrast_name, fdr_threshold=0.1):
    """Display detailed info for significant GCFs."""

    # Accept either 'qvalue' (current schema) or legacy 'fdr'.
    fdr_col = 'qvalue' if 'qvalue' in results_df.columns else 'fdr'
    sig = results_df[results_df[fdr_col] < fdr_threshold].copy()
    if len(sig) == 0:
        print(f"\n{contrast_name}: No significant GCFs (FDR < {fdr_threshold})")
        return
    
    print(f"\n{'='*80}")
    print(f"SIGNIFICANT GCFs: {contrast_name}")
    print(f"{'='*80}")
    print(f"Total significant: {len(sig)} GCFs")
    
    # Annotate with full details
    sig = annotate_gcf_results_detailed(sig, gcf_table)
    
    # Separate by direction
    enriched_case = sig[sig['beta'] > 0].sort_values('pvalue')
    enriched_ctrl = sig[sig['beta'] < 0].sort_values('pvalue')
    
    case_label = contrast_name.split('_')[0].upper()
    ctrl_label = contrast_name.split('_')[-2].upper() if '_vs_' in contrast_name else 'CONTROL'
    
    if len(enriched_case) > 0:
        print(f"\n--- Enriched in {case_label} (n={len(enriched_case)}) ---")
        display_cols = ['feature', 'beta', 'pvalue', 'fdr', 'gcf_category', 'gcf_class', 
                        'gcf_n_strains', 'mibig_ids', 'mibig_compounds']
        display_cols = [c for c in display_cols if c in enriched_case.columns]
        print(enriched_case[display_cols].head(10).to_string())
        
        # Summary by category
        cat_counts = enriched_case['gcf_category'].value_counts()
        print(f"\n  Category summary:")
        for cat, count in cat_counts.items():
            print(f"    {cat}: {count}")
    
    if len(enriched_ctrl) > 0:
        print(f"\n--- Enriched in {ctrl_label} (n={len(enriched_ctrl)}) ---")
        display_cols = ['feature', 'beta', 'pvalue', 'fdr', 'gcf_category', 'gcf_class',
                        'gcf_n_strains', 'gcf_n_bgcs']
        display_cols = [c for c in display_cols if c in enriched_ctrl.columns]
        print(enriched_ctrl[display_cols].head(10).to_string())
        
        cat_counts = enriched_ctrl['gcf_category'].value_counts()
        print(f"\n  Category summary:")
        for cat, count in cat_counts.items():
            print(f"    {cat}: {count}")

# =============================================================================
# LOAD AND BUILD GCF ANNOTATIONS
# =============================================================================

# print("Loading BGC/GCF annotations...")
# bgc_annot_df, gcf_mapping = load_bgc_annotations(SPECIES_DIR)

# if bgc_annot_df is not None:
#     print("\nBuilding MIBiG compound annotations...")
#     gcf_mibig = build_mibig_annotation_table(SPECIES_DIR)
#     
#     print("\nBuilding GCF annotation table...")
#     gcf_annotation_table = build_gcf_annotation_table(bgc_annot_df, gcf_mapping, gcf_mibig=gcf_mibig)
#     
#     print(f"\nGCF Summary for {SPECIES_DISPLAY}:")
#     print(f"  Total GCF families: {len(gcf_annotation_table)}")
#     print(f"  Total BGCs: {gcf_annotation_table['n_BGCs'].sum()}")
#     
    # Show category distribution
#     print(f"\n  GCFs by category:")
#     cat_dist = gcf_annotation_table['main_category'].value_counts()
#     for cat, count in cat_dist.items():
#         print(f"    {cat}: {count}")
#     
    # Show top GCFs by strain count
#     print(f"\n  Top 10 most common GCFs:")
#     top_gcfs = gcf_annotation_table.head(10)
#     print(top_gcfs[['GCF', 'n_strains', 'n_BGCs', 'main_category', 'main_class']].to_string())
#     
    # Save full table
#     gcf_table_path = f"{results_dir}/gcf_annotation_table.tsv"
#     gcf_annotation_table.to_csv(gcf_table_path, sep='\t', index=False)
#     print(f"\n  Saved GCF annotation table: {gcf_table_path}")
# else:
#     gcf_annotation_table = pd.DataFrame()

# =============================================================================
# INTERPRET GCF RESULTS
# =============================================================================

# if gcf_pav_results and not gcf_annotation_table.empty:
#     print("\n" + "="*80)
#     print("GCF PAV ASSOCIATION RESULTS - DETAILED")
#     print("="*80)
#     for contrast, results in gcf_pav_results.items():
#         display_significant_gcfs(results, gcf_annotation_table, contrast)

# if gcf_cnv_results and not gcf_annotation_table.empty:
#     print("\n" + "="*80)
#     print("GCF CNV ASSOCIATION RESULTS - DETAILED")
#     print("="*80)
#     for contrast, results in gcf_cnv_results.items():
#         display_significant_gcfs(results, gcf_annotation_table, contrast)

# =============================================================================
# BGC TYPE INTERPRETATION GUIDE
# =============================================================================

# print("\n" + "="*80)
# print("BGC TYPE INTERPRETATION GUIDE")
# print("="*80)
# print("""
# BGC Categories and their biological significance:
# 
# NRPS (Non-Ribosomal Peptide Synthetase):
#   - Peptide-based secondary metabolites
#   - Examples: siderophores (iron acquisition), antibiotics, toxins
#   - A. niger: includes production of malformins, tensyuic acid
# 
# PKS (Polyketide Synthase):
#   - Polyketide secondary metabolites  
#   - Examples: pigments (melanin), mycotoxins (ochratoxin, fumonisin)
#   - A. niger: funalenone, azanigerone, carbonarin
# 
# Terpene:
#   - Terpenoid compounds
#   - Examples: volatile compounds, hormones
#   - A. niger: includes kotanin production
# 
# NRPS-PKS hybrids:
#   - Complex bioactive compounds with both pathways
#   - Often highly bioactive molecules
#   - Examples: pseurotin, cytochalasin
# 
# RiPP (Ribosomally synthesized Post-translationally modified Peptides):
#   - Small modified peptides
# 
# Interpretation for isolation source associations:
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# • Clinical-enriched GCFs: 
#   - May encode virulence factors or toxins
#   - Siderophores for iron scavenging in host
#   - Immune evasion compounds
# 
# • Environmental-enriched GCFs:
#   - Soil competition factors
#   - Antimicrobials against other microbes
#   - UV/oxidative stress protection (melanins)
# 
# • Industrial-enriched GCFs:
#   - May have been selected during domestication
#   - Or lost due to reduced selection pressure
# """)

# =============================================================================
# SECTION 7: FUNCTIONAL CHARACTERISATION OF SIGNIFICANT HITS
# =============================================================================
# Per-species functional enrichment + BGC co-localisation.
# Results saved to NB0_Results/{species}/ for NB3 cross-species comparison.
# =============================================================================

import json as _json
import re
import os
import glob
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# print("=" * 70)
# print(f"SECTION 7: FUNCTIONAL CHARACTERISATION — A. {SPECIES_SHORT}")
# print("=" * 70)

# ── 7.0 Load data ────────────────────────────────────────────────────

# og_consensus_path = f"{MAIN_RESULTS}/{SPECIES_SHORT}_og_consensus.tsv"
# og_consensus = pd.read_csv(og_consensus_path, sep='\t')

# Find the pan-GWAS results for the primary contrast
# gwas_dir = f"{RESULTS_DIR}/pangwas_results"
# gwas_files = sorted(glob.glob(f"{gwas_dir}/pav_assoc_*.tsv"))

# if not gwas_files:
#     print("No pan-GWAS result files found. Skipping Section 7.")
# else:
#     print(f"Found {len(gwas_files)} PAV GWAS result files:")
#     for f in gwas_files:
#         print(f"  {os.path.basename(f)}")

# ── 7.1 Functional enrichment ────────────────────────────────────────

def run_enrichment(sig_df, bg_df, col, parse_mode='multi', min_count=3):
    """Fisher's exact enrichment test for each term in col."""
    def _get_terms(series, mode):
        terms = {}
        for idx, val in series.items():
            if pd.isna(val) or str(val).strip() == '':
                continue
            if mode == 'bool':
                if val:
                    terms.setdefault('True', set()).add(idx)
            elif mode == 'single':
                for ch in str(val).strip():
                    if ch.strip():
                        terms.setdefault(ch.strip(), set()).add(idx)
            else:
                for t in str(val).replace(';', ',').split(','):
                    t = t.strip()
                    if t:
                        terms.setdefault(t, set()).add(idx)
        return terms
    sig_terms = _get_terms(sig_df[col], parse_mode)
    bg_terms = _get_terms(bg_df[col], parse_mode)
    n_sig, n_bg = len(sig_df), len(bg_df)
    rows = []
    for term in sorted(set(list(sig_terms.keys()) + list(bg_terms.keys()))):
        a = len(sig_terms.get(term, set()))
        if a < min_count:
            continue
        c = len(bg_terms.get(term, set())) - a
        b = n_sig - a
        d = n_bg - n_sig - c
        if c < 0: c = 0
        if d < 0: d = 0
        odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        rows.append({'term': term, 'n_sig': a, 'n_bg': a + c,
                     'OR': odds, 'pvalue': pval})
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    _, df['fdr'], _, _ = multipletests(df['pvalue'], method='fdr_bh')
    df['log2OR'] = np.log2(df['OR'].replace(0, np.nan).replace(np.inf, np.nan))
    return df.sort_values('fdr')

# LAYERS = [
#     ('COG',          'COG_category',  'single'),
#     ('Pfam',         'PFAMs',         'multi'),
#     ('CAZy',         'CAZy',          'multi'),
#     ('KEGG_Pathway', 'KEGG_Pathway',  'multi'),
#     ('KEGG_KO',      'KEGG_ko',       'multi'),
#     ('GO',           'GOs',           'multi'),
#     ('InterPro',     'interpro_IPR',  'multi'),
#     ('Protease',     'is_protease',   'bool'),
#     ('Transporter',  'is_transporter','bool'),
#     ('Secreted',     'is_secreted',   'bool'),
# ]

# FDR_THRESH = 0.1

# all_enrichment_by_contrast = {}

# for gwas_file in gwas_files:
#     contrast = os.path.basename(gwas_file).replace('pav_assoc_', '').replace('.tsv', '')
#     gwas = pd.read_csv(gwas_file, sep='\t')
#     sig_ogs = set(gwas.loc[gwas['fdr'] < FDR_THRESH, 'feature'])
#     all_ogs = set(gwas['feature'])
# 
#     if len(sig_ogs) < 3:
#         print(f"\n{contrast}: {len(sig_ogs)} sig OGs — too few, skipping")
#         continue
# 
#     sig_df = og_consensus[og_consensus['Orthogroup'].isin(sig_ogs)].copy()
#     bg_df = og_consensus[og_consensus['Orthogroup'].isin(all_ogs)].copy()
#     beta_map = gwas.set_index('feature')['beta']
#     sig_df['beta'] = sig_df['Orthogroup'].map(beta_map)
# 
#     min_count = 2 if len(sig_ogs) < 30 else 3
# 
#     print(f"\n{'='*60}")
#     print(f"Contrast: {contrast} ({len(sig_ogs)} sig / {len(all_ogs)} tested)")
#     print(f"{'='*60}")
# 
#     contrast_results = {}
#     for layer_name, col, mode in LAYERS:
#         if col not in sig_df.columns:
#             continue
#         enr = run_enrichment(sig_df, bg_df, col, parse_mode=mode, min_count=min_count)
#         contrast_results[layer_name] = enr
#         n_sig_terms = (enr['fdr'] < 0.1).sum() if len(enr) > 0 else 0
#         if n_sig_terms > 0:
#             print(f"  {layer_name}: {n_sig_terms} significant")
#             for _, row in enr[enr['fdr'] < 0.1].iterrows():
#                 d = 'ENRICHED' if row['OR'] > 1 else 'DEPLETED'
#                 print(f"    {row['term']:45s} OR={row['OR']:.2f} FDR={row['fdr']:.4f} "
#                       f"({row['n_sig']}/{row['n_bg']}) {d}")
# 
#     all_enrichment_by_contrast[contrast] = contrast_results
# 
    # Save per-contrast enrichment
#     all_enr_rows = []
#     for layer, enr in contrast_results.items():
#         if len(enr) > 0:
#             e = enr.copy()
#             e['layer'] = layer
#             all_enr_rows.append(e)
#     if all_enr_rows:
#         out = pd.concat(all_enr_rows)
#         out_path = f"{RESULTS_DIR}/{SPECIES_SHORT}_{contrast}_functional_enrichment.csv"
#         out.to_csv(out_path, index=False)
#         print(f"  Saved: {out_path}")


# ── 7.2 BGC co-localisation ──────────────────────────────────────────

# print(f"\n{'='*70}")
# print(f"SECTION 7.2: BGC CO-LOCALISATION — A. {SPECIES_SHORT}")
# print(f"{'='*70}")

def parse_antismash_bgc_genes(as_json_path):
    """Parse antiSMASH JSON → gene_id → {products, mibig, region_id}."""
    with open(as_json_path) as f:
        data = _json.load(f)
    bgc_genes = {}
    for rec in data['records']:
        areas = rec.get('areas', [])
        if not areas:
            continue
        scaffold = rec['id']
        features = rec.get('features', [])
        modules = rec.get('modules', {})
        kcb = modules.get('antismash.modules.clusterblast', {})
        kcb_results = kcb.get('knowncluster', {}).get('results', [])
        for area_idx, area in enumerate(areas):
            a_start, a_end = area['start'], area['end']
            products = ','.join(sorted(area.get('products', [])))
            mibig_top = 'none'
            if kcb_results:
                for kcb_r in kcb_results:
                    if kcb_r.get('region_number', -1) == area_idx + 1:
                        ranking = kcb_r.get('ranking', [])
                        if ranking and len(ranking[0]) > 0:
                            mibig_top = f"{ranking[0][0].get('accession','')}: {ranking[0][0].get('description','')}"
                        break
            for feat in features:
                if feat.get('type') != 'CDS':
                    continue
                loc = feat.get('location', '')
                nums = re.findall(r'\d+', loc)
                if not nums:
                    continue
                cds_start, cds_end = int(nums[0]), int(nums[-1])
                if cds_start < a_end and cds_end > a_start:
                    quals = feat.get('qualifiers', {})
                    locus = quals.get('locus_tag', ['?'])
                    locus = locus[0] if isinstance(locus, list) else locus
                    bgc_genes[locus] = {'products': products, 'mibig': mibig_top}
    return bgc_genes

def simplify_bgc_class(products_str):
    ps = set(products_str.split(','))
    if 'NRPS' in ps or 'NRPS-like' in ps:
        return 'NRPS-PKS hybrid' if any('PKS' in p for p in ps) else 'NRPS'
    if any('PKS' in p for p in ps): return 'PKS'
    if 'terpene' in ps or 'terpene-precursor' in ps: return 'Terpene'
    if 'isocyanide' in ps or 'isocyanide-nrp' in ps: return 'Isocyanide'
    if any('siderophore' in p.lower() for p in ps): return 'Siderophore'
    if any('RiPP' in p for p in ps): return 'RiPP'
    if 'betalactone' in ps: return 'Betalactone'
    if 'indole' in ps: return 'Indole'
    return 'Other'

# Load OrthoFinder OG table
# og_tsv = glob.glob(f'{SPECIES_DIR}/orthofinder_output/Results_*/Orthogroups/Orthogroups.tsv')
# as_dir = f'{SPECIES_DIR}/antismash_output'

# if og_tsv and os.path.exists(as_dir):
#     og = pd.read_csv(og_tsv[0], sep='\t')
#     genome_cols = [c for c in og.columns if c != 'Orthogroup']
# 
    # Aggregate OG → BGC across all genomes
#     og_bgc_votes = defaultdict(Counter)
#     og_mibig_votes = defaultdict(Counter)
#     og_in_bgc = set()
#     og_genome_in_bgc = defaultdict(int)
#     og_genome_total = defaultdict(int)
#     n_processed = 0
# 
#     for col in genome_cols:
#         strain = col.replace('.proteins', '')
#         as_json = os.path.join(as_dir, strain, f'{strain}.json')
#         if not os.path.exists(as_json):
#             continue
#         g2o = {}
#         ogs_present = set()
#         for _, row in og.iterrows():
#             genes = str(row[col])
#             if genes == 'nan':
#                 continue
#             ogs_present.add(row['Orthogroup'])
#             for g in genes.split(', '):
#                 g2o[g.replace('-T1', '')] = row['Orthogroup']
#         bgc_genes = parse_antismash_bgc_genes(as_json)
#         n_processed += 1
#         ogs_in_bgc_this = set()
#         for gene, info in bgc_genes.items():
#             og_id = g2o.get(gene)
#             if og_id:
#                 og_bgc_votes[og_id][info['products']] += 1
#                 og_mibig_votes[og_id][info['mibig']] += 1
#                 og_in_bgc.add(og_id)
#                 ogs_in_bgc_this.add(og_id)
        # Track per primary contrast
#         primary_gwas = gwas_files[0] if gwas_files else None
#         if primary_gwas:
#             primary = pd.read_csv(primary_gwas, sep='\t')
#             primary_sig = set(primary.loc[primary['fdr'] < 0.1, 'feature'])
#             for og_id in ogs_present & primary_sig:
#                 og_genome_total[og_id] += 1
#                 if og_id in ogs_in_bgc_this:
#                     og_genome_in_bgc[og_id] += 1
# 
#     print(f"  Processed {n_processed} genomes")
#     print(f"  OGs in BGC regions: {len(og_in_bgc)}")
# 
    # For primary contrast: build detail table
#     if gwas_files:
#         primary = pd.read_csv(gwas_files[0], sep='\t')
#         primary_sig = set(primary.loc[primary['fdr'] < 0.1, 'feature'])
#         primary_tested = set(primary['feature'])
#         primary_beta = primary.set_index('feature')['beta'].to_dict()
#         contrast_name = os.path.basename(gwas_files[0]).replace('pav_assoc_','').replace('.tsv','')
# 
#         sig_in_bgc = {og for og in primary_sig if og in og_in_bgc}
#         bg_in_bgc = {og for og in primary_tested if og in og_in_bgc}
# 
#         a = len(sig_in_bgc)
#         b = len(primary_sig) - a
#         c = len(bg_in_bgc) - a
#         d = len(primary_tested) - len(primary_sig) - c
#         overall_or, overall_p = fisher_exact([[a, b], [c, d]])
#         print(f"\n  {contrast_name}: {a}/{len(primary_sig)} sig OGs in BGCs "
#               f"({a/len(primary_sig)*100:.1f}%) vs {len(bg_in_bgc)}/{len(primary_tested)} bg "
#               f"({len(bg_in_bgc)/len(primary_tested)*100:.1f}%)")
#         print(f"  Overall BGC enrichment: OR={overall_or:.2f}, p={overall_p:.4f}")
# 
        # Build detail table
#         detail_rows = []
#         for og_id in sorted(sig_in_bgc):
#             beta = primary_beta.get(og_id, 0)
#             n_in = og_genome_in_bgc.get(og_id, 0)
#             n_tot = og_genome_total.get(og_id, 0)
#             consistency = n_in / n_tot * 100 if n_tot > 0 else 0
#             top_product = og_bgc_votes[og_id].most_common(1)[0] if og_id in og_bgc_votes else ('?', 0)
#             top_mibig = og_mibig_votes[og_id].most_common(1)[0] if og_id in og_mibig_votes else ('none', 0)
#             mibig_str = top_mibig[0]
#             if mibig_str == 'none' and len(og_mibig_votes.get(og_id, {})) > 1:
#                 second = og_mibig_votes[og_id].most_common(2)
#                 if len(second) > 1:
#                     mibig_str = second[1][0]
#             desc = og_consensus.loc[og_consensus['Orthogroup'] == og_id, 'Description'].iloc[0] \
#                    if og_id in og_consensus['Orthogroup'].values else ''
#             pfam = og_consensus.loc[og_consensus['Orthogroup'] == og_id, 'PFAMs'].iloc[0] \
#                    if og_id in og_consensus['Orthogroup'].values else ''
#             detail_rows.append({
#                 'OG': og_id, 'beta': beta,
#                 'direction': 'enriched' if beta > 0 else 'depleted',
#                 'BGC_class': simplify_bgc_class(top_product[0]),
#                 'BGC_raw': top_product[0],
#                 'MIBiG': mibig_str,
#                 'n_genomes_in_bgc': n_in, 'n_genomes_total': n_tot,
#                 'consistency_pct': consistency,
#                 'Description': desc, 'PFAMs': pfam,
#             })
#         detail_df = pd.DataFrame(detail_rows)
# 
#         if len(detail_df) > 0:
#             print(f"\n  Significant OGs in BGC regions:")
#             for _, row in detail_df.sort_values('consistency_pct', ascending=False).iterrows():
#                 arrow = '\u2191' if row['direction'] == 'enriched' else '\u2193'
#                 print(f"    {arrow} {row['OG']} ({row['BGC_class']}): "
#                       f"{str(row['MIBiG'])[:60]} — {row['consistency_pct']:.0f}% consistency")
# 
#             detail_df.to_csv(f"{RESULTS_DIR}/{SPECIES_SHORT}_sig_ogs_in_bgc_detail.csv", index=False)
#             print(f"\n  Saved: {RESULTS_DIR}/{SPECIES_SHORT}_sig_ogs_in_bgc_detail.csv")
# 
        # Per-class enrichment
#         og_to_class = {}
#         for og_id, votes in og_bgc_votes.items():
#             og_to_class[og_id] = simplify_bgc_class(votes.most_common(1)[0][0])
# 
#         class_rows = []
#         for bgc_class in sorted(set(og_to_class.values())):
#             n_sig_cls = sum(1 for og in primary_sig if og_to_class.get(og) == bgc_class)
#             n_bg_cls = sum(1 for og in primary_tested if og_to_class.get(og) == bgc_class)
#             if n_sig_cls < 1: continue
#             a2, b2 = n_sig_cls, len(primary_sig) - n_sig_cls
#             c2 = n_bg_cls - n_sig_cls
#             d2 = len(primary_tested) - len(primary_sig) - c2
#             if c2 < 0: c2 = 0
#             if d2 < 0: d2 = 0
#             or_val, pval = fisher_exact([[a2, b2], [c2, d2]])
#             class_rows.append({'class': bgc_class, 'n_sig': a2, 'n_bg': n_bg_cls,
#                                'OR': or_val, 'pvalue': pval})
#         class_df = pd.DataFrame(class_rows)
#         if len(class_df) > 0:
#             _, class_df['fdr'], _, _ = multipletests(class_df['pvalue'], method='fdr_bh')
#             class_df.to_csv(f"{RESULTS_DIR}/{SPECIES_SHORT}_bgc_class_enrichment.csv", index=False)
# 
        # Save MIBiG mapping for cross-species comparison
#         mibig_rows = []
#         for og_id in primary_sig:
#             votes = og_mibig_votes.get(og_id, {})
#             non_none = {k: v for k, v in votes.items() if k != 'none'}
#             if non_none:
#                 top = max(non_none, key=non_none.get)
#                 mibig_rows.append({
#                     'OG': og_id, 'MIBiG': top,
#                     'beta': primary_beta.get(og_id, 0),
#                     'species': SPECIES_SHORT
#                 })
#         if mibig_rows:
#             pd.DataFrame(mibig_rows).to_csv(
#                 f"{RESULTS_DIR}/{SPECIES_SHORT}_sig_og_mibig.csv", index=False)
# 
        # Save GCF mapping
#         bs_dir = f'{SPECIES_DIR}/bigscape_output'
#         gcf_map_files = glob.glob(f'{bs_dir}/output_files/*_c0.*/*/*.tsv')
#         if gcf_map_files:
#             gcf_dfs = []
#             for f in gcf_map_files:
#                 d = pd.read_csv(f, sep='\t')
#                 if 'Family' in d.columns:
#                     gcf_dfs.append(d)
#             if gcf_dfs:
#                 gcf_combined = pd.concat(gcf_dfs, ignore_index=True)
#                 gcf_combined = gcf_combined[gcf_combined['GBK'].str.startswith(('GCA','GCF'), na=False)]
#                 gcf_map = dict(zip(gcf_combined['GBK'], gcf_combined['Family']))
# 
#                 og_gcf_votes_local = defaultdict(Counter)
#                 for col_name in genome_cols:
#                     strain = col_name.replace('.proteins', '')
#                     as_json = os.path.join(as_dir, strain, f'{strain}.json')
#                     if not os.path.exists(as_json): continue
#                     g2o_local = {}
#                     for _, row in og.iterrows():
#                         genes = str(row[col_name])
#                         if genes == 'nan': continue
#                         for g in genes.split(', '):
#                             g2o_local[g.replace('-T1', '')] = row['Orthogroup']
# 
#                     with open(as_json) as fh:
#                         data = _json.load(fh)
#                     for rec in data['records']:
#                         areas = rec.get('areas', [])
#                         if not areas: continue
#                         scaffold = rec['id']
#                         features = rec.get('features', [])
#                         for area_idx, area in enumerate(areas):
#                             a_s, a_e = area['start'], area['end']
#                             region_gbk = f"{strain}_{scaffold}.region{area_idx+1:03d}"
#                             gcf_id = gcf_map.get(region_gbk)
#                             if gcf_id is None:
#                                 for k in gcf_map:
#                                     if strain in k and scaffold in k and f"region{area_idx+1:03d}" in k:
#                                         gcf_id = gcf_map[k]
#                                         break
#                             if gcf_id:
#                                 for feat in features:
#                                     if feat.get('type') != 'CDS': continue
#                                     loc = feat.get('location', '')
#                                     nums = re.findall(r'\d+', loc)
#                                     if not nums: continue
#                                     cs, ce = int(nums[0]), int(nums[-1])
#                                     if cs < a_e and ce > a_s:
#                                         quals = feat.get('qualifiers', {})
#                                         locus = quals.get('locus_tag', ['?'])
#                                         locus = locus[0] if isinstance(locus, list) else locus
#                                         og_id = g2o_local.get(locus)
#                                         if og_id:
#                                             og_gcf_votes_local[og_id][gcf_id] += 1
# 
#                 gcf_rows = []
#                 for og_id in primary_sig:
#                     votes = og_gcf_votes_local.get(og_id, {})
#                     if votes:
#                         top_gcf = votes.most_common(1)[0][0]
#                         gcf_rows.append({
#                             'OG': og_id, 'GCF': top_gcf,
#                             'beta': primary_beta.get(og_id, 0),
#                             'species': SPECIES_SHORT
#                         })
#                 if gcf_rows:
#                     pd.DataFrame(gcf_rows).to_csv(
#                         f"{RESULTS_DIR}/{SPECIES_SHORT}_sig_og_gcf.csv", index=False)
#                     print(f"  Saved: {RESULTS_DIR}/{SPECIES_SHORT}_sig_og_gcf.csv")
# 
# else:
#     print("  OrthoFinder or antiSMASH output not found — skipping BGC analysis")

# print(f"\n{'='*70}")
# print(f"Section 7 complete for A. {SPECIES_SHORT}")
# print(f"{'='*70}")


# ── 7.3 Functional annotation overlap with other species (record) ────
# For the current species' significant OGs, list which Pfam/COG/CAZy/GO
# terms they carry. This is saved so that NB3 can compute the
# set-intersection overlap between species' hit annotation sets.

# if gwas_files:
#     print(f"\n{'='*70}")
#     print(f"SECTION 7.3: FUNCTIONAL ANNOTATION SETS FOR CROSS-SPECIES OVERLAP")
#     print(f"{'='*70}")
#     
#     primary = pd.read_csv(gwas_files[0], sep='\t')
#     contrast_name = os.path.basename(gwas_files[0]).replace('pav_assoc_','').replace('.tsv','')
#     sig_ogs_set = set(primary.loc[primary['fdr'] < 0.1, 'feature'])
#     sig_annot = og_consensus[og_consensus['Orthogroup'].isin(sig_ogs_set)]
#     
#     annotation_sets = {}
#     for layer_name, col, mode in [('COG', 'COG_category', 'single'),
#                                     ('Pfam', 'PFAMs', 'multi'),
#                                     ('CAZy', 'CAZy', 'multi'),
#                                     ('GO', 'GOs', 'multi'),
#                                     ('KEGG_Pathway', 'KEGG_Pathway', 'multi'),
#                                     ('KEGG_KO', 'KEGG_ko', 'multi'),
#                                     ('InterPro', 'interpro_IPR', 'multi')]:
#         if col not in sig_annot.columns:
#             continue
#         terms = set()
#         for val in sig_annot[col].dropna():
#             if mode == 'single':
#                 for ch in str(val).strip():
#                     if ch.strip():
#                         terms.add(ch.strip())
#             else:
#                 for t in str(val).replace(';', ',').split(','):
#                     t = t.strip()
#                     if t:
#                         terms.add(t)
#         annotation_sets[layer_name] = terms
#         print(f"  {layer_name}: {len(terms)} unique terms in {len(sig_ogs_set)} sig OGs")
#     
    # Save as CSV for NB3 cross-species overlap
#     rows = []
#     for layer, terms in annotation_sets.items():
#         for t in sorted(terms):
#             rows.append({'layer': layer, 'term': t, 'species': SPECIES_SHORT, 'contrast': contrast_name})
#     if rows:
#         out_path = f"{RESULTS_DIR}/{SPECIES_SHORT}_{contrast_name}_annotation_sets.csv"
#         pd.DataFrame(rows).to_csv(out_path, index=False)
#         print(f"\n  Saved: {out_path}")


# =============================================================================
# NB2 REFACTORED FUNCTIONS (notebook cell extractions)
# =============================================================================

def load_all_species_gwas_data(species_list, nb0_results, nb1_results,
                               annotation_layers):
    """
    Load per-species data: PAV, CNV, kinship, PCs, BGC/GCF matrices,
    OG consensus annotations, and functional count matrices.

    Parameters
    ----------
    species_list : list of str
        Species identifiers.
    nb0_results : Path
        Path to NB0_Results directory.
    nb1_results : Path
        Path to NB1_Results directory.
    annotation_layers : list of str
        Annotation layer names (e.g. ['COG', 'Pfam', ...]).

    Returns
    -------
    dict
        {species: {key: pd.DataFrame, ...}}
    """
    from pathlib import Path
    species_data = {}

    for sp in species_list:
        sp_dir = Path(nb1_results) / sp
        sp_data = {}

        # Per-species phenotype data (may be in NB0 root or per-species subdir)
        for pheno_path in [Path(nb0_results) / sp / 'phenotype_data.tsv',
                            Path(nb0_results) / f'{sp}_phenotype_data.tsv',
                            Path(nb0_results) / 'phenotype_classified_for_gwas.csv']:
            if pheno_path.exists():
                sep = ',' if pheno_path.suffix == '.csv' else '\t'
                df = pd.read_csv(pheno_path, sep=sep)
                if 'Species' in df.columns:
                    df = df[df['Species'].astype(str).str.contains(sp, case=False, na=False)]
                # Filter ANI-excluded samples by row
                excluded = ANI_EXCLUDED.get(sp, set())
                if excluded:
                    for col in ('sample_id', 'Assembly Accession', 'accession'):
                        if col in df.columns:
                            df = df[~df[col].astype(str).apply(
                                lambda s, exc=excluded: any(a in s for a in exc))]
                            break
                sp_data['phenotype'] = df
                break

        # Standard matrices — named like {sp}_{type}.tsv in NB1_Results/{sp}/
        # Each is loaded then filtered to drop ANI-excluded genomes at runtime.
        matrix_keys = [
            ('pav',     f'{sp}_pav.tsv',     'cols'),
            ('cnv',     f'{sp}_cnv.tsv',     'cols'),
            ('kinship', f'{sp}_kinship.tsv', 'both'),
            ('pcs',     f'{sp}_snp_pcs.tsv', 'auto'),
            ('bgc_pav', f'{sp}_bgc_pav.tsv', 'auto'),
            ('gcf_pav', f'{sp}_gcf_pav.tsv', 'auto'),
            ('gcf_cnv', f'{sp}_gcf_cnv.tsv', 'auto'),
        ]
        for key, fname, axis in matrix_keys:
            fpath = sp_dir / fname
            if fpath.exists():
                df = pd.read_csv(fpath, sep='\t', index_col=0)
                df = filter_ani_excluded(df, sp, axis=axis)
                sp_data[key] = df

        # OG consensus annotations
        og_path = sp_dir / f'{sp}_og_consensus.tsv'
        if og_path.exists():
            sp_data['og_consensus'] = pd.read_csv(og_path, sep='\t')

        # Functional count matrices (per annotation layer)
        for layer in annotation_layers:
            func_path = sp_dir / f'{sp}_functional_{layer.lower()}_counts.tsv'
            if func_path.exists():
                df = pd.read_csv(func_path, sep='\t', index_col=0)
                df = filter_ani_excluded(df, sp, axis='auto')
                sp_data[f'functional_{layer.lower()}'] = df

        species_data[sp] = sp_data
        loaded_keys = list(sp_data.keys())
        print(f'{sp}: loaded {len(loaded_keys)} data objects -- {loaded_keys}')

    return species_data


def create_one_vs_rest_phenotypes(pheno_df, pheno_col=None, id_col=None,
                                   exclude=('Unknown', 'unknown', '', 'Lab'),
                                   min_cases=3, min_controls=3):
    """Create one-vs-rest binary phenotype DataFrames for each category.

    Returns
    -------
    dict
        ``{contrast_name: pd.DataFrame}`` where each DataFrame has columns
        ``sample_id`` and ``phenotype`` (0/1).
    """
    # Identify phenotype column
    if pheno_col is None:
        for candidate in ['isolation_class', 'Phenotype', 'phenotype',
                          'Niche', 'niche', 'PhenotypeClassification']:
            if candidate in pheno_df.columns:
                pheno_col = candidate
                break
    if pheno_col is None or pheno_col not in pheno_df.columns:
        return {}

    # Identify sample ID column
    if id_col is None:
        for c in ['sample_id', 'Assembly Accession', 'accession']:
            if c in pheno_df.columns:
                id_col = c
                break
        if id_col is None:
            id_col = pheno_df.columns[0]

    # Filter out excluded categories
    known = pheno_df[~pheno_df[pheno_col].isin(exclude)].copy()
    known = known[known[pheno_col].notna()]

    contrasts = {}
    for target_cat in known[pheno_col].unique():
        contrast_name = str(target_cat).lower().replace('-', '_').replace(' ', '_') + '_vs_rest'
        is_case = (known[pheno_col] == target_cat).astype(int)
        n_case = int(is_case.sum())
        n_control = int((~is_case.astype(bool)).sum())
        if n_case < min_cases or n_control < min_controls:
            continue
        contrasts[contrast_name] = pd.DataFrame({
            'sample_id': known[id_col].values,
            'phenotype': is_case.values,
        })
    return contrasts


def build_all_phenotype_contrasts(species_list, species_data, nb2_results, force=False):
    """
    Build one-vs-rest phenotype contrasts per species and save them.

    Cached to ``{nb2_results}/{sp}/phenotypes/pheno_*.tsv``.

    Returns
    -------
    dict
        {species: {contrast_name: pd.Series}}
    """
    from pathlib import Path
    species_contrasts = {}

    for sp in species_list:
        sp_out = Path(nb2_results) / sp / 'phenotypes'
        sp_out.mkdir(parents=True, exist_ok=True)

        # --- Load from cache ---
        cached_files = list(sp_out.glob('pheno_*.tsv'))
        if not force and cached_files:
            contrasts = {}
            for fp in cached_files:
                cname = fp.stem.replace('pheno_', '')
                df = pd.read_csv(fp, sep='\t')
                # Normalise to {'sample_id', 'phenotype'} schema
                if 'sample_id' not in df.columns:
                    df = df.rename(columns={df.columns[0]: 'sample_id'})
                if 'phenotype' not in df.columns:
                    df = df.rename(columns={df.columns[-1]: 'phenotype'})
                contrasts[cname] = df[['sample_id', 'phenotype']]
            species_contrasts[sp] = contrasts
            print(f'{sp}: loaded {len(contrasts)} cached contrasts')
            continue

        if 'phenotype' not in species_data[sp]:
            print(f'{sp}: no phenotype data, skipping.')
            continue

        pheno_df = species_data[sp]['phenotype']
        contrasts = create_one_vs_rest_phenotypes(pheno_df)
        species_contrasts[sp] = contrasts

        for cname, cdf in contrasts.items():
            out_path = sp_out / f'pheno_{cname}.tsv'
            cdf.to_csv(out_path, sep='\t', index=False)

        # Summarise
        contrast_summary = {
            c: f'{int(v["phenotype"].sum())} cases / {int((~v["phenotype"].astype(bool)).sum())} controls'
            for c, v in contrasts.items()
        }
        print(f'\n{sp} -- {len(contrasts)} contrasts:')
        for c, s in contrast_summary.items():
            print(f'  {c}: {s}')

    return species_contrasts


def run_all_associations(species_list, species_contrasts, species_data,
                         nb2_results, feature_layers, annotation_layers,
                         fdr_threshold=0.05, n_pcs=0, force=False):
    """
    Run LMM association testing for all species x contrast x feature layer
    combinations. Applies FDR correction and saves per-test result files.

    Cached to ``{nb2_results}/{sp}/pangwas_results/{layer}_assoc_{contrast}.tsv``.
    Pass ``force=True`` to rerun.

    Returns
    -------
    dict
        {(species, contrast, layer): pd.DataFrame}
    """
    from pathlib import Path
    all_results = {}

    for sp in species_list:
        if sp not in species_contrasts or not species_contrasts[sp]:
            print(f'\n{sp}: no contrasts defined, skipping.')
            continue

        sp_results_dir = Path(nb2_results) / sp / 'pangwas_results'
        sp_results_dir.mkdir(parents=True, exist_ok=True)

        kinship = species_data[sp].get('kinship')
        pcs = species_data[sp].get('pcs')

        if kinship is None:
            print(f'\n{sp}: no kinship matrix, skipping association tests.')
            continue

        print(f'\n{"=" * 60}')
        print(f'{sp}: running association tests')
        print(f'{"=" * 60}')

        for contrast_name, phenotype in species_contrasts[sp].items():
            print(f'\n  Contrast: {contrast_name}')

            # Determine which feature layers to test.
            # PAV / CNV are written as (features x samples) i.e. OG-rows by
            # strain-columns. BGC PAV / GCF PAV / GCF CNV are written as
            # (samples x features) i.e. strain-rows by BGC/GCF-columns, so
            # they must be transposed before LMM testing (which iterates
            # rows as features).
            transpose_layers = {'bgc_pav', 'gcf_pav', 'gcf_cnv'}
            layers_to_test = []
            for layer in feature_layers:
                if layer in species_data[sp]:
                    fm = species_data[sp][layer]
                    if layer in transpose_layers:
                        fm = fm.T
                    layers_to_test.append((layer, fm))

            # Also add functional count layers
            for layer in annotation_layers:
                key = f'functional_{layer.lower()}'
                if key in species_data[sp]:
                    layers_to_test.append((key, species_data[sp][key]))

            for layer_name, feature_matrix in layers_to_test:
                out_path = sp_results_dir / f'{layer_name}_assoc_{contrast_name}.tsv'

                # --- Load from cache ---
                if not force and out_path.exists():
                    results = pd.read_csv(out_path, sep='\t', index_col=0)
                    all_results[(sp, contrast_name, layer_name)] = results
                    n_sig = int(results.get('significant', pd.Series(dtype=bool)).sum()) \
                            if 'significant' in results.columns else 0
                    print(f'    Layer: {layer_name} ({feature_matrix.shape[1]} features) '
                          f'... cached ({n_sig} significant)')
                    continue

                print(f'    Layer: {layer_name} ({feature_matrix.shape[1]} features) ... ',
                      end='')

                try:
                    results = run_association_test(
                        feature_matrix=feature_matrix,
                        phenotype_df=phenotype,
                        pcs_df=pcs,
                        grm_df=kinship,
                        n_pcs=n_pcs,
                    )

                    if len(results) > 0 and 'pvalue' in results.columns:
                        reject, qvals, _, _ = multipletests(
                            results['pvalue'].dropna(), method='fdr_bh')
                        results.loc[results['pvalue'].notna(), 'qvalue'] = qvals
                        results['significant'] = results['qvalue'] < fdr_threshold
                        n_sig = results['significant'].sum()
                    else:
                        n_sig = 0

                    results.to_csv(out_path, sep='\t')
                    all_results[(sp, contrast_name, layer_name)] = results
                    print(f'{n_sig} significant (FDR < {fdr_threshold})')

                except Exception as e:
                    print(f'FAILED: {e}')
                    continue

    return all_results


def run_power_analysis_grid(species_list, species_contrasts, nb2_results,
                            fdr_threshold=0.05, force=False):
    """
    Run power analysis across species/contrasts for a range of odds ratios.

    Parameters
    ----------
    species_list : list of str
    species_contrasts : dict
    nb2_results : Path
    fdr_threshold : float

    Returns
    -------
    pd.DataFrame
        Power analysis results.
    """
    from pathlib import Path
    out_path = Path(nb2_results) / 'power_analysis_results.tsv'

    # --- Load from cache ---
    if not force and out_path.exists():
        power_df = pd.read_csv(out_path, sep='\t')
        print(f'Loaded cached power analysis ({len(power_df)} rows) from {out_path}')
        return power_df

    power_rows = []

    for sp in species_list:
        if sp not in species_contrasts:
            continue
        for contrast_name, phenotype in species_contrasts[sp].items():
            ph = phenotype['phenotype'] if isinstance(phenotype, pd.DataFrame) else phenotype
            n_cases = int(ph.sum())
            n_controls = int((~ph.astype(bool)).sum())

            # Power for a range of odds ratios at gene_freq=0.2 (typical for
            # variable pangenome genes; 0.5 is the worst-case frequency for
            # Fisher's test on small samples and would overestimate the
            # minimum detectable OR).
            for target_or in [2.0, 3.0, 5.0, 10.0]:
                power = power_fisher(
                    n_case=n_cases,
                    n_control=n_controls,
                    gene_freq=0.2,
                    odds_ratio=target_or,
                    alpha=fdr_threshold,
                )
                power_rows.append({
                    'species': sp,
                    'contrast': contrast_name,
                    'n_cases': n_cases,
                    'n_controls': n_controls,
                    'ratio': f'{n_cases}:{n_controls}',
                    'target_OR': target_or,
                    'power': round(power, 4),
                })

            # Minimum detectable odds ratio at 80% power — take best-case across
            # gene_freq range (0.1 / 0.2 / 0.3) since the best frequency
            # reflects the true detectability ceiling.
            min_or_vals = []
            for gf in [0.1, 0.2, 0.3]:
                mo = find_min_detectable_or(
                    n_case=n_cases, n_control=n_controls, gene_freq=gf,
                    alpha=fdr_threshold, target_power=0.80,
                )
                min_or_vals.append(mo)
            min_or = min(min_or_vals)
            power_rows.append({
                'species': sp,
                'contrast': contrast_name,
                'n_cases': n_cases,
                'n_controls': n_controls,
                'ratio': f'{n_cases}:{n_controls}',
                'target_OR': 'min_OR@80%power',
                'power': round(min_or, 2) if min_or != float('inf') else min_or,
            })

    power_df = pd.DataFrame(power_rows)
    out_path = Path(nb2_results) / 'power_analysis_results.tsv'
    power_df.to_csv(out_path, sep='\t', index=False)
    print(f'Saved: {out_path}')

    return power_df


def run_functional_enrichment_all(species_list, species_data, all_results,
                                  nb2_results, annotation_layers,
                                  fdr_threshold=0.05, min_n_bg=3, force=False):
    """
    Run functional enrichment (Fisher's exact test) of significant PAV/CNV
    orthogroups against each annotation layer.

    Cached to ``{nb2_results}/{sp}/{sp}_functional_enrichment.csv``.

    Returns
    -------
    dict
        {species: pd.DataFrame} of enrichment results.
    """
    from pathlib import Path
    from collections import defaultdict
    enrichment_results = {}

    for sp in species_list:
        sp_dir = Path(nb2_results) / sp
        sp_dir.mkdir(parents=True, exist_ok=True)
        out_path = sp_dir / f'{sp}_functional_enrichment.csv'

        # --- Load from cache ---
        if not force and out_path.exists():
            enrichment_results[sp] = pd.read_csv(out_path)
            print(f'{sp}: loaded cached enrichment ({len(enrichment_results[sp])} rows)')
            continue

        if 'og_consensus' not in species_data[sp]:
            print(f'{sp}: no OG consensus annotations, skipping enrichment.')
            continue

        og_annot = species_data[sp]['og_consensus']
        # Index og_annot by Orthogroup ID so we iterate by OG, not row number.
        og_annot_by_og = og_annot.set_index('Orthogroup', drop=False)

        # Build the per-contrast significant-OG sets and a single shared
        # background of all tested OGs (the union of PAV/CNV across contrasts).
        # The 'feature' column carries the OG ID; integer index is post-sort
        # row position and must NOT be used as an OG identifier.
        per_contrast_sig = {}
        all_ogs = set()
        for (s, contrast, layer), res in all_results.items():
            if s != sp or layer not in ('pav', 'cnv'):
                continue
            if 'feature' not in res.columns:
                continue
            all_ogs.update(res['feature'].astype(str).tolist())
            if 'significant' in res.columns:
                sig_set = per_contrast_sig.setdefault(contrast, set())
                sig_set.update(res.loc[res['significant'], 'feature']
                               .astype(str).tolist())

        per_contrast_sig = {c: s for c, s in per_contrast_sig.items() if s}
        if not per_contrast_sig:
            print(f'{sp}: no significant orthogroups, skipping enrichment.')
            continue

        print(f'\n{sp}: enrichment for {len(per_contrast_sig)} contrast(s):')
        for c, s in per_contrast_sig.items():
            print(f'  {c}: {len(s)} significant OGs')

        # Run enrichment per (contrast, annotation layer), with FDR applied
        # within each (contrast, layer) block (the Theme convention).
        enrichment_rows = []
        for contrast_name, sig_ogs in per_contrast_sig.items():
            for layer in annotation_layers:
                if layer not in og_annot_by_og.columns:
                    continue

                annot_series = og_annot_by_og[layer].dropna()

                # Expand multi-term annotations (semicolon-separated)
                term_to_ogs = defaultdict(set)
                for og, terms_str in annot_series.items():
                    og = str(og)
                    if og not in all_ogs:
                        continue
                    for term in str(terms_str).split(';'):
                        term = term.strip()
                        if term and term != 'nan':
                            term_to_ogs[term].add(og)

                # Drop uninformative terms (COG-S/R, Pfam DUFs, GO roots, etc.)
                # before any Fisher tests / FDR correction. Standard practice in
                # functional-enrichment work; see funpan_utils.is_informative_term.
                from funpan_utils import is_informative_term
                term_to_ogs = {t: ogs for t, ogs in term_to_ogs.items()
                               if is_informative_term(layer, t)}

                # Fisher's exact test per term, one-sided (enrichment).
                # Drop terms appearing in < min_n_bg background OGs to avoid
                # singleton/doubleton-driven FDR inflation (standard practice;
                # term_to_ogs at n_term==1 produces extreme but unstable OR).
                n_total = len(all_ogs)
                n_sig = len(sig_ogs)
                layer_rows = []
                for term, term_ogs in term_to_ogs.items():
                    n_term = len(term_ogs)
                    if n_term < min_n_bg:
                        continue
                    n_sig_and_term = len(sig_ogs & term_ogs)

                    table = [
                        [n_sig_and_term, n_sig - n_sig_and_term],
                        [n_term - n_sig_and_term,
                         n_total - n_sig - n_term + n_sig_and_term],
                    ]
                    odds_ratio, pvalue = stats.fisher_exact(
                        table, alternative='greater')

                    layer_rows.append({
                        'contrast': contrast_name,
                        'annotation_layer': layer,
                        'term': term,
                        'n_sig_with_term': n_sig_and_term,
                        'n_sig_total': n_sig,
                        'n_bg_with_term': n_term,
                        'n_bg_total': n_total,
                        'odds_ratio': odds_ratio,
                        'pvalue': pvalue,
                    })

                if not layer_rows:
                    continue

                # FDR correction within this (contrast, layer) block
                ldf = pd.DataFrame(layer_rows)
                ldf['qvalue'] = multipletests(ldf['pvalue'], method='fdr_bh')[1]
                ldf['significant'] = ldf['qvalue'] < fdr_threshold
                enrichment_rows.append(ldf)

        if not enrichment_rows:
            print(f'  No terms tested.')
            continue

        enr_df = pd.concat(enrichment_rows, ignore_index=True)
        enr_df = enr_df.sort_values(['contrast', 'annotation_layer', 'pvalue'])

        # Save
        out_path = sp_dir / f'{sp}_functional_enrichment.csv'
        enr_df.to_csv(out_path, index=False)
        enrichment_results[sp] = enr_df

        n_sig_terms = int(enr_df['significant'].sum())
        print(f'  {len(enr_df)} terms tested across all (contrast, layer) blocks; '
              f'{n_sig_terms} significant (FDR < {fdr_threshold} per block)')
        if n_sig_terms > 0:
            print(f'  Significant terms:')
            print(enr_df[enr_df['significant']].head(20).to_string(index=False))

    return enrichment_results


def _build_bgc_og_mapping(species, genus='Aspergillus'):
    """Build a BGC-region-to-orthogroup mapping from antiSMASH JSON + OrthoFinder.

    For each antiSMASH region, identifies CDS features whose coordinates fall
    within the region, then maps those gene IDs to orthogroups via the
    OrthoFinder Orthogroups.tsv.

    Returns
    -------
    pd.DataFrame with columns: orthogroup, gene_id, accession, bgc_id,
    scaffold, product, mibig, region_start, region_end
    """
    import json as _json, re, glob
    species_dir = f'/datadrive/Species/{genus}/{species}'
    as_base = os.path.join(species_dir, 'antismash_output')
    og_tsvs = glob.glob(os.path.join(species_dir, 'orthofinder_output',
                                      '*', 'Orthogroups', 'Orthogroups.tsv'))
    if not og_tsvs or not os.path.isdir(as_base):
        return pd.DataFrame()

    # Load OrthoFinder — build gene_id -> orthogroup lookup
    og_df = pd.read_csv(sorted(og_tsvs)[-1], sep='\t')
    gene_to_og = {}
    # Build gene_id -> orthogroup with aliases. antiSMASH JSON qualifiers and
    # OrthoFinder Orthogroups.tsv use overlapping but distinct ID conventions:
    #   antiSMASH JSON:    'FUN_002899' (no transcript)  or  'ncbi_FUN_000318-T1'
    #   OrthoFinder TSV:   'FUN_000263-T1'
    # We register every gene under (a) its raw form, (b) the form without any
    # trailing -T<digits> transcript suffix, (c) both forms with a 'ncbi_'
    # prefix prepended. Lookup by the antiSMASH gid then succeeds regardless
    # of which convention the JSON uses.
    _t_re = re.compile(r'-T\d+$')
    for _, row in og_df.iterrows():
        og_id = row['Orthogroup']
        for col in og_df.columns[1:]:
            val = row[col]
            if pd.isna(val):
                continue
            for gene in str(val).split(', '):
                gene = gene.strip()
                if not gene:
                    continue
                gene_to_og[gene] = og_id
                base = _t_re.sub('', gene)
                if base != gene:
                    gene_to_og[base] = og_id
                gene_to_og['ncbi_' + gene] = og_id
                if base != gene:
                    gene_to_og['ncbi_' + base] = og_id

    rows = []
    for as_json in glob.glob(os.path.join(as_base, '*', '*.json')):
        acc_m = re.search(r'(GC[AF]_\d+\.\d+)', as_json)
        accession = acc_m.group(1) if acc_m else os.path.basename(as_json)
        try:
            with open(as_json) as f:
                data = _json.load(f)
        except Exception:
            continue

        for rec in data.get('records', []):
            areas = rec.get('areas', [])
            if not areas:
                continue
            scaffold = rec.get('id', '?')
            features = rec.get('features', [])

            for area_idx, area in enumerate(areas):
                a_start, a_end = area.get('start', 0), area.get('end', 0)
                products = ','.join(sorted(area.get('products', [])))
                bgc_id = f'{accession}__{scaffold}__region{area_idx + 1}'

                for feat in features:
                    if feat.get('type') != 'CDS':
                        continue
                    loc_nums = re.findall(r'\d+', feat.get('location', ''))
                    if not loc_nums:
                        continue
                    cds_start, cds_end = int(loc_nums[0]), int(loc_nums[-1])
                    if cds_end < a_start or cds_start > a_end:
                        continue

                    qualifiers = feat.get('qualifiers', {})
                    gene_ids = (qualifiers.get('locus_tag', []) +
                                qualifiers.get('gene', []) +
                                qualifiers.get('protein_id', []))
                    for gid in gene_ids:
                        og_id = gene_to_og.get(gid)
                        if og_id:
                            rows.append({
                                'orthogroup': og_id,
                                'gene_id': gid,
                                'accession': accession,
                                'bgc_id': bgc_id,
                                'scaffold': scaffold,
                                'product': products,
                                'region_start': a_start,
                                'region_end': a_end,
                            })

    return pd.DataFrame(rows).drop_duplicates()


def run_bgc_colocalisation(species_list, all_results, nb1_results,
                           nb2_results, force=False):
    """
    Map significant orthogroups onto antiSMASH-predicted BGC regions.

    Cached to ``{nb2_results}/{sp}/{sp}_sig_ogs_in_bgc_detail.csv``.

    Returns
    -------
    dict
        {species: pd.DataFrame} of BGC co-localisation details.
    """
    from pathlib import Path
    coloc_results = {}

    for sp in species_list:
        sp_dir = Path(nb2_results) / sp
        out_path = sp_dir / f'{sp}_sig_ogs_in_bgc_detail.csv'

        # --- Load from cache ---
        if not force and out_path.exists():
            coloc_results[sp] = pd.read_csv(out_path)
            print(f'{sp}: loaded cached BGC co-localisation ({len(coloc_results[sp])} rows)')
            continue

        # Collect significant OG IDs from PAV/CNV associations.
        # Use the 'feature' column (OG ID), not the integer row index.
        sig_ogs = set()
        for (s, contrast, layer), res in all_results.items():
            if s != sp or layer not in ('pav', 'cnv'):
                continue
            if 'significant' in res.columns and 'feature' in res.columns:
                sig_ogs.update(res.loc[res['significant'], 'feature']
                               .astype(str).tolist())

        if not sig_ogs:
            continue

        # Load BGC-to-OG mapping — build on the fly if not cached
        bgc_og_map_path = Path(nb1_results) / sp / 'bgc_og_mapping.tsv'
        if not bgc_og_map_path.exists():
            print(f'{sp}: building BGC-OG mapping from antiSMASH output + OrthoFinder...')
            try:
                bgc_og_map = _build_bgc_og_mapping(sp)
                if bgc_og_map.empty:
                    print(f'{sp}: no BGC data found, skipping co-localisation.')
                    continue
                bgc_og_map_path.parent.mkdir(parents=True, exist_ok=True)
                bgc_og_map.to_csv(bgc_og_map_path, sep='\t', index=False)
                print(f'  Saved BGC-OG mapping to {bgc_og_map_path}')
            except Exception as e:
                print(f'{sp}: failed to build BGC-OG mapping ({e}), skipping.')
                continue
        else:
            bgc_og_map = pd.read_csv(bgc_og_map_path, sep='\t')

        # Filter to significant OGs
        sig_in_bgc = bgc_og_map[bgc_og_map['orthogroup'].isin(sig_ogs)].copy()

        if sig_in_bgc.empty:
            print(f'{sp}: no significant OGs map to BGC regions.')
            continue

        # Save detailed mapping
        sp_dir.mkdir(parents=True, exist_ok=True)
        out_path = sp_dir / f'{sp}_sig_ogs_in_bgc_detail.csv'
        sig_in_bgc.to_csv(out_path, index=False)

        n_ogs_in_bgc = sig_in_bgc['orthogroup'].nunique()
        n_bgcs = (sig_in_bgc['bgc_id'].nunique()
                  if 'bgc_id' in sig_in_bgc.columns else 'N/A')
        print(f'{sp}: {n_ogs_in_bgc} significant OGs fall within '
              f'{n_bgcs} BGC regions')
        print(sig_in_bgc.head(10).to_string())

        coloc_results[sp] = sig_in_bgc

    return coloc_results


def compile_gwas_summary(species_list, all_results, enrichment_results):
    """
    Compile a cross-species summary of Pan-GWAS results.

    Parameters
    ----------
    species_list : list of str
    all_results : dict
        {(species, contrast, layer): pd.DataFrame}
    enrichment_results : dict
        {species: pd.DataFrame}

    Returns
    -------
    pd.DataFrame
        Cross-species summary table.
    """
    cross_species_rows = []

    for sp in species_list:
        row = {'species': sp}

        n_sig_pav = 0
        n_sig_cnv = 0
        contrasts_tested = set()
        for (s, contrast, layer), res in all_results.items():
            if s != sp:
                continue
            contrasts_tested.add(contrast)
            if 'significant' not in res.columns:
                continue
            n = int(res['significant'].sum())
            if layer == 'pav':
                n_sig_pav += n
            elif layer == 'cnv':
                n_sig_cnv += n

        row['n_contrasts'] = len(contrasts_tested)
        row['n_significant_PAV'] = n_sig_pav
        row['n_significant_CNV'] = n_sig_cnv

        # Top enriched terms — prefer FDR-significant; fall back to top by p-value
        if sp in enrichment_results and not enrichment_results[sp].empty:
            df_enr = enrichment_results[sp]
            sig_terms = df_enr[df_enr['significant'] == True] if 'significant' in df_enr.columns else df_enr.iloc[0:0]

            if len(sig_terms) > 0:
                top = sig_terms.nsmallest(5, 'qvalue')
                row['top_enriched (FDR<0.05)'] = '; '.join(top['term'].astype(str).tolist())
                row['top_p_values'] = '; '.join(f'{p:.1e}' for p in top['pvalue'])
            else:
                # No FDR-sig — show top 5 by raw p-value, flagged
                top = df_enr.nsmallest(5, 'pvalue') if 'pvalue' in df_enr.columns else df_enr.head(0)
                if len(top) > 0:
                    row['top_enriched (FDR<0.05)'] = 'none survive FDR'
                    row['top_p_values'] = '; '.join(
                        f'{t} (p={p:.1e})' for t, p in zip(top['term'].astype(str), top['pvalue'])
                    )
                else:
                    row['top_enriched (FDR<0.05)'] = 'N/A'
                    row['top_p_values'] = ''
        else:
            row['top_enriched (FDR<0.05)'] = 'N/A (no significant hits)'
            row['top_p_values'] = ''

        cross_species_rows.append(row)

    return pd.DataFrame(cross_species_rows)
