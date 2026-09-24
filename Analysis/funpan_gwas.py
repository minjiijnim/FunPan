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


# =============================================================================
# ASSOCIATION TESTING (PER-SPECIES)
# =============================================================================


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================


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
        Uses LMM with EMMA/EMMAX approach, eigendecomposes the kinship matrix,
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
    # FALLBACK PATH (no kinship, logistic regression or Fisher's exact)
    # =========================================================================
    else:
        print(f"\n  No kinship matrix, falling back to logistic/Fisher")
        
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


# print("Diagnostic plotting functions loaded.")

# =============================================================================
# SECTION 6.4: CAZy FAMILY EXPANSION PAN-GWAS
# =============================================================================
# Test if phenotype groups differ in CAZyme repertoire SIZE (not individual families)
# Hypothesis: Industrial/pathogenic strains may have expanded specific enzyme classes
#
# Aggregate per-genome counts:
#   - Total CAZymes per genome
#   - Count per CAZy class (GH, GT, PL, CE, AA, CBM)

import re
from scipy.stats import mannwhitneyu
from collections import defaultdict


from pathlib import Path
from collections import Counter, defaultdict


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

        # Standard matrices, named like {sp}_{type}.tsv in NB1_Results/{sp}/
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

            # Minimum detectable odds ratio at 80% power, take best-case across
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
    species_dir = str(species_path(genus, species))
    as_base = os.path.join(species_dir, 'antismash_output')
    og_tsvs = glob.glob(os.path.join(species_dir, 'orthofinder_output',
                                      '*', 'Orthogroups', 'Orthogroups.tsv'))
    if not og_tsvs or not os.path.isdir(as_base):
        return pd.DataFrame()

    # Load OrthoFinder, build gene_id -> orthogroup lookup
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

        # Load BGC-to-OG mapping, build on the fly if not cached
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

        # Top enriched terms, prefer FDR-significant; fall back to top by p-value
        if sp in enrichment_results and not enrichment_results[sp].empty:
            df_enr = enrichment_results[sp]
            sig_terms = df_enr[df_enr['significant'] == True] if 'significant' in df_enr.columns else df_enr.iloc[0:0]

            if len(sig_terms) > 0:
                top = sig_terms.nsmallest(5, 'qvalue')
                row['top_enriched (FDR<0.05)'] = '; '.join(top['term'].astype(str).tolist())
                row['top_p_values'] = '; '.join(f'{p:.1e}' for p in top['pvalue'])
            else:
                # No FDR-sig, show top 5 by raw p-value, flagged
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


# =============================================================================
# NB2 FIGURES, TABLES AND FIRST-RUN SUPPORT
# =============================================================================

def _label_feature_focus(name):
    """Display-only shortening of a feature name (the data itself is untouched).

    BGC features carry an assembly alias between the accession and the scaffold,
    and the alias varies -- ASM221148v2, CSR3, and others. It is redundant with
    the GCA/GCF accession right in front of it and only eats label width, so ANY
    alias is dropped for plotting:
        GCA_002211485.2_ASM221148v2_scaffold_22.region007.gbk_region_7
          -> GCA_002211485.2_scaffold_22.region007.gbk_region_7
        GCA_026261755.1_CSR3_scaffold_8.region009.gbk_region_9
          -> GCA_026261755.1_scaffold_8.region009.gbk_region_9
    Names without an accession + scaffold (PAV / CNV orthogroups) pass through.
    """
    return re.sub(r'^(GC[AF]_\d+\.\d+)_.+?(?=_(?:scaffold|contig))', r'\1',
                  str(name), flags=re.IGNORECASE)

def _safe_load_assoc_focus(path):
    if not path.exists():
        return None, 'missing'
    if path.stat().st_size < 10:
        return None, 'empty'
    try:
        df = pd.read_csv(path, sep='\t')
    except Exception:
        return None, 'empty'
    if df.empty or 'beta' not in df.columns or 'pvalue' not in df.columns:
        return None, 'empty'
    df = df.dropna(subset=['beta', 'pvalue']).copy()
    if df.empty:
        return None, 'empty'
    return df, 'ok'


def plot_focus_volcano_grid(species_list, nb2_results, out_path=None,
                            fontsize=30, dpi=800,
                            fig_w_per_layer=10.0, fig_h_per_row=10.0,
                            n_top_labels=5, q_threshold=0.05, y_max_clip=50,
                            layers=('pav', 'cnv', 'bgc_pav'),
                            exclude_contrasts=('environmental_vs_rest',),
                            jitter_x=0.06, jitter_y=0.06, jitter_seed=0,
                            y_min_half_range=1.0, gap_trigger_factor=1.25,
                            upper_height_ratio=0.30, overflow_padding=0.05,
                            species_colors=None):
    """Volcano matrix of the phenotype-of-interest contrasts.

    Rows are (species, contrast) pairs found in
    ``{nb2_results}/{sp}/pangwas_results/{layer}_assoc_{contrast}.tsv``; columns
    are ``layers``. Contrasts in ``exclude_contrasts`` are skipped.

    - Jitter spreads stacked points and is visual only.
    - The x-axis is symmetric around beta = 0.
    - The BH cutoff line is centred on the y-axis. When significant points sit
      far above that range, the panel is split into a lower and a smaller upper
      sub-axis, each with accurate -log10(p) ticks, and a wave marks the gap.
    - ``gap_trigger_factor`` sets when the split happens (y_max above factor x
      the centred upper limit); ``upper_height_ratio`` is the share of panel
      height given to the upper sub-axis.

    Returns the Figure, or None when no contrast is found.
    """
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D

    SPECIES = list(species_list)
    NB2_RESULTS = Path(nb2_results)
    FS_FOCUS = fontsize
    DPI_FOCUS = dpi
    FIG_W_PER_LAYER_FOCUS = fig_w_per_layer
    FIG_H_PER_ROW_FOCUS = fig_h_per_row
    N_TOP_LABELS_FOCUS = n_top_labels
    Q_THRESHOLD_FOCUS = q_threshold
    Y_MAX_CLIP_FOCUS = y_max_clip
    LAYERS_FOCUS = list(layers)
    EXCLUDE_CONTRASTS = list(exclude_contrasts)
    JITTER_X = jitter_x
    JITTER_Y = jitter_y
    JITTER_SEED = jitter_seed
    Y_MIN_HALF_RANGE = y_min_half_range
    GAP_TRIGGER_FACTOR = gap_trigger_factor
    UPPER_HEIGHT_RATIO = upper_height_ratio
    OVERFLOW_PADDING = overflow_padding
    fig = None

    # ----- Pre-compute panel data and decide broken vs single axis -----
    panel_specs = []   # list of dicts: sp, contrast, df, status, sig_mask, cutoff_y,
                       # p_cutoff, y_data_max, ylim_lo, ylim_lower_top, broken,
                       # upper_lo, upper_hi
    row_keys_focus = []
    for sp in SPECIES:
        sp_dir = NB2_RESULTS / sp / 'pangwas_results'
        if not sp_dir.exists():
            continue
        contrasts_seen = set()
        for layer in LAYERS_FOCUS:
            for f in sorted(sp_dir.glob(f'{layer}_assoc_*.tsv')):
                c = f.stem.replace(f'{layer}_assoc_', '')
                if c in EXCLUDE_CONTRASTS:
                    continue
                contrasts_seen.add(c)
        for c in sorted(contrasts_seen):
            row_keys_focus.append((sp, c))

    if not row_keys_focus:
        print('No phenotype-of-interest contrasts found.')
    else:
        nrows = len(row_keys_focus)
        ncols = len(LAYERS_FOCUS)
        fig = plt.figure(figsize=(FIG_W_PER_LAYER_FOCUS * ncols,
                                  FIG_H_PER_ROW_FOCUS * nrows))
        outer_gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                                      hspace=0.2, wspace=0.1)
        SP_COLORS = species_colors if species_colors is not None else SPECIES_COLORS
        rng_master = np.random.default_rng(JITTER_SEED)

        for ri, (sp, contrast) in enumerate(row_keys_focus):
            color = SP_COLORS.get(sp, '#444444')
            for ci, layer in enumerate(LAYERS_FOCUS):
                spec = outer_gs[ri, ci]
                path = NB2_RESULTS / sp / 'pangwas_results' / f'{layer}_assoc_{contrast}.tsv'
                df, status = _safe_load_assoc_focus(path)

                # Empty / missing panel
                if status != 'ok':
                    ax = fig.add_subplot(spec)
                    msg = 'no TSV' if status == 'missing' else 'no features tested'
                    ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                            ha='center', va='center',
                            fontsize=FS_FOCUS - 2, color='#666')
                    ax.set_xticks([]); ax.set_yticks([])
                    if ri == 0:
                        ax.set_title(layer.upper(), fontsize=FS_FOCUS + 2,
                                     fontweight='bold', color='#222')
                    if ci == 0:
                        ax.text(-0.27, 0.5,
                                f"A. {sp}\n{contrast.replace('_', ' ')}",
                                transform=ax.transAxes,
                                rotation=90, ha='center', va='center',
                                fontsize=FS_FOCUS + 1, fontweight='bold', color=color)
                    continue

                df['pvalue_clip'] = df['pvalue'].clip(lower=1e-300)
                df['neglog10p']   = -np.log10(df['pvalue_clip'])
                if Y_MAX_CLIP_FOCUS is not None:
                    df['neglog10p'] = df['neglog10p'].clip(upper=Y_MAX_CLIP_FOCUS)
                if 'significant' in df.columns:
                    sig_mask = df['significant'].astype(bool)
                elif 'qvalue' in df.columns:
                    sig_mask = df['qvalue'].fillna(1.0) < Q_THRESHOLD_FOCUS
                else:
                    sig_mask = pd.Series(False, index=df.index)

                if sig_mask.any():
                    p_cutoff = float(df.loc[sig_mask, 'pvalue'].max())
                    cutoff_y = -np.log10(max(p_cutoff, 1e-300))
                    if Y_MAX_CLIP_FOCUS is not None:
                        cutoff_y = min(cutoff_y, Y_MAX_CLIP_FOCUS)
                else:
                    p_cutoff = None; cutoff_y = None

                y_data_max = float(df['neglog10p'].max())

                # Decide if we need a broken axis (real gap between regular cloud
                # and a separate cluster of extreme outliers).
                broken = False
                if cutoff_y is not None:
                    half = max(cutoff_y, Y_MIN_HALF_RANGE)
                    ylim_lo = max(0.0, cutoff_y - half)
                    lower_top = cutoff_y + half  # cutoff centered top
                    if y_data_max > lower_top * GAP_TRIGGER_FACTOR:
                        broken = True
                        # Define the overflow band exactly around the high points.
                        overflow_vals = df.loc[df['neglog10p'] > lower_top, 'neglog10p']
                        upper_lo = float(overflow_vals.min())
                        upper_hi = float(overflow_vals.max())
                        span = max(upper_hi - upper_lo, 1.0)
                        upper_lo -= span * OVERFLOW_PADDING
                        upper_hi += span * OVERFLOW_PADDING
                else:
                    half = max(Y_MIN_HALF_RANGE, y_data_max / 2 if y_data_max > 0 else Y_MIN_HALF_RANGE)
                    ylim_lo = 0.0
                    lower_top = max(y_data_max * 1.10, Y_MIN_HALF_RANGE * 2)

                # ----- Build axes for this panel -----
                if broken:
                    inner = gridspec.GridSpecFromSubplotSpec(
                        2, 1, subplot_spec=spec,
                        height_ratios=[UPPER_HEIGHT_RATIO, 1.0 - UPPER_HEIGHT_RATIO],
                        hspace=0.07)
                    ax_lower = fig.add_subplot(inner[1])
                    ax_upper = fig.add_subplot(inner[0], sharex=ax_lower)
                else:
                    ax_lower = fig.add_subplot(spec)
                    ax_upper = None

                # Column header on the top row, on whichever axis is uppermost
                if ri == 0:
                    top_ax = ax_upper if ax_upper is not None else ax_lower
                    top_ax.set_title(layer.upper(), fontsize=FS_FOCUS + 2,
                                     fontweight='bold', color='#222', pad=8)
                # Row label on the left column, on the lower axis
                if ci == 0:
                    ax_lower.text(-0.27, 0.5,
                                  f"A. {sp}\n{contrast.replace('_', ' ')}",
                                  transform=ax_lower.transAxes,
                                  rotation=90, ha='center', va='center',
                                  fontsize=FS_FOCUS + 1, fontweight='bold', color=color)

                # ----- Jitter (visual only; raw values unchanged) -----
                n = len(df)
                jx = rng_master.uniform(-JITTER_X, JITTER_X, size=n) if JITTER_X else np.zeros(n)
                jy = rng_master.uniform(-JITTER_Y, JITTER_Y, size=n) if JITTER_Y else np.zeros(n)
                beta_plot = df['beta'].values + jx
                nlp_plot  = df['neglog10p'].values + jy

                # Plot on both sub-axes; matplotlib's clip will handle ranges.
                for axx in ([ax_lower] if not broken else [ax_lower, ax_upper]):
                    axx.scatter(beta_plot[~sig_mask.values], nlp_plot[~sig_mask.values],
                                s=12, color='#9e9e9e', alpha=0.7, edgecolor='none')
                    axx.scatter(beta_plot[sig_mask.values], nlp_plot[sig_mask.values],
                                s=30, color=color, alpha=0.78,
                                edgecolor='white', linewidth=0.4)
                    if cutoff_y is not None:
                        axx.axhline(cutoff_y, color='#3f3f3f', linestyle='--',
                                    linewidth=1.2, alpha=0.9)
                    axx.axvline(0, color='#7a7a7a', linestyle=':',
                                linewidth=1.0, alpha=0.85)
                    axx.grid(True, linestyle=':', alpha=0.55, color='#8f8f8f')
                    axx.tick_params(axis='both', labelsize=FS_FOCUS - 3)

                # X-limits: symmetric around beta = 0, shared across both sub-axes
                x_max_abs = float(np.nanmax(np.abs(df['beta'].values)))
                x_max_abs = max(x_max_abs + abs(JITTER_X), 0.5)
                xlim = (-x_max_abs * 1.10, x_max_abs * 1.10)
                ax_lower.set_xlim(xlim)
                if ax_upper is not None:
                    ax_upper.set_xlim(xlim)

                # Y-limits
                ax_lower.set_ylim(ylim_lo, lower_top)
                if ax_upper is not None:
                    ax_upper.set_ylim(upper_lo, upper_hi)

                # Hide the spines at the broken edges and draw the wave there
                if ax_upper is not None:
                    ax_lower.spines['top'].set_visible(False)
                    ax_upper.spines['bottom'].set_visible(False)
                    ax_upper.tick_params(labeltop=False, top=False, bottom=False)
                    ax_upper.set_xticklabels([])

                    # Draw TWO PARALLEL waves in figure coords across the gap.
                    # Both waves use the same +sin phase (same direction) so peaks
                    # of one align with peaks of the other -- they're parallel.
                    # The vertical offset keeps them visibly separated.
                    import matplotlib.lines as _mlines
                    bbox_up = ax_upper.get_position()
                    bbox_lo = ax_lower.get_position()
                    gap_lo = bbox_lo.y1
                    gap_hi = bbox_up.y0
                    gap_y  = (gap_lo + gap_hi) / 2.0
                    gap_h  = max(gap_hi - gap_lo, 1e-6)
                    x_left  = max(bbox_up.x0, bbox_lo.x0)
                    x_right = min(bbox_up.x1, bbox_lo.x1)
                    _xs   = np.linspace(x_left, x_right, 240)
                    _phase = np.sin(12 * np.pi * (_xs - x_left) / (x_right - x_left))
                    _off  = gap_h * 0.12           # separation between the two parallel waves
                    _amp  = gap_h * 0.14           # shared amplitude
                    _y_top = (gap_y + _off) + _amp * _phase           # upper wave
                    _y_bot = (gap_y - _off) + _amp * _phase           # lower wave, SAME phase = parallel
                    for _ys in (_y_top, _y_bot):
                        fig.add_artist(_mlines.Line2D(_xs, _ys, color='white',
                                                      linewidth=4.0, solid_capstyle='round',
                                                      zorder=10, transform=fig.transFigure))
                        fig.add_artist(_mlines.Line2D(_xs, _ys, color="#0000004D",
                                                      linewidth=1, solid_capstyle='round',
                                                      zorder=11, transform=fig.transFigure))

                # Top-hit labels (use actual beta / -log10p; matplotlib will place them
                # on whichever sub-axis the point falls in, automatically clipped)
                if sig_mask.any() and N_TOP_LABELS_FOCUS > 0:
                    top = (df[sig_mask]
                           .assign(absbeta=df.loc[sig_mask, 'beta'].abs())
                           .sort_values('absbeta', ascending=False)
                           .head(N_TOP_LABELS_FOCUS))
                    STEP_PTS = 27         # initial stagger between pairs (pts)
                    MIN_SEP_PTS = 16     # min vertical separation between any two labels (pts)
                    _placed_y_fig = []    # fig-y positions of labels placed in THIS panel
                    _ppx = fig.dpi / 72.0
                    _fh = fig.get_figheight()
                    _min_sep_figy = MIN_SEP_PTS / (_fh * 72.0)
                    for _i, (_, r) in enumerate(top.iterrows()):
                        target_ax = ax_lower
                        if ax_upper is not None and r['neglog10p'] >= upper_lo:
                            target_ax = ax_upper
                        _pair = _i // 2
                        _side = _i % 2
                        if r['beta'] >= 0:
                            _xoff = -3; _ha = 'right'
                        else:
                            _xoff =  3; _ha = 'left'
                        if _side == 0:
                            _yoff =  6 + _pair * STEP_PTS
                            _va = 'bottom'
                            _dir = +1
                        else:
                            _yoff = -6 - _pair * STEP_PTS
                            _va = 'top'
                            _dir = -1
                        _disp_dot = target_ax.transData.transform((r['beta'], r['neglog10p']))
                        _disp = (_disp_dot[0] + _xoff * _ppx, _disp_dot[1] + _yoff * _ppx)
                        _bb = target_ax.get_position()
                        _figxy = fig.transFigure.inverted().transform(_disp)
                        _mx, _my = 0.004, 0.005
                        # Flip side once if initially out of bounds
                        if _figxy[1] > _bb.y1 - _my:
                            _yoff = -abs(_yoff); _va = 'top'; _dir = -1
                            _disp = (_disp_dot[0] + _xoff*_ppx, _disp_dot[1] + _yoff*_ppx)
                            _figxy = fig.transFigure.inverted().transform(_disp)
                        elif _figxy[1] < _bb.y0 + _my:
                            _yoff = abs(_yoff); _va = 'bottom'; _dir = +1
                            _disp = (_disp_dot[0] + _xoff*_ppx, _disp_dot[1] + _yoff*_ppx)
                            _figxy = fig.transFigure.inverted().transform(_disp)
                        # Greedy collision: push along _dir until far enough from
                        # every label already placed on this panel
                        _y = _figxy[1]
                        for _attempt in range(20):
                            _collide = any(abs(_y - yp) < _min_sep_figy for yp in _placed_y_fig)
                            if not _collide:
                                break
                            _y += _dir * _min_sep_figy
                            # If we've wandered out of the panel along _dir, flip
                            if _y > _bb.y1 - _my or _y < _bb.y0 + _my:
                                _dir = -_dir
                                _y = _figxy[1] + _dir * _min_sep_figy
                                _va = 'bottom' if _dir > 0 else 'top'
                        _figx = min(max(_figxy[0], _bb.x0 + _mx), _bb.x1 - _mx)
                        _figy = min(max(_y, _bb.y0 + _my), _bb.y1 - _my)
                        _placed_y_fig.append(_figy)
                        fig.text(_figx, _figy, _label_feature_focus(r['feature']),
                                 ha=_ha, va=_va,
                                 fontsize=FS_FOCUS - 4, color=color, alpha=0.95,
                                 fontweight='bold', zorder=200)

                n_total = len(df); n_sig = int(sig_mask.sum())
                note = f'n={n_total:,}  •  sig={n_sig:,}'
                if p_cutoff is not None:
                    note += f'\nBH p≤{p_cutoff:.1e}'
                # Place the note on the lower axis (always present)
                ax_lower.text(0.02, 0.02, note, transform=ax_lower.transAxes,
                              ha='left', va='bottom', fontsize=FS_FOCUS - 3,
                              color='#1f1f1f',
                              bbox=dict(boxstyle='round,pad=0.25', fc='white',
                                        ec='none', alpha=0.85))

                # Axis labels: x on the bottom (lower axis) of bottom-most row,
                # y label on whichever lower axis is in column 0
                if ri == nrows - 1:
                    ax_lower.set_xlabel(r'$\beta$', fontsize=FS_FOCUS - 1)
                if ci == 0:
                    ax_lower.set_ylabel(r'$-\log_{10}(p)$', fontsize=FS_FOCUS - 1)
                    if ax_upper is not None:
                        ax_upper.set_ylabel(r'$-\log_{10}(p)$', fontsize=FS_FOCUS - 1)

        jitter_note = (f"  •  jitter ±{JITTER_X:g}β / ±{JITTER_Y:g} -log10(p)"
                       if (JITTER_X or JITTER_Y) else '')
        fig.suptitle("Phenotype-of-interest volcano matrix",
                     fontsize=FS_FOCUS + 1, fontweight='bold', y=0.995)

        if out_path:
            fig.savefig(out_path, dpi=DPI_FOCUS, bbox_inches='tight')
            print(f'Saved: {out_path}')
        print(f"Rows (species x phenotype-of-interest): {nrows}   "
              f"Columns (layers): {ncols}")

    return fig


def association_summary_table(all_results):
    """Tested and significant counts, plus genomic inflation, per test.

    Returns
    -------
    pd.DataFrame
        One row per (species, contrast, feature layer), sorted.
    """
    rows = []
    for (sp, contrast, layer), res in all_results.items():
        n_tested = len(res)
        n_sig = int(res.get('significant', pd.Series(dtype=bool)).sum()) if 'significant' in res.columns else 0
        lambda_gc = calculate_genomic_inflation(res['pvalue'].dropna().values) \
            if 'pvalue' in res.columns and len(res['pvalue'].dropna()) > 0 else np.nan
        rows.append({
            'species': sp,
            'contrast': contrast,
            'feature_layer': layer,
            'n_tested': n_tested,
            'n_significant': n_sig,
            'lambda_gc': round(lambda_gc, 3) if not np.isnan(lambda_gc) else np.nan,
        })
    return (pd.DataFrame(rows)
            .sort_values(['species', 'contrast', 'feature_layer'])
            .reset_index(drop=True))


def plot_significance_heatmap(all_results, out_path=None, figsize=(6, 3.5),
                              dpi=800, cmap='YlOrRd'):
    """Heatmap of significant associations per species and feature layer.

    Counts are summed over contrasts.

    Returns
    -------
    fig, heatmap_df : Figure (None if there are no results) and the counts.
    """
    import seaborn as sns

    sig_counts = {}
    for (sp, contrast, layer), res in all_results.items():
        key = (sp, layer)
        n = int(res['significant'].sum()) if 'significant' in res.columns else 0
        sig_counts[key] = sig_counts.get(key, 0) + n
    if not sig_counts:
        return None, pd.DataFrame()

    heatmap_df = pd.Series(sig_counts).unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(heatmap_df, annot=True, fmt='d', cmap=cmap, ax=ax, linewidths=0.5)
    ax.set_title('Number of Significant Associations per Feature Layer', fontweight='bold')
    ax.set_xlabel('Feature Layer')
    ax.set_ylabel('Species')
    plt.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')
    return fig, heatmap_df


def compute_power_curves(species_list, species_contrasts, nb2_results,
                         fdr_threshold=0.05, or_range=None, gene_freq=0.2,
                         n_sim=2000, force=False):
    """Fisher's exact test power against odds ratio for every contrast.

    Cached to ``{nb2_results}/power_curves.tsv``; pass ``force=True`` to rerun.

    Returns
    -------
    pd.DataFrame
        Columns species, contrast, n_case, n_control, odds_ratio, power.
    """
    cache = Path(nb2_results) / 'power_curves.tsv'
    if or_range is None:
        or_range = np.arange(1.5, 15.5, 0.5)

    if not force and cache.exists():
        curves_df = pd.read_csv(cache, sep='\t')
        print(f'Loaded cached power curves from {cache}')
        return curves_df

    plot_configs = []
    for sp in species_list:
        if sp in species_contrasts:
            for cname, pheno in species_contrasts[sp].items():
                ph = pheno['phenotype'] if hasattr(pheno, 'columns') else pheno
                n_c = int(ph.sum())
                n_ctrl = int((~ph.astype(bool)).sum())
                plot_configs.append((sp, cname, n_c, n_ctrl))

    print('Computing power curves...')
    rows = []
    for sp, cname, n_c, n_ctrl in plot_configs:
        for or_val in or_range:
            p = power_fisher(n_case=n_c, n_control=n_ctrl,
                             gene_freq=gene_freq, odds_ratio=float(or_val),
                             alpha=fdr_threshold, n_sim=n_sim)
            rows.append({'species': sp, 'contrast': cname,
                         'n_case': n_c, 'n_control': n_ctrl,
                         'odds_ratio': or_val, 'power': p})
    curves_df = pd.DataFrame(rows)
    curves_df.to_csv(cache, sep='\t', index=False)
    print(f'Saved power curves to {cache}')
    return curves_df


def plot_power_curves(curves_df, out_path=None, figsize=(8, 5), dpi=150):
    """Power against odds ratio, one line per species and contrast."""
    fig, ax = plt.subplots(figsize=figsize)
    for (sp, cname), grp in curves_df.groupby(['species', 'contrast']):
        n_c = grp['n_case'].iloc[0]; n_ctrl = grp['n_control'].iloc[0]
        ax.plot(grp['odds_ratio'], grp['power'],
                label=f'{sp} ({cname}, {n_c}:{n_ctrl})', linewidth=1.5)
    ax.axhline(0.8, color='grey', linestyle='--', alpha=0.6, label='80% power')
    ax.set_xlabel('Odds Ratio')
    ax.set_ylabel('Statistical Power')
    ax.set_title("Power Analysis: Fisher's Exact Test per Species/Contrast")
    ax.legend(fontsize=7, loc='lower right')
    ax.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')
    return fig


def plot_underpowered_contrasts(nb2_results, out_path=None, fontsize=14, dpi=800,
                                figsize=(6.0, 5.0),
                                exclude_contrasts=('environmental_vs_rest',),
                                underpowered_threshold=10.0):
    """Minimum detectable odds ratio at 80% power per contrast.

    Reads ``{nb2_results}/power_analysis_results.tsv``. Bars above
    ``underpowered_threshold`` are grey, the rest green.

    Returns
    -------
    fig, table : Figure and the plotted rows, or (None, None) if there is
    nothing to plot.
    """
    FS_UP = fontsize
    path = Path(nb2_results) / 'power_analysis_results.tsv'
    if not path.exists():
        print(f'No power-analysis TSV at {path}; run the power analysis first.')
        return None, None

    pa = pd.read_csv(path, sep='\t')
    mde = pa[pa['target_OR'] == 'min_OR@80%power'].copy()
    mde['min_OR'] = pd.to_numeric(mde['power'], errors='coerce')
    mde = mde.dropna(subset=['min_OR'])
    mde = mde[~mde['contrast'].isin(list(exclude_contrasts))].copy()
    if mde.empty:
        print('No phenotype-of-interest contrasts left after filtering.')
        return None, None

    mde['label'] = ('A. ' + mde['species'] + '\n'
                    + mde['contrast'].str.replace('_', ' ', regex=False)
                    + '\n(ratio ' + mde['ratio'] + ')')
    mde = mde.sort_values('min_OR')
    colors = ['#888888' if v > underpowered_threshold else '#2ca02c'
              for v in mde['min_OR']]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.barh(mde['label'], mde['min_OR'], color=colors,
                   edgecolor='white', linewidth=0.6)
    xmax = max(mde['min_OR'].max(), underpowered_threshold) * 1.05
    for b, v in zip(bars, mde['min_OR']):
        ax.text(v + xmax * 0.01, b.get_y() + b.get_height() / 2,
                f'{v:.1f}', va='center', ha='left', fontsize=FS_UP - 2,
                color='#888888' if v > underpowered_threshold else '#2ca02c',
                fontweight='bold')
    ax.axvline(underpowered_threshold, color='grey', linestyle='--', linewidth=1.0,
               label=f'OR = {underpowered_threshold:g} (underpowered above)')
    ax.set_xlim(0, xmax * 1.18)
    ax.set_xlabel('Min. detectable OR at 80% power', fontsize=FS_UP)
    ax.set_title('Power Analysis Results', fontsize=FS_UP + 2, fontweight='bold')
    ax.tick_params(axis='x', labelsize=FS_UP - 1)
    ax.tick_params(axis='y', labelsize=FS_UP - 2)
    ax.legend(frameon=True, fontsize=FS_UP - 2, loc='lower right', framealpha=0.95)
    ax.grid(True, axis='x', linestyle=':', alpha=0.5)
    plt.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {out_path}')

    table = mde[['species', 'contrast', 'ratio', 'min_OR']].reset_index(drop=True)
    table['underpowered'] = table['min_OR'] > underpowered_threshold
    return fig, table


def nb2_preflight(species_list, nb0_results, nb1_results, nb2_results, verbose=True):
    """Check NB2's inputs exist and create the output directories.

    Looks for NB0's phenotype table and, per species, the NB1 matrices the
    association tests read. Nothing is computed or written beyond directories.

    Returns
    -------
    pd.DataFrame
        One row per species with a boolean per required NB1 input.
    """
    nb1_results = Path(nb1_results)
    nb2_results = Path(nb2_results)
    nb2_results.mkdir(parents=True, exist_ok=True)

    needed = ['pav', 'cnv', 'kinship', 'bgc_pav', 'gcf_pav', 'gcf_cnv', 'og_consensus']
    rows = []
    for sp in species_list:
        (nb2_results / sp / 'pangwas_results').mkdir(parents=True, exist_ok=True)
        row = {'species': sp}
        for key in needed:
            row[key] = (nb1_results / sp / f'{sp}_{key}.tsv').exists()
        rows.append(row)
    status = pd.DataFrame(rows).set_index('species')

    pheno = Path(nb0_results) / 'phenotype_classified_for_gwas.csv'
    if verbose:
        print('NB2 inputs')
        print('=' * 60)
        print(f'  NB0 phenotype table: {"found" if pheno.exists() else "MISSING, run NB0 first"}')
        print(f'  -> {pheno}')
        print(f'  NB1 matrices under: {nb1_results}')
        print()
        print(status.replace({True: 'ok', False: 'MISSING'}).to_string())
        missing = int((~status).sum().sum())
        print()
        if missing or not pheno.exists():
            print(f'{missing} missing NB1 input(s); run NB1 before the sections that need them.')
        else:
            print('All inputs present.')
    return status


def print_nb2_summary(all_results, species_list, nb2_results):
    """Print every NB2 output with its full path, and the significant-hit counts.

    Returns
    -------
    counts : pd.DataFrame
        Significant associations per (species, contrast) and feature layer.
    files : pd.DataFrame
        Every expected output with its full path and whether it exists.
    """
    import glob

    nb2_results = Path(nb2_results)
    summary = association_summary_table(all_results)
    counts = (summary.pivot_table(index=['species', 'contrast'],
                                  columns='feature_layer',
                                  values='n_significant', aggfunc='sum')
              .fillna(0).astype(int))
    counts.columns.name = None

    top_level = ['power_analysis_results.tsv', 'power_curves.tsv',
                 'power_analysis_curves.png', 'underpowered_contrasts_nonenv.png',
                 'volcano_grid_focus.png', 'summary_heatmap.png']
    per_species = ['{sp}/phenotypes/pheno_*.tsv',
                   '{sp}/pangwas_results/*_assoc_*.tsv',
                   '{sp}/{sp}_functional_enrichment.csv',
                   '{sp}/{sp}_sig_ogs_in_bgc_detail.csv']

    rows = []
    for name in top_level:
        path = nb2_results / name
        rows.append({'species': 'all', 'path': str(path), 'n': int(path.exists()),
                     'exists': path.exists()})
    for sp in species_list:
        for pattern in per_species:
            path = nb2_results / pattern.format(sp=sp)
            if '*' in path.name:
                n = len(glob.glob(str(path)))
                rows.append({'species': sp, 'path': str(path), 'n': n, 'exists': n > 0})
            else:
                rows.append({'species': sp, 'path': str(path), 'n': int(path.exists()),
                             'exists': path.exists()})
    files = pd.DataFrame(rows)

    print(f'{"=" * 70}')
    print('NB2 Pan-GWAS complete.')
    print(f'All outputs saved to: {nb2_results}')
    for _, r in files.iterrows():
        mark = ' ' if r['exists'] else '!'
        extra = f"  ({r['n']} files)" if '*' in r['path'] else ''
        print(f"  {mark} {r['path']}{extra}")
    n_missing = int((~files['exists']).sum())
    if n_missing:
        print(f'\n  ! = missing ({n_missing} path(s))')
    print('=' * 70)
    print('\nSignificant associations per contrast and feature layer:')

    return counts, files
