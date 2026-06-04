#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_phylo.py
===============
Phylogenetic signal analysis, ancestral niche reconstruction, gain/loss events,
Mash clustering, and rare genome characterization.

Extracted from NB0_Analysis.ipynb (cells 48, 62, 72) and NB4_Analysis.ipynb.
"""

from funpan_utils import *
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter
from scipy import stats
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# try:
#     from Bio import SeqIO, Phylo
# except ImportError:
#     SeqIO = None
#     Phylo = None



# --- From NB0 Cell 48 ---
# =============================================================================
# SECTION 9: PHYLOGENETIC TREE WITH PHENOTYPE COLORING
# =============================================================================
# Load IQ-TREE results and visualize the phylogeny with strains colored by phenotype
# Trees were built from Gubbins recombination-filtered SNP alignments using IQ-TREE

from Bio import Phylo
import re
import matplotlib.patches as mpatches

# =============================================================================
# CONFIGURATION - Change these to adjust visualization
# =============================================================================
# TRANSFORM = 'sqrt'  # Options: 'none', 'log', 'sqrt', or a float for power (e.g., 0.3)
# =============================================================================

def parse_iqtree_log(iqtree_file):
    """Parse IQ-TREE log file to extract key statistics."""
    stats = {}
    if not os.path.exists(iqtree_file):
        return stats
    
    with open(iqtree_file, 'r') as f:
        content = f.read()
    
    patterns = {
        'n_sequences': r'Input data: (\d+) sequences',
        'n_sites': r'Input data: \d+ sequences with (\d+) nucleotide sites',
        'parsimony_informative': r'Number of parsimony informative sites: (\d+)',
        'best_model': r'Best-fit model according to BIC: ([^\n]+)',
        'log_likelihood': r'Log-likelihood of the tree: ([-\d.]+)',
        'tree_length': r'Total tree length \(sum of branch lengths\): ([\d.]+)',
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            val = match.group(1).strip()
            try:
                stats[key] = int(val) if '.' not in val else float(val)
            except:
                stats[key] = val
    return stats


def clean_tip_labels(tree):
    """Clean up long genome file names to just accession IDs."""
    for clade in tree.get_terminals():
        match = re.search(r'(GC[AF]_\d+\.\d+)', clade.name)
        if match:
            clade.name = match.group(1)
    return tree


def apply_transform(depths, transform):
    """Apply transformation to depths."""
    min_depth = min(depths.values())
    shifted = {k: v - min_depth for k, v in depths.items()}
    
    if transform == 'none':
        return shifted
    
    max_depth = max(shifted.values())
    if max_depth == 0:
        return shifted
    
    if transform == 'log':
        transformed = {k: np.log1p(v) for k, v in shifted.items()}
    elif transform == 'sqrt':
        transformed = {k: np.sqrt(v) for k, v in shifted.items()}
    elif isinstance(transform, (int, float)):
        transformed = {k: np.power(v / max_depth, transform) * max_depth for k, v in shifted.items()}
    else:
        transformed = shifted
    
    return transformed


# =============================================================================
# PHENOTYPE COLOR SCHEMES
# =============================================================================

PHENOTYPE_COLORS = {
    'Human-pathogenic': '#e74c3c',      # Red
    'Clinical': '#e74c3c',               # Red
    'Plant-pathogenic': '#f1c40f',       # Yellow
    'Animal-pathogenic': '#f39c12',      # Orange
    'Environmental': '#27ae60',          # Green
    'Industrial': '#9b59b6',             # Purple
    'Industrial-trait': '#9b59b6',       # Purple
    'Lab': '#95a5a6',                    # Gray
    'Unknown': '#bdc3c7',                # Light gray
}

def get_color_for_phenotype(phenotype):
    if pd.isna(phenotype):
        return PHENOTYPE_COLORS['Unknown']
    return PHENOTYPE_COLORS.get(phenotype, PHENOTYPE_COLORS['Unknown'])


def get_phenotype_for_sample(sample_id, df_meta):
    """Get phenotype for a sample by matching accession."""
    norm_id = normalize_sample_id(sample_id)
    for idx, row in df_meta.iterrows():
        if 'Assembly Accession' in row.index:
            if normalize_sample_id(str(row['Assembly Accession'])) == norm_id:
                return row.get('Phenotype', 'Unknown')
        if 'sample_id' in row.index:
            if normalize_sample_id(str(row['sample_id'])) == norm_id:
                return row.get('Phenotype', 'Unknown')
    return 'Unknown'


# =============================================================================
# PHYLOGRAM PLOTTING WITH PHENOTYPE COLORING
# =============================================================================

def plot_phylogram(tree, ax, df_meta=None, title="", transform='none'):
    """Plot phylogram with optional branch length transformation and phenotype coloring."""
    terminals = tree.get_terminals()
    n_tips = len(terminals)
    
    tip_y = {t.name: i for i, t in enumerate(terminals)}
    
    raw_depths = tree.depths()
    depths = apply_transform(raw_depths, transform)
    max_depth = max(depths.values()) if depths.values() else 1
    
    def get_y(clade):
        if clade.is_terminal():
            return tip_y[clade.name]
        return np.mean([get_y(c) for c in clade.clades])
    
    # Get phenotype colors if metadata provided
    tip_colors = {}
    tip_phenotypes = []
    if df_meta is not None:
        for t in terminals:
            phenotype = get_phenotype_for_sample(t.name, df_meta)
            tip_phenotypes.append(phenotype)
            tip_colors[t.name] = get_color_for_phenotype(phenotype)
    else:
        for t in terminals:
            tip_phenotypes.append('Unknown')
            tip_colors[t.name] = 'steelblue'
    
    # Draw branches
    for clade in tree.find_clades():
        x = depths[clade]
        y = get_y(clade)
        
        if not clade.is_terminal() and clade.clades:
            children_y = [get_y(c) for c in clade.clades]
            ax.plot([x, x], [min(children_y), max(children_y)], color='#333333', lw=0.5)
        
        path = tree.get_path(clade)
        if len(path) >= 2:
            parent = path[-2]
            px = depths[parent]
            ax.plot([px, x], [y, y], color='#333333', lw=0.5)
    
    # Draw tip points and labels
    for t in terminals:
        tx = depths[t]
        ty = tip_y[t.name]
        color = tip_colors[t.name]
        ax.plot(tx, ty, 'o', ms=4, color=color, alpha=0.85, markeredgecolor='white', markeredgewidth=0.3)
        ax.text(tx + max_depth*0.01, ty, t.name, fontsize=6, va='center', color=color, alpha=0.9)
    
    ax.set_xlim(-max_depth*0.02, max_depth*1.15)
    ax.set_ylim(-1, n_tips)
    
    if transform == 'none':
        xlabel = "Substitutions per site"
    elif transform == 'log':
        xlabel = "Substitutions per site (log scale)"
    elif transform == 'sqrt':
        xlabel = "Substitutions per site (sqrt scale)"
    else:
        xlabel = f"Substitutions per site (^{transform})"
    ax.set_xlabel(xlabel, fontsize=9)
    
    ax.set_yticks([])
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return tip_phenotypes


# =============================================================================
# LOAD AND PLOT IQ-TREE RESULTS
# =============================================================================

# print("="*80)
# print("SECTION 9: PHYLOGENETIC TREE WITH PHENOTYPE COLORING")
# print("="*80)
# print(f"\nSpecies: {SPECIES_DISPLAY}")
# print("Input: Recombination-filtered SNP alignments from Gubbins")
# print("Method: IQ-TREE with ModelFinder + UFBoot (1000 replicates)")
# print(f"Transform: {TRANSFORM}")

# iqtree_dir = f"{SPECIES_DIR}/iqtree_output_snv"
# tree_file = os.path.join(iqtree_dir, "gubbins_tree.treefile")
# contree_file = os.path.join(iqtree_dir, "gubbins_tree.contree")
# iqtree_log = os.path.join(iqtree_dir, "gubbins_tree.iqtree")

# if not os.path.exists(tree_file):
#     print(f"\n[ERROR] Tree file not found: {tree_file}")
#     print("Please run IQ-TREE analysis first (run_iqtree_gubbins.sh)")
#     iqtree_tree = None
# else:
#     stats = parse_iqtree_log(iqtree_log)
#     
#     tree_to_load = contree_file if os.path.exists(contree_file) else tree_file
#     print(f"\nLoading tree from: {tree_to_load}")
#     iqtree_tree = Phylo.read(tree_to_load, "newick")
#     
    # Remove .ref duplicate tips (parsnp artifact) BEFORE cleaning names
#     ref_tips = [t for t in iqtree_tree.get_terminals() if t.name and t.name.endswith('.ref')]
#     for tip in ref_tips:
#         iqtree_tree.prune(tip)
#         print(f"  Pruned .ref duplicate: {tip.name}")
#     
#     iqtree_tree = clean_tip_labels(iqtree_tree)
#     if SPECIES_SHORT in ANI_EXCLUDED:
#         for exc_acc in ANI_EXCLUDED[SPECIES_SHORT]:
#             exc_base = re.match(r'(GC[AF]_\d+)', exc_acc).group(1)
#             for tip in list(iqtree_tree.get_terminals()):
#                 if exc_base in tip.name:
#                     iqtree_tree.prune(tip)
#                     print(f"  Pruned ANI-excluded: {tip.name}")
#     
#     terminals = iqtree_tree.get_terminals()
#     n_tips = len(terminals)
#     print(f"Tree contains {n_tips} tips (genomes)")
#     
#     df_meta = df_analysis if 'df_analysis' in dir() else df_classified
#     
#     print(f"\nIQ-TREE Statistics:")
#     print(f"  Sequences: {stats.get('n_sequences', 'N/A')}")
#     print(f"  SNP sites: {stats.get('n_sites', 'N/A'):,}" if isinstance(stats.get('n_sites'), int) else f"  SNP sites: {stats.get('n_sites', 'N/A')}")
#     print(f"  Parsimony informative: {stats.get('parsimony_informative', 'N/A'):,}" if isinstance(stats.get('parsimony_informative'), int) else f"  Parsimony informative: {stats.get('parsimony_informative', 'N/A')}")
#     print(f"  Best model: {stats.get('best_model', 'N/A')}")
#     print(f"  Log-L: {stats.get('log_likelihood', 'N/A')}")
#     
    # Plot phylogram
#     print("\n" + "-"*70)
#     print(f"Phylogram (transform={TRANSFORM})")
#     print("-"*70)
#     
#     fig_height = max(6, n_tips * 0.18)
#     fig, ax = plt.subplots(figsize=(10, fig_height))
#     
#     title = f"A. {SPECIES_SHORT} (n={n_tips}) | Model: {stats.get('best_model', 'N/A')}"
#     tip_phenotypes = plot_phylogram(iqtree_tree, ax, df_meta=df_meta, title=title, transform=TRANSFORM)
#     
    # Add legend
#     phenotype_counts = pd.Series(tip_phenotypes).value_counts()
#     legend_handles = []
#     for phenotype in phenotype_counts.index:
#         color = get_color_for_phenotype(phenotype)
#         count = phenotype_counts[phenotype]
#         legend_handles.append(mpatches.Patch(color=color, label=f'{phenotype} (n={count})'))
#     
#     ax.legend(handles=legend_handles, loc='lower right', fontsize=8, framealpha=0.9, title='Phenotype')
#     
#     plt.tight_layout()
#     output_path = f"{RESULTS_DIR}/phylogeny_phenotypes.png"
#     plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
#     print(f"Saved: {output_path}")
#     
    # Summary
#     print("\n" + "-"*70)
#     print("Phenotype Distribution in Tree")
#     print("-"*70)
#     print(f"\n{'Phenotype':<25} {'Count':>8} {'Percent':>10}")
#     print("-"*45)
#     for phenotype, count in phenotype_counts.items():
#         pct = 100 * count / n_tips
#         print(f"{phenotype:<25} {count:>8} {pct:>9.1f}%")
#     print(f"\n{'Total':<25} {n_tips:>8} {100.0:>9.1f}%")

# print("\n" + "="*80)
# print("TREE LOADING COMPLETE")
# print("="*80)


# --- From NB0 Cell 62 ---
# =============================================================================
# SECTION 10.1: MASH SETUP & DISTANCE COMPUTATION (PER-SPECIES)
# =============================================================================

import subprocess
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform, pdist
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples

# Configuration -- Mash defaults. Override via env var $MASH if your binary
# lives elsewhere; otherwise we auto-locate at module load time.
def _find_mash_binary():
    """Resolve the Mash executable: $MASH > PATH > known conda envs."""
    cand = os.environ.get('MASH')
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    on_path = shutil.which('mash')
    if on_path:
        return on_path
    for env in ('funpan', 'pipeline_test3', 'pipeline_test2'):
        p = f'/home/user/anaconda3/envs/{env}/bin/mash'
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return 'mash'  # let subprocess raise a useful error if it's still missing

MASH_BIN = _find_mash_binary()
KMER_SIZE = 21
SKETCH_SIZE = 1000

# Phenotype color map (reuse from earlier in notebook)
# PHENOTYPE_COLORS = {
#     'Human-pathogenic': '#e74c3c', 'Animal-pathogenic': '#d35400',
#     'Plant-pathogenic': '#27ae60', 'Industrial-trait': '#3498db',
#     'Environmental': '#9b59b6', 'Lab': '#95a5a6', 'Unknown': '#bdc3c7'
# }

# Load classified metadata for phenotype annotations
# meta_path = '/datadrive/Analysis/NB0_Results/qc_passed_classified_all.csv'
# df_meta_mash = pd.read_csv(meta_path)
# acc_to_phenotype = dict(zip(
#     df_meta_mash['Assembly Accession'].str.extract(r'(GC[AF]_\d+\.\d+)')[0],
#     df_meta_mash['Phenotype']
# ))

# ---------------------------------------------------------------------------
# Helper: Find elbow point using second derivative of inertia curve
# ---------------------------------------------------------------------------
def find_elbow(k_range, inertias):
    """Find elbow point using the maximum second derivative of inertia curve."""
    if len(k_range) < 3:
        return k_range[0]
    inertias = np.array(inertias)
    inertias_norm = (inertias - inertias.min()) / (inertias.max() - inertias.min() + 1e-10)
    d2 = np.diff(inertias_norm, n=2)
    elbow_idx = np.argmax(d2) + 1
    return k_range[elbow_idx]

# ---------------------------------------------------------------------------
# Helper: Robust k selection using elbow + silhouette with edge-case handling
# ---------------------------------------------------------------------------
def choose_optimal_k(k_range, inertias, sil_scores, n_samples):
    """
    Choose optimal k with robust handling of edge cases:
    1. If both methods agree → use that
    2. If silhouette is flat (range < 0.10) → trust elbow (no k is clearly better)
    3. If silhouette peak at high k (avg cluster < 5) → artifact, trust elbow
    4. If silhouette monotonically increasing in upper half → artifact, trust elbow
    5. If silhouette peak is near elbow (within ±2) → use silhouette (refines elbow)
    6. Otherwise → pick best silhouette in the neighborhood of the elbow
    """
    k_range = list(k_range)
    sil_scores = np.array(sil_scores)
    
    elbow_k = find_elbow(k_range, inertias)
    best_sil_k = k_range[np.argmax(sil_scores)]
    sil_range = sil_scores.max() - sil_scores.min()
    
    # Case 1: Both agree
    if elbow_k == best_sil_k:
        reason = "elbow and silhouette agree"
        return elbow_k, elbow_k, best_sil_k, reason
    
    # Case 2: Silhouette essentially flat
    if sil_range < 0.10:
        reason = f"silhouette flat (range={sil_range:.3f}); trusting elbow"
        return elbow_k, elbow_k, best_sil_k, reason
    
    # Case 3: Silhouette peak at high k → small-cluster artifact
    if best_sil_k > n_samples // 5:
        reason = (f"silhouette peak k={best_sil_k} too high for n={n_samples} "
                  f"(avg cluster < 5); trusting elbow")
        return elbow_k, elbow_k, best_sil_k, reason
    
    # Case 4: Silhouette monotonically increasing in upper half
    upper_half = sil_scores[len(sil_scores)//2:]
    if len(upper_half) >= 3 and all(upper_half[i] <= upper_half[i+1]
                                     for i in range(len(upper_half)-1)):
        reason = "silhouette monotonically increasing at high k; trusting elbow"
        return elbow_k, elbow_k, best_sil_k, reason
    
    # Case 5: Silhouette peak close to elbow
    if abs(best_sil_k - elbow_k) <= 2:
        reason = (f"silhouette peak k={best_sil_k} close to elbow k={elbow_k}; "
                  f"using silhouette")
        return best_sil_k, elbow_k, best_sil_k, reason
    
    # Case 6: Real disagreement — best silhouette near elbow neighborhood
    neighborhood = range(max(elbow_k - 1, k_range[0]),
                         min(elbow_k + 3, k_range[-1] + 1))
    neighborhood_sils = {k: sil_scores[k_range.index(k)]
                         for k in neighborhood if k in k_range}
    best_neighbor_k = max(neighborhood_sils, key=neighborhood_sils.get)
    reason = (f"elbow k={elbow_k}, sil peak k={best_sil_k} (far apart); "
              f"best sil near elbow is k={best_neighbor_k}")
    return best_neighbor_k, elbow_k, best_sil_k, reason

# ---------------------------------------------------------------------------
# Helper: Run Mash sketch + dist for a list of genome files
# ---------------------------------------------------------------------------
def run_mash_distances(genome_files, output_dir, label):
    """Run mash sketch and mash dist, return parsed distance matrix."""
    os.makedirs(output_dir, exist_ok=True)
    sketch_prefix = os.path.join(output_dir, f'{label}_sketch')
    dist_file = os.path.join(output_dir, f'{label}_distances.tsv')
    
    # Check if distance file already exists
    if os.path.exists(dist_file):
        print(f"  Loading pre-computed Mash distances from {dist_file}")
        return parse_mash_distances(dist_file)
    
    # Write genome file list
    file_list_path = os.path.join(output_dir, f'{label}_genome_list.txt')
    with open(file_list_path, 'w') as f:
        for gf in genome_files:
            f.write(gf + '\n')
    
    # Mash sketch
    print(f"  Running mash sketch (k={KMER_SIZE}, s={SKETCH_SIZE}) on {len(genome_files)} genomes...")
    cmd_sketch = [
        MASH_BIN, 'sketch', '-k', str(KMER_SIZE), '-s', str(SKETCH_SIZE),
        '-o', sketch_prefix, '-l', file_list_path
    ]
    result = subprocess.run(cmd_sketch, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR in mash sketch: {result.stderr}")
        return None
    
    # Mash dist
    print(f"  Running mash dist (all-vs-all)...")
    sketch_file = sketch_prefix + '.msh'
    cmd_dist = [MASH_BIN, 'dist', sketch_file, sketch_file]
    with open(dist_file, 'w') as f:
        result = subprocess.run(cmd_dist, stdout=f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"  ERROR in mash dist: {result.stderr}")
        return None
    
    print(f"  Mash distances saved to {dist_file}")
    return parse_mash_distances(dist_file)

def parse_mash_distances(dist_file):
    """Parse mash dist tabular output into a square distance matrix DataFrame."""
    acc_re = re.compile(r'(GC[AF]_\d+\.\d+)')
    
    rows = []
    with open(dist_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                ref_acc = acc_re.search(parts[0])
                qry_acc = acc_re.search(parts[1])
                if ref_acc and qry_acc:
                    rows.append({
                        'ref': ref_acc.group(1),
                        'query': qry_acc.group(1),
                        'distance': float(parts[2])
                    })
    
    df = pd.DataFrame(rows)
    dist_matrix = df.pivot(index='ref', columns='query', values='distance')
    dist_matrix = dist_matrix.fillna(0)
    
    # Ensure symmetric
    all_accs = sorted(set(dist_matrix.index) | set(dist_matrix.columns))
    dist_matrix = dist_matrix.reindex(index=all_accs, columns=all_accs, fill_value=0)
    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    # numpy 2.x reindex can return a read-only view; force a writable copy
    arr = np.array(dist_matrix.values, copy=True)
    np.fill_diagonal(arr, 0)
    dist_matrix = pd.DataFrame(arr, index=dist_matrix.index, columns=dist_matrix.columns)

    return dist_matrix

# ---------------------------------------------------------------------------
# Run Mash for each species
# ---------------------------------------------------------------------------
# mash_results = {}

# for species in SPECIES_LIST_MASH:
#     print(f"\n{'='*60}")
#     print(f"  {SPECIES_DISPLAY_MASH[species]} - Mash Distance Computation")
#     print(f"{'='*60}")
#     
#     genome_dir = f'/datadrive/Species/{GENUS}/{species}/filtered_genome'
    # Build ANI exclusion patterns for this species
#     _exc_bases = set()
#     if species in ANI_EXCLUDED:
#         for exc_acc in ANI_EXCLUDED[species]:
#             _exc_bases.add(re.match(r'(GC[AF]_\d+)', exc_acc).group(1))
#     genome_files = sorted([
#         str(p) for p in Path(genome_dir).glob('*.fna')
#         if '_renamed' not in p.name
#         and not any(exc in p.name for exc in _exc_bases)
#     ])
#     print(f"  Found {len(genome_files)} genome files")
#     
#     sp_output_dir = os.path.join(MASH_OUTPUT_BASE, species)
#     dist_matrix = run_mash_distances(genome_files, sp_output_dir, species)
#     
#     if dist_matrix is not None:
#         mash_results[species] = dist_matrix
#         print(f"  Distance matrix shape: {dist_matrix.shape}")
#     else:
#         print(f"  FAILED - skipping {species}")

# print(f"\n{'='*60}")
# print("Mash distance computation complete!")
# for sp, dm in mash_results.items():
#     print(f"  {SPECIES_DISPLAY_MASH[sp]}: {dm.shape[0]} genomes")
# print(f"{'='*60}")


# --- From NB0 Cell 72 ---
# =============================================================================
# CELL: Issue 7 — Heap's Law by Phenotype Subset (Fumigatus Clinical vs Environmental)
# Tests whether pangenome closure is a sampling artefact of single ecological context
# =============================================================================

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# print("=" * 70)
# print("HEAP'S LAW BY PHENOTYPE SUBSET — A. fumigatus")
# print("=" * 70)
# print("Question: Is pangenome closure robust within a single ecological stratum?")

# pav = pd.read_csv('/datadrive/Analysis/NB1_Results/fumigatus/fumigatus_pav.tsv', sep='\t', index_col=0)
# burden = pd.read_csv('/datadrive/Analysis/NB4_Results/rare_gene_burden_all_genomes.csv')
# fum_b = burden[burden['Species'] == 'A. fumigatus']

# clinical_genomes = set(fum_b[fum_b['Phenotype'] == 'Human-pathogenic']['Genome'])
# env_genomes = set(fum_b[fum_b['Phenotype'] == 'Environmental']['Genome'])

# pav_c = pav[[c for c in pav.columns if c in clinical_genomes]]
# pav_e = pav[[c for c in pav.columns if c in env_genomes]]

def heaps_law(n, kappa, gamma):
    return kappa * (n ** gamma)

def fit_heaps(pav_sub, n_reps=100, seed=42):
    rng = np.random.default_rng(seed)
    n_g = pav_sub.shape[1]
    steps = list(range(1, n_g + 1))
    means = []
    for n in steps:
        sizes = [(pav_sub.iloc[:, rng.choice(n_g, size=n, replace=False)].sum(axis=1) > 0).sum()
                 for _ in range(n_reps)]
        means.append(np.mean(sizes))
    popt, _ = curve_fit(heaps_law, steps, means, p0=[5000, 0.2], maxfev=5000, bounds=([0, 0], [np.inf, 1]))
    return popt[1], means

# print(f"\nClinical only (n={pav_c.shape[1]})...")
# g_c, _ = fit_heaps(pav_c)
# print(f"  gamma = {g_c:.4f}")

# print(f"Environmental only (n={pav_e.shape[1]})...")
# g_e, _ = fit_heaps(pav_e)
# print(f"  gamma = {g_e:.4f}")

# print(f"All strains (n={pav.shape[1]})...")
# g_all, _ = fit_heaps(pav)
# print(f"  gamma = {g_all:.4f}")

# print("\nCONCLUSION:")
# print(f"  Clinical-only gamma ({g_c:.4f}) and environmental-only gamma ({g_e:.4f})")
# print(f"  are both << 0.5 (open threshold), confirming pangenome closure")
# print(f"  within each ecological stratum independently.")

# results = pd.DataFrame({'subset': ['clinical_only', 'environmental_only', 'all_strains'],
#                          'n_genomes': [pav_c.shape[1], pav_e.shape[1], pav.shape[1]],
#                          'heaps_gamma': [g_c, g_e, g_all]})
# results.to_csv('/datadrive/Analysis/NB0_Results/fumigatus_heaps_phenotype_subset.tsv', sep='\t', index=False)
# print("\nSaved to NB0_Results/fumigatus_heaps_phenotype_subset.tsv")



# --- From NB4 Cell 1 ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.signal import find_peaks
from statsmodels.stats.multitest import multipletests
from collections import Counter
import warnings, os, sys, re

warnings.filterwarnings('ignore')
sns.set_style('whitegrid')
plt.rcParams.update({'font.size': 11, 'figure.dpi': 120})

# sys.path.insert(0, '/datadrive/Analysis')

def get_accession(col):
    """Extract accession (e.g. 'GCA_001599455.1') from OrthoFinder column names
    like 'GCA_001599455.1_JCM_10253_assembly_v001.proteins'."""
    m = re.match(r'(GC[AF]_\d+\.\d+)', col)
    return m.group(1) if m else col

# ---------------------------------------------------------------------------
# Shared configuration (mirrors NB1 / NB0 / NB3 conventions)
# ---------------------------------------------------------------------------
# GENUS = 'Aspergillus'
# SPECIES_LIST = ['fumigatus', 'flavus', 'niger', 'oryzae']
# BASE_PATH = f"/datadrive/Species/{GENUS}"

# Combined OrthoFinder output (same as NB3)
# ORTHOFINDER_DIR = os.path.join(BASE_PATH, 'all_combined', 'orthofinder_output', 'Results_Dec14')
# GENECOUNT_PATH = os.path.join(ORTHOFINDER_DIR, 'Orthogroups', 'Orthogroups.GeneCount.tsv')
# ORTHOGROUPS_PATH = os.path.join(ORTHOFINDER_DIR, 'Orthogroups', 'Orthogroups.tsv')

# Upstream theme results
# NB0_RESULTS = '/datadrive/Analysis/NB0_Results'
# NB3_RESULTS = '/datadrive/Analysis/NB3_Results'
# NB1_RESULTS = '/datadrive/Analysis/NB1_Results'

# NB4 output directory
# RESULTS_DIR = '/datadrive/Analysis/NB4_Results'
# os.makedirs(RESULTS_DIR, exist_ok=True)

# Per-species config — paths built dynamically from BASE_PATH / NB1_RESULTS
# SPECIES_CONFIG = {
#     sp: {
#         'label': f'A. {sp}',
#         'annot_path': os.path.join(BASE_PATH, sp, f'{sp}_ortho_annot_long.tsv'),
#         'pav_path': os.path.join(NB1_RESULTS, sp, f'{sp}_pav.tsv'),
#         'og_consensus_path': os.path.join(NB1_RESULTS, sp, f'{sp}_og_consensus.tsv'),
#         'species_filter': f'A. {sp}',
#     }
#     for sp in SPECIES_LIST
# }

# Load phenotype labels (produced by NB0)
# pheno_df = pd.read_csv(os.path.join(NB0_RESULTS, 'phenotype_classified_for_gwas.csv'))
# pheno_map = dict(zip(pheno_df['Assembly Accession'], pheno_df['Phenotype']))
# species_map = dict(zip(pheno_df['Assembly Accession'], pheno_df['Species']))
# print(f"Loaded phenotype labels for {len(pheno_map)} genomes")
# print(pheno_df.groupby('Species')['Phenotype'].value_counts().to_string())

# Verify all paths exist
# print("\n--- Path verification ---")
# for sp, cfg in SPECIES_CONFIG.items():
#     print(f"  {cfg['label']} annotation:    {os.path.exists(cfg['annot_path'])}")
#     print(f"  {cfg['label']} PAV matrix:    {os.path.exists(cfg['pav_path'])}")
#     print(f"  {cfg['label']} OG consensus:  {os.path.exists(cfg['og_consensus_path'])}")
# print(f"  OrthoFinder GeneCount:  {os.path.exists(GENECOUNT_PATH)}")
# print(f"  OrthoFinder Orthogroups: {os.path.exists(ORTHOGROUPS_PATH)}")
# print(f"  NB0 phenotype data:  {os.path.exists(os.path.join(NB0_RESULTS, 'phenotype_classified_for_gwas.csv'))}")


# --- From NB4 Cell 3 ---
def load_species_data(species_key):
    """Load pre-built PAV + OG consensus (with Pangenome_Class) from NB1_Results,
    plus the long annotation table for protein-level queries."""
    cfg = SPECIES_CONFIG[species_key]
    print(f"\n{'='*60}")
    print(f"Loading {cfg['label']}...")

    # Load long annotation table (for protein-level lookups)
    long_annot = pd.read_csv(cfg['annot_path'], sep='\t', low_memory=False)
    print(f"  Long annotation: {len(long_annot):,} rows, {long_annot['Orthogroup'].nunique():,} OGs, "
          f"{long_annot['Assembly Accession'].nunique()} genomes")

    # Load pre-built PAV matrix from NB1_Results, drop any ANI-excluded
    # genomes that may still be present in the on-disk file.
    pav = pd.read_csv(cfg['pav_path'], sep='\t', index_col=0)
    sp_short = cfg.get('species') or cfg['label'].replace('A. ', '').strip()
    try:
        from funpan_utils import filter_ani_excluded as _fae
        pav = _fae(pav, sp_short, axis='cols')
    except Exception:
        pass
    print(f"  PAV matrix (from NB1_Results, ANI-filtered): "
          f"{pav.shape[0]:,} OGs x {pav.shape[1]} genomes")

    # Load pre-classified OG consensus table from NB1_Results
    og_consensus = pd.read_csv(cfg['og_consensus_path'], sep='\t')
    n_core = (og_consensus['Pangenome_Class'] == 'Core').sum()
    n_acc = (og_consensus['Pangenome_Class'] == 'Accessory').sum()
    n_rare = (og_consensus['Pangenome_Class'] == 'Rare').sum()
    n_total = len(og_consensus)
    print(f"  OG consensus (from NB1_Results): {n_total} OGs")
    print(f"    Core: {n_core} ({n_core/n_total:.1%}), Accessory: {n_acc} ({n_acc/n_total:.1%}), "
          f"Rare: {n_rare} ({n_rare/n_total:.1%})")

    # Recover thresholds from the PAV + classification
    # Core threshold = min presence count among Core OGs
    # Rare threshold = max presence count among Rare OGs
    gene_freq = pav.sum(axis=1)
    core_ogs = set(og_consensus.loc[og_consensus['Pangenome_Class'] == 'Core', 'Orthogroup'])
    rare_ogs = set(og_consensus.loc[og_consensus['Pangenome_Class'] == 'Rare', 'Orthogroup'])
    core_freqs = gene_freq[gene_freq.index.isin(core_ogs)]
    rare_freqs = gene_freq[gene_freq.index.isin(rare_ogs)]
    core_n = int(core_freqs.min()) if len(core_freqs) > 0 else pav.shape[1]
    rare_n = int(rare_freqs.max()) if len(rare_freqs) > 0 else 1
    n_strains = pav.shape[1]
    print(f"  Thresholds (recovered): core >= {core_n}/{n_strains}, rare <= {rare_n}/{n_strains}")

    return {
        'long_annot': long_annot,
        'pav': pav,
        'og_consensus': og_consensus,
        'core_n': core_n,
        'rare_n': rare_n,
        'n_strains': n_strains,
    }

# Load all species
# species_data = {}
# for sp in SPECIES_LIST:
#     species_data[sp] = load_species_data(sp)


# --- From NB4 Cell 10 ---
def annotation_rate_by_class(og_consensus, annotation_col, classes=('Core', 'Accessory', 'Rare')):
    """Compute fraction of OGs with a given annotation, per pangenome class."""
    results = {}
    for cls in classes:
        sub = og_consensus[og_consensus['Pangenome_Class'] == cls]
        if len(sub) == 0:
            results[cls] = {'n': 0, 'annotated': 0, 'rate': 0}
            continue
        has_annot = sub[annotation_col].apply(
            lambda x: bool(x) and str(x).strip() not in ('', '-', 'nan')
        ).sum()
        results[cls] = {'n': len(sub), 'annotated': int(has_annot), 'rate': has_annot / len(sub)}
    return results


def cog_distribution_by_class(og_consensus, classes=('Core', 'Accessory', 'Rare')):
    """Get COG category distribution per pangenome class."""
    COG_NAMES = {
        'C': 'Energy', 'D': 'Cell cycle', 'E': 'Amino acid', 'F': 'Nucleotide',
        'G': 'Carbohydrate', 'H': 'Coenzyme', 'I': 'Lipid', 'J': 'Translation',
        'K': 'Transcription', 'L': 'Replication', 'M': 'Cell wall', 'N': 'Cell motility',
        'O': 'PTM/chaperone', 'P': 'Inorganic ion', 'Q': 'Secondary metabolites',
        'S': 'Unknown function', 'T': 'Signal transduction', 'U': 'Trafficking/secretion',
        'V': 'Defense', 'W': 'Extracellular', 'X': 'Mobilome', 'Z': 'Cytoskeleton',
    }
    results = {}
    for cls in classes:
        sub = og_consensus[og_consensus['Pangenome_Class'] == cls]
        cog_counts = Counter()
        for cog_str in sub['COG_category'].dropna():
            for c in str(cog_str):
                if c in COG_NAMES:
                    cog_counts[c] += 1
        total = sum(cog_counts.values())
        results[cls] = {k: v / total if total > 0 else 0 for k, v in cog_counts.items()}
        results[cls]['_total'] = total
    return results, COG_NAMES


# Compute annotation rates across species
# print("=== Annotation rates by pangenome class ===\n")
# annot_cols = ['PFAMs', 'CAZy', 'GOs', 'Description']
# all_rates = []

# for sp in ['fumigatus', 'flavus', 'niger', 'oryzae']:
#     ogc = species_data[sp]['og_consensus']
#     label = SPECIES_CONFIG[sp]['label']
#     print(f"\n{label}:")
#     for col in annot_cols:
#         if col in ogc.columns:
#             rates = annotation_rate_by_class(ogc, col)
#             for cls, info in rates.items():
#                 all_rates.append({
#                     'Species': label, 'Annotation': col, 
#                     'Class': cls, 'Rate': info['rate'], 
#                     'n': info['n'], 'annotated': info['annotated']
#                 })
#                 print(f"  {col} - {cls}: {info['annotated']}/{info['n']} ({info['rate']:.1%})")

# rates_df = pd.DataFrame(all_rates)
# rates_df.to_csv(f'{RESULTS_DIR}/annotation_rates_by_class.csv', index=False)


# --- From NB4 Cell 25 ---
# ── C1 sensitivity: fixed absolute rare-genome thresholds ─────────────────────────────────────────────
# Current definition: rare = present in ≤5% of conspecific genomes (species-specific threshold)
#   fumigatus n=89  → ≤4 genomes
#   flavus    n=71  → ≤3 genomes
#   niger     n=30  → ≤1 genome
#   oryzae    n=33  → ≤1 genome
#
# Sensitivity check: reassign rare/accessory at fixed thresholds n≤1, n≤2, n≤3
# Key questions:
#   (a) Does A. oryzae rare count remain stable across thresholds?
#   (b) Does A. fumigatus rare-burden phenotype association survive at stricter thresholds?

import pandas as pd
import numpy as np
from scipy import stats
import os, re as _re

# RESULTS_DIR_SENS = '/datadrive/Analysis/NB4_Results'
# NB1_RESULTS = '/datadrive/Analysis/NB1_Results'
# NB0_RESULTS = '/datadrive/Analysis/NB0_Results'

# SPECIES_CONFIG_LOCAL = {
#     'fumigatus': {'label': 'A. fumigatus', 'n': 89, 'current_rare_n': 4, 'pheno_group': 'pathogenic'},
#     'flavus':    {'label': 'A. flavus',    'n': 71, 'current_rare_n': 3, 'pheno_group': 'pathogenic'},
#     'niger':     {'label': 'A. niger',     'n': 30, 'current_rare_n': 1, 'pheno_group': 'industrial'},
#     'oryzae':    {'label': 'A. oryzae',    'n': 33, 'current_rare_n': 1, 'pheno_group': 'industrial'},
# }

# Load phenotype map
# Phenotype column values: 'Human-pathogenic', 'Plant-pathogenic', 'Animal-pathogenic',
#                          'Environmental', 'Industrial-trait', 'Lab'
# pheno_df_local = pd.read_csv(os.path.join(NB0_RESULTS, 'phenotype_classified_for_gwas.csv'))
# pheno_map_local = dict(zip(pheno_df_local['Assembly Accession'], pheno_df_local['Phenotype']))

def get_acc_local(col):
    m = _re.match(r'(GC[AF]_\d+\.\d+)', col)
    return m.group(1) if m else col

# thresholds = [1, 2, 3]
# rows = []
# burden_rows = []

# for sp, cfg in SPECIES_CONFIG_LOCAL.items():
#     og_path = os.path.join(NB1_RESULTS, sp, f'{sp}_og_consensus.tsv')
#     pav_path = os.path.join(NB1_RESULTS, sp, f'{sp}_pav.tsv')
#     ogc = pd.read_csv(og_path, sep='\t')
#     pav = pd.read_csv(pav_path, sep='\t', index_col=0)
# 
#     n_total = len(ogc)
#     n_genomes = pav.shape[1]
# 
#     print(f"\n{'='*60}")
#     print(f"{cfg['label']}  (n_genomes={n_genomes}, n_OGs={n_total})")
#     print(f"  Current rare threshold: \u2264{cfg['current_rare_n']} genomes (\u22645% of {n_genomes})")
# 
#     current_rare_count = (ogc['Pangenome_Class'] == 'Rare').sum()
#     rows.append({'Species': cfg['label'], 'Threshold_type': 'Current (≤5%)',
#                  'Threshold_n': cfg['current_rare_n'], 'n_Rare': current_rare_count,
#                  'Pct_Rare': current_rare_count/n_total*100})
# 
#     for thresh in thresholds:
#         rare_ogs_fixed = set(ogc.loc[ogc['n_genomes'] <= thresh, 'Orthogroup'])
#         n_rare = len(rare_ogs_fixed)
#         rows.append({'Species': cfg['label'], 'Threshold_type': f'Fixed n\u2264{thresh}',
#                      'Threshold_n': thresh, 'n_Rare': n_rare,
#                      'Pct_Rare': n_rare/n_total*100})
#         print(f"  Fixed n\u2264{thresh}: {n_rare} rare OGs ({n_rare/n_total:.1%})")
# 
#     if sp in ('fumigatus', 'flavus'):
#         pheno_map_sp = {col: pheno_map_local.get(get_acc_local(col), 'Unknown') for col in pav.columns}
# 
#         for thresh in thresholds:
#             rare_fixed = set(ogc.loc[ogc['n_genomes'] <= thresh, 'Orthogroup'])
#             if len(rare_fixed) < 2:
#                 burden_rows.append({'Species': cfg['label'], 'Threshold': f'n\u2264{thresh}',
#                                     'n_rare_OGs': len(rare_fixed), 'Group1': 'NA', 'Group2': 'NA',
#                                     'n_g1': 0, 'n_g2': 0,
#                                     'mean_g1': np.nan, 'mean_g2': np.nan, 'U_stat': np.nan, 'p_value': np.nan})
#                 continue
# 
#             rare_pav = pav.loc[pav.index.isin(rare_fixed)]
#             burden = rare_pav.sum(axis=0)
# 
            # pathogenic = any *-pathogenic label; environmental = Environmental
#             g1_burden = burden[[c for c in pav.columns if 'pathogenic' in pheno_map_sp.get(c, '').lower()]].values
#             g2_burden = burden[[c for c in pav.columns if pheno_map_sp.get(c, '') == 'Environmental']].values
# 
#             if len(g1_burden) >= 2 and len(g2_burden) >= 2:
#                 u_stat, p_val = stats.mannwhitneyu(g1_burden, g2_burden, alternative='two-sided')
#                 direction = 'higher in pathogenic' if np.mean(g1_burden) > np.mean(g2_burden) else 'lower in pathogenic'
#                 print(f"  Burden test n\u2264{thresh}: pathogenic mean={np.mean(g1_burden):.1f} (n={len(g1_burden)}) "
#                       f"vs env mean={np.mean(g2_burden):.1f} (n={len(g2_burden)}), "
#                       f"U={u_stat:.0f}, P={p_val:.2e} [{direction}]")
#                 burden_rows.append({'Species': cfg['label'], 'Threshold': f'n\u2264{thresh}',
#                                     'n_rare_OGs': len(rare_fixed),
#                                     'Group1': 'Pathogenic', 'Group2': 'Environmental',
#                                     'n_g1': len(g1_burden), 'n_g2': len(g2_burden),
#                                     'mean_g1': np.mean(g1_burden), 'mean_g2': np.mean(g2_burden),
#                                     'U_stat': u_stat, 'p_value': p_val})
#             else:
#                 print(f"  Burden test n\u2264{thresh}: insufficient phenotype-labelled genomes "
#                       f"(pathogenic n={len(g1_burden)}, env n={len(g2_burden)})")
#                 burden_rows.append({'Species': cfg['label'], 'Threshold': f'n\u2264{thresh}',
#                                     'n_rare_OGs': len(rare_fixed),
#                                     'Group1': 'Pathogenic', 'Group2': 'Environmental',
#                                     'n_g1': len(g1_burden), 'n_g2': len(g2_burden),
#                                     'mean_g1': np.nan, 'mean_g2': np.nan,
#                                     'U_stat': np.nan, 'p_value': np.nan})

# print("\n\n=== SENSITIVITY SUMMARY TABLE ===")
# summary_df = pd.DataFrame(rows)
# summary_pivot = summary_df.pivot_table(index='Species', columns='Threshold_type', values='n_Rare', aggfunc='first')
# print(summary_pivot.to_string())

# print("\n=== PHENOTYPE BURDEN ASSOCIATION STABILITY ===")
# burden_df = pd.DataFrame(burden_rows)
# if len(burden_df) > 0:
#     print(burden_df[['Species','Threshold','n_rare_OGs','mean_g1','mean_g2','p_value']].to_string(index=False))

# summary_df.to_csv(os.path.join(RESULTS_DIR_SENS, 'rare_threshold_sensitivity.csv'), index=False)
# burden_df.to_csv(os.path.join(RESULTS_DIR_SENS, 'rare_threshold_burden_sensitivity.csv'), index=False)
# print("\nSaved: rare_threshold_sensitivity.csv, rare_threshold_burden_sensitivity.csv")



# --- From NB4 Cell 27 ---
# =========================================================================
# 9.1  Gene length distribution by pangenome class
# =========================================================================
# Parse protein lengths from OrthoFinder Orthogroup_Sequences/*.fa
# Each OG file contains all member proteins → compute median length per OG

from pathlib import Path
from Bio import SeqIO
import warnings, os, sys, re
warnings.filterwarnings('ignore')

# SPECIES_DEEP = ['fumigatus', 'flavus', 'niger', 'oryzae']

# RESULTS_DIR = '/datadrive/Analysis/NB4_Results'
# os.makedirs(RESULTS_DIR, exist_ok=True)

def og_protein_lengths(og_seq_dir, og_list):
    """Return dict {OG: [len1, len2, ...]} for requested OGs."""
    lengths = {}
    for og in og_list:
        fa = Path(og_seq_dir) / f"{og}.fa"
        if fa.exists():
            lengths[og] = [len(rec.seq) for rec in SeqIO.parse(fa, "fasta")]
    return lengths

# print("Computing protein lengths per orthogroup ...")
# length_data = {}   # species → DataFrame(Orthogroup, median_len, Pangenome_Class)

# for sp in SPECIES_DEEP:
#     d = species_data[sp]
#     ogc = d['og_consensus']
#     og_dir = f"/datadrive/Species/Aspergillus/{sp}/orthofinder_output/"
#     res_dir = sorted(Path(og_dir).glob("Results_*/Orthogroup_Sequences"))[0]
# 
#     all_ogs = ogc['Orthogroup'].tolist()
#     lens = og_protein_lengths(str(res_dir), all_ogs)
# 
#     rows = []
#     for _, r in ogc.iterrows():
#         og = r['Orthogroup']
#         if og in lens and lens[og]:
#             rows.append({
#                 'Orthogroup': og,
#                 'median_len': np.median(lens[og]),
#                 'mean_len': np.mean(lens[og]),
#                 'min_len': min(lens[og]),
#                 'n_proteins': len(lens[og]),
#                 'Pangenome_Class': r['Pangenome_Class']
#             })
#     length_data[sp] = pd.DataFrame(rows)
#     print(f"  {sp}: {len(rows)} OGs with length data")

# --- Visualization ---
# fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=False)
# axes = axes.flatten()
# class_order = ['Core', 'Accessory', 'Rare']
# class_colors = {'Core': 'forestgreen', 'Accessory': 'steelblue', 'Rare': 'coral'}

# for idx, sp in enumerate(SPECIES_DEEP):
#     ax = axes[idx]
#     df = length_data[sp]
#     data_to_plot = [df.loc[df['Pangenome_Class'] == c, 'median_len'].values for c in class_order]
# 
    # Clip outliers for visualization: cap at 99th percentile across all classes
#     all_vals = np.concatenate([d for d in data_to_plot if len(d) > 0])
#     upper = np.percentile(all_vals, 99)
# 
    # Skip empty classes for violin (oryzae may have very few rare)
#     non_empty = [len(d) > 0 for d in data_to_plot]
#     if all(non_empty):
#         vp = ax.violinplot(data_to_plot, positions=[1, 2, 3], showmedians=True, showextrema=False)
#         for i, body in enumerate(vp['bodies']):
#             body.set_facecolor(list(class_colors.values())[i])
#             body.set_alpha(0.6)
#         vp['cmedians'].set_color('black')
#     else:
#         bp = ax.boxplot(data_to_plot, positions=[1, 2, 3], patch_artist=True, widths=0.5)
#         for i, patch in enumerate(bp['boxes']):
#             patch.set_facecolor(list(class_colors.values())[i])
#             patch.set_alpha(0.6)
# 
    # Add median labels
#     for i, c in enumerate(class_order):
#         if len(data_to_plot[i]) > 0:
#             med = np.median(data_to_plot[i])
#             n = len(data_to_plot[i])
#             ax.text(i + 1, med + 20, f"med={med:.0f}\nn={n}", ha='center', va='bottom', fontsize=8)
# 
#     ax.set_xticks([1, 2, 3])
#     ax.set_xticklabels(class_order, fontsize=10)
#     ax.set_ylabel('Median protein length (aa)' if idx % 2 == 0 else '', fontsize=11)
#     ax.set_title(f'A. {sp}', fontsize=12, fontweight='bold')
#     ax.set_ylim(0, upper * 1.15)
#     ax.grid(axis='y', linestyle=':', alpha=0.4)

# plt.suptitle('Protein Length Distribution by Pangenome Class', fontsize=14, fontweight='bold', y=1.02)
# plt.tight_layout()
# plt.savefig(f'{RESULTS_DIR}/rare_gene_length_distribution.png', dpi=150, bbox_inches='tight')

# Statistical test: are rare genes shorter?
from scipy.stats import mannwhitneyu
# print("\nMann-Whitney U test (Rare vs Core median protein length):")
# for sp in SPECIES_DEEP:
#     df = length_data[sp]
#     core_lens = df.loc[df['Pangenome_Class'] == 'Core', 'median_len']
#     rare_lens = df.loc[df['Pangenome_Class'] == 'Rare', 'median_len']
#     if len(rare_lens) >= 2 and len(core_lens) >= 2:
#         stat, pval = mannwhitneyu(rare_lens, core_lens, alternative='less')
#         print(f"  {sp}: Rare median={rare_lens.median():.0f}aa, Core median={core_lens.median():.0f}aa, "
#               f"p={pval:.2e} ({'*' if pval < 0.05 else 'ns'})")
#     else:
#         print(f"  {sp}: Rare n={len(rare_lens)} (too few for test), "
#               f"Core median={core_lens.median():.0f}aa")


# --- From NB4 Cell 29 ---
# =========================================================================
# 9.3  DIAMOND BLASTp: Rare vs Core/Accessory to detect fragments
# =========================================================================
import subprocess, tempfile, shutil

# CONDA_PREFIX = '/home/user/anaconda3/envs/pipeline_test2/bin'
# DIAMOND = f'{CONDA_PREFIX}/diamond'

def get_representative_seqs(og_seq_dir, og_list):
    """For each OG, return the longest protein sequence as representative."""
    reps = {}
    for og in og_list:
        fa = Path(og_seq_dir) / f"{og}.fa"
        if fa.exists():
            best = None
            for rec in SeqIO.parse(fa, "fasta"):
                if best is None or len(rec.seq) > len(best.seq):
                    best = rec
            if best:
                best.id = og
                best.description = og
                reps[og] = best
    return reps

# print("Running DIAMOND BLASTp: Rare OGs vs Core+Accessory OGs ...")
# print("=" * 80)

# diamond_results = {}

# for sp in SPECIES_DEEP:
#     d = species_data[sp]
#     ogc = d['og_consensus']
#     og_dir = f"/datadrive/Species/Aspergillus/{sp}/orthofinder_output/"
#     res_dir = str(sorted(Path(og_dir).glob("Results_*/Orthogroup_Sequences"))[0])
# 
#     rare_ogs = ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'].tolist()
#     nonrare_ogs = ogc.loc[ogc['Pangenome_Class'].isin(['Core', 'Accessory']), 'Orthogroup'].tolist()
# 
#     if len(rare_ogs) == 0:
#         print(f"\n{sp}: 0 rare OGs — skipping DIAMOND")
#         diamond_results[sp] = pd.DataFrame(columns=['Orthogroup', 'blast_class', 'pident', 'qcovhsp', 'scovhsp', 'best_target'])
#         continue
# 
#     print(f"\n{sp}: {len(rare_ogs)} rare queries vs {len(nonrare_ogs)} core+accessory targets")
# 
#     rare_reps = get_representative_seqs(res_dir, rare_ogs)
#     nonrare_reps = get_representative_seqs(res_dir, nonrare_ogs)
#     print(f"  Representatives: {len(rare_reps)} rare, {len(nonrare_reps)} non-rare")
# 
#     tmpdir = tempfile.mkdtemp(prefix=f"diamond_{sp}_")
#     query_fa = f"{tmpdir}/rare_query.fa"
#     db_fa = f"{tmpdir}/nonrare_db.fa"
#     db_path = f"{tmpdir}/nonrare_db"
#     out_tsv = f"{tmpdir}/diamond_out.tsv"
# 
#     SeqIO.write(rare_reps.values(), query_fa, "fasta")
#     SeqIO.write(nonrare_reps.values(), db_fa, "fasta")
# 
#     subprocess.run([DIAMOND, 'makedb', '--in', db_fa, '-d', db_path],
#                    capture_output=True, check=True)
#     subprocess.run([
#         DIAMOND, 'blastp',
#         '-q', query_fa, '-d', db_path, '-o', out_tsv,
#         '--sensitive', '--max-target-seqs', '5',
#         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
#         'qcovhsp', 'scovhsp', 'evalue', 'bitscore',
#         '--threads', '8'
#     ], capture_output=True, check=True)
# 
#     cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
#             'qcovhsp', 'scovhsp', 'evalue', 'bitscore']
#     if os.path.getsize(out_tsv) > 0:
#         hits = pd.read_csv(out_tsv, sep='\t', names=cols)
#         best_hits = hits.sort_values('bitscore', ascending=False).drop_duplicates('qseqid', keep='first')
#     else:
#         best_hits = pd.DataFrame(columns=cols)
# 
#     classifications = []
#     for og in rare_ogs:
#         if og in rare_reps:
#             row = best_hits[best_hits['qseqid'] == og]
#             if row.empty:
#                 classifications.append({'Orthogroup': og, 'blast_class': 'No hit', 'pident': 0,
#                                         'qcovhsp': 0, 'scovhsp': 0, 'best_target': ''})
#             else:
#                 r = row.iloc[0]
#                 qcov, scov, pident = r['qcovhsp'], r['scovhsp'], r['pident']
#                 if qcov >= 70 and scov < 50 and pident >= 30:
#                     bclass = 'Fragment of longer gene'
#                 elif qcov >= 70 and scov >= 70 and pident >= 30:
#                     bclass = 'Full-length homolog'
#                 else:
#                     bclass = 'No significant similarity'
#                 classifications.append({
#                     'Orthogroup': og, 'blast_class': bclass, 'pident': pident,
#                     'qcovhsp': qcov, 'scovhsp': scov, 'best_target': r['sseqid']
#                 })
# 
#     diamond_results[sp] = pd.DataFrame(classifications)
#     shutil.rmtree(tmpdir)
# 
#     counts = diamond_results[sp]['blast_class'].value_counts()
#     print(f"  Results:")
#     for cls, cnt in counts.items():
#         print(f"    {cls}: {cnt} ({100*cnt/len(classifications):.1f}%)")

# --- Visualization ---
# fig, axes = plt.subplots(2, 2, figsize=(14, 10))
# axes = axes.flatten()
# blast_class_colors = {
#     'Full-length homolog': 'steelblue',
#     'Fragment of longer gene': 'coral',
#     'No significant similarity': 'lightgrey',
#     'No hit': "#8D8585"
# }
# blast_class_order = ['Full-length homolog', 'Fragment of longer gene',
#                      'No significant similarity', 'No hit']

# for idx, sp in enumerate(SPECIES_DEEP):
#     ax = axes[idx]
#     df = diamond_results[sp]
#     if df.empty:
#         ax.text(0.5, 0.5, f'A. {sp}\nNo rare OGs', ha='center', va='center',
#                 fontsize=12, transform=ax.transAxes)
#         ax.set_title(f'A. {sp}', fontsize=11, fontweight='bold')
#         ax.axis('off')
#         continue
#     counts = df['blast_class'].value_counts()
#     vals = [counts.get(c, 0) for c in blast_class_order]
#     colors = [blast_class_colors[c] for c in blast_class_order]
#     wedges, _ = ax.pie(vals, colors=colors, startangle=90,
#                        wedgeprops=dict(linewidth=0.5, edgecolor='white'))
#     total = sum(vals)
#     for i, (v, w) in enumerate(zip(vals, wedges)):
#         if v > 0:
#             ang = (w.theta1 + w.theta2) / 2
#             x = 0.65 * np.cos(np.radians(ang))
#             y = 0.65 * np.sin(np.radians(ang))
#             ax.text(x, y, f"{100*v/total:.0f}%", ha='center', va='center', fontsize=9, fontweight='bold')
#     ax.set_title(f'A. {sp}\n(n={total} rare OGs)', fontsize=11, fontweight='bold')

# Shared legend
from matplotlib.patches import Patch
# legend_patches = [Patch(facecolor=blast_class_colors[c], label=c) for c in blast_class_order]
# fig.legend(handles=legend_patches, fontsize=9, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.02))
# plt.suptitle('DIAMOND BLASTp Classification of Rare Orthogroups', fontsize=14, fontweight='bold', y=1.02)
# plt.tight_layout()
# plt.savefig(f'{RESULTS_DIR}/rare_diamond_blast_classification.png', dpi=150, bbox_inches='tight')


# --- From NB4 Cell 34 ---
# =========================================================================
# 9.5c  NCBI BLASTp: Taxonomic origin of xenolog candidates
# =========================================================================
# Results are saved incrementally to a TSV — safe to re-run after a crash.
# Already-completed OGs are skipped on re-run.

import requests, time, xml.etree.ElementTree as ET

NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
# MAX_HITS = 5
# RESULTS_FILE = f'{RESULTS_DIR}/rare_xenolog_ncbi_blast_results.tsv'
# MAX_RETRIES = 3

def ncbi_remote_blast(sequence, program='blastp', database='nr', max_hits=5, max_wait=300):
    """Submit a sequence to NCBI BLAST and return parsed top hits."""
    params = {
        'CMD': 'Put', 'PROGRAM': program, 'DATABASE': database,
        'QUERY': str(sequence), 'HITLIST_SIZE': str(max_hits), 'FORMAT_TYPE': 'XML',
    }
    resp = requests.post(NCBI_BLAST_URL, data=params, timeout=60)
    resp.raise_for_status()

    rid = None
    for line in resp.text.split('\n'):
        if 'RID = ' in line:
            rid = line.split('RID = ')[1].strip()
            break
    if not rid:
        return None

    elapsed = 0
    while elapsed < max_wait:
        time.sleep(15)
        elapsed += 15
        try:
            check = requests.get(NCBI_BLAST_URL,
                                 params={'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'},
                                 timeout=120)
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError,
                requests.exceptions.Timeout) as e:
            print(f"      Network error polling RID {rid}: {e} — retrying ...")
            continue

        if 'Status=WAITING' in check.text:
            continue
        if 'Status=FAILED' in check.text:
            return None
        if '<?xml' in check.text:
            try:
                root = ET.fromstring(check.text)
                hits = []
                for hit in root.iter('Hit'):
                    hit_def = hit.find('Hit_def').text if hit.find('Hit_def') is not None else ''
                    hit_acc = hit.find('Hit_accession').text if hit.find('Hit_accession') is not None else ''
                    hsp = hit.find('.//Hsp')
                    if hsp is not None:
                        identity = float(hsp.find('Hsp_identity').text)
                        align_len = float(hsp.find('Hsp_align-len').text)
                        pident = 100 * identity / align_len if align_len > 0 else 0
                        evalue = hsp.find('Hsp_evalue').text
                    else:
                        pident, evalue = 0, 'N/A'
                    organism = ''
                    if '[' in hit_def and ']' in hit_def:
                        organism = hit_def[hit_def.rfind('[') + 1:hit_def.rfind(']')]
                    hits.append({
                        'accession': hit_acc, 'description': hit_def[:120],
                        'organism': organism, 'pident': pident, 'evalue': evalue,
                    })
                    if len(hits) >= max_hits:
                        break
                return hits
            except ET.ParseError:
                return None
    return None


# Collect xenolog candidate sequences
# all_candidates = []
# for sp in SPECIES_DEEP:
#     xc = xenolog_candidates.get(sp)
#     if xc is None or xc.empty:
#         continue
#     og_dir = f"/datadrive/Species/Aspergillus/{sp}/orthofinder_output/"
#     res_dir = str(sorted(Path(og_dir).glob("Results_*/Orthogroup_Sequences"))[0])
#     reps = get_representative_seqs(res_dir, xc['Orthogroup'].tolist())
#     for _, row in xc.iterrows():
#         og = row['Orthogroup']
#         if og in reps:
#             all_candidates.append({
#                 'species': sp, 'Orthogroup': og,
#                 'sequence': str(reps[og].seq),
#                 'GC_deviation': row['GC_deviation'],
#                 'description': row.get('Description', ''),
#             })

# Load existing results (resume support)
# completed_ogs = set()
# if os.path.exists(RESULTS_FILE):
#     existing = pd.read_csv(RESULTS_FILE, sep='\t')
#     completed_ogs = set(existing['Orthogroup'].unique())
#     print(f"Resuming: {len(completed_ogs)} OGs already completed in {RESULTS_FILE}")

# remaining = [c for c in all_candidates if c['Orthogroup'] not in completed_ogs]
# n_total = len(all_candidates)
# n_remaining = len(remaining)

# print(f"Total xenolog candidates: {n_total}")
# print(f"Already done: {n_total - n_remaining}, remaining: {n_remaining}")
# if n_remaining > 0:
#     print(f"Estimated time: ~{n_remaining * 45 // 60} - {n_remaining * 75 // 60} minutes")
# print("=" * 80)

# Write header if file doesn't exist
# if not os.path.exists(RESULTS_FILE):
#     with open(RESULTS_FILE, 'w') as f:
#         f.write('\t'.join(['species', 'Orthogroup', 'GC_deviation', 'hit_rank',
#                            'accession', 'organism', 'description', 'pident',
#                            'evalue', 'is_aspergillus', 'is_fungal']) + '\n')

# for i, cand in enumerate(remaining):
#     print(f"\n  [{i+1}/{n_remaining}] {cand['species']} / {cand['Orthogroup']} "
#           f"(GC Δ={cand['GC_deviation']:+.1f}%) ...")
# 
#     hits = None
#     for attempt in range(MAX_RETRIES):
#         try:
#             hits = ncbi_remote_blast(cand['sequence'], max_hits=MAX_HITS)
#             break
#         except Exception as e:
#             print(f"      Attempt {attempt+1}/{MAX_RETRIES} failed: {e}")
#             if attempt < MAX_RETRIES - 1:
#                 time.sleep(30)
# 
#     rows_to_write = []
#     if hits is None or len(hits) == 0:
#         print(f"    → No BLAST results returned")
#         rows_to_write.append({
#             'species': cand['species'], 'Orthogroup': cand['Orthogroup'],
#             'GC_deviation': cand['GC_deviation'],
#             'hit_rank': 1, 'accession': '', 'organism': 'No hit',
#             'description': '', 'pident': 0, 'evalue': '',
#             'is_aspergillus': False, 'is_fungal': False,
#         })
#     else:
#         for rank, h in enumerate(hits, 1):
#             org_lower = h['organism'].lower()
#             is_asp = 'aspergillus' in org_lower
#             is_fungal = any(t in org_lower for t in [
#                 'aspergillus', 'penicillium', 'fusarium', 'neurospora',
#                 'saccharomyces', 'candida', 'trichoderma', 'botrytis',
#                 'magnaporthe', 'ustilago', 'cryptococcus', 'talaromyces',
#                 'cladosporium', 'alternaria', 'colletotrichum', 'mycosphaerella',
#                 'sclerotinia', 'rhizopus', 'mucor', 'ascomycet', 'basidiomycet',
#                 'eurotiomycet', 'sordariomycet', 'dothideomycet', 'fungi',
#             ])
#             rows_to_write.append({
#                 'species': cand['species'], 'Orthogroup': cand['Orthogroup'],
#                 'GC_deviation': cand['GC_deviation'],
#                 'hit_rank': rank, 'accession': h['accession'],
#                 'organism': h['organism'], 'description': h['description'],
#                 'pident': h['pident'], 'evalue': h['evalue'],
#                 'is_aspergillus': is_asp, 'is_fungal': is_fungal,
#             })
#             print(f"    Hit {rank}: {h['organism']} | {h['pident']:.1f}% ident | "
#                   f"{'Aspergillus' if is_asp else ('Fungal' if is_fungal else '** NON-FUNGAL **')}")
# 
    # Save incrementally
#     with open(RESULTS_FILE, 'a') as f:
#         for r in rows_to_write:
#             f.write('\t'.join(str(r[c]) for c in ['species', 'Orthogroup', 'GC_deviation',
#                     'hit_rank', 'accession', 'organism', 'description', 'pident',
#                     'evalue', 'is_aspergillus', 'is_fungal']) + '\n')
# 
#     time.sleep(3)

# --- Load all results (including previous runs) and visualize ---
# print(f"\nResults saved to: {RESULTS_FILE}")
# blast_tax_df = pd.read_csv(RESULTS_FILE, sep='\t')
# print(f"Total rows: {len(blast_tax_df)}, OGs completed: {blast_tax_df['Orthogroup'].nunique()}")

# --- Summary ---
# print("\n" + "=" * 80)
# print("TAXONOMIC ORIGIN SUMMARY (top hit per xenolog candidate)")
# print("=" * 80)

# top_hits = blast_tax_df[blast_tax_df['hit_rank'] == 1]
# for sp in SPECIES_DEEP:
#     sp_hits = top_hits[top_hits['species'] == sp]
#     if sp_hits.empty:
#         continue
#     n = len(sp_hits)
#     n_asp = sp_hits['is_aspergillus'].sum()
#     n_fungal = sp_hits['is_fungal'].sum() - n_asp
#     n_nonfungal = n - sp_hits['is_fungal'].sum()
#     print(f"\n  {sp} ({n} candidates):")
#     print(f"    Closest hit is Aspergillus: {n_asp} ({100*n_asp/n:.0f}%)")
#     print(f"    Closest hit is other fungi: {n_fungal} ({100*n_fungal/n:.0f}%)")
#     print(f"    Closest hit is non-fungal:  {n_nonfungal} ({100*n_nonfungal/n:.0f}%) ← strongest HGT evidence")

# --- Visualization ---
# fig, axes = plt.subplots(1, len(SPECIES_DEEP), figsize=(4 * len(SPECIES_DEEP), 5))
# if len(SPECIES_DEEP) == 1:
#     axes = [axes]

# tax_colors = {'Aspergillus': 'forestgreen', 'Other fungi': 'steelblue',
#               'Non-fungal': '#e74c3c', 'No hit': 'lightgray'}
# tax_order = ['Aspergillus', 'Other fungi', 'Non-fungal', 'No hit']

# for idx, sp in enumerate(SPECIES_DEEP):
#     ax = axes[idx]
#     sp_hits = top_hits[top_hits['species'] == sp]
#     if sp_hits.empty:
#         ax.text(0.5, 0.5, f'A. {sp}\nNo candidates', ha='center', va='center',
#                 fontsize=11, transform=ax.transAxes)
#         ax.set_title(f'A. {sp}', fontweight='bold')
#         ax.axis('off')
#         continue
# 
#     tax_labels = []
#     for _, r in sp_hits.iterrows():
#         if r['organism'] == 'No hit':
#             tax_labels.append('No hit')
#         elif r['is_aspergillus']:
#             tax_labels.append('Aspergillus')
#         elif r['is_fungal']:
#             tax_labels.append('Other fungi')
#         else:
#             tax_labels.append('Non-fungal')
# 
#     counts = pd.Series(tax_labels).value_counts()
#     vals = [counts.get(c, 0) for c in tax_order]
#     nonzero = [(v, c, tax_colors[c]) for v, c in zip(vals, tax_order) if v > 0]
# 
#     if nonzero:
#         wedges, _ = ax.pie([x[0] for x in nonzero],
#                            colors=[x[2] for x in nonzero], startangle=90,
#                            wedgeprops=dict(linewidth=0.5, edgecolor='white'))
#         total = sum(v for v, _, _ in nonzero)
#         for j, (v, c, _) in enumerate(nonzero):
#             ang = (wedges[j].theta1 + wedges[j].theta2) / 2
#             x_pos = 0.6 * np.cos(np.radians(ang))
#             y_pos = 0.6 * np.sin(np.radians(ang))
#             ax.text(x_pos, y_pos, f"{v}\n({100*v/total:.0f}%)",
#                     ha='center', va='center', fontsize=9, fontweight='bold')
# 
#     ax.set_title(f'A. {sp}\n(n={len(sp_hits)} candidates)', fontsize=11, fontweight='bold')

from matplotlib.patches import Patch
# legend_patches = [Patch(facecolor=tax_colors[c], label=c) for c in tax_order]
# fig.legend(handles=legend_patches, fontsize=9, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.05))
# plt.suptitle('Taxonomic Origin of Xenolog Candidates (NCBI nr top hit)',
#              fontsize=13, fontweight='bold', y=1.03)
# plt.tight_layout()
# plt.savefig(f'{RESULTS_DIR}/rare_xenolog_ncbi_taxonomy.png', dpi=150, bbox_inches='tight')

# --- Detail table: non-fungal hits ---
# nonfungal = top_hits[~top_hits['is_fungal'] & (top_hits['organism'] != 'No hit')]
# if not nonfungal.empty:
#     print("\n" + "=" * 80)
#     print("NON-FUNGAL TOP HITS (strongest HGT candidates)")
#     print("=" * 80)
#     for _, r in nonfungal.iterrows():
#         print(f"\n  {r['species']} / {r['Orthogroup']}:")
#         print(f"    Organism: {r['organism']}")
#         print(f"    Identity: {r['pident']:.1f}%, E-value: {r['evalue']}")
#         print(f"    GC deviation: {r['GC_deviation']:+.1f}%")
#         print(f"    Description: {r['description']}")
# else:
#     print("\nNo non-fungal top hits found — all xenolog candidates have fungal best matches.")



# --- Functions extracted from NB4_RareGenome notebook ---

def get_accession(col):
    """Extract GCA/GCF accession from OrthoFinder column name."""
    m = re.match(r'(GC[AF]_\d+\.\d+)', col)
    return m.group(1) if m else col


def load_species_data(species_key):
    """Load PAV matrix + OG consensus (with Pangenome_Class) from NB1_Results."""
    cfg = SPECIES_CONFIG[species_key]
    print(f"\n{'='*60}")
    print(f"Loading {cfg['label']}...")

    # Long annotation table
    long_annot = pd.read_csv(cfg['annot_path'], sep='\t', low_memory=False)
    print(f"  Long annotation: {len(long_annot):,} rows, {long_annot['Orthogroup'].nunique():,} OGs, "
          f"{long_annot['Assembly Accession'].nunique()} genomes")

    # PAV matrix, dropping ANI-excluded genomes
    pav = pd.read_csv(cfg['pav_path'], sep='\t', index_col=0)
    sp_short = cfg.get('species') or cfg['label'].replace('A. ', '').strip()
    try:
        from funpan_utils import filter_ani_excluded as _fae
        pav = _fae(pav, sp_short, axis='cols')
    except Exception:
        pass
    print(f"  PAV matrix (ANI-filtered): {pav.shape[0]:,} OGs x {pav.shape[1]} genomes")

    # OG consensus with Pangenome_Class
    og_consensus = pd.read_csv(cfg['og_consensus_path'], sep='\t')
    n_core = (og_consensus['Pangenome_Class'] == 'Core').sum()
    n_acc = (og_consensus['Pangenome_Class'] == 'Accessory').sum()
    n_rare = (og_consensus['Pangenome_Class'] == 'Rare').sum()
    n_total = len(og_consensus)
    print(f"  OG consensus: {n_total} OGs")
    print(f"    Core: {n_core} ({n_core/n_total:.1%}), Accessory: {n_acc} ({n_acc/n_total:.1%}), "
          f"Rare: {n_rare} ({n_rare/n_total:.1%})")

    # Recover thresholds
    gene_freq = pav.sum(axis=1)
    core_ogs = set(og_consensus.loc[og_consensus['Pangenome_Class'] == 'Core', 'Orthogroup'])
    rare_ogs = set(og_consensus.loc[og_consensus['Pangenome_Class'] == 'Rare', 'Orthogroup'])
    core_freqs = gene_freq[gene_freq.index.isin(core_ogs)]
    rare_freqs = gene_freq[gene_freq.index.isin(rare_ogs)]
    core_n = int(core_freqs.min()) if len(core_freqs) > 0 else pav.shape[1]
    rare_n = int(rare_freqs.max()) if len(rare_freqs) > 0 else 1
    n_strains = pav.shape[1]
    print(f"  Thresholds: core >= {core_n}/{n_strains}, rare <= {rare_n}/{n_strains}")

    return {
        'long_annot': long_annot,
        'pav': pav,
        'og_consensus': og_consensus,
        'core_n': core_n,
        'rare_n': rare_n,
        'n_strains': n_strains,
    }


# Load all species
# species_data = {}
# for sp in SPECIES_LIST:
#     species_data[sp] = load_species_data(sp)

# Annotation rates by pangenome class
def annotation_rate_by_class(og_consensus, annotation_col, classes=('Core', 'Accessory', 'Rare')):
    """Compute fraction of OGs with a given annotation, per pangenome class."""
    results = {}
    for cls in classes:
        sub = og_consensus[og_consensus['Pangenome_Class'] == cls]
        if len(sub) == 0:
            results[cls] = {'n': 0, 'annotated': 0, 'rate': 0}
            continue
        has_annot = sub[annotation_col].apply(
            lambda x: bool(x) and str(x).strip() not in ('', '-', 'nan')
        ).sum()
        results[cls] = {'n': len(sub), 'annotated': int(has_annot), 'rate': has_annot / len(sub)}
    return results


# print('=== Annotation rates by pangenome class ===\n')
# annot_cols = ['PFAMs', 'CAZy', 'GOs', 'Description']
# all_rates = []

# for sp in SPECIES_LIST:
#     ogc = species_data[sp]['og_consensus']
#     label = SPECIES_CONFIG[sp]['label']
#     print(f'\n{label}:')
#     for col in annot_cols:
#         if col in ogc.columns:
#             rates = annotation_rate_by_class(ogc, col)
#             for cls, info in rates.items():
#                 all_rates.append({
#                     'Species': label, 'Annotation': col,
#                     'Class': cls, 'Rate': info['rate'],
#                     'n': info['n'], 'annotated': info['annotated']
#                 })
#                 print(f"  {col} - {cls}: {info['annotated']}/{info['n']} ({info['rate']:.1%})")

# rates_df = pd.DataFrame(all_rates)
# rates_df.to_csv(f'{RESULTS_DIR}/annotation_rates_by_class.csv', index=False)
# print(f'\nSaved: {RESULTS_DIR}/annotation_rates_by_class.csv')

# 6.1 Protein length distribution by pangenome class

# SPECIES_DEEP = ['fumigatus', 'flavus', 'niger', 'oryzae']

def og_protein_lengths(og_seq_dir, og_list):
    """Return dict {OG: [len1, len2, ...]} for requested OGs."""
    lengths = {}
    for og in og_list:
        fa = Path(og_seq_dir) / f'{og}.fa'
        if fa.exists():
            lengths[og] = [len(rec.seq) for rec in SeqIO.parse(fa, 'fasta')]
    return lengths


# print('Computing protein lengths per orthogroup ...')
# length_data = {}  # species -> DataFrame

# for sp in SPECIES_DEEP:
#     d = species_data[sp]
#     ogc = d['og_consensus']
#     og_dir = f'/datadrive/Species/Aspergillus/{sp}/orthofinder_output/'
#     res_dirs = sorted(Path(og_dir).glob('Results_*/Orthogroup_Sequences'))
#     if not res_dirs:
#         print(f'  {sp}: Orthogroup_Sequences not found, skipping')
#         continue
#     res_dir = res_dirs[0]
# 
#     all_ogs = ogc['Orthogroup'].tolist()
#     lens = og_protein_lengths(str(res_dir), all_ogs)
# 
#     rows = []
#     for _, r in ogc.iterrows():
#         og = r['Orthogroup']
#         if og in lens and lens[og]:
#             rows.append({
#                 'Orthogroup': og,
#                 'median_len': np.median(lens[og]),
#                 'mean_len': np.mean(lens[og]),
#                 'min_len': min(lens[og]),
#                 'n_proteins': len(lens[og]),
#                 'Pangenome_Class': r['Pangenome_Class'],
#             })
#     length_data[sp] = pd.DataFrame(rows)
#     print(f'  {sp}: {len(rows)} OGs with length data')

# 6.2 DIAMOND BLASTp: Rare OGs vs Core+Accessory

# CONDA_PREFIX = '/home/user/anaconda3/envs/pipeline_test2/bin'
# DIAMOND = f'{CONDA_PREFIX}/diamond'


def get_representative_seqs(og_seq_dir, og_list):
    """For each OG, return the longest protein sequence as representative."""
    reps = {}
    for og in og_list:
        fa = Path(og_seq_dir) / f'{og}.fa'
        if fa.exists():
            best = None
            for rec in SeqIO.parse(fa, 'fasta'):
                if best is None or len(rec.seq) > len(best.seq):
                    best = rec
            if best:
                best.id = og
                best.description = og
                reps[og] = best
    return reps


# print('Running DIAMOND BLASTp: Rare OGs vs Core+Accessory OGs ...')
# print('=' * 80)

# diamond_results = {}

# for sp in SPECIES_DEEP:
#     d = species_data[sp]
#     ogc = d['og_consensus']
#     og_dir = f'/datadrive/Species/Aspergillus/{sp}/orthofinder_output/'
#     res_dirs = sorted(Path(og_dir).glob('Results_*/Orthogroup_Sequences'))
#     if not res_dirs:
#         print(f'\n{sp}: Orthogroup_Sequences not found, skipping')
#         diamond_results[sp] = pd.DataFrame()
#         continue
#     res_dir = str(res_dirs[0])
# 
#     rare_ogs = ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'].tolist()
#     nonrare_ogs = ogc.loc[ogc['Pangenome_Class'].isin(['Core', 'Accessory']), 'Orthogroup'].tolist()
# 
#     if len(rare_ogs) == 0:
#         print(f'\n{sp}: 0 rare OGs -- skipping DIAMOND')
#         diamond_results[sp] = pd.DataFrame(
#             columns=['Orthogroup', 'blast_class', 'pident', 'qcovhsp', 'scovhsp', 'best_target'])
#         continue
# 
#     print(f'\n{sp}: {len(rare_ogs)} rare queries vs {len(nonrare_ogs)} core+accessory targets')
# 
#     rare_reps = get_representative_seqs(res_dir, rare_ogs)
#     nonrare_reps = get_representative_seqs(res_dir, nonrare_ogs)
#     print(f'  Representatives: {len(rare_reps)} rare, {len(nonrare_reps)} non-rare')
# 
#     tmpdir = tempfile.mkdtemp(prefix=f'diamond_{sp}_')
#     query_fa = f'{tmpdir}/rare_query.fa'
#     db_fa = f'{tmpdir}/nonrare_db.fa'
#     db_path = f'{tmpdir}/nonrare_db'
#     out_tsv = f'{tmpdir}/diamond_out.tsv'
# 
#     SeqIO.write(rare_reps.values(), query_fa, 'fasta')
#     SeqIO.write(nonrare_reps.values(), db_fa, 'fasta')
# 
#     subprocess.run([DIAMOND, 'makedb', '--in', db_fa, '-d', db_path],
#                    capture_output=True, check=True)
#     subprocess.run([
#         DIAMOND, 'blastp',
#         '-q', query_fa, '-d', db_path, '-o', out_tsv,
#         '--sensitive', '--max-target-seqs', '5',
#         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
#         'qcovhsp', 'scovhsp', 'evalue', 'bitscore',
#         '--threads', '8'
#     ], capture_output=True, check=True)
# 
#     cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
#             'qcovhsp', 'scovhsp', 'evalue', 'bitscore']
#     if os.path.getsize(out_tsv) > 0:
#         hits = pd.read_csv(out_tsv, sep='\t', names=cols)
#         best_hits = hits.sort_values('bitscore', ascending=False).drop_duplicates('qseqid', keep='first')
#     else:
#         best_hits = pd.DataFrame(columns=cols)
# 
#     classifications = []
#     for og in rare_ogs:
#         if og in rare_reps:
#             row = best_hits[best_hits['qseqid'] == og]
#             if row.empty:
#                 classifications.append({'Orthogroup': og, 'blast_class': 'No hit',
#                                         'pident': 0, 'qcovhsp': 0, 'scovhsp': 0, 'best_target': ''})
#             else:
#                 r = row.iloc[0]
#                 qcov, scov, pident = r['qcovhsp'], r['scovhsp'], r['pident']
#                 if qcov >= 70 and scov < 50 and pident >= 30:
#                     bclass = 'Fragment of longer gene'
#                 elif qcov >= 70 and scov >= 70 and pident >= 30:
#                     bclass = 'Full-length homolog'
#                 else:
#                     bclass = 'No significant similarity'
#                 classifications.append({
#                     'Orthogroup': og, 'blast_class': bclass, 'pident': pident,
#                     'qcovhsp': qcov, 'scovhsp': scov, 'best_target': r['sseqid']
#                 })
# 
#     diamond_results[sp] = pd.DataFrame(classifications)
#     shutil.rmtree(tmpdir)
# 
#     counts = diamond_results[sp]['blast_class'].value_counts()
#     print(f'  Results:')
#     for cls, cnt in counts.items():
#         print(f'    {cls}: {cnt} ({100*cnt/len(classifications):.1f}%)')

# 7.1 Identify xenolog candidates based on GC content deviation

def compute_gc_content(sequence):
    """Compute GC content of a protein-coding DNA or protein sequence."""
    seq = str(sequence).upper()
    gc = seq.count('G') + seq.count('C')
    total = len(seq)
    return 100.0 * gc / total if total > 0 else 0


# xenolog_candidates = {}
GC_THRESHOLD = 2.0  # standard deviations from mean

# for sp in SPECIES_DEEP:
#     d = species_data[sp]
#     ogc = d['og_consensus']
#     rare_ogc = ogc[ogc['Pangenome_Class'] == 'Rare'].copy()
# 
#     if len(rare_ogc) < 5:
#         print(f'{sp}: too few rare OGs ({len(rare_ogc)}) for xenolog detection')
#         xenolog_candidates[sp] = pd.DataFrame()
#         continue
# 
#     og_dir = f'/datadrive/Species/Aspergillus/{sp}/orthofinder_output/'
#     res_dirs = sorted(Path(og_dir).glob('Results_*/Orthogroup_Sequences'))
#     if not res_dirs:
#         print(f'{sp}: Orthogroup_Sequences not found')
#         xenolog_candidates[sp] = pd.DataFrame()
#         continue
#     res_dir = str(res_dirs[0])
# 
    # Compute GC of all OG representative sequences
#     reps = get_representative_seqs(res_dir, ogc['Orthogroup'].tolist())
#     gc_data = {og: compute_gc_content(reps[og].seq) for og in reps}
# 
#     gc_series = pd.Series(gc_data)
#     gc_mean = gc_series.mean()
#     gc_std = gc_series.std()
# 
#     rare_gc = gc_series[gc_series.index.isin(rare_ogc['Orthogroup'])]
#     rare_ogc = rare_ogc.copy()
#     rare_ogc['GC_content'] = rare_ogc['Orthogroup'].map(gc_data)
#     rare_ogc['GC_deviation'] = (rare_ogc['GC_content'] - gc_mean) / gc_std
# 
#     candidates = rare_ogc[rare_ogc['GC_deviation'].abs() > GC_THRESHOLD].copy()
#     xenolog_candidates[sp] = candidates
# 
#     print(f'{sp}: {len(candidates)} xenolog candidates '
#           f'(|GC deviation| > {GC_THRESHOLD} SD, mean GC={gc_mean:.1f}%)')

# 7.2 NCBI remote BLASTp for taxonomic origin

NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
# MAX_HITS = 5
# RESULTS_FILE = f'{RESULTS_DIR}/rare_xenolog_ncbi_blast_results.tsv'
# MAX_RETRIES = 3


def ncbi_remote_blast(sequence, program='blastp', database='nr', max_hits=5, max_wait=300):
    """Submit a sequence to NCBI BLAST and return parsed top hits."""
    params = {
        'CMD': 'Put', 'PROGRAM': program, 'DATABASE': database,
        'QUERY': str(sequence), 'HITLIST_SIZE': str(max_hits), 'FORMAT_TYPE': 'XML',
    }
    resp = requests.post(NCBI_BLAST_URL, data=params, timeout=60)
    resp.raise_for_status()

    rid = None
    for line in resp.text.split('\n'):
        if 'RID = ' in line:
            rid = line.split('RID = ')[1].strip()
            break
    if not rid:
        return None

    elapsed = 0
    while elapsed < max_wait:
        time.sleep(15)
        elapsed += 15
        try:
            check = requests.get(NCBI_BLAST_URL,
                                 params={'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'},
                                 timeout=120)
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError,
                requests.exceptions.Timeout) as e:
            print(f'      Network error polling RID {rid}: {e} -- retrying ...')
            continue

        if 'Status=WAITING' in check.text:
            continue
        if 'Status=FAILED' in check.text:
            return None
        if '<?xml' in check.text:
            try:
                root = ET.fromstring(check.text)
                hits = []
                for hit in root.iter('Hit'):
                    hit_def = hit.find('Hit_def').text if hit.find('Hit_def') is not None else ''
                    hit_acc = hit.find('Hit_accession').text if hit.find('Hit_accession') is not None else ''
                    hsp = hit.find('.//Hsp')
                    if hsp is not None:
                        identity = float(hsp.find('Hsp_identity').text)
                        align_len = float(hsp.find('Hsp_align-len').text)
                        pident = 100 * identity / align_len if align_len > 0 else 0
                        evalue = hsp.find('Hsp_evalue').text
                    else:
                        pident, evalue = 0, 'N/A'
                    organism = ''
                    if '[' in hit_def and ']' in hit_def:
                        organism = hit_def[hit_def.rfind('[') + 1:hit_def.rfind(']')]
                    hits.append({
                        'accession': hit_acc, 'description': hit_def[:120],
                        'organism': organism, 'pident': pident, 'evalue': evalue,
                    })
                    if len(hits) >= max_hits:
                        break
                return hits
            except ET.ParseError:
                return None
    return None


# Collect xenolog candidate sequences
# all_candidates = []
# for sp in SPECIES_DEEP:
#     xc = xenolog_candidates.get(sp)
#     if xc is None or xc.empty:
#         continue
#     og_dir = f'/datadrive/Species/Aspergillus/{sp}/orthofinder_output/'
#     res_dir = str(sorted(Path(og_dir).glob('Results_*/Orthogroup_Sequences'))[0])
#     reps = get_representative_seqs(res_dir, xc['Orthogroup'].tolist())
#     for _, row in xc.iterrows():
#         og = row['Orthogroup']
#         if og in reps:
#             all_candidates.append({
#                 'species': sp, 'Orthogroup': og,
#                 'sequence': str(reps[og].seq),
#                 'GC_deviation': row['GC_deviation'],
#                 'description': row.get('Description', ''),
#             })

# Load existing results (resume support)
# completed_ogs = set()
# if os.path.exists(RESULTS_FILE):
#     existing = pd.read_csv(RESULTS_FILE, sep='\t')
#     completed_ogs = set(existing['Orthogroup'].unique())
#     print(f'Resuming: {len(completed_ogs)} OGs already completed in {RESULTS_FILE}')

# remaining = [c for c in all_candidates if c['Orthogroup'] not in completed_ogs]
# n_total = len(all_candidates)
# n_remaining = len(remaining)

# print(f'Total xenolog candidates: {n_total}')
# print(f'Already done: {n_total - n_remaining}, remaining: {n_remaining}')
# if n_remaining > 0:
#     print(f'Estimated time: ~{n_remaining * 45 // 60} - {n_remaining * 75 // 60} minutes')
# print('=' * 80)

# Write header if file doesn't exist
# if not os.path.exists(RESULTS_FILE):
#     with open(RESULTS_FILE, 'w') as f:
#         f.write('\t'.join(['species', 'Orthogroup', 'GC_deviation', 'hit_rank',
#                            'accession', 'organism', 'description', 'pident',
#                            'evalue', 'is_aspergillus', 'is_fungal']) + '\n')

# for i, cand in enumerate(remaining):
#     print(f"\n  [{i+1}/{n_remaining}] {cand['species']} / {cand['Orthogroup']} "
#           f"(GC delta={cand['GC_deviation']:+.1f}%) ...")
# 
#     hits = None
#     for attempt in range(MAX_RETRIES):
#         try:
#             hits = ncbi_remote_blast(cand['sequence'], max_hits=MAX_HITS)
#             break
#         except Exception as e:
#             print(f'      Attempt {attempt+1}/{MAX_RETRIES} failed: {e}')
#             if attempt < MAX_RETRIES - 1:
#                 time.sleep(30)
# 
#     rows_to_write = []
#     if hits is None or len(hits) == 0:
#         print(f'    -> No BLAST results returned')
#         rows_to_write.append({
#             'species': cand['species'], 'Orthogroup': cand['Orthogroup'],
#             'GC_deviation': cand['GC_deviation'],
#             'hit_rank': 1, 'accession': '', 'organism': 'No hit',
#             'description': '', 'pident': 0, 'evalue': '',
#             'is_aspergillus': False, 'is_fungal': False,
#         })
#     else:
#         for rank, h in enumerate(hits, 1):
#             org_lower = h['organism'].lower()
#             is_asp = 'aspergillus' in org_lower
#             is_fungal = any(t in org_lower for t in [
#                 'aspergillus', 'penicillium', 'fusarium', 'neurospora',
#                 'saccharomyces', 'candida', 'trichoderma', 'botrytis',
#                 'magnaporthe', 'ustilago', 'cryptococcus', 'talaromyces',
#                 'cladosporium', 'alternaria', 'colletotrichum', 'mycosphaerella',
#                 'sclerotinia', 'rhizopus', 'mucor', 'ascomycet', 'basidiomycet',
#                 'eurotiomycet', 'sordariomycet', 'dothideomycet', 'fungi',
#             ])
#             rows_to_write.append({
#                 'species': cand['species'], 'Orthogroup': cand['Orthogroup'],
#                 'GC_deviation': cand['GC_deviation'],
#                 'hit_rank': rank, 'accession': h['accession'],
#                 'organism': h['organism'], 'description': h['description'],
#                 'pident': h['pident'], 'evalue': h['evalue'],
#                 'is_aspergillus': is_asp, 'is_fungal': is_fungal,
#             })
#             print(f"    Hit {rank}: {h['organism']} | {h['pident']:.1f}% ident | "
#                   f"{'Aspergillus' if is_asp else ('Fungal' if is_fungal else '** NON-FUNGAL **')}")
# 
    # Save incrementally
#     with open(RESULTS_FILE, 'a') as f:
#         for r in rows_to_write:
#             f.write('\t'.join(str(r[c]) for c in ['species', 'Orthogroup', 'GC_deviation',
#                     'hit_rank', 'accession', 'organism', 'description', 'pident',
#                     'evalue', 'is_aspergillus', 'is_fungal']) + '\n')
# 
#     time.sleep(3)

# print(f'\nResults saved to: {RESULTS_FILE}')
# if os.path.exists(RESULTS_FILE):
#     blast_tax_df = pd.read_csv(RESULTS_FILE, sep='\t')
#     print(f'Total rows: {len(blast_tax_df)}, OGs completed: {blast_tax_df["Orthogroup"].nunique()}')

# Helper functions for phylogenetic analysis

def normalize_sample_id(sample_id):
    """Normalize sample ID by removing version suffix."""
    if pd.isna(sample_id):
        return ''
    s = str(sample_id).strip()
    if re.match(r'^GC[AF]_\d+\.\d+$', s):
        return s.rsplit('.', 1)[0]
    return s


def parsimony_score(tree, tip_states):
    """Fitch parsimony: count minimum character-state changes on tree."""
    state_sets = {}
    for tip in tree.get_terminals():
        name = tip.name
        state_sets[id(tip)] = {tip_states.get(name, 'Unknown')}

    def _fitch_down(clade):
        if clade.is_terminal():
            return state_sets[id(clade)]
        child_sets = [_fitch_down(c) for c in clade.clades]
        inter = child_sets[0]
        for cs in child_sets[1:]:
            inter = inter & cs
        if inter:
            state_sets[id(clade)] = inter
        else:
            state_sets[id(clade)] = child_sets[0] | child_sets[1]
        return state_sets[id(clade)]

    _fitch_down(tree.root)

    changes = 0
    def _fitch_up(clade, parent_state):
        nonlocal changes
        my_set = state_sets[id(clade)]
        if parent_state in my_set:
            chosen = parent_state
        else:
            chosen = list(my_set)[0]
            changes += 1
        for c in clade.clades:
            _fitch_up(c, chosen)

    root_state = list(state_sets[id(tree.root)])[0]
    for c in tree.root.clades:
        _fitch_up(c, root_state)
    return changes


def fritz_purvis_d(tree, tip_states, n_permutations=1000):
    """Fritz & Purvis D statistic for a binary trait."""
    observed = parsimony_score(tree, tip_states)

    tips = [t.name for t in tree.get_terminals()]
    values = [tip_states.get(t, 0) for t in tips]
    n1 = sum(values)
    n0 = len(values) - n1

    if n1 == 0 or n0 == 0:
        return {'D': np.nan, 'p_brownian': np.nan, 'p_random': np.nan, 'observed_changes': observed}

    # Random permutation distribution
    random_scores = []
    rng = np.random.default_rng(42)
    for _ in range(n_permutations):
        perm_vals = list(values)
        rng.shuffle(perm_vals)
        perm_states = dict(zip(tips, perm_vals))
        random_scores.append(parsimony_score(tree, perm_states))

    random_mean = np.mean(random_scores)

    # Brownian expectation
    brownian_expected = n1 * n0 / len(values) * 0.5

    denom = random_mean - brownian_expected
    if abs(denom) < 1e-10:
        D = 0
    else:
        D = (observed - brownian_expected) / denom

    p_random = np.mean([s <= observed for s in random_scores])
    p_brownian = np.mean([s >= observed for s in random_scores])

    return {'D': D, 'p_brownian': p_brownian, 'p_random': p_random,
            'observed_changes': observed, 'random_mean': random_mean}


def fitch_gain_loss(tree, pav_matrix, og_list, phenotype_map):
    """Fitch parsimony for gene gain/loss on each branch."""
    tips = [t.name for t in tree.get_terminals()]
    events = []

    for og in og_list:
        if og not in pav_matrix.index:
            continue
        tip_states = {}
        for tip_name in tips:
            # Match tip name to PAV columns
            matched_col = None
            for col in pav_matrix.columns:
                if tip_name in col or get_accession(col) == tip_name:
                    matched_col = col
                    break
            if matched_col is not None:
                tip_states[tip_name] = int(pav_matrix.loc[og, matched_col])
            else:
                tip_states[tip_name] = 0

        n_present = sum(tip_states.values())
        if n_present == 0 or n_present == len(tips):
            continue

        changes = parsimony_score(tree, tip_states)
        events.append({
            'Orthogroup': og,
            'n_present': n_present,
            'n_absent': len(tips) - n_present,
            'parsimony_changes': changes,
        })

    return pd.DataFrame(events)


# print('Phylogenetic helper functions defined')



# --- Functions extracted from NB4_RareGenome notebook (final pass) ---


def compute_rare_genome_summary(species_list, species_data, species_config, results_dir=None):
    """Compute rare genome summary table across all species.

    Parameters
    ----------
    species_list : list of str
    species_data : dict
        ``{species: {og_consensus, pav, n_strains, core_n, rare_n, ...}}``.
    species_config : dict
        ``{species: {label, ...}}``.
    results_dir : str or None
        If provided, save CSV to this directory.

    Returns
    -------
    pd.DataFrame
        Summary table with one row per species.
    """
    rows = []
    for sp in species_list:
        d = species_data[sp]
        ogc = d['og_consensus']
        pav = d['pav']

        n_core = (ogc['Pangenome_Class'] == 'Core').sum()
        n_acc = (ogc['Pangenome_Class'] == 'Accessory').sum()
        n_rare = (ogc['Pangenome_Class'] == 'Rare').sum()
        n_total = len(ogc)

        rare_ogs = set(ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'])
        rare_pav = pav.loc[pav.index.isin(rare_ogs)]
        mean_rare_per_genome = rare_pav.sum(axis=0).mean() if not rare_pav.empty else 0

        rare_ogc = ogc[ogc['Pangenome_Class'] == 'Rare']
        pfam_rate = rare_ogc['PFAMs'].apply(
            lambda x: bool(x) and str(x).strip() not in ('', '-', 'nan')
        ).mean() if 'PFAMs' in rare_ogc.columns and len(rare_ogc) > 0 else 0
        cazy_rate = rare_ogc['CAZy'].apply(
            lambda x: bool(x) and str(x).strip() not in ('', '-', 'nan')
        ).mean() if 'CAZy' in rare_ogc.columns and len(rare_ogc) > 0 else 0

        rows.append({
            'Species': species_config[sp]['label'],
            'Genomes': d['n_strains'],
            'Total OGs': n_total,
            'Core': n_core,
            'Accessory': n_acc,
            'Rare': n_rare,
            'Mean rare genes/genome': round(mean_rare_per_genome, 1),
            'Rare w/ Pfam': f'{pfam_rate:.1%}',
            'Rare w/ CAZy': f'{cazy_rate:.1%}',
            'Core threshold': d['core_n'],
            'Rare threshold': d['rare_n'],
        })

    summary_df = pd.DataFrame(rows)

    if results_dir is not None:
        summary_df.to_csv(f'{results_dir}/rare_genome_summary.csv', index=False)
        print(f'Saved: {results_dir}/rare_genome_summary.csv')

    return summary_df


def run_rare_cog_enrichment(species_list, species_data, species_config,
                            results_dir=None):
    """Fisher's exact test for COG category enrichment in rare vs non-rare OGs.

    Parameters
    ----------
    species_list : list of str
    species_data : dict
    species_config : dict
    results_dir : str or None

    Returns
    -------
    pd.DataFrame
        COG enrichment results with q-values.
    """
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    from collections import Counter
    import funpan_utils as fpu

    all_cog_results = []

    for sp in species_list:
        ogc = species_data[sp]['og_consensus']
        label = species_config[sp]['label']

        if 'COG_category' not in ogc.columns:
            print(f'{label}: COG_category column not found, skipping')
            continue

        rare_ogs = ogc[ogc['Pangenome_Class'] == 'Rare']
        nonrare_ogs = ogc[ogc['Pangenome_Class'].isin(['Core', 'Accessory'])]

        if len(rare_ogs) < 2:
            print(f'{label}: only {len(rare_ogs)} rare OGs, skipping enrichment')
            continue

        def count_cogs(subset):
            counts = Counter()
            for cog_str in subset['COG_category'].dropna():
                for c in str(cog_str):
                    if c in fpu.COG_DESCRIPTIONS:
                        counts[c] += 1
            return counts

        rare_cogs = count_cogs(rare_ogs)
        nonrare_cogs = count_cogs(nonrare_ogs)

        n_rare_total = sum(rare_cogs.values())
        n_nonrare_total = sum(nonrare_cogs.values())

        for cog in sorted(set(list(rare_cogs.keys()) + list(nonrare_cogs.keys()))):
            # Skip COG-S (Function unknown) and COG-R (General function
            # prediction only) -- standard practice in functional enrichment.
            if not fpu.is_informative_term('COG_category', cog):
                continue
            a = rare_cogs.get(cog, 0)
            b = nonrare_cogs.get(cog, 0)
            c = n_rare_total - a
            d = n_nonrare_total - b

            if a + b < 3:
                continue

            odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            all_cog_results.append({
                'Species': label,
                'COG': cog,
                'COG_description': fpu.COG_DESCRIPTIONS.get(cog, ''),
                'Rare_count': a,
                'NonRare_count': b,
                'Rare_freq': a / n_rare_total if n_rare_total > 0 else 0,
                'NonRare_freq': b / n_nonrare_total if n_nonrare_total > 0 else 0,
                'Odds_Ratio': odds,
                'p_value': pval,
            })

    cog_enrich_df = pd.DataFrame(all_cog_results)
    if len(cog_enrich_df) > 0:
        cog_enrich_df['q_value'] = multipletests(
            cog_enrich_df['p_value'], method='fdr_bh')[1]
        cog_enrich_df['Significant'] = cog_enrich_df['q_value'] < 0.05
        cog_enrich_df = cog_enrich_df.sort_values(['Species', 'q_value'])

    if results_dir is not None:
        cog_enrich_df.to_csv(
            f'{results_dir}/cog_enrichment_rare_vs_nonrare.csv', index=False)
        print(f'Saved: {results_dir}/cog_enrichment_rare_vs_nonrare.csv')

    sig = cog_enrich_df[cog_enrich_df['Significant']] if len(cog_enrich_df) > 0 else pd.DataFrame()
    print(f'\nSignificant enrichments (FDR < 0.05): {len(sig)}')

    return cog_enrich_df


def compute_rare_burden_by_phenotype(species_list, species_data, species_config,
                                     pheno_map, get_accession_fn, results_dir=None,
                                     exclude=('Unknown', 'unknown', '', 'Lab', None),
                                     min_cases=3, min_controls=3):
    """Per-species, one-vs-rest Mann-Whitney U test of rare gene burden.

    Mirrors the contrast convention used in NB2 (Pan-GWAS): for every phenotype
    category present in a species, compare burden of "case" strains (with that
    phenotype) against the rest of the *known* strains (i.e. all other strains
    with a non-excluded phenotype label).

    Parameters
    ----------
    species_list : list of str
    species_data : dict
    species_config : dict
    pheno_map : dict
        ``{accession: phenotype_label}``.
    get_accession_fn : callable
        Function to extract accession from PAV column name.
    results_dir : str or None
    exclude : tuple
        Phenotype labels to drop entirely (treated as "unknown").
    min_cases, min_controls : int
        Skip a contrast if either group has fewer samples.

    Returns
    -------
    pd.DataFrame
        Burden test results, one row per (species, phenotype) contrast.
        Columns: Species, Group1, Group2, n_g1, n_g2, mean_g1, mean_g2,
        U_stat, p_value, direction.
    """
    from scipy.stats import mannwhitneyu

    burden_results = []

    for sp in species_list:
        d = species_data[sp]
        pav = d['pav']
        ogc = d['og_consensus']
        label = species_config[sp]['label']

        rare_ogs = set(ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'])
        rare_pav = pav.loc[pav.index.isin(rare_ogs)]

        if rare_pav.empty or len(rare_ogs) < 2:
            print(f'{label}: insufficient rare OGs ({len(rare_ogs)}), skipping')
            continue

        burden = rare_pav.sum(axis=0)
        pheno_map_sp = {
            col: pheno_map.get(get_accession_fn(col), 'Unknown')
            for col in pav.columns
        }

        # Restrict to "known" strains (drop excluded labels) so controls are meaningful
        known_cols = [c for c in pav.columns if pheno_map_sp.get(c) not in exclude]
        if len(known_cols) < (min_cases + min_controls):
            print(f'{label}: only {len(known_cols)} known strains, skipping')
            continue

        # Distinct phenotype categories present
        cats = sorted({pheno_map_sp[c] for c in known_cols})

        print(f"\n{'='*60}")
        print(f'{label}: rare gene burden, one-vs-rest')

        for target in cats:
            case_cols = [c for c in known_cols if pheno_map_sp[c] == target]
            ctrl_cols = [c for c in known_cols if pheno_map_sp[c] != target]
            n_case, n_ctrl = len(case_cols), len(ctrl_cols)
            if n_case < min_cases or n_ctrl < min_controls:
                print(f'  [skip] {target}: n_case={n_case}, n_ctrl={n_ctrl}')
                continue

            g1_burden = burden[case_cols].values
            g2_burden = burden[ctrl_cols].values
            u_stat, p_val = mannwhitneyu(g1_burden, g2_burden, alternative='two-sided')
            direction = 'higher' if np.mean(g1_burden) > np.mean(g2_burden) else 'lower'
            print(f'  {target} (n={n_case}, mean={np.mean(g1_burden):.1f}) vs '
                  f'rest (n={n_ctrl}, mean={np.mean(g2_burden):.1f})  '
                  f'U={u_stat:.0f}, P={p_val:.2e}  [{target} {direction}]')
            burden_results.append({
                'Species': label, 'Group1': target, 'Group2': 'Rest',
                'n_g1': n_case, 'n_g2': n_ctrl,
                'mean_g1': float(np.mean(g1_burden)), 'mean_g2': float(np.mean(g2_burden)),
                'U_stat': float(u_stat), 'p_value': float(p_val),
                'direction': direction,
            })

    burden_df = pd.DataFrame(burden_results)
    if len(burden_df) > 0 and results_dir is not None:
        burden_df.to_csv(f'{results_dir}/rare_gene_burden_by_phenotype.csv', index=False)
        print(f'\nSaved: {results_dir}/rare_gene_burden_by_phenotype.csv')

    return burden_df


def build_full_genome_burden_table(species_list, species_data, species_config,
                                    pheno_map, get_accession_fn, results_dir=None):
    """Build per-genome burden table (core / accessory / rare gene counts).

    Parameters
    ----------
    species_list : list of str
    species_data : dict
    species_config : dict
    pheno_map : dict
    get_accession_fn : callable
    results_dir : str or None

    Returns
    -------
    pd.DataFrame
        One row per genome.
    """
    all_genome_burden = []

    for sp in species_list:
        d = species_data[sp]
        pav = d['pav']
        ogc = d['og_consensus']
        label = species_config[sp]['label']

        rare_ogs = set(ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'])
        acc_ogs = set(ogc.loc[ogc['Pangenome_Class'] == 'Accessory', 'Orthogroup'])
        core_ogs = set(ogc.loc[ogc['Pangenome_Class'] == 'Core', 'Orthogroup'])

        rare_pav = pav.loc[pav.index.isin(rare_ogs)]
        acc_pav = pav.loc[pav.index.isin(acc_ogs)]
        core_pav = pav.loc[pav.index.isin(core_ogs)]

        for genome_col in pav.columns:
            acc = get_accession_fn(genome_col)
            phenotype = pheno_map.get(acc, 'Unknown')
            all_genome_burden.append({
                'Species': label,
                'Genome': acc,
                'Phenotype': phenotype,
                'Core_genes': int(core_pav[genome_col].sum()) if not core_pav.empty else 0,
                'Accessory_genes': int(acc_pav[genome_col].sum()) if not acc_pav.empty else 0,
                'Rare_genes': int(rare_pav[genome_col].sum()) if not rare_pav.empty else 0,
                'Total_genes': int(pav[genome_col].sum()),
            })

    burden_all = pd.DataFrame(all_genome_burden)
    print(f'Total genomes with phenotype labels: {len(burden_all)}')
    print(burden_all.groupby(['Species', 'Phenotype'])['Rare_genes'].agg(
        ['mean', 'std', 'count']).round(1))

    if results_dir is not None:
        burden_all.to_csv(f'{results_dir}/rare_gene_burden_all_genomes.csv', index=False)
        print(f'\nSaved: {results_dir}/rare_gene_burden_all_genomes.csv')

    return burden_all


def _find_diamond_binary():
    """Locate the DIAMOND executable.

    Tries (in order): ``$DIAMOND`` env var, ``diamond`` on PATH, then known
    conda env locations. Returns the path or raises FileNotFoundError.
    """
    cand = os.environ.get('DIAMOND')
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    on_path = shutil.which('diamond')
    if on_path:
        return on_path
    for env in ('funpan', 'pipeline_test3', 'pipeline_test2'):
        p = f'/home/user/anaconda3/envs/{env}/bin/diamond'
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    raise FileNotFoundError(
        'DIAMOND binary not found. Install via `conda install -c bioconda diamond` '
        'or set the $DIAMOND environment variable to the executable path.'
    )


def find_gc_outlier_rare_ogs(species_list, species_data, base_path,
                              gc_threshold=2.0):
    """Find rare OGs with anomalous GC content -- candidate xenologs.

    For each species:
      1. Compute G+C frequency of each OG's representative protein sequence
         (longest member). Note this is *amino-acid* G+C, a coarse proxy for
         codon GC and therefore for HGT signal.
      2. Compute the species-wide mean and SD across *all* OGs.
      3. Return rare OGs whose deviation exceeds ``gc_threshold`` SDs.

    Returns
    -------
    dict
        ``{species: DataFrame}`` with columns Orthogroup, GC_content,
        GC_deviation, Description (if available).
    """
    candidates = {}
    for sp in species_list:
        if sp not in species_data or 'og_consensus' not in species_data[sp]:
            candidates[sp] = pd.DataFrame()
            continue
        ogc = species_data[sp]['og_consensus']
        rare_ogc = ogc[ogc['Pangenome_Class'] == 'Rare'].copy()
        if rare_ogc.empty:
            candidates[sp] = pd.DataFrame()
            continue

        og_seq_root = Path(base_path) / sp / 'orthofinder_output'
        res_dirs = sorted(og_seq_root.glob('Results_*/Orthogroup_Sequences'))
        if not res_dirs:
            print(f'  {sp}: Orthogroup_Sequences not found, skipping GC filter')
            candidates[sp] = pd.DataFrame()
            continue
        og_seq_dir = str(res_dirs[-1])

        all_ogs = ogc['Orthogroup'].tolist()
        reps = get_representative_seqs(og_seq_dir, all_ogs)
        gc_data = {og: compute_gc_content(reps[og].seq) for og in reps}
        gc_series = pd.Series(gc_data, dtype=float)
        gc_mean = float(gc_series.mean())
        gc_std = float(gc_series.std())

        rare_ogc['GC_content'] = rare_ogc['Orthogroup'].map(gc_data)
        rare_ogc['GC_deviation'] = (rare_ogc['GC_content'] - gc_mean) / gc_std
        keep_cols = ['Orthogroup', 'GC_content', 'GC_deviation']
        if 'Description' in rare_ogc.columns:
            keep_cols.append('Description')
        outliers = rare_ogc.loc[
            rare_ogc['GC_deviation'].abs() > gc_threshold, keep_cols
        ].copy()
        candidates[sp] = outliers
        print(f'  {sp}: {len(outliers)} GC outliers '
              f'(|dev|>{gc_threshold} SD; mean GC={gc_mean:.1f}%, n_rare={len(rare_ogc)})')
    return candidates


def run_ncbi_blast_xenologs(species_list, species_data, species_config,
                             base_path, results_dir,
                             n_per_species=20, sleep_s=2.0, force=False,
                             gc_filter=False, gc_threshold=2.0):
    """Pick rare-OG representatives and submit each one to NCBI remote BLAST.

    Resumable: an existing TSV at
    ``{results_dir}/rare_xenolog_ncbi_blast_results.tsv`` is loaded, and only
    new OGs are queried (unless ``force=True``).

    Parameters
    ----------
    species_list : list of str
    species_data : dict
        Per-species ``{'og_consensus': DataFrame}``.
    species_config : dict
    base_path : str
        Path under which ``{species}/orthofinder_output/Results_*/Orthogroup_Sequences``
        lives.
    results_dir : str
    n_per_species : int
        Number of rare OGs to query per species.
    sleep_s : float
        Pause between requests to be polite to NCBI.
    force : bool
        If True, re-BLAST every OG even if already cached.

    Returns
    -------
    pd.DataFrame
        One row per BLAST hit (multiple per OG).
    """
    out_fp = Path(results_dir) / 'rare_xenolog_ncbi_blast_results.tsv'

    # Always back up an existing non-empty cache before any destructive op
    if out_fp.exists() and out_fp.stat().st_size > 50:  # >50 bytes = has data
        bak_fp = out_fp.with_suffix(out_fp.suffix + '.bak')
        shutil.copy2(out_fp, bak_fp)
        print(f'Backed up existing cache -> {bak_fp}')

    if force:
        completed = set()
        existing = pd.DataFrame()
    elif out_fp.exists() and out_fp.stat().st_size > 0:
        try:
            existing = pd.read_csv(out_fp, sep='\t')
        except pd.errors.EmptyDataError:
            existing = pd.DataFrame()
        completed = (set(zip(existing['species'].astype(str),
                             existing['orthogroup'].astype(str)))
                     if 'species' in existing.columns else set())
        print(f'Resuming: {len(completed)} (species, OG) pairs already in {out_fp}')
    else:
        completed = set()
        existing = pd.DataFrame()

    if gc_filter:
        print(f'GC-deviation filter active (|GC_dev| > {gc_threshold} SD); '
              f'BLASTing all candidates (n_per_species ignored)')
        gc_candidates = find_gc_outlier_rare_ogs(
            species_list, species_data, base_path, gc_threshold=gc_threshold,
        )
    else:
        gc_candidates = None

    # First pass: build the to-do queue across all species so we can print
    # an overall progress index (e.g. "[12/47] flavus OG0001234").
    queue = []  # list of (sp, og, og_seq_dir)
    for sp in species_list:
        if sp not in species_data or 'og_consensus' not in species_data[sp]:
            continue
        ogc = species_data[sp]['og_consensus']
        if gc_filter:
            cand_df = gc_candidates.get(sp, pd.DataFrame())
            rare_ogs = cand_df['Orthogroup'].tolist() if not cand_df.empty else []
        else:
            rare_ogs = ogc.loc[
                ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'
            ].head(n_per_species).tolist()
        if not rare_ogs:
            continue

        og_seq_root = Path(base_path) / sp / 'orthofinder_output'
        res_dirs = sorted(og_seq_root.glob('Results_*/Orthogroup_Sequences'))
        if not res_dirs:
            print(f'  {sp}: Orthogroup_Sequences not found under {og_seq_root}, skipping',
                  flush=True)
            continue
        og_seq_dir = str(res_dirs[-1])

        todo = [og for og in rare_ogs if (sp, og) not in completed]
        skipped = len(rare_ogs) - len(todo)
        print(f'  {sp}: queued {len(todo)} new rare OGs (skipping {skipped} cached)',
              flush=True)
        for og in todo:
            queue.append((sp, og, og_seq_dir))

    total = len(queue)
    if total == 0:
        print('Nothing to BLAST.', flush=True)
        out_fp.parent.mkdir(parents=True, exist_ok=True)
        existing.to_csv(out_fp, sep='\t', index=False)
        return existing

    print(f'\nStarting BLAST run -- {total} OGs to query (~{total * 30 // 60} min estimate)\n',
          flush=True)

    # Cache reps per og_seq_dir to avoid re-parsing the directory each time
    reps_cache = {}
    new_rows = []
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    run_start = time.time()

    for idx, (sp, og, og_seq_dir) in enumerate(queue, start=1):
        if og_seq_dir not in reps_cache:
            # Pre-load all OGs for this species' dir at first encounter
            sp_todo_ogs = [_og for _sp, _og, _d in queue if _d == og_seq_dir]
            reps_cache[og_seq_dir] = get_representative_seqs(og_seq_dir, sp_todo_ogs)
        reps = reps_cache[og_seq_dir]
        rec = reps.get(og)

        og_start = time.time()
        prefix = f'[{idx:>3}/{total}] {sp} {og}'

        if rec is None:
            # Mark as completed with a sentinel so we don't retry every run
            new_rows.append({
                'species': sp, 'orthogroup': og, 'accession': 'NO_REPRESENTATIVE',
                'description': 'no representative sequence found',
                'organism': '', 'pident': float('nan'), 'evalue': float('nan'),
            })
            print(f'{prefix}  -> no representative sequence, marking done', flush=True)
            continue

        print(f'{prefix}  submitting (len={len(rec.seq)}aa) ...', flush=True)
        failed = False
        try:
            hits = ncbi_remote_blast(str(rec.seq), max_hits=3) or []
        except Exception as exc:
            print(f'{prefix}  FAILED ({exc})', flush=True)
            hits = []
            failed = True

        elapsed = time.time() - og_start
        run_elapsed = time.time() - run_start
        avg_per_og = run_elapsed / idx
        eta = avg_per_og * (total - idx)
        top_org = hits[0].get('organism', '?') if hits else 'no hits'
        print(f'{prefix}  -> {len(hits)} hits, top="{top_org}", '
              f'{elapsed:.1f}s  |  ETA {eta/60:.1f} min',
              flush=True)

        if hits:
            for h in hits:
                new_rows.append({'species': sp, 'orthogroup': og, **h})
        else:
            # Always record at least one placeholder so the OG is in the
            # `completed` set on the next run and we don't keep re-querying it.
            sentinel_acc = 'BLAST_FAILED' if failed else 'NO_HIT'
            sentinel_desc = ('BLAST submission failed' if failed
                             else 'no NCBI hits returned')
            new_rows.append({
                'species': sp, 'orthogroup': og,
                'accession': sentinel_acc, 'description': sentinel_desc,
                'organism': '', 'pident': float('nan'), 'evalue': float('nan'),
            })

        # Incremental save: every OG, so a kernel kill loses nothing
        new_df_so_far = pd.DataFrame(new_rows)
        snapshot = (pd.concat([existing, new_df_so_far], ignore_index=True)
                    if not new_df_so_far.empty else existing)
        snapshot.to_csv(out_fp, sep='\t', index=False)

        time.sleep(sleep_s)

    new_df = pd.DataFrame(new_rows)
    blast_df = pd.concat([existing, new_df], ignore_index=True) if not new_df.empty else existing
    blast_df.to_csv(out_fp, sep='\t', index=False)

    total_min = (time.time() - run_start) / 60
    print(f'\nDone in {total_min:.1f} min. Saved {len(blast_df)} total BLAST hits '
          f'({len(new_df)} new) to {out_fp}', flush=True)
    return blast_df


def run_diamond_rare_vs_nonrare(species_list, species_data, species_config,
                                base_path, results_dir,
                                threads=8, force=False):
    """Run DIAMOND blastp of rare-OG representatives vs Core+Accessory representatives.

    For each species:
      1. Pick the longest protein in each OG as its representative.
      2. blastp the rare reps (queries) against Core+Accessory reps (database).
      3. Classify the best hit per query as one of:
         * ``Full-length homolog``       qcov>=70, scov>=70, pident>=30
         * ``Fragment of longer gene``   qcov>=70, scov<50,  pident>=30
         * ``No significant similarity`` else (a hit was found but weak)
         * ``No hit``                    no DIAMOND hit at all

    Cached per species at ``{results_dir}/{sp}/diamond_rare_classification.tsv``.
    Pass ``force=True`` to ignore the cache and re-run.

    Parameters
    ----------
    species_list : list of str
    species_data : dict
        Must contain ``og_consensus`` per species with ``Pangenome_Class`` column.
    species_config : dict
    base_path : str
        Path to the species directory (e.g. ``/datadrive/Species/Aspergillus``)
        whose subdirectories contain ``orthofinder_output/Results_*/Orthogroup_Sequences``.
    results_dir : str
        NB4 results directory; one subfolder per species is created.
    threads : int
    force : bool

    Returns
    -------
    dict
        ``{species: DataFrame}`` with columns
        ``Orthogroup, blast_class, pident, qcovhsp, scovhsp, best_target``.
    """
    results = {}
    diamond_bin = None  # lazy: only required if any species actually needs to run

    for sp in species_list:
        label = species_config.get(sp, {}).get('label', sp)
        sp_dir = Path(results_dir) / sp
        sp_dir.mkdir(parents=True, exist_ok=True)
        cache_fp = sp_dir / 'diamond_rare_classification.tsv'

        if not force and cache_fp.exists():
            df = pd.read_csv(cache_fp, sep='\t')
            results[sp] = df
            counts = df['blast_class'].value_counts() if 'blast_class' in df else {}
            print(f'{label}: loaded cached DIAMOND results ({len(df)} rare OGs) — '
                  + ', '.join(f'{k}: {v}' for k, v in counts.items()))
            continue

        if sp not in species_data or 'og_consensus' not in species_data[sp]:
            print(f'{label}: no og_consensus, skipping')
            results[sp] = pd.DataFrame()
            continue

        ogc = species_data[sp]['og_consensus']
        rare_ogs = ogc.loc[ogc['Pangenome_Class'] == 'Rare', 'Orthogroup'].tolist()
        nonrare_ogs = ogc.loc[ogc['Pangenome_Class'].isin(['Core', 'Accessory']),
                              'Orthogroup'].tolist()

        empty_cols = ['Orthogroup', 'blast_class', 'pident', 'qcovhsp',
                      'scovhsp', 'best_target']
        if not rare_ogs:
            df = pd.DataFrame(columns=empty_cols)
            df.to_csv(cache_fp, sep='\t', index=False)
            results[sp] = df
            print(f'{label}: 0 rare OGs — wrote empty cache')
            continue

        og_seq_root = Path(base_path) / sp / 'orthofinder_output'
        res_dirs = sorted(og_seq_root.glob('Results_*/Orthogroup_Sequences'))
        if not res_dirs:
            print(f'{label}: Orthogroup_Sequences not found under {og_seq_root}, skipping')
            results[sp] = pd.DataFrame(columns=empty_cols)
            continue
        og_seq_dir = str(res_dirs[-1])

        if diamond_bin is None:
            diamond_bin = _find_diamond_binary()
            print(f'Using DIAMOND at {diamond_bin}')

        print(f'\n{label}: {len(rare_ogs)} rare queries vs {len(nonrare_ogs)} core+accessory targets')
        rare_reps = get_representative_seqs(og_seq_dir, rare_ogs)
        nonrare_reps = get_representative_seqs(og_seq_dir, nonrare_ogs)
        print(f'  Representatives: {len(rare_reps)} rare, {len(nonrare_reps)} non-rare')

        if not rare_reps or not nonrare_reps:
            df = pd.DataFrame(columns=empty_cols)
            df.to_csv(cache_fp, sep='\t', index=False)
            results[sp] = df
            print(f'{label}: missing representative sequences, wrote empty cache')
            continue

        tmpdir = tempfile.mkdtemp(prefix=f'diamond_{sp}_')
        try:
            query_fa = f'{tmpdir}/rare_query.fa'
            db_fa = f'{tmpdir}/nonrare_db.fa'
            db_path = f'{tmpdir}/nonrare_db'
            out_tsv = f'{tmpdir}/diamond_out.tsv'

            SeqIO.write(list(rare_reps.values()), query_fa, 'fasta')
            SeqIO.write(list(nonrare_reps.values()), db_fa, 'fasta')

            subprocess.run([diamond_bin, 'makedb', '--in', db_fa, '-d', db_path],
                           capture_output=True, check=True)
            subprocess.run([
                diamond_bin, 'blastp',
                '-q', query_fa, '-d', db_path, '-o', out_tsv,
                '--sensitive', '--max-target-seqs', '5',
                '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
                'qlen', 'slen', 'qcovhsp', 'scovhsp', 'evalue', 'bitscore',
                '--threads', str(threads),
            ], capture_output=True, check=True)

            cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                    'qcovhsp', 'scovhsp', 'evalue', 'bitscore']
            if os.path.getsize(out_tsv) > 0:
                hits = pd.read_csv(out_tsv, sep='\t', names=cols)
                best_hits = (hits.sort_values('bitscore', ascending=False)
                                  .drop_duplicates('qseqid', keep='first'))
            else:
                best_hits = pd.DataFrame(columns=cols)

            classifications = []
            for og in rare_ogs:
                if og not in rare_reps:
                    continue
                row = best_hits[best_hits['qseqid'] == og]
                if row.empty:
                    classifications.append({
                        'Orthogroup': og, 'blast_class': 'No hit',
                        'pident': 0.0, 'qcovhsp': 0.0, 'scovhsp': 0.0,
                        'best_target': '',
                    })
                else:
                    r = row.iloc[0]
                    qcov, scov, pident = r['qcovhsp'], r['scovhsp'], r['pident']
                    if qcov >= 70 and scov >= 70 and pident >= 30:
                        bclass = 'Full-length homolog'
                    elif qcov >= 70 and scov < 50 and pident >= 30:
                        bclass = 'Fragment of longer gene'
                    else:
                        bclass = 'No significant similarity'
                    classifications.append({
                        'Orthogroup': og, 'blast_class': bclass,
                        'pident': float(pident), 'qcovhsp': float(qcov),
                        'scovhsp': float(scov), 'best_target': str(r['sseqid']),
                    })
            df = pd.DataFrame(classifications)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

        df.to_csv(cache_fp, sep='\t', index=False)
        results[sp] = df

        counts = df['blast_class'].value_counts() if not df.empty else {}
        for cls, cnt in counts.items():
            print(f'    {cls}: {cnt} ({100 * cnt / len(df):.1f}%)')
        print(f'  Cached: {cache_fp}')

    return results


def prepare_rare_characterisation_inputs(species_list, species_data, species_config,
                                         species_display, base_path, results_dir,
                                         min_rare=50, run_diamond=True,
                                         diamond_threads=8, force_diamond=False,
                                         annot_cols=('PFAMs', 'CAZy', 'KEGG_ko', 'GOs', 'Description')):
    """Build the four inputs required by :func:`plot_rare_characterisation_figure`.

    Wraps the setup that previously lived inline in NB4 so the notebook stays a
    one-line call. Returns a dict with keys:

    * ``species_deep`` — list of species with at least ``min_rare`` rare OGs
    * ``length_data`` — ``{species: DataFrame}`` of (Orthogroup, median_len, Pangenome_Class)
    * ``rates_df`` — long-format annotation-rate table across Core/Accessory/Rare
    * ``diamond_results`` — ``{species: DataFrame}`` from
      :func:`run_diamond_rare_vs_nonrare`. Empty dict if ``run_diamond=False``.

    Parameters
    ----------
    species_list : list of str
    species_data : dict
    species_config : dict
    species_display : dict
        ``{species_key: display_label}``.
    base_path : str
        Path under which ``{species}/orthofinder_output/Results_*`` lives.
    results_dir : str
    min_rare : int
        Threshold (in rare OG count) for inclusion in ``species_deep``.
    run_diamond : bool
        If True, run / load DIAMOND BLASTp results. If False, ``diamond_results``
        is returned as ``{}`` (Panel B will render as "No data" placeholders).
    diamond_threads : int
    force_diamond : bool
        Re-run DIAMOND even if cached results exist.
    annot_cols : iterable of str
        Annotation columns to compute Core/Accessory/Rare coverage for.
    """
    species_deep = [
        sp for sp in species_list
        if sp in species_data
        and 'og_consensus' in species_data[sp]
        and (species_data[sp]['og_consensus']['Pangenome_Class'] == 'Rare').sum() > min_rare
    ]
    print(f'Species with >{min_rare} rare OGs (deep characterisation): {species_deep}')

    # 1. Per-OG protein length, restricted to species_deep
    length_data = {}
    for sp in species_deep:
        ogc = species_data[sp]['og_consensus']
        og_seq_root = Path(base_path) / sp / 'orthofinder_output'
        res_dirs = sorted(og_seq_root.glob('Results_*/Orthogroup_Sequences'))
        if not res_dirs:
            print(f'  {sp}: no OG sequences found, skipping length analysis')
            continue
        og_seq_dir = str(res_dirs[-1])
        all_ogs = ogc['Orthogroup'].tolist()
        lens = og_protein_lengths(og_seq_dir, all_ogs)
        rows = []
        for _, r in ogc.iterrows():
            og = r['Orthogroup']
            if og in lens and lens[og]:
                rows.append({
                    'Orthogroup': og,
                    'median_len': float(np.median(lens[og])),
                    'Pangenome_Class': r['Pangenome_Class'],
                })
        length_data[sp] = pd.DataFrame(rows)
        print(f'  {sp}: {len(rows)} OGs with length data')

    # 2. Annotation rates per class for every species in species_list (full table)
    rates_rows = []
    for sp in species_list:
        if sp not in species_data or 'og_consensus' not in species_data[sp]:
            continue
        ogc = species_data[sp]['og_consensus']
        for col in annot_cols:
            if col not in ogc.columns:
                continue
            for cls, info in annotation_rate_by_class(ogc, col).items():
                rates_rows.append({
                    'Species': species_display.get(sp, sp),
                    'Annotation': col, 'Class': cls,
                    'Rate': float(info['rate']),
                    'n': int(info['n']),
                    'annotated': int(info['annotated']),
                })
    rates_df = pd.DataFrame(rates_rows)
    print(f'Annotation rates table: {rates_df.shape}')

    # 3. DIAMOND classifications (cached)
    if run_diamond:
        diamond_results = run_diamond_rare_vs_nonrare(
            species_deep, species_data, species_config, base_path, results_dir,
            threads=diamond_threads, force=force_diamond,
        )
    else:
        diamond_results = {}
        print('DIAMOND disabled (run_diamond=False)')

    return {
        'species_deep': species_deep,
        'length_data': length_data,
        'rates_df': rates_df,
        'diamond_results': diamond_results,
    }


_KINGDOM_PATTERNS = [
    ('Aspergillus (self)', ['aspergillus']),
    ('Other fungi',        ['penicillium', 'fusarium', 'candida', 'saccharomyces',
                            'cryptococcus', 'neurospora', 'magnaporthe', 'colletotrichum',
                            'botrytis', 'trichoderma', 'sphaerospora', 'sordariomycetes',
                            'pezizomycotina', 'rhizopus', 'phaeoacremonium', 'mucor',
                            'schizosaccharomyces', 'agaricomycetes', 'metarhizium',
                            'beauveria', 'mycena', 'fungus', 'fungi']),
    ('Bacteria',           ['bacterium', 'bacteroidota', 'escherichia', 'bacillus',
                            'pseudomonas', 'staphylococcus', 'streptomyces',
                            'mycobacterium', 'salmonella', 'lactobacillus',
                            'clostridium', 'rhizobium', 'sphingomonas',
                            'stenotrophomonas', 'firmicutes', 'proteobacteria',
                            'actinobacteria', 'cyanobacteria',
                            'thermodesulfobacteriota']),
    ('Protist',            ['perkinsus', 'plasmodium', 'toxoplasma', 'trypanosoma',
                            'leishmania', 'tetrahymena', 'paramecium', 'phytophthora',
                            'amoeba', 'dictyostelium', 'naegleria', 'giardia',
                            'protist', 'protozoa']),
    ('Plant',              ['arabidopsis', 'oryza sativa', 'plant', 'viridiplantae',
                            'algae']),
    ('Animal',             ['drosophila', 'mouse', 'human', 'homo sapiens',
                            'caenorhabditis', 'rhipicephalus', 'penaeus',
                            'metazoa', 'animal', 'mammal']),
    ('Virus',              ['virus', 'viridae']),
]


def classify_organism_kingdom(organism):
    """Heuristic kingdom assignment from a BLAST hit's organism string.

    Returns one of: ``Aspergillus (self)``, ``Other fungi``, ``Bacteria``,
    ``Protist``, ``Plant``, ``Animal``, ``Virus``, ``Unknown``.
    """
    if not isinstance(organism, str) or not organism.strip():
        return 'Unknown'
    o = organism.lower()
    for kingdom, keywords in _KINGDOM_PATTERNS:
        if any(k in o for k in keywords):
            return kingdom
    return 'Unknown'


def summarize_xenolog_results(blast_df, evalue_cutoff=1e-5):
    """Reduce BLAST hits to one row per OG, classify by kingdom, flag confident hits.

    Parameters
    ----------
    blast_df : pd.DataFrame
        Long-format BLAST output (one row per hit) with columns
        ``species, orthogroup, organism, pident, evalue``.
    evalue_cutoff : float
        Hits with e-value above this are flagged as ``Weak/no hit`` regardless
        of the organism.

    Returns
    -------
    pd.DataFrame
        One row per OG with columns species, orthogroup, organism, pident,
        evalue, kingdom, kingdom_filtered.
    """
    if blast_df is None or blast_df.empty:
        return pd.DataFrame()

    # Best hit per OG = lowest e-value (after coercing to numeric)
    bd = blast_df.copy()
    bd['evalue'] = pd.to_numeric(bd['evalue'], errors='coerce')
    bd = bd.dropna(subset=['evalue'])
    top = (bd.sort_values(['orthogroup', 'evalue'])
              .drop_duplicates('orthogroup', keep='first')
              .reset_index(drop=True))
    top['kingdom'] = top['organism'].apply(classify_organism_kingdom)
    top['kingdom_filtered'] = np.where(top['evalue'] > evalue_cutoff,
                                       'Weak / no hit', top['kingdom'])
    return top


def plot_xenolog_taxonomy(blast_df, species_config, results_dir=None,
                           evalue_cutoff=1e-5):
    """Visualise xenolog BLAST results: kingdom counts + per-OG hit-quality scatter.

    Saves
    -----
    {results_dir}/rare_xenolog_taxonomy.png / .pdf
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    top = summarize_xenolog_results(blast_df, evalue_cutoff=evalue_cutoff)
    if top.empty:
        print('No BLAST hits to plot.')
        return None

    kingdom_order = ['Aspergillus (self)', 'Other fungi', 'Bacteria',
                     'Protist', 'Plant', 'Animal', 'Virus',
                     'Unknown', 'Weak / no hit']
    kingdom_colors = {
        'Aspergillus (self)': '#7f7f7f',
        'Other fungi':        '#bcbd22',
        'Bacteria':           '#1f77b4',
        'Protist':            '#9467bd',
        'Plant':              '#2ca02c',
        'Animal':             '#e377c2',
        'Virus':              '#17becf',
        'Unknown':            '#cccccc',
        'Weak / no hit':      '#dddddd',
    }

    species_display = {sp: cfg.get('label', sp) for sp, cfg in species_config.items()}
    species_keys = sorted(top['species'].unique())

    fig = plt.figure(figsize=(15, 9))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.1], width_ratios=[1, 1.4],
                            hspace=0.4, wspace=0.3)

    # --- Panel A: stacked bars of kingdom counts per species (filtered by e-val) ---
    ax_a = fig.add_subplot(gs[0, 0])
    pivot = (top.groupby(['species', 'kingdom_filtered']).size()
                 .unstack(fill_value=0)
                 .reindex(index=species_keys, columns=kingdom_order, fill_value=0))
    bottom = np.zeros(len(species_keys))
    for k in kingdom_order:
        vals = pivot[k].values
        if vals.sum() == 0:
            continue
        ax_a.bar([species_display[s] for s in species_keys], vals, bottom=bottom,
                  color=kingdom_colors[k], label=k, edgecolor='white', linewidth=0.8)
        bottom += vals
    ax_a.set_ylabel('Number of rare OGs')
    ax_a.set_title(f'A. BLAST top-hit kingdom (e-value $\\leq$ {evalue_cutoff:g})',
                   fontweight='bold', loc='left')
    ax_a.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left',
                title='Kingdom of top hit')

    # --- Panel B: per-OG hit-quality scatter ---
    ax_b = fig.add_subplot(gs[0, 1])
    safe_e = top['evalue'].clip(lower=1e-200)
    top_plot = top.copy()
    top_plot['neg_log_e'] = -np.log10(safe_e)
    for k in kingdom_order:
        sub = top_plot[top_plot['kingdom_filtered'] == k]
        if sub.empty:
            continue
        is_strong = ~sub['kingdom_filtered'].isin(['Aspergillus (self)',
                                                     'Other fungi',
                                                     'Weak / no hit',
                                                     'Unknown'])
        ax_b.scatter(sub['pident'], sub['neg_log_e'], color=kingdom_colors[k],
                      s=80 if is_strong.any() else 35,
                      edgecolor='black' if is_strong.any() else 'none',
                      linewidth=0.7, alpha=0.85, label=k)
    ax_b.axhline(-np.log10(evalue_cutoff), color='red', linestyle='--', lw=1, alpha=0.6)
    ax_b.text(ax_b.get_xlim()[0] + 1, -np.log10(evalue_cutoff) + 1.5,
              f'e-value cutoff ({evalue_cutoff:g})', color='red', fontsize=8)
    ax_b.set_xlabel('% identity')
    ax_b.set_ylabel('-log10(e-value)')
    ax_b.set_title('B. Hit quality of best BLAST match per rare OG',
                   fontweight='bold', loc='left')
    ax_b.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left',
                title='Kingdom')

    # --- Panel C: confident non-Aspergillus / non-fungal hits, annotated ---
    ax_c = fig.add_subplot(gs[1, :])
    confident = top[
        ~top['kingdom_filtered'].isin(['Aspergillus (self)', 'Other fungi',
                                        'Weak / no hit', 'Unknown'])
    ].sort_values('evalue').copy()
    if confident.empty:
        ax_c.text(0.5, 0.5,
                  f'No confident non-fungal top hits (e $\\leq$ {evalue_cutoff:g})',
                  ha='center', va='center', fontsize=12, transform=ax_c.transAxes)
        ax_c.axis('off')
    else:
        confident['neg_log_e'] = -np.log10(confident['evalue'].clip(lower=1e-200))
        labels = confident.apply(
            lambda r: f"{r['orthogroup']}  →  {r['organism']}  "
                      f"({species_display[r['species']]})", axis=1)
        y_pos = np.arange(len(confident))
        bar_colors = [kingdom_colors[k] for k in confident['kingdom']]
        ax_c.barh(y_pos, confident['neg_log_e'], color=bar_colors,
                   edgecolor='black', linewidth=0.5)
        ax_c.set_yticks(y_pos)
        ax_c.set_yticklabels(labels, fontsize=8)
        ax_c.invert_yaxis()
        ax_c.set_xlabel('-log10(e-value)  (higher = more significant)')
        ax_c.set_title(f'C. Confident non-fungal top hits — candidate xenologs '
                       f'(e $\\leq$ {evalue_cutoff:g})',
                       fontweight='bold', loc='left')
        for i, r in enumerate(confident.itertuples()):
            ax_c.text(r.neg_log_e + 0.3, i,
                      f'pident={r.pident:.0f}%', va='center', fontsize=7)

    plt.suptitle('NCBI BLASTp results for GC-outlier rare OGs',
                  fontsize=13, fontweight='bold', y=0.995)

    if results_dir is not None:
        from pathlib import Path
        out = Path(results_dir)
        out.mkdir(parents=True, exist_ok=True)
        plt.savefig(out / 'rare_xenolog_taxonomy.png', dpi=180, bbox_inches='tight')
        plt.savefig(out / 'rare_xenolog_taxonomy.pdf', bbox_inches='tight')
        print(f'Saved: {out / "rare_xenolog_taxonomy.png"}')

    return top


DIAMOND_CLASS_INTERPRETATION = {
    'Full-length homolog': {
        'rule':     'qcov >= 70%, scov >= 70%, pident >= 30%',
        'meaning':  'Rare OG aligns end-to-end with a Core/Accessory protein',
        'biology':  'Likely paralog or in-genome duplicate that ended up rare due to copy-number variation. Not biologically novel.',
        'rare_genome_status': 'NOT a real rare gene — paralog of a common gene',
        'color':    'steelblue',
    },
    'Fragment of longer gene': {
        'rule':     'qcov >= 70%, scov < 50%, pident >= 30%',
        'meaning':  'Rare OG covers only a small portion of a longer Core/Accessory protein',
        'biology':  'Typically a pseudogene fragment, an assembly truncation, or a domain-only homolog of a longer protein. Often a degradation artefact, not a gain.',
        'rare_genome_status': 'NOT a real rare gene — partial / degraded fragment',
        'color':    'coral',
    },
    'No significant similarity': {
        'rule':     'has a hit but qcov, scov or pident below thresholds',
        'meaning':  'Faint similarity to a Core/Accessory gene but not classifiable',
        'biology':  'Highly diverged paralog, ambiguous case. May or may not be a real rare gene.',
        'rare_genome_status': 'AMBIGUOUS — possible diverged copy',
        'color':    'lightgrey',
    },
    'No hit': {
        'rule':     'no DIAMOND hit at all against Core+Accessory pool',
        'meaning':  'No homolog in the species\' Core or Accessory genome',
        'biology':  'Genuine novel/orphan rare gene. Candidates for strain-specific innovations, horizontally transferred genes, or lineage-specific evolution. THIS is the biologically interesting set.',
        'rare_genome_status': 'TRUE rare gene — novel / orphan / candidate xenolog',
        'color':    '#8D8585',
    },
}


def summarize_diamond_classifications(diamond_results, species_config,
                                       results_dir=None):
    """Build per-species count + percentage table of DIAMOND BLAST classifications.

    Parameters
    ----------
    diamond_results : dict
        ``{species: DataFrame}`` from :func:`run_diamond_rare_vs_nonrare`.
    species_config : dict
    results_dir : str or None
        If given, writes ``{results_dir}/diamond_classification_summary.csv``.

    Returns
    -------
    pd.DataFrame
        One row per species; columns are the four BLAST classes (counts) plus
        ``Total``, ``Pct_TrueRare``, ``Pct_Paralog``, ``Pct_Fragment``.
    """
    rows = []
    classes = list(DIAMOND_CLASS_INTERPRETATION.keys())
    for sp, df in diamond_results.items():
        label = species_config.get(sp, {}).get('label', sp)
        if df is None or df.empty or 'blast_class' not in df.columns:
            rows.append({'Species': label, **{c: 0 for c in classes},
                         'Total': 0, 'Pct_TrueRare': float('nan'),
                         'Pct_Paralog': float('nan'),
                         'Pct_Fragment': float('nan')})
            continue
        counts = df['blast_class'].value_counts()
        total = int(counts.sum())
        row = {'Species': label, 'Total': total}
        for c in classes:
            row[c] = int(counts.get(c, 0))
        row['Pct_TrueRare']  = round(100 * row['No hit'] / total, 1) if total else 0
        row['Pct_Paralog']   = round(100 * row['Full-length homolog'] / total, 1) if total else 0
        row['Pct_Fragment']  = round(100 * row['Fragment of longer gene'] / total, 1) if total else 0
        rows.append(row)
    summary = pd.DataFrame(rows)
    if results_dir is not None:
        out_fp = Path(results_dir) / 'diamond_classification_summary.csv'
        summary.to_csv(out_fp, index=False)
        print(f'Saved: {out_fp}')
    return summary


def plot_diamond_classification_breakdown(diamond_results, species_config,
                                           results_dir=None):
    """Detailed per-species DIAMOND breakdown with biological interpretation.

    Renders a stacked-bar chart of the four classes per species, plus an inset
    legend that explains what each class means in terms of rare-genome biology.

    Parameters
    ----------
    diamond_results : dict
        ``{species: DataFrame}`` from :func:`run_diamond_rare_vs_nonrare`.
    species_config : dict
    results_dir : str or None
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches

    classes = list(DIAMOND_CLASS_INTERPRETATION.keys())
    colors = {c: DIAMOND_CLASS_INTERPRETATION[c]['color'] for c in classes}

    summary = summarize_diamond_classifications(diamond_results, species_config,
                                                 results_dir=results_dir)
    if summary.empty:
        print('No DIAMOND results to plot.')
        return None

    import textwrap
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.1, 1.3], wspace=0.35,
                            left=0.06, right=0.97, top=0.92, bottom=0.10)

    # ---- Left panel: per-species stacked count bar with absolute and % labels ----
    ax_left = fig.add_subplot(gs[0, 0])
    species_labels = summary['Species'].tolist()
    x = np.arange(len(species_labels))
    bottom = np.zeros(len(species_labels))
    for c in classes:
        vals = summary[c].values
        ax_left.bar(x, vals, bottom=bottom,
                    color=colors[c], edgecolor='white', linewidth=0.6, label=c)
        for xi, (v, b, tot) in enumerate(zip(vals, bottom, summary['Total'].values)):
            if v > 0:
                pct = 100 * v / tot if tot else 0
                ax_left.text(xi, b + v / 2, f'{v}\n({pct:.0f}%)',
                             ha='center', va='center', fontsize=8, color='white'
                             if c in ('Full-length homolog', 'No hit') else 'black')
        bottom += vals
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(species_labels, fontsize=9)
    ax_left.set_ylabel('Number of rare OGs')
    ax_left.set_title('DIAMOND BLASTp classification of rare OGs',
                      fontweight='bold', loc='left')
    ax_left.legend(fontsize=8, loc='upper right', title='BLAST class', frameon=True)

    # ---- Right panel: biological interpretation card layout ----
    # Use a sub-gridspec with one row per class so each entry has its own
    # bounded box and lines can't bleed into the next entry.
    ax_right_title = fig.add_subplot(gs[0, 1])
    ax_right_title.axis('off')
    ax_right_title.set_title('What each class means biologically',
                              fontweight='bold', loc='left', pad=10)

    n_classes = len(classes)
    inner = gridspec.GridSpecFromSubplotSpec(
        n_classes, 1, subplot_spec=gs[0, 1], hspace=0.45,
    )
    for i, c in enumerate(classes):
        ax_card = fig.add_subplot(inner[i, 0])
        ax_card.set_xlim(0, 1)
        ax_card.set_ylim(0, 1)
        ax_card.axis('off')

        info = DIAMOND_CLASS_INTERPRETATION[c]

        # Left strip: colour swatch identifies the class
        ax_card.add_patch(mpatches.Rectangle((0.0, 0.05), 0.05, 0.9,
                                              color=info['color'], lw=0))
        # Title + rule (top of card)
        ax_card.text(0.07, 0.93, c, fontsize=11, fontweight='bold', va='top')
        ax_card.text(0.07, 0.74,
                     textwrap.fill(f"Rule: {info['rule']}", width=70),
                     fontsize=8, style='italic', va='top', color='dimgrey')

        # Biology body (wrapped for readability)
        ax_card.text(0.07, 0.55,
                     textwrap.fill(info['biology'], width=70),
                     fontsize=8.5, va='top')

        # Status line (bottom of card)
        ax_card.text(0.07, 0.13,
                     f"→ {info['rare_genome_status']}",
                     fontsize=9, fontweight='bold', color=info['color'], va='top')
        # Card border
        ax_card.add_patch(mpatches.Rectangle((0, 0), 1, 1, fill=False,
                                              edgecolor='lightgrey', lw=0.7))

    if results_dir is not None:
        out_fp = Path(results_dir) / 'diamond_classification_breakdown.png'
        plt.savefig(out_fp, dpi=180, bbox_inches='tight')
        plt.savefig(str(out_fp).replace('.png', '.pdf'), bbox_inches='tight')
        print(f'Saved: {out_fp}')
    return summary


def plot_rare_characterisation_figure(species_deep, length_data, diamond_results,
                                       burden_all, rates_df, species_config,
                                       results_dir=None):
    """Create the multi-panel rare genome characterisation figure.

    Parameters
    ----------
    species_deep : list of str
        Species with enough rare OGs for deep characterisation.
    length_data : dict
        ``{species: DataFrame}`` with median_len, Pangenome_Class columns.
    diamond_results : dict
        ``{species: DataFrame}`` with blast_class column.
    burden_all : pd.DataFrame
        Full genome burden table (from :func:`build_full_genome_burden_table`).
    rates_df : pd.DataFrame
        Annotation rates (from annotation_rate_by_class calls).
    species_config : dict
    results_dir : str or None
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import seaborn as sns

    class_order = ['Core', 'Accessory', 'Rare']
    class_colors = {'Core': 'forestgreen', 'Accessory': 'steelblue', 'Rare': 'coral'}

    fig = plt.figure(figsize=(16, 14), constrained_layout=False)
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           height_ratios=[1, 1, 1.2],
                           hspace=0.45, wspace=0.3,
                           left=0.07, right=0.97, top=0.95, bottom=0.06)

    # ---- Panel A: Gene length by pangenome class (violin) ----
    ax_a = fig.add_subplot(gs[0, 0])
    plot_species = [sp for sp in species_deep if sp in length_data]

    all_len_data = []
    for sp in plot_species:
        df = length_data[sp]
        df_copy = df[['median_len', 'Pangenome_Class']].copy()
        df_copy['Species'] = f'A. {sp}'
        all_len_data.append(df_copy)

    if all_len_data:
        all_len_df = pd.concat(all_len_data)
        upper = np.percentile(all_len_df['median_len'].dropna(), 99)
        all_len_df = all_len_df[all_len_df['median_len'] <= upper]
        sns.boxplot(data=all_len_df, x='Species', y='median_len', hue='Pangenome_Class',
                    hue_order=class_order, palette=class_colors, ax=ax_a,
                    fliersize=1, linewidth=0.8)
        ax_a.set_ylabel('Median protein length (aa)')
        ax_a.set_xlabel('')
    ax_a.set_title('A. Protein length by pangenome class', fontweight='bold', loc='left')
    ax_a.legend(fontsize=8, title='Class')

    # ---- Panel B: DIAMOND BLAST classification (pie charts) ----
    # Nested 2x2 grid inside gs[0, 1] keeps pies in-bounds and avoids
    # collision with Panel C below.
    blast_class_colors = {
        'Full-length homolog': 'steelblue',
        'Fragment of longer gene': 'coral',
        'No significant similarity': 'lightgrey',
        'No hit': '#8D8585',
    }
    blast_class_order = ['Full-length homolog', 'Fragment of longer gene',
                         'No significant similarity', 'No hit']

    # Reserve a strip on the right of the row for a single shared legend
    gs_b_outer = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs[0, 1], width_ratios=[3, 2], wspace=0.05,
    )
    gs_b = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=gs_b_outer[0, 0], hspace=0.45, wspace=0.15,
    )
    ax_b_title = fig.add_subplot(gs[0, 1])
    ax_b_title.axis('off')
    ax_b_title.set_title('B. DIAMOND BLAST classification of rare OGs',
                         fontweight='bold', loc='left', pad=12)

    for i, sp in enumerate(species_deep[:4]):
        sub_ax = fig.add_subplot(gs_b[i // 2, i % 2])
        sub_ax.set_aspect('equal')
        df = diamond_results.get(sp, pd.DataFrame())
        if df.empty or 'blast_class' not in df.columns:
            sub_ax.text(0.5, 0.5, f'A. {sp}\nNo data', ha='center', va='center',
                        fontsize=9, transform=sub_ax.transAxes)
            sub_ax.axis('off')
            continue
        counts = df['blast_class'].value_counts()
        vals = [counts.get(c, 0) for c in blast_class_order]
        colors = [blast_class_colors[c] for c in blast_class_order]
        nonzero = [(v, c, col) for v, c, col in zip(vals, blast_class_order, colors) if v > 0]
        if nonzero:
            sub_ax.pie([x[0] for x in nonzero], colors=[x[2] for x in nonzero],
                       startangle=90, wedgeprops=dict(linewidth=0.5, edgecolor='white'))
        sub_ax.set_title(f'A. {sp} (n={sum(vals)})', fontsize=8, fontweight='bold')

    # Shared legend for Panel B with one-line interpretations
    ax_b_legend = fig.add_subplot(gs_b_outer[0, 1])
    ax_b_legend.axis('off')
    legend_meanings = {
        'Full-length homolog':       'paralog / duplicate',
        'Fragment of longer gene':   'pseudogene / truncation',
        'No significant similarity': 'highly diverged',
        'No hit':                    'novel / orphan / HGT',
    }
    handles = [mpatches.Patch(color=blast_class_colors[c],
                               label=f'{c}\n  ({legend_meanings[c]})')
               for c in blast_class_order]
    ax_b_legend.legend(handles=handles, loc='center left', fontsize=7.5,
                       frameon=False, handlelength=1.3, handleheight=1.3,
                       borderpad=0.4, labelspacing=1.0)

    # ---- Panel C: Rare gene burden by phenotype (box) ----
    # Constrain to ~60% of row width so the boxplot doesn't stretch flat.
    gs_c = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[1, :], wspace=0)
    ax_c = fig.add_subplot(gs_c[0, 1:4])
    burden_plot = burden_all[
        burden_all['Species'].isin(['A. fumigatus', 'A. flavus', 'A. niger'])].copy()
    burden_plot['Species_Phenotype'] = (
        burden_plot['Species'] + '\n' +
        burden_plot['Phenotype'].str.replace('-', '\n', regex=False)
    )

    phenotype_colors = {
        'Human-pathogenic': '#e74c3c', 'Animal-pathogenic': '#f39c12',
        'Plant-pathogenic': '#f1c40f', 'Environmental': '#27ae60',
        'Industrial-trait': '#9b59b6', 'Lab': '#95a5a6', 'Unknown': '#bdc3c7',
    }
    sns.boxplot(data=burden_plot, x='Species', y='Rare_genes', hue='Phenotype',
                palette=phenotype_colors, ax=ax_c, fliersize=2, linewidth=0.8)
    ax_c.set_ylabel('Rare genes per genome')
    ax_c.set_xlabel('')
    ax_c.set_title('C. Rare gene burden by phenotype', fontweight='bold', loc='left')
    ax_c.legend(fontsize=8, title='Phenotype', bbox_to_anchor=(1.02, 1), loc='upper left')

    # ---- Panel D: Accessory vs Rare scatter (niger + oryzae) ----
    ax_d = fig.add_subplot(gs[2, 0])
    for sp, color, marker in [('niger', '#00A087', 'o'), ('oryzae', '#F39B7F', 's')]:
        sub = burden_all[burden_all['Species'] == species_config[sp]['label']]
        ax_d.scatter(sub['Accessory_genes'], sub['Rare_genes'], c=color,
                     marker=marker, s=60, alpha=0.7, edgecolors='black', linewidth=0.5,
                     label=species_config[sp]['label'])
    ax_d.set_xlabel('Accessory genes per genome')
    ax_d.set_ylabel('Rare genes per genome')
    ax_d.set_title('D. Accessory vs Rare gene count', fontweight='bold', loc='left')
    ax_d.legend(fontsize=9)

    # ---- Panel E: Annotation rate by class ----
    ax_e = fig.add_subplot(gs[2, 1])
    pfam_rates = rates_df[rates_df['Annotation'] == 'PFAMs'].copy()
    if not pfam_rates.empty:
        x_pos = np.arange(len(pfam_rates['Species'].unique()))
        species_names = pfam_rates['Species'].unique()
        width = 0.25
        for i, cls in enumerate(class_order):
            sub = pfam_rates[pfam_rates['Class'] == cls].set_index('Species').reindex(species_names)
            ax_e.bar(x_pos + i * width, sub['Rate'].values, width=width,
                     label=cls, color=class_colors[cls], alpha=0.8,
                     edgecolor='white', linewidth=0.5)
        ax_e.set_xticks(x_pos + width)
        ax_e.set_xticklabels(species_names, fontsize=9)
    ax_e.set_ylabel('Pfam annotation rate')
    ax_e.set_title('E. Pfam annotation rate by pangenome class', fontweight='bold', loc='left')
    ax_e.legend(fontsize=9)

    if results_dir is not None:
        plt.savefig(f'{results_dir}/rare_genome_characterisation_figure.png',
                    dpi=200, bbox_inches='tight')
        plt.savefig(f'{results_dir}/rare_genome_characterisation_figure.pdf',
                    bbox_inches='tight')
    print('Saved main characterisation figure')


def run_phylogenetic_analysis_all_species(species_list, species_data, species_config,
                                           pheno_map, pheno_df, results_dir,
                                           genus='Aspergillus',
                                           ani_excluded=None, transform='sqrt'):
    """Run phylogenetic signal, gain/loss, and ancestral niche reconstruction.

    Parameters
    ----------
    species_list : list of str
    species_data : dict
    species_config : dict
    pheno_map : dict
    pheno_df : pd.DataFrame
    results_dir : str
    genus : str
    ani_excluded : dict or None
    transform : str
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    if ani_excluded is None:
        ani_excluded = {}

    for sp in species_list:
        label = species_config[sp]['label']
        sp_results_dir = os.path.join(results_dir, sp)
        os.makedirs(sp_results_dir, exist_ok=True)

        print(f"\n{'='*80}")
        print(f'{label} -- Phylogenetic Analysis')
        print(f"{'='*80}")

        # Load IQ-TREE tree
        tree_candidates = [
            f'/datadrive/Species/{genus}/{sp}/iqtree_output_snv/gubbins_tree.contree',
            f'/datadrive/Species/{genus}/{sp}/iqtree_output_snv/gubbins_tree.treefile',
        ]
        tree_file = None
        for tc in tree_candidates:
            if os.path.exists(tc):
                tree_file = tc
                break

        if tree_file is None:
            print(f'  [SKIP] No IQ-TREE tree file found for {sp}')
            continue

        print(f'  Loading tree from: {tree_file}')
        tree = Phylo.read(tree_file, 'newick')

        # Clean tip labels
        ref_tips = [t for t in tree.get_terminals()
                    if t.name and t.name.endswith('.ref')]
        for tip in ref_tips:
            tree.prune(tip)
        tree = clean_tip_labels(tree)

        # Prune ANI-excluded tips
        if sp in ani_excluded:
            for exc_acc in ani_excluded[sp]:
                exc_base = re.match(r'(GC[AF]_\d+)', exc_acc).group(1)
                for tip in list(tree.get_terminals()):
                    if exc_base in tip.name:
                        tree.prune(tip)
                        print(f'  Pruned ANI-excluded: {tip.name}')

        terminals = tree.get_terminals()
        n_tips = len(terminals)
        print(f'  Tree: {n_tips} tips')

        # ---- Phylogenetic signal ----
        print(f'\n  --- Fritz & Purvis D (phylogenetic signal) ---')

        tip_phenotypes = {}
        for t in terminals:
            tip_phenotypes[t.name] = pheno_map.get(t.name, 'Unknown')

        present_phenotypes = set(tip_phenotypes.values()) - {'Unknown'}
        print(f'  Phenotypes in tree: {present_phenotypes}')

        d_results = []
        for pheno in sorted(present_phenotypes):
            binary_states = {t: (1 if tip_phenotypes[t] == pheno else 0)
                             for t in tip_phenotypes}
            n1 = sum(binary_states.values())
            if n1 < 3 or n1 > n_tips - 3:
                print(f'    {pheno}: n={n1}, too few for D test')
                continue

            d_stat = fritz_purvis_d(tree, binary_states)
            d_results.append({'Phenotype': pheno, **d_stat, 'n': n1})
            sig_str = '*' if d_stat['p_random'] < 0.05 else 'ns'
            print(f"    {pheno}: D={d_stat['D']:.3f}, "
                  f"p(random)={d_stat['p_random']:.3f} {sig_str}")

        # Plot phylogenetic signal
        if d_results:
            meta_df = pheno_df[['Assembly Accession', 'Phenotype']].copy()
            fig, (ax1, ax2) = plt.subplots(1, 2,
                                            figsize=(14, max(5, n_tips * 0.15)))

            plot_phylogram(tree, ax1, df_meta=meta_df,
                           title=f'{label} Phylogram', transform=transform)

            d_df = pd.DataFrame(d_results)
            colors = [get_color_for_phenotype(p) for p in d_df['Phenotype']]
            ax2.barh(d_df['Phenotype'], d_df['D'], color=colors,
                     edgecolor='black', linewidth=0.5)
            ax2.axvline(x=0, color='green', linestyle='--', alpha=0.7,
                        label='D=0 (Brownian)')
            ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.7,
                        label='D=1 (Random)')
            ax2.set_xlabel('Fritz & Purvis D')
            ax2.set_title('Phylogenetic Signal per Phenotype')
            ax2.legend(fontsize=8)

            plt.tight_layout()
            output_path = os.path.join(sp_results_dir, 'phylogenetic_signal.png')
            plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
            print(f'  Saved: {output_path}')

        # ---- Gene gain/loss ----
        print(f'\n  --- Gene gain/loss events ---')
        pav = species_data[sp]['pav']
        ogc = species_data[sp]['og_consensus']

        variable_ogs = ogc.loc[
            ogc['Pangenome_Class'].isin(['Accessory', 'Rare']), 'Orthogroup'
        ].tolist()
        events_df = fitch_gain_loss(tree, pav, variable_ogs, pheno_map)

        if not events_df.empty:
            events_path = os.path.join(sp_results_dir, 'gain_loss_events.tsv')
            events_df.to_csv(events_path, sep='\t', index=False)
            print(f'  Gain/loss events: {len(events_df)} OGs analysed')
            print(f'  Mean parsimony changes: {events_df["parsimony_changes"].mean():.2f}')
            print(f'  Saved: {events_path}')

            summary = events_df.describe()
            summary_path = os.path.join(sp_results_dir, 'gain_loss_summary.tsv')
            summary.to_csv(summary_path, sep='\t')
            print(f'  Saved: {summary_path}')

        # ---- Ancestral niche reconstruction ----
        print(f'\n  --- Ancestral niche reconstruction ---')

        state_sets = {}

        def fitch_reconstruct(clade):
            if clade.is_terminal():
                pheno = tip_phenotypes.get(clade.name, 'Unknown')
                state_sets[id(clade)] = {pheno}
                return {pheno}
            child_sets = [fitch_reconstruct(c) for c in clade.clades]
            inter = child_sets[0]
            for cs in child_sets[1:]:
                inter = inter & cs
            if inter:
                state_sets[id(clade)] = inter
            else:
                state_sets[id(clade)] = set().union(*child_sets)
            return state_sets[id(clade)]

        fitch_reconstruct(tree.root)

        # Plot ancestral reconstruction
        fig_height = max(6, n_tips * 0.18)
        fig, ax = plt.subplots(figsize=(12, fig_height))

        phenotype_colors_local = {
            'Human-pathogenic': '#e74c3c', 'Animal-pathogenic': '#f39c12',
            'Plant-pathogenic': '#f1c40f', 'Environmental': '#27ae60',
            'Industrial-trait': '#9b59b6', 'Lab': '#95a5a6', 'Unknown': '#bdc3c7',
        }

        meta_df = pheno_df[['Assembly Accession', 'Phenotype']].copy()
        plot_phylogram(tree, ax, df_meta=meta_df,
                       title=f'{label} -- Ancestral Niche Reconstruction',
                       transform=transform)

        depths = tree.depths()
        depths = apply_transform(depths, transform)

        for clade in tree.find_clades():
            if not clade.is_terminal() and id(clade) in state_sets:
                states = state_sets[id(clade)]
                if len(states) == 1:
                    state = list(states)[0]
                    color = phenotype_colors_local.get(state, '#bdc3c7')
                    x = depths.get(clade, 0)
                    children_y = []
                    for t in clade.get_terminals():
                        for i_idx, term in enumerate(tree.get_terminals()):
                            if term.name == t.name:
                                children_y.append(i_idx)
                                break
                    if children_y:
                        y = np.mean(children_y)
                        ax.plot(x, y, 's', ms=6, color=color, alpha=0.6,
                                markeredgecolor='black', markeredgewidth=0.3)

        legend_handles = [mpatches.Patch(color=c, label=p)
                          for p, c in phenotype_colors_local.items()
                          if p in set(tip_phenotypes.values())]
        ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
                  title='Phenotype')

        plt.tight_layout()
        output_path = os.path.join(sp_results_dir,
                                   'ancestral_niche_reconstruction.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {output_path}')


def run_mash_per_species(species_list_mash, species_display_mash,
                          genus, ani_excluded, results_dir):
    """Compute per-species Mash distance matrices.

    Parameters
    ----------
    species_list_mash : list of str
    species_display_mash : dict
    genus : str
    ani_excluded : dict
    results_dir : str

    Returns
    -------
    dict
        ``{species: pd.DataFrame}`` distance matrices.
    """
    from pathlib import Path
    mash_output_base = os.path.join(results_dir, 'mash_clustering')
    os.makedirs(mash_output_base, exist_ok=True)

    mash_results = {}

    for species in species_list_mash:
        print(f"\n{'='*60}")
        print(f'  {species_display_mash[species]} - Mash Distance Computation')
        print(f"{'='*60}")

        genome_dir = f'/datadrive/Species/{genus}/{species}/filtered_genome'

        _exc_bases = set()
        if species in ani_excluded:
            for exc_acc in ani_excluded[species]:
                _exc_bases.add(re.match(r'(GC[AF]_\d+)', exc_acc).group(1))

        genome_files = sorted([
            str(p) for p in Path(genome_dir).glob('*.fna')
            if '_renamed' not in p.name
            and not any(exc in p.name for exc in _exc_bases)
        ])
        print(f'  Found {len(genome_files)} genome files')

        if len(genome_files) < 3:
            print(f'  Too few genomes for clustering, skipping')
            continue

        sp_output_dir = os.path.join(mash_output_base, species)
        dist_matrix = run_mash_distances(genome_files, sp_output_dir, species)

        if dist_matrix is not None:
            mash_results[species] = dist_matrix
            print(f'  Distance matrix shape: {dist_matrix.shape}')

            dm_path = os.path.join(sp_output_dir,
                                   f'{species}_mash_distance_matrix.tsv')
            dist_matrix.to_csv(dm_path, sep='\t')
        else:
            print(f'  FAILED - skipping {species}')

    print(f"\n{'='*60}")
    print('Per-species Mash distance computation complete!')
    for sp_name, dm in mash_results.items():
        print(f'  {species_display_mash[sp_name]}: {dm.shape[0]} genomes')

    return mash_results


def run_mash_kmeans(species_list_mash, species_display_mash,
                     mash_results, acc_to_phenotype, results_dir):
    """K-means clustering on per-species Mash distance matrices.

    Parameters
    ----------
    species_list_mash : list of str
    species_display_mash : dict
    mash_results : dict
        ``{species: pd.DataFrame}`` distance matrices.
    acc_to_phenotype : dict
    results_dir : str
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    mash_output_base = os.path.join(results_dir, 'mash_clustering')

    for species in species_list_mash:
        if species not in mash_results:
            continue

        dist_mat = mash_results[species]
        sp_dir = os.path.join(mash_output_base, species)
        os.makedirs(sp_dir, exist_ok=True)

        print(f"\n{'='*60}")
        print(f'  {species_display_mash[species]} - Mash Clustering')
        print(f"{'='*60}")

        n_samples = dist_mat.shape[0]
        max_k = min(n_samples - 1, 15)
        k_range = list(range(2, max_k + 1))

        if len(k_range) < 2:
            print(f'  Too few samples (n={n_samples}) for clustering')
            continue

        current_accs = list(dist_mat.index)
        dist_array = dist_mat.values

        inertias = []
        sil_scores = []

        for k in k_range:
            km = KMeans(n_clusters=k, random_state=42, n_init=10).fit(dist_array)
            inertias.append(km.inertia_)
            if k < n_samples:
                sil = (silhouette_score(dist_array, km.labels_, metric='precomputed')
                       if np.allclose(dist_array, dist_array.T)
                       else silhouette_score(dist_array, km.labels_))
                sil_scores.append(sil)
            else:
                sil_scores.append(0)

        optimal_k, elbow_k, best_sil_k, reason = choose_optimal_k(
            k_range, inertias, sil_scores, n_samples)
        print(f'  Optimal k={optimal_k} ({reason})')

        km_final = KMeans(n_clusters=optimal_k, random_state=42, n_init=10).fit(dist_array)
        labels = km_final.labels_

        cluster_df = pd.DataFrame({
            'accession': current_accs,
            'mash_cluster': [f'M{l+1}' for l in labels],
            'phenotype': [acc_to_phenotype.get(a, 'Unknown') for a in current_accs],
        })
        cluster_df.to_csv(os.path.join(sp_dir, f'{species}_primary_clusters.csv'),
                          index=False)

        print(f'  Cluster sizes:')
        for cl, count in cluster_df['mash_cluster'].value_counts().sort_index().items():
            print(f'    {cl}: {count}')

        # Heatmap
        cluster_names = sorted(cluster_df['mash_cluster'].unique())
        cluster_palette = dict(zip(cluster_names,
                                   sns.color_palette('Set2', len(cluster_names))))
        pheno_palette = PHENOTYPE_COLORS

        row_colors = pd.DataFrame({
            'Cluster': cluster_df.set_index('accession')['mash_cluster'].map(
                cluster_palette),
            'Phenotype': cluster_df.set_index('accession')['phenotype'].map(
                pheno_palette),
        })

        g = sns.clustermap(
            dist_mat, method='ward', row_colors=row_colors, col_colors=row_colors,
            cmap='viridis_r', figsize=(10, 10), linewidths=0,
            xticklabels=False, yticklabels=False
        )
        g.fig.suptitle(
            f'{species_display_mash[species]} Mash Distances (k={optimal_k} clusters)',
            y=1.02, fontsize=13, fontweight='bold')
        g.fig.savefig(os.path.join(sp_dir, f'{species}_mash_heatmap.pdf'),
                      bbox_inches='tight', dpi=150)
        print(f'  Saved heatmap: {sp_dir}/{species}_mash_heatmap.pdf')


def run_mash_all_species(species_list_mash, species_display_mash,
                          genus, ani_excluded, acc_to_phenotype, results_dir):
    """All-species combined Mash clustering.

    Parameters
    ----------
    species_list_mash : list of str
    species_display_mash : dict
    genus : str
    ani_excluded : dict
    acc_to_phenotype : dict
    results_dir : str

    Returns
    -------
    pd.DataFrame or None
        Combined distance matrix, or None on failure.
    """
    from pathlib import Path
    from matplotlib.patches import Patch
    import matplotlib.pyplot as plt
    import seaborn as sns

    mash_output_base = os.path.join(results_dir, 'mash_clustering')
    all_sp_dir = os.path.join(mash_output_base, 'all_species')
    os.makedirs(all_sp_dir, exist_ok=True)

    all_genome_files = []
    genome_to_species = {}

    for species in species_list_mash:
        genome_dir = f'/datadrive/Species/{genus}/{species}/filtered_genome'

        _exc_bases = set()
        if species in ani_excluded:
            for exc_acc in ani_excluded[species]:
                _exc_bases.add(re.match(r'(GC[AF]_\d+)', exc_acc).group(1))

        genome_files = sorted([
            str(p) for p in Path(genome_dir).glob('*.fna')
            if '_renamed' not in p.name
            and not any(exc in p.name for exc in _exc_bases)
        ])

        for gf in genome_files:
            acc_match = re.search(r'(GC[AF]_\d+\.\d+)', gf)
            if acc_match:
                genome_to_species[acc_match.group(1)] = species_display_mash[species]
        all_genome_files.extend(genome_files)

    print(f'Total genomes across all species: {len(all_genome_files)}')

    all_dist_matrix = run_mash_distances(all_genome_files, all_sp_dir, 'all_species')

    if all_dist_matrix is not None:
        print(f'Combined distance matrix: {all_dist_matrix.shape}')

        all_dist_matrix.to_csv(
            os.path.join(all_sp_dir, 'all_species_mash_distance_matrix.csv'))

        # Species side-bar colours = the project-wide species palette
        # (funpan_utils.SPECIES_COLORS_DISPLAY).
        sp_color_map = {
            'A. fumigatus': '#e74c3c', 'A. flavus': '#ffd500',
            'A. niger': '#1a1aff', 'A. oryzae': '#27ae60',
        }
        row_species = pd.Series(
            {acc: genome_to_species.get(acc, 'Unknown')
             for acc in all_dist_matrix.index},
            name='Species'
        )
        row_colors_all = row_species.map(sp_color_map).fillna('#CCCCCC')

        g = sns.clustermap(
            all_dist_matrix, method='ward',
            row_colors=row_colors_all, col_colors=row_colors_all,
            cmap='viridis_r', figsize=(12, 12), linewidths=0,
            xticklabels=False, yticklabels=False
        )
        g.ax_cbar.tick_params(labelsize=15)
        g.ax_cbar.set_ylabel('Mash distance', fontsize=16)
        g.fig.suptitle(
            f'All Aspergillus Species Mash Distances '
            f'(n={all_dist_matrix.shape[0]} genomes)',
            y=1.02, fontsize=20, fontweight='bold')

        legend_patches = [Patch(color=c, label=s) for s, c in sp_color_map.items()]
        g.fig.legend(handles=legend_patches, loc='lower right',
                     bbox_to_anchor=(0.95, 0.05), fontsize=15,
                     title='Species', title_fontsize=16)

        g.fig.savefig(os.path.join(all_sp_dir, 'all_species_mash_heatmap.pdf'),
                      bbox_inches='tight', dpi=150)
        print(f'Saved: {all_sp_dir}/all_species_mash_heatmap.pdf')

        species_df = pd.DataFrame({
            'accession': all_dist_matrix.index,
            'species': [genome_to_species.get(a, 'Unknown')
                        for a in all_dist_matrix.index],
            'phenotype': [acc_to_phenotype.get(a, 'Unknown')
                          for a in all_dist_matrix.index],
        })
        species_df.to_csv(
            os.path.join(all_sp_dir, 'all_species_genome_metadata.csv'), index=False)
        print(f'Saved: {all_sp_dir}/all_species_genome_metadata.csv')
    else:
        print('Combined Mash analysis failed')

    return all_dist_matrix


def plot_rare_ogs_characterization(species_list, species_data, species_config,
                                    results_dir=None, top_cog_n=8):
    """Multi-panel visualization characterising rare orthogroups across species.

    Per species, panels show:
      (a) annotation coverage rates (Description / Pfam / CAZy / KEGG / GO / EC)
      (b) top COG functional categories
      (c) protein-count distribution (how many protein members per rare OG)

    Parameters
    ----------
    species_list : list of str
    species_data : dict
        ``{sp: {'og_consensus': DataFrame}}``.
    species_config : dict
        ``{sp: {'label': str}}``.
    results_dir : str or None
        If provided, save figure to this directory.
    top_cog_n : int
        Number of top COG categories to show per species.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from collections import Counter

    n_sp = len(species_list)
    fig, axes = plt.subplots(n_sp, 3, figsize=(15, 3.5 * n_sp))
    if n_sp == 1:
        axes = axes.reshape(1, -1)

    annot_cols = ['Description', 'PFAMs', 'CAZy', 'KEGG_ko', 'GOs', 'EC']
    # Fixed COG -> colour mapping so the same category always gets the same colour
    # across species (uses the canonical 26-letter ordering from COG_DESCRIPTIONS).
    # tab20 + tab20b gives 40 distinct hues — enough for the 26 COG categories
    # without the tab20/tab10 collision (e.g. K and P both landing on orange).
    _cog_keys = list(COG_DESCRIPTIONS.keys())
    _palette = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors)
    COG_COLORS = {c: _palette[i] for i, c in enumerate(_cog_keys)}
    # Semantically meaningful swaps so the most common categories get
    # intuitive colours (and all 26 stay unique):
    #   Q (secondary metabolites) -> bright red   (swap with Y)
    #   G (carbohydrate metabolism) -> bright green (swap with B)
    #   E (amino acid metabolism)   -> strong blue  (swap with J)
    COG_COLORS['Y'], COG_COLORS['Q'] = COG_COLORS['Q'], COG_COLORS['Y']
    COG_COLORS['B'], COG_COLORS['G'] = COG_COLORS['G'], COG_COLORS['B']
    COG_COLORS['J'], COG_COLORS['E'] = COG_COLORS['E'], COG_COLORS['J']

    for i, sp in enumerate(species_list):
        if sp not in species_data or 'og_consensus' not in species_data[sp]:
            for ax in axes[i]:
                ax.text(0.5, 0.5, f'No data for {sp}', ha='center', va='center')
                ax.axis('off')
            continue

        ogc = species_data[sp]['og_consensus']
        rare = ogc[ogc['Pangenome_Class'] == 'Rare'].copy()
        label = species_config.get(sp, {}).get('label', sp)
        n_rare = len(rare)

        # --- Panel A: annotation coverage rates -----------------------------
        ax = axes[i, 0]
        if n_rare == 0:
            ax.text(0.5, 0.5, f'{label}\nNo rare OGs', ha='center', va='center',
                    transform=ax.transAxes)
            ax.axis('off')
        else:
            rates = []
            for col in annot_cols:
                if col in rare.columns:
                    has_annot = rare[col].apply(
                        lambda x: bool(x) and str(x).strip() not in ('', '-', 'nan')
                    ).mean()
                    rates.append((col, has_annot * 100))
            if rates:
                names, vals = zip(*rates)
                bars = ax.barh(range(len(names)), vals, color='steelblue', alpha=0.8)
                ax.set_yticks(range(len(names)))
                ax.set_yticklabels(names, fontsize=9)
                ax.set_xlim(0, 100)
                ax.set_xlabel('% of rare OGs with annotation')
                ax.set_title(f'{label} — annotation coverage (n={n_rare})',
                             fontsize=10, fontweight='bold')
                for bar, v in zip(bars, vals):
                    ax.text(v + 1, bar.get_y() + bar.get_height() / 2,
                            f'{v:.0f}%', va='center', fontsize=8)
                ax.grid(axis='x', linestyle=':', alpha=0.4)

        # --- Panel B: top COG categories ------------------------------------
        ax = axes[i, 1]
        if n_rare == 0 or 'COG_category' not in rare.columns:
            ax.text(0.5, 0.5, 'No COG data', ha='center', va='center',
                    transform=ax.transAxes)
            ax.axis('off')
        else:
            cogs = rare['COG_category'].dropna().astype(str)
            # Each entry can be multi-letter like "JKL"; expand
            counter = Counter()
            for s in cogs:
                for ch in s:
                    if ch.isalpha() and ch != 'S':  # exclude 'S' (function unknown)
                        counter[ch] += 1
            top = counter.most_common(top_cog_n)
            if top:
                cogs_lbl = [c for c, _ in top]
                cogs_val = [v for _, v in top]
                cogs_full = [f"{c}: {COG_DESCRIPTIONS.get(c, '?')[:30]}"
                             for c in cogs_lbl]
                # Consistent colour per COG category across all species
                colors = [COG_COLORS.get(c, '#888888') for c in cogs_lbl]
                bars = ax.barh(range(len(cogs_lbl)), cogs_val, color=colors)
                ax.set_yticks(range(len(cogs_lbl)))
                ax.set_yticklabels(cogs_full, fontsize=8)
                ax.invert_yaxis()
                ax.set_xlabel('Number of rare OGs')
                ax.set_title(f'{label} — top COG categories',
                             fontsize=10, fontweight='bold')
                for bar, v in zip(bars, cogs_val):
                    ax.text(v + 0.1, bar.get_y() + bar.get_height() / 2,
                            f'{v}', va='center', fontsize=8)
                ax.grid(axis='x', linestyle=':', alpha=0.4)
            else:
                ax.text(0.5, 0.5, 'No COG hits', ha='center', va='center',
                        transform=ax.transAxes)
                ax.axis('off')

        # --- Panel C: protein-count distribution (how big each rare OG is) ---
        ax = axes[i, 2]
        if n_rare == 0 or 'n_proteins' not in rare.columns:
            ax.text(0.5, 0.5, 'No size data', ha='center', va='center',
                    transform=ax.transAxes)
            ax.axis('off')
        else:
            sizes = rare['n_proteins'].dropna().astype(int)
            if len(sizes):
                ax.hist(sizes, bins=min(30, max(5, sizes.max())),
                        color='coral', edgecolor='white', alpha=0.85)
                ax.set_xlabel('Proteins per rare OG')
                ax.set_ylabel('Count')
                ax.set_title(f'{label} — rare OG size distribution\n'
                             f'(median = {int(sizes.median())}, max = {int(sizes.max())})',
                             fontsize=10, fontweight='bold')
                ax.grid(axis='y', linestyle=':', alpha=0.4)

    plt.suptitle('Rare orthogroup characterisation per species',
                 fontsize=13, fontweight='bold', y=1.001)
    plt.tight_layout()

    if results_dir is not None:
        out_path = os.path.join(results_dir, 'rare_ogs_characterization.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f'Saved: {out_path}')
