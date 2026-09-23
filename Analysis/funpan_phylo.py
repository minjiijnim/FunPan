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


from Bio import Phylo
import re
import matplotlib.patches as mpatches

# =============================================================================
# CONFIGURATION - Change these to adjust visualization
# =============================================================================
# TRANSFORM = 'sqrt'  # Options: 'none', 'log', 'sqrt', or a float for power (e.g., 0.3)
# =============================================================================


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
        p = os.path.expanduser(f'~/anaconda3/envs/{env}/bin/mash')
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return 'mash'  # let subprocess raise a useful error if it's still missing

MASH_BIN = _find_mash_binary()
KMER_SIZE = 21
SKETCH_SIZE = 1000


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
    
    # Case 6: Real disagreement: best silhouette near elbow neighborhood
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


import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')


def heaps_law(n, kappa, gamma):
    return kappa * (n ** gamma)


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


import pandas as pd
import numpy as np
from scipy import stats
import os, re as _re


from pathlib import Path
from Bio import SeqIO
import warnings, os, sys, re
warnings.filterwarnings('ignore')


from scipy.stats import mannwhitneyu


import subprocess, tempfile, shutil


from matplotlib.patches import Patch


import requests, time, xml.etree.ElementTree as ET

NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'


from matplotlib.patches import Patch


def get_accession(col):
    """Extract GCA/GCF accession from OrthoFinder column name."""
    m = re.match(r'(GC[AF]_\d+\.\d+)', col)
    return m.group(1) if m else col


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


def og_protein_lengths(og_seq_dir, og_list):
    """Return dict {OG: [len1, len2, ...]} for requested OGs."""
    lengths = {}
    for og in og_list:
        fa = Path(og_seq_dir) / f'{og}.fa'
        if fa.exists():
            lengths[og] = [len(rec.seq) for rec in SeqIO.parse(fa, 'fasta')]
    return lengths


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


def compute_gc_content(sequence):
    """Compute GC content of a protein-coding DNA or protein sequence."""
    seq = str(sequence).upper()
    gc = seq.count('G') + seq.count('C')
    total = len(seq)
    return 100.0 * gc / total if total > 0 else 0


# xenolog_candidates = {}
GC_THRESHOLD = 2.0  # standard deviations from mean


NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'


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
        p = os.path.expanduser(f'~/anaconda3/envs/{env}/bin/diamond')
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    raise FileNotFoundError(
        'DIAMOND binary not found. Install via `conda install -c bioconda diamond` '
        'or set the $DIAMOND environment variable to the executable path.'
    )


def find_gc_outlier_rare_ogs(species_list, species_data, base_path,
                              gc_threshold=2.0, results_dir=None, force=False):
    """Find rare OGs with anomalous GC content -- candidate xenologs.

    For each species:
      1. Compute G+C frequency of each OG's representative protein sequence
         (longest member). Note this is *amino-acid* G+C, a coarse proxy for
         codon GC and therefore for HGT signal.
      2. Compute the species-wide mean and SD across *all* OGs.
      3. Return rare OGs whose deviation exceeds ``gc_threshold`` SDs.

    Per-OG GC content is cached in
    ``{results_dir}/{species}/og_gc_content.tsv`` when ``results_dir`` is
    given. Cached values are reused and only orthogroups missing from the
    cache are read from the sequence files. ``force=True`` recomputes every
    orthogroup and rewrites the cache. The threshold is applied after loading,
    so changing ``gc_threshold`` does not need a recompute.

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

        all_ogs = ogc['Orthogroup'].astype(str).tolist()

        cache_fp = (Path(results_dir) / sp / 'og_gc_content.tsv'
                    if results_dir is not None else None)
        gc_data = {}
        if cache_fp is not None and cache_fp.exists() and not force:
            cached = pd.read_csv(cache_fp, sep='\t')
            gc_data = dict(zip(cached['Orthogroup'].astype(str),
                               cached['GC_content'].astype(float)))
        missing = [og for og in all_ogs if og not in gc_data]

        if missing:
            og_seq_root = Path(base_path) / sp / 'orthofinder_output'
            res_dirs = sorted(og_seq_root.glob('Results_*/Orthogroup_Sequences'))
            if not res_dirs:
                print(f'  {sp}: Orthogroup_Sequences not found, skipping GC filter')
                candidates[sp] = pd.DataFrame()
                continue
            og_seq_dir = str(res_dirs[-1])
            print(f'  {sp}: reading GC content for {len(missing)} orthogroups '
                  f'({len(gc_data)} from cache)')
            reps = get_representative_seqs(og_seq_dir, missing)
            gc_data.update({og: compute_gc_content(reps[og].seq) for og in reps})
            if cache_fp is not None:
                cache_fp.parent.mkdir(parents=True, exist_ok=True)
                (pd.DataFrame({'Orthogroup': list(gc_data),
                               'GC_content': list(gc_data.values())})
                 .to_csv(cache_fp, sep='\t', index=False))
                print(f'  {sp}: cached GC content -> {cache_fp}')
        else:
            print(f'  {sp}: GC content loaded from {cache_fp}')

        gc_series = pd.Series({og: gc_data[og] for og in all_ogs if og in gc_data},
                              dtype=float)
        gc_mean = float(gc_series.mean())
        gc_std = float(gc_series.std())

        rare_ogc['GC_content'] = rare_ogc['Orthogroup'].astype(str).map(gc_data)
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
            results_dir=results_dir,
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
            print(f'{label}: loaded cached DIAMOND results ({len(df)} rare OGs), '
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
            print(f'{label}: 0 rare OGs, wrote empty cache')
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

    Builds the inputs the NB4 characterisation figures need, so the notebook stays a
    one-line call. Returns a dict with keys:

    * ``species_deep``: list of species with at least ``min_rare`` rare OGs
    * ``length_data``: ``{species: DataFrame}`` of (Orthogroup, median_len, Pangenome_Class)
    * ``rates_df``: long-format annotation-rate table across Core/Accessory/Rare
    * ``diamond_results``: ``{species: DataFrame}`` from
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
    {results_dir}/rare_xenolog_taxonomy.png
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
        ax_c.set_title(f'C. Confident non-fungal top hits: candidate xenologs '
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
        print(f'Saved: {out / "rare_xenolog_taxonomy.png"}')

    return top


DIAMOND_CLASS_INTERPRETATION = {
    'Full-length homolog': {
        'rule':     'qcov >= 70%, scov >= 70%, pident >= 30%',
        'meaning':  'Rare OG aligns end-to-end with a Core/Accessory protein',
        'biology':  'Likely paralog or in-genome duplicate that ended up rare due to copy-number variation. Not biologically novel.',
        'rare_genome_status': 'NOT a real rare gene: paralog of a common gene',
        'color':    'steelblue',
    },
    'Fragment of longer gene': {
        'rule':     'qcov >= 70%, scov < 50%, pident >= 30%',
        'meaning':  'Rare OG covers only a small portion of a longer Core/Accessory protein',
        'biology':  'Typically a pseudogene fragment, an assembly truncation, or a domain-only homolog of a longer protein. Often a degradation artefact, not a gain.',
        'rare_genome_status': 'NOT a real rare gene: partial / degraded fragment',
        'color':    'coral',
    },
    'No significant similarity': {
        'rule':     'has a hit but qcov, scov or pident below thresholds',
        'meaning':  'Faint similarity to a Core/Accessory gene but not classifiable',
        'biology':  'Highly diverged paralog, ambiguous case. May or may not be a real rare gene.',
        'rare_genome_status': 'AMBIGUOUS: possible diverged copy',
        'color':    'lightgrey',
    },
    'No hit': {
        'rule':     'no DIAMOND hit at all against Core+Accessory pool',
        'meaning':  'No homolog in the species\' Core or Accessory genome',
        'biology':  'Genuine novel/orphan rare gene. Candidates for strain-specific innovations, horizontally transferred genes, or lineage-specific evolution. THIS is the biologically interesting set.',
        'rare_genome_status': 'TRUE rare gene: novel / orphan / candidate xenolog',
        'color':    '#8D8585',
    },
}


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
            'A. fumigatus': '#e74c3c', 'A. flavus': '#27ae60',
            'A. niger': '#1a1aff', 'A. oryzae': '#ffd500',
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


# =============================================================================
# NB4 ANALYSIS AND FIGURE FUNCTIONS
# =============================================================================

# Phenotype palette shared with the NB0 phenotype-distribution figure.
NB4_PHENO_COLORS = {
    'Human-pathogenic': '#e74c3c', 'Plant-pathogenic': '#e67e22',
    'Animal-pathogenic': '#7d3c98', 'Environmental': '#27ae60',
    'Industrial-trait': '#3498db',
}


def plot_xenolog_neighbourhood(base_path, results_dir, flank=8):
    """Gene neighbourhood of each confident xenolog.

    Reads gene to orthogroup assignments per genome column from the
    OrthoFinder Orthogroups tables, colours genes by pangenome class and
    labels them with their eggNOG annotation.
    """
    FLANK = flank
    _BASE = str(base_path)
    _OUT = str(results_dir)

    import os, glob, re
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Patch

    THRESH = {'fumigatus': (85, 6), 'flavus': (68, 6)}
    CLS_COLOR = {'Core':'#2c6fbb','Accessory':'#e08b2d','Rare':'#9aa0a6','XENOLOG':'#c0392b'}

    XENOLOGS = [
     {'label':'OG0014462','sp':'flavus','acc':'GCA_052815365.1','contig':'scaffold_69',
      'xen':['FUN_012569','FUN_012574'],'foreign':'S. maltophilia Zot toxin (Bacteria, 98.4%)'},
     {'label':'OG0010936','sp':'fumigatus','acc':'GCA_049901915.1','contig':'scaffold_31',
      'xen':['FUN_006532','FUN_006533'],'foreign':'S. maltophilia DUF3299 (Bacteria, 65%)'},
    ]

    def _load_run(sp):
        ogf=glob.glob(f'{_BASE}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups.tsv')[0]
        unf=glob.glob(f'{_BASE}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups_UnassignedGenes.tsv')[0]
        with open(ogf) as fh:
            header=fh.readline().rstrip('\n').split('\t'); rows=[l.rstrip('\n').split('\t') for l in fh]
        with open(unf) as fh:
            fh.readline(); urows=[l.rstrip('\n').split('\t') for l in fh]
        og_ng={r[0]:sum(1 for x in r[1:] if x.strip()) for r in rows}
        for r in urows: og_ng[r[0]]=1
        return header, rows, urows, og_ng

    def _gene2og(header, rows, urows, acc):
        c=next((i for i,h in enumerate(header) if acc.split('.')[0] in h), None)
        g={}
        for r in rows+urows:
            if c<len(r) and r[c].strip():
                for p in r[c].split(','): g[p.strip()]=r[0]
        return g

    def _gff(sp, acc):
        path=glob.glob(f'{_BASE}/{sp}/funannotate_output/{acc}_*_genomic_renamed/predict_results/{acc}_*.gff3')[0]
        d={}
        for ln in open(path):
            if ln.startswith('#'): continue
            f=ln.rstrip('\n').split('\t')
            if len(f)<9 or f[2]!='gene': continue
            d[f[8].split('ID=')[1].split(';')[0]]=(f[0],int(f[3]),int(f[4]),f[6])
        return d

    def _eggnog_labels(sp, acc):
        """protein base -> concise functional label from per-genome eggNOG annotations."""
        hits=glob.glob(f'{_BASE}/{sp}/eggnog_output/{acc}_*.emapper.annotations')
        lab={}
        if not hits: return lab
        rows=[l.rstrip('\n').split('\t') for l in open(hits[0]) if not l.startswith('##')]
        hdr=[h.lstrip('#') for h in rows[0]]; idx={h:i for i,h in enumerate(hdr)}
        def col(r,n):
            i=idx.get(n); return (r[i] if i is not None and i<len(r) else '').strip()
        PFMAP={'Zot':'Zot toxin','DUF3299':'DUF3299 (bacterial)','MFS_1':'MFS transporter',
               'Aldo_ket_red':'aldo-keto reductase','ACR_tran':'RND efflux pump','MacB_PCD':'MacB transporter',
               'FtsX':'FtsX/MacB transporter','LysR_substrate':'LysR regulator','HTH_1':'LysR regulator',
               'Bac_luciferase':'F420 luciferase','adh_short_C2':'NAD epimerase','DUF2796':'DUF2796',
               'Mrr_cat':'restriction endonuclease','KfrA_N':'plasmid rep (KfrA)','DUF3693':'phage protein',
               'Phage_DNA_bind':'phage ssDNA-binding'}
        for r in rows[1:]:
            base=r[0].split('-')[0]
            nm=col(r,'Preferred_name'); pf=col(r,'PFAMs'); desc=col(r,'Description')
            dom=pf.split(',')[0] if pf and pf!='-' else ''
            if nm and nm!='-': s=nm
            elif dom in PFMAP: s=PFMAP[dom]
            elif desc and desc not in ('-','') and not re.match(r'(Belongs to|Catalyzes|COG\d)',desc):
                s=' '.join(re.sub(r'\s*\(.*?\)','',desc).split()[:3])
            elif dom: s=dom
            else: s=''
            if s: lab[base]=s[:26]
        return lab

    def _cls(sp, ng):
        core,rare=THRESH[sp]
        return 'Core' if ng>=core else ('Rare' if ng<=rare else 'Accessory')

    tracks=[]
    for X in XENOLOGS:
        header,rows,urows,og_ng=_load_run(X['sp'])
        g2o=_gene2og(header,rows,urows,X['acc']); genes=_gff(X['sp'],X['acc']); egg=_eggnog_labels(X['sp'],X['acc'])
        on=sorted([g for g in genes if genes[g][0]==X['contig']], key=lambda g:genes[g][1])
        idx=[on.index(x) for x in X['xen'] if x in on]
        lo,hi=max(0,min(idx)-FLANK),min(len(on),max(idx)+FLANK+1)
        row=[]
        for g in on[lo:hi]:
            c,s,e,st=genes[g]; og=g2o.get(g+'-T1',g2o.get(g,'NA')); ng=og_ng.get(og,0); isx=g in X['xen']
            name=egg.get(g,'') or 'hypothetical'
            row.append({'gene':g,'start':s,'end':e,'strand':st,'og':og,'ng':ng,'name':name,
                        'cls':('XENOLOG' if isx else _cls(X['sp'],ng)),'is_xen':isx})
        tracks.append({**X,'genes':row,'ncontig':len(on),'clen':max(genes[g][2] for g in on)})

    fig,axes=plt.subplots(len(tracks),1,figsize=(14,3.9*len(tracks)+0.8))
    if len(tracks)==1: axes=[axes]
    def _arrow(ax,x0,x1,strand,y,color,h,head):
        L=abs(x1-x0); head=min(head,L*0.9)
        pts=([(x0,y-h),(x1-head,y-h),(x1,y),(x1-head,y+h),(x0,y+h)] if strand=='+'
             else [(x1,y-h),(x0+head,y-h),(x0,y),(x0+head,y+h),(x1,y+h)])
        ax.add_patch(Polygon(pts,closed=True,facecolor=color,edgecolor='black',lw=0.7,zorder=3))

    for ax,t in zip(axes,tracks):
        xen=[g for g in t['genes'] if g['is_xen']]
        orig=sum((g['start']+g['end'])/2 for g in xen)/len(xen)
        xs=[(g['start']-orig)/1000 for g in t['genes']]+[(g['end']-orig)/1000 for g in t['genes']]
        span=max(xs)-min(xs); head=span*0.012
        ax.plot([min(xs),max(xs)],[0,0],color='#cccccc',lw=1.0,zorder=1); done=False
        for i,g in enumerate(t['genes']):
            x0=(g['start']-orig)/1000; x1=(g['end']-orig)/1000; xm=(x0+x1)/2
            _arrow(ax,x0,x1,g['strand'],0,CLS_COLOR[g['cls']],0.24,head)
            # functional NAME label, alternating above/below to avoid overlap
            above = (i%2==0)
            yl = 0.34 if above else -0.34
            ax.annotate(g['name'],(xm,yl),ha='left',va=('bottom' if above else 'top'),
                        fontsize=6.6,rotation=32,color=('#111' if g['is_xen'] else '#333'),
                        fontweight=('bold' if g['is_xen'] else 'normal'))
            if g['cls']=='Core':
                ax.annotate(f"core {g['ng']}/{ {'fumigatus':89,'flavus':70}[t['sp']] }",(xm,0),(xm,-0.02),
                            ha='center',va='center',fontsize=0,alpha=0)
            if g['is_xen'] and not done:
                xxs=[(gg['start']-orig)/1000 for gg in xen]+[(gg['end']-orig)/1000 for gg in xen]
                ax.annotate('xenolog',(sum(xxs)/len(xxs),0.95),ha='center',va='bottom',
                            fontsize=9,color=CLS_COLOR['XENOLOG'],fontweight='bold'); done=True
        ax.set_title(f"{t['label']}  |  {t['acc']}  {t['contig']} ({t['clen']/1000:.0f} kb)  |  best foreign hit: {t['foreign']}",
                     fontsize=10,fontweight='bold')
        ax.set_xlim(min(xs)-span*0.03,max(xs)+span*0.22); ax.set_ylim(-1.5,1.4); ax.axis('off')
        bar=round(span*0.2/5)*5 or 5
        ax.plot([min(xs),min(xs)+bar],[-1.32,-1.32],color='k',lw=1.4)
        ax.annotate(f'{bar:.0f} kb',(min(xs)+bar/2,-1.42),ha='center',va='top',fontsize=7.5)

    _leg=[Patch(facecolor=CLS_COLOR[k],edgecolor='k',label=l) for k,l in
          [('Core','Core (all strains)'),('Accessory','Accessory'),('Rare','Rare / strain-specific'),
           ('XENOLOG','Xenolog (foreign best hit)')]]
    fig.legend(handles=_leg,loc='upper center',ncol=4,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
    fig.suptitle('Genomic neighbourhood of the confident xenologs (genes labelled by function, coloured by pangenome class)',
                 fontsize=11.5,fontweight='bold',y=1.045)
    plt.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(_OUT,'xenolog_neighbourhood.png'),dpi=200,bbox_inches='tight')
    for t in tracks:
        print(f"\n{t['label']} {t['acc']} {t['contig']}:")
        for g in t['genes']:
            print(f"   {g['gene']} [{g['cls']:9s} {g['ng']:>2}] {g['name']}")

    return fig


def plot_diamond_classification_pies(nb4_results, species_list=None, fs_dp=15,
                                     dpi_dp=400, fig_w_dp=15.0, fig_h_dp=4.6,
                                     nontruly_alpha=1.0):
    """DIAMOND classification of rare orthogroups, one pie per species.

    Counts the classes directly from the per-species
    ``{nb4_results}/{species}/diamond_rare_classification.tsv`` tables.
    """
    FS_DP = fs_dp
    DPI_DP = dpi_dp
    FIG_W_DP = fig_w_dp
    FIG_H_DP = fig_h_dp
    NONTRULY_ALPHA = nontruly_alpha
    NB4_RESULTS = Path(nb4_results)
    species_list = list(species_list) if species_list is not None else list(SPECIES_LIST)
    fig = None

    import pandas as pd, numpy as np, matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    classes = ['Full-length homolog', 'Fragment of longer gene',
               'No significant similarity', 'No hit']
    short   = ['Full-length homolog (paralog)', 'Fragment of longer gene',
               'No significant similarity', 'No hit (truly rare)']
    colors  = ["#313030", '#595959', '#969696', "#464ee8"]  # black -> light grey ramp

    # Class counts per species, read from the cached DIAMOND tables.
    _rows = []
    for _sp in species_list:
        _fp = NB4_RESULTS / _sp / 'diamond_rare_classification.tsv'
        if not _fp.exists():
            print(f'No DIAMOND table for {_sp} at {_fp}.')
            continue
        _vc = pd.read_csv(_fp, sep='\t')['blast_class'].value_counts()
        _row = {'Species': f'A. {_sp}', 'Total': int(_vc.sum())}
        _row.update({c: int(_vc.get(c, 0)) for c in classes})
        _rows.append(_row)

    if not _rows:
        print('No DIAMOND tables found; run the classification section first.')
    else:
        d = pd.DataFrame(_rows)

        fig, axes = plt.subplots(1, len(d), figsize=(FIG_W_DP, FIG_H_DP),
                                 squeeze=False, gridspec_kw=dict(wspace=-0.03))
        for ax, (_, row) in zip(axes.ravel(), d.iterrows()):
            sizes  = [int(row[c]) for c in classes]
            wedges, _t, autotxts = ax.pie(
                sizes, labels=None, colors=colors,
                autopct='%1.1f%%',
                startangle=90,
                wedgeprops=dict(linewidth=1.0, edgecolor='white'),
                textprops=dict(fontsize=FS_DP - 2))
            # De-emphasise the 3 non-truly-rare wedges; keep "No hit" solid.
            for wi, wedge in enumerate(wedges):
                wedge.set_alpha(1.0 if wi == 3 else NONTRULY_ALPHA)
            # % text inside each wedge: white on dark slices, black on light slices
            for tt, col, n in zip(autotxts, colors, sizes):
                lum = int(col.lstrip('#')[0:2], 16)   # grey ramp -> R==G==B
                tt.set_color('white' if lum < 128 else 'black')
                tt.set_fontsize(FS_DP - 2); tt.set_fontweight('bold')
                tt.set_text(f'{tt.get_text()}\nn = {n:,}')   # move n= inside the wedge
            ax.set_aspect('equal')
            ax.set_title(f"{row['Species']}\n"
                         f"Total = {int(row['Total']):,} rare OGs",
                         fontsize=FS_DP, fontweight='bold', pad=3)
        # One shared legend at the bottom-right of the figure
        handles = [Patch(facecolor=c, edgecolor='white', label=s)
                   for c, s in zip(colors, short)]
        fig.legend(handles=handles, title='Rare-OG class',
                   loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=4,
                   fontsize=FS_DP - 2, title_fontsize=FS_DP - 2, frameon=True)
        fig.suptitle('DIAMOND Classification of Rare Orthogroups',
                     fontsize=FS_DP + 2, fontweight='bold', y=0.99)
        plt.tight_layout(rect=[0, 0.14, 1, 0.92])
        out = NB4_RESULTS / 'diamond_classification_pies.png'
        fig.savefig(out, dpi=DPI_DP, bbox_inches='tight')
        print(f'Saved -> {out}  ')

    return fig


def build_pooled_covariances(nb1_results, nb4_results, tree_path, species=['fumigatus', 'flavus', 'niger', 'oryzae']):
    """Assemble the covariance matrices used by the pooled burden models.

    Places the per-species kinship matrices block-diagonally and builds the
    combined species-tree covariance. Writes pooled_kinship_nb2.tsv and
    combined_tree_vcv.tsv.
    """
    SPECIES = species
    NB1_RESULTS = Path(nb1_results)
    NB4_RESULTS = Path(nb4_results)
    TREE = str(tree_path)


    import numpy as np, pandas as pd, re, warnings
    from Bio import Phylo
    from scipy.linalg import block_diag
    import funpan_utils as fpu
    warnings.filterwarnings('ignore')

    _ACC = re.compile(r'(GC[AF]_\d+\.\d+)')
    meta = pd.read_csv(NB4_RESULTS / 'rare_gene_burden_all_genomes.csv')

    keep, Ks = [], []
    for sp in SPECIES:
        K = pd.read_csv(NB1_RESULTS / sp / f'{sp}_kinship.tsv', sep='\t', index_col=0)
        K.index = [(_ACC.search(s).group(1) if _ACC.search(s) else s) for s in K.index]; K.columns = K.index
        m = meta[(meta.Species == f'A. {sp}') & (meta.Phenotype != 'Unknown')]
        g = [x for x in m.Genome if x in set(K.index)]
        keep += g
        Kn = K.loc[g, g].values
        Ks.append(Kn / np.mean(np.diag(Kn)))
    Kbd = block_diag(*Ks)
    pd.DataFrame(Kbd, index=keep, columns=keep).to_csv(NB4_RESULTS / 'pooled_kinship_nb2.tsv', sep='\t')

    # --- combined species-tree covariance ---
    T = Phylo.read(TREE, 'newick'); tip = {}
    for t in T.get_terminals():
        mm = _ACC.search(t.name or '')
        if mm: tip[mm.group(1)] = t
    paths = {}
    for g in keep:
        d, acc = 0.0, []
        for cl in T.get_path(tip[g]):
            d += (cl.branch_length or 0.0); acc.append((id(cl), d))
        paths[g] = acc
    n = len(keep); C = np.zeros((n, n)); seen = {g: {c: d for c, d in paths[g]} for g in keep}
    for i, gi in enumerate(keep):
        for j in range(i, n):
            sj = seen[keep[j]]; sh = 0.0
            for c, d in paths[gi]:
                if c in sj: sh = d
                else: break
            C[i, j] = C[j, i] = sh
        C[i, i] = paths[gi][-1][1] if paths[gi] else 0.0
    C = C / np.mean(np.diag(C))
    pd.DataFrame(C, index=keep, columns=keep).to_csv(NB4_RESULTS / 'combined_tree_vcv.tsv', sep='\t')

    spv = np.array([next(s for s in SPECIES if g in set(meta.loc[meta.Species == f'A. {s}', 'Genome'])) for g in keep])
    off = ~np.eye(n, dtype=bool)
    print(f'Pooled covariance for the burden models = block-diagonal NB2 per-species kinship, {n} x {n}')
    print(f'Saved -> {NB4_RESULTS / "pooled_kinship_nb2.tsv"}')
    print(f'Tree VCV (phylANOVA + diagnostic only) -> {NB4_RESULTS / "combined_tree_vcv.tsv"}\n')
    print('Within-species INDEPENDENT variation -- why the NB2 GRM is used to fit and the tree is not:')
    print(f"{'species':13s} {'n':>4s} {'NB2 GRM':>22s} {'combined protein tree':>24s}")
    for sp in SPECIES:
        m_ = (spv == sp)[:, None] & (spv == sp)[None, :] & off
        i_ = np.flatnonzero(spv == sp)
        kb = Kbd[np.ix_(i_, i_)]; ct = C[np.ix_(i_, i_)]
        kshare = 100 * (np.mean(np.diag(kb)) - kb[~np.eye(len(i_), dtype=bool)].mean()) / np.mean(np.diag(kb))
        cshare = 100 * (np.mean(np.diag(ct)) - ct[~np.eye(len(i_), dtype=bool)].mean()) / np.mean(np.diag(ct))
        print(f'  A. {sp:11s} {len(i_):4d} {kshare:21.1f}% {cshare:23.1f}%')
    print('\nThe NB2 GRM retains essentially all within-species variation; the protein tree retains ~1%.')
    print('Between-species entries are 0 in the block-diagonal matrix by construction -- handled by')
    print('species fixed effects in every model, not by the covariance.')

    return {'kinship': pd.DataFrame(Kbd, index=keep, columns=keep),
            'tree_vcv': pd.DataFrame(C, index=keep, columns=keep)}


def run_burden_kinship_control(nb1_results, nb4_results, npc=5):
    """Burden contrasts corrected for population structure.

    Reports the naive p-value, the SNP-PC-corrected p-value, the kinship
    mixed-model p-value and the burden heritability per contrast.
    """
    NPC = npc
    NB1_RESULTS = Path(nb1_results)
    NB4_RESULTS = Path(nb4_results)

    import numpy as np, pandas as pd, re, warnings
    from scipy.optimize import minimize_scalar
    from scipy import stats
    import statsmodels.formula.api as smf
    import funpan_utils as fpu
    warnings.filterwarnings('ignore')

    _VCV = pd.read_csv(NB4_RESULTS / 'pooled_kinship_nb2.tsv', sep='\t', index_col=0)
    _excl = set().union(*fpu.ANI_EXCLUDED.values())
    _ACC = re.compile(r'(GC[AF]_\d+\.\d+)')
    meta = pd.read_csv(NB4_RESULTS / 'rare_gene_burden_all_genomes.csv')

    def _truly(sp):
        pav = pd.read_csv(NB1_RESULTS / sp / f'{sp}_pav.tsv', sep='\t'); pav['og'] = pav['Orthogroup'].astype(str)
        gcols = [c for c in pav.columns if c not in ('Orthogroup', 'og') and not any(a in c for a in _excl)]
        dc = pd.read_csv(NB4_RESULTS / sp / 'diamond_rare_classification.tsv', sep='\t')
        nohit = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'].astype(str))
        b = (pav[pav['og'].isin(nohit)][gcols] > 0).sum(0)
        return {g: int(v) for g, v in b.items()}

    def _load(sp, case):
        tr = _truly(sp)
        K = pd.read_csv(NB1_RESULTS / sp / f'{sp}_kinship.tsv', sep='\t', index_col=0)
        K.index = [(_ACC.search(s).group(1) if _ACC.search(s) else s) for s in K.index]; K.columns = K.index
        pcs = pd.read_csv(NB1_RESULTS / sp / f'{sp}_snp_pcs.tsv', sep='\t')
        pcs['acc'] = pcs['sample_id'].apply(lambda s: _ACC.search(s).group(1) if _ACC.search(s) else None)
        m = meta[meta.Species == f'A. {sp}'][['Genome', 'Phenotype', 'Species']].copy(); m['burden'] = m['Genome'].map(tr)
        m = m[(m.Phenotype != 'Unknown') & (m.Genome.isin(K.index))].merge(
            pcs[['acc'] + [f'PC{i}' for i in range(1, NPC + 1)]], left_on='Genome', right_on='acc', how='inner')
        m['is_case'] = (m.Phenotype == case).astype(float)
        Kn = K.loc[m.Genome, m.Genome].values; Kn = Kn / np.mean(np.diag(Kn))
        return m.reset_index(drop=True), Kn

    def _grm_lmm(y, X, K):
        n = len(y); lam, U = np.linalg.eigh(K); lam = np.maximum(lam, 1e-6)
        Uty, UtX = U.T @ y, U.T @ X
        def nll(ld):
            d = lam + np.exp(ld); XtWX = (UtX.T * (1 / d)) @ UtX; beta = np.linalg.solve(XtWX, (UtX.T * (1 / d)) @ Uty)
            r = Uty - UtX @ beta; s2 = np.sum(r * r / d) / n
            return 0.5 * (n * np.log(2 * np.pi * s2) + np.sum(np.log(d)) + n)
        res = minimize_scalar(nll, bounds=(-8, 8), method='bounded'); delta = np.exp(res.x); d = lam + delta
        XtWX = (UtX.T * (1 / d)) @ UtX; beta = np.linalg.solve(XtWX, (UtX.T * (1 / d)) @ Uty)
        r = Uty - UtX @ beta; s2 = np.sum(r * r / d) / n; se = np.sqrt(np.diag(s2 * np.linalg.inv(XtWX)))
        p = 2 * stats.t.sf(np.abs(beta / se), n - X.shape[1])
        return beta[1], p[1], 1 / (1 + delta)

    def _pc_block(M, pair):
        """Per-species PC covariates, zero-padded across species (for the combined GRM)."""
        cols = []
        for sp in pair:
            mask = (M.Species == f'A. {sp}').astype(float).values
            for i in range(1, NPC + 1):
                cols.append(np.nan_to_num(M[f'PC{i}'].values) * mask)
        return cols

    rows = []
    # --- per species ---
    for sp, case in [('fumigatus', 'Human-pathogenic'), ('flavus', 'Human-pathogenic'), ('niger', 'Industrial-trait')]:
        m, K = _load(sp, case)
        y = np.log1p(m['burden'].values.astype(float)); one = np.ones(len(y))
        pcc = [m[f'PC{i}'].values for i in range(1, NPC + 1)]
        naive = smf.ols("np.log1p(burden) ~ is_case", m).fit().pvalues['is_case']
        pc = smf.ols("np.log1p(burden) ~ is_case + " + "+".join(f'PC{i}' for i in range(1, NPC + 1)), m).fit().pvalues['is_case']
        _, pk, h2 = _grm_lmm(y, np.column_stack([one, m.is_case.values]), K)
        _, pkpc, _ = _grm_lmm(y, np.column_stack([one, m.is_case.values] + pcc), K)
        rows.append({'contrast': f'A. {sp} | {case} vs rest', 'n': len(m), 'naive_p': naive,
                     'PC_only_p': pc, 'kinship_only_p': pk, 'kinship_plus_PC_p': pkpc, 'burden_h2': round(h2, 2)})
    # --- combined (NB2 per-species kinship + species fixed effect) ---
    for pair, case, sp2 in [(('fumigatus', 'flavus'), 'Human-pathogenic', 'flavus'),
                            (('niger', 'oryzae'), 'Industrial-trait', 'oryzae')]:
        ms, Ks = [], []
        for sp in pair:
            m, K = _load(sp, case); ms.append(m); Ks.append(K)
        M = pd.concat(ms, ignore_index=True); Kc = _VCV.loc[M.Genome, M.Genome].values
        M['is_sp2'] = (M.Species == f'A. {sp2}').astype(float)
        y = np.log1p(M['burden'].values.astype(float)); base = [np.ones(len(M)), M.is_case.values, M.is_sp2.values]
        naive = smf.mixedlm("np.log1p(burden) ~ is_case", M, groups=M['Species']).fit(reml=False).pvalues['is_case']
        _, pk, h2 = _grm_lmm(y, np.column_stack(base), Kc)
        _, pkpc, _ = _grm_lmm(y, np.column_stack(base + _pc_block(M, pair)), Kc)
        rows.append({'contrast': f'{"+".join(pair)} (combined) | {case} vs rest', 'n': len(M), 'naive_p': naive,
                     'PC_only_p': np.nan, 'kinship_only_p': pk, 'kinship_plus_PC_p': pkpc, 'burden_h2': round(h2, 2)})

    res = pd.DataFrame(rows)
    res.to_csv(NB4_RESULTS / 'rare_burden_kinship_control.csv', index=False)
    disp = res.copy()
    for c in ('naive_p', 'PC_only_p', 'kinship_only_p', 'kinship_plus_PC_p'):
        disp[c] = disp[c].apply(lambda v: f'{v:.1e}' if pd.notna(v) else '-')
    print('Does the truly-rare burden survive correction for phylogenetic relatedness?')
    print('(burden_h2 = fraction of burden variance explained by kinship;')
    print(' kinship_only = GRM random effect; kinship_plus_PC = GRM + top-5 SNP PCs as fixed covariates)\n')
    print(disp[['contrast', 'n', 'naive_p', 'PC_only_p', 'kinship_only_p', 'kinship_plus_PC_p', 'burden_h2']].to_string(index=False))
    print(f"\nSaved -> {NB4_RESULTS / 'rare_burden_kinship_control.csv'}")

    return res


def plot_burden_four_panel(nb1_results, nb4_results, qc_table, fs=23, dpi=400, fig_w=13.2, fig_h=13.9, c_rest='#8a97a3', bg_pair1='#eaf1f8', bg_pair2='#e9f4ea', species=['fumigatus', 'flavus', 'niger', 'oryzae']):
    """The four burden contrasts in two pairs on a shared axis.
    """
    FS = fs
    DPI = dpi
    FIG_W = fig_w
    FIG_H = fig_h
    C_REST = c_rest
    BG_PAIR1 = bg_pair1
    BG_PAIR2 = bg_pair2
    SPECIES = species
    NB1_RESULTS = Path(nb1_results)
    NB4_RESULTS = Path(nb4_results)
    QC_TABLE = str(qc_table)

    C_PHENO = {'Human-pathogenic': '#e74c3c', 'Plant-pathogenic': '#e67e22',
               'Animal-pathogenic': '#7d3c98', 'Environmental': '#27ae60',
               'Industrial-trait': '#3498db'}

    import numpy as np, pandas as pd, re, warnings, zlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, FancyBboxPatch
    from scipy.optimize import minimize_scalar
    from scipy import stats
    import funpan_utils as fpu
    warnings.filterwarnings('ignore')

    _excl = set().union(*fpu.ANI_EXCLUDED.values()); _ACC = re.compile(r'(GC[AF]_\d+\.\d+)')
    meta = pd.read_csv(NB4_RESULTS / 'rare_gene_burden_all_genomes.csv')


    def _truly(sp):
        pav = pd.read_csv(NB1_RESULTS / sp / f'{sp}_pav.tsv', sep='\t'); pav['og'] = pav['Orthogroup'].astype(str)
        gcols = [c for c in pav.columns if c not in ('Orthogroup', 'og') and not any(a in c for a in _excl)]
        dc = pd.read_csv(NB4_RESULTS / sp / 'diamond_rare_classification.tsv', sep='\t')
        nohit = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'].astype(str))
        return {g: int(v) for g, v in (pav[pav['og'].isin(nohit)][gcols] > 0).sum(0).items()}


    bm = {}
    for sp in SPECIES:
        bm.update(_truly(sp))
    KIN = pd.read_csv(NB4_RESULTS / 'pooled_kinship_nb2.tsv', sep='\t', index_col=0)
    TREE = pd.read_csv(NB4_RESULTS / 'combined_tree_vcv.tsv', sep='\t', index_col=0)
    M = meta.set_index('Genome').loc[list(KIN.index)].reset_index()
    M.columns = ['Genome'] + list(M.columns[1:])
    M['burden'] = M['Genome'].map(bm)
    qc = pd.read_csv(QC_TABLE); qc['Genome'] = qc['Assembly Accession']
    M = M.merge(qc[['Genome', 'total_len']], on='Genome', how='left')
    M['Mb'] = M.total_len / 1e6; M['rate'] = M.burden / M.Mb
    M = M.reset_index(drop=True)


    def _fit(sub, cov, species_fe):
        """log1p(burden) - log(Mb) ~ case [+ species FE], with the given covariance."""
        n = len(sub); Ks = cov.loc[sub.Genome, sub.Genome].values
        lam, U = np.linalg.eigh(Ks); lam = np.maximum(lam, 1e-8)
        y = np.log1p(sub.burden.values.astype(float)) - np.log(sub.Mb.values)
        cols = [np.ones(n), sub['case'].values]
        if species_fe:
            pres = [s for s in SPECIES if (sub.Species == f'A. {s}').any()]
            cols += [(sub.Species == f'A. {s}').astype(float).values for s in pres[1:]]
        X = np.column_stack(cols); Uty, UtX = U.T @ y, U.T @ X
        def nll(ld):
            d = lam + np.exp(ld); A = (UtX.T * (1/d)) @ UtX; b = np.linalg.solve(A, (UtX.T * (1/d)) @ Uty)
            r = Uty - UtX @ b
            return 0.5 * (n * np.log(2 * np.pi * np.sum(r*r/d)/n) + np.sum(np.log(d)) + n)
        dl = np.exp(minimize_scalar(nll, bounds=(-10, 10), method='bounded').x); d = lam + dl
        A = (UtX.T * (1/d)) @ UtX; b = np.linalg.solve(A, (UtX.T * (1/d)) @ Uty); r = Uty - UtX @ b
        s2 = np.sum(r*r/d) / n; se = np.sqrt(np.diag(s2 * np.linalg.inv(A)))
        return float(np.exp(b[1])), float(2 * stats.t.sf(abs(b[1] / se[1]), n - X.shape[1]))


    # Species names are set in italic via mathtext: \mathbfit inside the bold panel titles and
    # pair headings, \mathit inside the regular-weight tick labels. _PLAIN() strips the markup
    # again so the CSV keeps readable text.
    def _BI(s): return r'$\mathbfit{' + s.replace(' ', r'\ ') + '}$'     # bold italic


    def _IT(s): return r'$\mathit{' + s.replace(' ', r'\ ') + '}$'       # regular italic


    def _PLAIN(s):
        return re.sub(r'\$\\math(?:bfit|it)\{(.*?)\}\$',
                      lambda m: m.group(1).replace('\\ ', ' '), s)


    SPP = M.Species
    PANELS = [
     ('A', f"Human-pathogenic vs Rest\n({_BI('A. fumigatus')} + {_BI('A. flavus')})",
      SPP.isin(['A. fumigatus', 'A. flavus']) & (M.Phenotype == 'Human-pathogenic'),
      SPP.isin(['A. fumigatus', 'A. flavus']) & (M.Phenotype != 'Human-pathogenic'),
      'Human-\npathogenic', 'Rest', C_PHENO['Human-pathogenic'], C_REST, KIN, True),
     ('B', f"Industrial vs Rest\n({_BI('A. niger')} + {_BI('A. oryzae')})",
      SPP.isin(['A. niger', 'A. oryzae']) & (M.Phenotype == 'Industrial-trait'),
      SPP.isin(['A. niger', 'A. oryzae']) & (M.Phenotype != 'Industrial-trait'),
      'Industrial', 'Rest', C_PHENO['Industrial-trait'], C_REST, KIN, True),
     ('C', f"{_BI('A. oryzae')} Industrial vs\n{_BI('A. flavus')} Human-pathogenic",
      (SPP == 'A. oryzae') & (M.Phenotype == 'Industrial-trait'),
      (SPP == 'A. flavus') & (M.Phenotype == 'Human-pathogenic'),
      f"{_IT('A. oryzae')}\nIndustrial", f"{_IT('A. flavus')}\nHuman-pathogenic",
      C_PHENO['Industrial-trait'], C_PHENO['Human-pathogenic'], KIN, False),
     ('D', f"{_BI('A. oryzae')} Industrial vs\n{_BI('A. flavus')} Plant-pathogenic",
      (SPP == 'A. oryzae') & (M.Phenotype == 'Industrial-trait'),
      (SPP == 'A. flavus') & (M.Phenotype == 'Plant-pathogenic'),
      f"{_IT('A. oryzae')}\nIndustrial", f"{_IT('A. flavus')}\nPlant-pathogenic",
      C_PHENO['Industrial-trait'], C_PHENO['Plant-pathogenic'], KIN, False),
    ]
    res = []
    for tag, ttl, a, b, la, lb, ca, cb, cov, fe in PANELS:
        A, B = M[a.values], M[b.values]
        sub = pd.concat([A.assign(case=1.0), B.assign(case=0.0)], ignore_index=True)
        f, p = _fit(sub, cov, fe)
        res.append(dict(panel=tag, contrast=_PLAIN(ttl).replace('\n', ' '), n_case=len(A), n_rest=len(B),
                        mean_case=A.rate.mean(), mean_rest=B.rate.mean(), fold=f, p=p,
                        covariance='NB2 kinship + species FE' if fe else 'NB2 kinship, no species FE'))
    R = pd.DataFrame(res)
    R.to_csv(NB4_RESULTS / 'burden_manuscript_four_panel.csv', index=False)


    def _stars(p): return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'


    def _jit(k, key):
        return np.random.RandomState(zlib.crc32(key.encode()) % (2 ** 31)).uniform(-0.10, 0.10, k)


    gmax = M.rate.max()
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = fig.add_gridspec(2, 2, hspace=1.05, wspace=0.08,
                          top=0.915, bottom=0.145, left=0.130, right=0.980)
    axes = [fig.add_subplot(gs[0, 0])]
    axes += [fig.add_subplot(gs[r, c], sharey=axes[0]) for r, c in [(0, 1), (1, 0), (1, 1)]]
    for ax, (tag, ttl, a, b, la, lb, ca, cb, cov, fe), (_, r) in zip(axes, PANELS, R.iterrows()):
        A, B = M[a.values], M[b.values]
        groups, cols, labs = [A, B], [ca, cb], [la, lb]
        bp = ax.boxplot([g['rate'].values for g in groups], positions=[0, 0.72], widths=0.34,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color='black', linewidth=1.4),
                        whiskerprops=dict(color='#444'), capprops=dict(color='#444'))
        for patch, c in zip(bp['boxes'], cols):
            patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor('black')
        for x, g, c, key in zip([0, 0.72], groups, cols, [tag + la, tag + lb]):
            ax.scatter(x + _jit(len(g), key), g['rate'].values, s=26, color=c, alpha=0.9,
                       edgecolor='black', linewidth=0.3, zorder=3)
            ax.scatter([x], [g['rate'].mean()], marker='D', s=68, color='white',
                       edgecolor='black', linewidth=1.3, zorder=5)
        trans = ax.get_xaxis_transform(); by = 0.79
        ax.plot([0, 0, 0.72, 0.72], [by, by + 0.02, by + 0.02, by], transform=trans, color='black',
                linewidth=1.2, clip_on=False)
        sig = r.p < 0.05
        arrow = '\u25b2' if r.fold > 1 else '\u25bc'
        ax.text(0.36, by + 0.04, f"{_stars(r.p)}   {r.fold:.2f}\u00d7 {arrow}\np = {r.p:.1e}", transform=trans,
                ha='center', va='bottom', fontsize=FS - 5, fontweight='bold',
                color='black' if sig else '#666',
                bbox=dict(boxstyle='round,pad=0.30', facecolor='#ffe066' if sig else '#ffffff',
                          edgecolor='#c9a227' if sig else '#c8ccd2', linewidth=1.0))
        ax.set_xticks([0, 0.72])
        ax.set_xticklabels([f"{labs[0]}\nn={len(A)}\nx\u0304={A.rate.mean():.2f}",
                            f"{labs[1]}\nn={len(B)}\nx\u0304={B.rate.mean():.2f}"], fontsize=FS - 7)
        ax.tick_params(axis='y', labelsize=FS - 7)
        ax.set_title(ttl, fontsize=FS - 5, fontweight='bold', pad=10)
        ax.set_facecolor('none')                       # the pair block behind shows through
        ax.set_xlim(-0.45, 1.17)
        ax.set_ylim(-0.05, gmax * 1.34)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)
        ax.spines['bottom'].set_bounds(-0.32, 1.04)      # keep the axis line inside the pair block
    for ax in (axes[1], axes[3]):
        ax.tick_params(labelleft=False); ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0)

    # --- one large coloured block behind each pair ---------------------------------------
    # Sized from the true drawn extent of each row (get_tightbbox: titles, tick labels,
    # annotation boxes and all), then padded equally on every side, so nothing bleeds out
    # and the plots sit centred inside their block.
    for ax in (axes[0], axes[2]):
        ax.set_ylabel('Truly-rare orthogroups\nper Mb', fontsize=FS - 5)
    fig.canvas.draw()
    _rend = fig.canvas.get_renderer()
    _inv = fig.transFigure.inverted()
    _bb = [axes[k].get_tightbbox(_rend).transformed(_inv) for k in range(4)]
    PAD_X, PAD_Y = 0.022, 0.020
    NOTE_H = 0.034                                  # room reserved for the methods note
    BX0 = min(b.x0 for b in _bb) - PAD_X
    BX1 = max(b.x1 for b in _bb) + PAD_X


    def _block(i, j, colour, title, note):
        y0 = min(_bb[i].y0, _bb[j].y0) - PAD_Y - NOTE_H
        y1 = max(_bb[i].y1, _bb[j].y1) + PAD_Y
        fig.add_artist(FancyBboxPatch((BX0, y0), BX1 - BX0, y1 - y0,
                                      boxstyle='round,pad=0.004,rounding_size=0.014',
                                      transform=fig.transFigure, facecolor=colour,
                                      edgecolor='none', linewidth=0, zorder=0))
        fig.text((BX0 + BX1) / 2, y1 + 0.012, title, ha='center', va='bottom',
                 fontsize=FS - 2, fontweight='bold', color='#222')
        fig.text((BX0 + BX1) / 2, y0 + 0.012, note, ha='center', va='bottom',
                 fontsize=FS - 8, color='#3f4b56')
        return y0


    NOTE1 = ('Log-transformed, size-normalised and kinship-corrected; species differences removed')
    NOTE2 = ('Log-transformed, size-normalised and kinship-corrected')

    _block(0, 1, BG_PAIR1, 'PAIR 1 : lineage-matched phenotype contrasts', NOTE1)
    _BLK2_Y0 = _block(2, 3, BG_PAIR2, 'PAIR 2 : domestication: ' + _BI('A. oryzae')
                        + ' vs wild progenitor ' + _BI('A. flavus'), NOTE2)

    handles = [Patch(facecolor=C_PHENO[p], alpha=0.38, edgecolor='black', label=p)
               for p in ['Human-pathogenic', 'Industrial-trait', 'Plant-pathogenic']] + \
              [Patch(facecolor=C_REST, alpha=0.38, edgecolor='black', label='Rest'),
               plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='w', markeredgecolor='black',
                          markersize=9, label='group mean')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, _BLK2_Y0 - 0.014), ncol=len(handles),
               frameon=False, fontsize=FS - 7)

    out = NB4_RESULTS / 'burden_manuscript_four_panel.png'
    print(f'Saved -> {out} \n')
    show = R.copy()
    for c in ('mean_case', 'mean_rest'):
        show[c] = show[c].map('{:.2f}'.format)
    show['fold'] = show['fold'].map('{:.2f}x'.format); show['p'] = show['p'].map('{:.2e}'.format)
    print(show.to_string(index=False))
    print(f"\nSaved -> {NB4_RESULTS / 'burden_manuscript_four_panel.csv'}")

    return fig, R


def plot_compartment_length_pfam_orf(nb1_results, nb4_results, base_path, edge_bp=100, fs=13, dpi=400, fig_w=11.4, fig_h=4.4, compartments=['Core', 'Accessory', 'Rare', 'Truly-rare'], sp_keys=['fumigatus', 'flavus', 'niger', 'oryzae']):
    """Protein length, Pfam coverage and ORF completeness per compartment.

    Per-OG lengths and per-gene ORF flags are cached under nb4_results.
    """
    EDGE_BP = edge_bp
    FS = fs
    DPI = dpi
    FIG_W = fig_w
    FIG_H = fig_h
    COMPARTMENTS = compartments
    SP_KEYS = sp_keys
    NB1_RESULTS = Path(nb1_results)
    NB4_RESULTS = Path(nb4_results)
    BASE = str(base_path)

    COMP_COLOR = {'Core': '#2ca02c', 'Accessory': "#fdbd4d", 'Rare': "#db4c3c", 'Truly-rare': "#1E6BD6"}
    SP_MARK = {'fumigatus': 'o', 'flavus': 's', 'niger': '^', 'oryzae': 'D'}
    SP_OFF = {'fumigatus': -0.185, 'flavus': -0.062, 'niger': 0.062, 'oryzae': 0.185}

    import glob, os, re, warnings
    import numpy as np, pandas as pd, matplotlib.pyplot as plt
    import funpan_utils as fpu
    warnings.filterwarnings('ignore')

    _excl = set().union(*fpu.ANI_EXCLUDED.values()); _ACC = re.compile(r'GC[AF]_\d+\.\d+')
    LEN_CACHE = NB4_RESULTS / 'og_length_by_class.csv'
    ORF_CACHE = NB4_RESULTS / 'orf_completeness_by_class.csv'
    STOPS = {'TAA', 'TAG', 'TGA'}

    def _priv(sp):
        priv = set()
        for f in glob.glob(f'{BASE}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups.tsv') + \
                 glob.glob(f'{BASE}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups_UnassignedGenes.tsv'):
            with open(f) as fh:
                hdr = fh.readline().rstrip('\n').split('\t'); ex = {i for i, h in enumerate(hdr) if any(a in h for a in _excl)}
                if not ex: continue
                for l in fh:
                    r = l.rstrip('\n').split('\t'); p = [i for i in range(1, len(r)) if r[i].strip()]
                    if p and all(i in ex for i in p): priv.add(r[0])
        return priv

    def _og_median_len(fa):
        lens = []; cur = 0
        with open(fa) as fh:
            for line in fh:
                if line[:1] == '>':
                    if cur: lens.append(cur); cur = 0
                else: cur += len(line.strip())
        if cur: lens.append(cur)
        return float(np.median(lens)) if lens else np.nan

    # ---------- per-OG protein length (cached) ----------
    if not LEN_CACHE.exists():
        print('Parsing Orthogroup_Sequences for per-OG protein lengths (first run, ~2-3 min)...')
        lrows = []
        for sp in SP_KEYS:
            seqdir = glob.glob(f'{BASE}/{sp}/orthofinder_output/*/Orthogroup_Sequences')[0]
            ogc = pd.read_csv(NB1_RESULTS / sp / f'{sp}_og_consensus.tsv', sep='\t')
            cls = dict(zip(ogc['Orthogroup'].astype(str), ogc['Pangenome_Class']))
            dc = pd.read_csv(NB4_RESULTS / sp / 'diamond_rare_classification.tsv', sep='\t')
            nohit = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'].astype(str)); priv = _priv(sp)
            for fn in os.scandir(seqdir):
                og = fn.name[:-3]
                if og in priv or cls.get(og) is None: continue
                lrows.append({'Species': sp, 'og': og, 'Class': cls[og], 'TrulyRare': og in nohit,
                              'median_len': _og_median_len(fn.path)})
        pd.DataFrame(lrows).to_csv(LEN_CACHE, index=False)
        print(f'Saved per-OG lengths -> {LEN_CACHE.name}')

    # ---------- intrinsic ORF completeness, per species x compartment (cached) ----------
    def _cds_flags(fp):
        d = {}; k = None; buf = []
        with open(fp) as fh:
            for l in fh:
                if l[0] == '>':
                    if k:
                        s = ''.join(buf).upper(); d[k] = (s[:3] == 'ATG', s[-3:] in STOPS, len(s) % 3 == 0)
                    k = l[1:].split()[0]; buf = []
                else: buf.append(l.strip())
        if k:
            s = ''.join(buf).upper(); d[k] = (s[:3] == 'ATG', s[-3:] in STOPS, len(s) % 3 == 0)
        return d

    def _gff_genes(fp):
        d = {}
        with open(fp) as fh:
            for l in fh:
                if l[0] == '#': continue
                f = l.rstrip('\n').split('\t')
                if len(f) < 9 or f[2] != 'gene': continue
                m = re.search(r'ID=([^;]+)', f[8])
                if m: d[m.group(1)] = (f[0], int(f[3]), int(f[4]))
        return d

    def _fa_lens(fp):
        d = {}; k = None; n = 0
        with open(fp) as fh:
            for l in fh:
                if l[0] == '>':
                    if k: d[k] = n
                    k = l[1:].split()[0]; n = 0
                else: n += len(l.strip())
        if k: d[k] = n
        return d

    if ORF_CACHE.exists():
        orf_sp = pd.read_csv(ORF_CACHE)
        print(f'Loaded cached ORF completeness from {ORF_CACHE.name}')
    else:
        print('Scanning Funannotate CDS/GFF3 for intrinsic ORF completeness (first run, ~3 min)...')
        from collections import defaultdict, Counter
        rows = []
        for sp in SP_KEYS:
            ogc = pd.read_csv(NB1_RESULTS / sp / f'{sp}_og_consensus.tsv', sep='\t')
            dc = pd.read_csv(NB4_RESULTS / sp / 'diamond_rare_classification.tsv', sep='\t')
            nohit = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'].astype(str))
            cls = dict(zip(ogc['Orthogroup'].astype(str), ogc['Pangenome_Class'])); priv = _priv(sp)
            OGt = glob.glob(f'{BASE}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups.tsv')[0]
            OG = pd.read_csv(OGt, sep='\t')
            acols = [c for c in OG.columns if c != 'Orthogroup' and not any(a in c for a in _excl)]
            g2c = defaultdict(dict); C = defaultdict(Counter)
            for _, r in OG.iterrows():
                og = str(r['Orthogroup']); k = cls.get(og)
                if og in priv or k not in ('Core', 'Accessory', 'Rare'): continue
                comps = [k] if k != 'Rare' else (['Rare', 'Truly-rare'] if og in nohit else ['Rare'])
                for c in acols:
                    cell = r[c]
                    if isinstance(cell, str) and cell.strip():
                        for g in cell.replace(' ', '').split(','):
                            if g: g2c[c][g] = comps
            for c in acols:
                asm = c.replace('.proteins', '')
                P = f'{BASE}/{sp}/funannotate_output/{asm}_genomic_renamed/predict_results'
                cds = f'{P}/{asm}.cds-transcripts.fa'
                if not os.path.exists(cds): continue
                CD = _cds_flags(cds)
                GF = _gff_genes(f'{P}/{asm}.gff3') if os.path.exists(f'{P}/{asm}.gff3') else {}
                SL = _fa_lens(f'{P}/{asm}.scaffolds.fa') if os.path.exists(f'{P}/{asm}.scaffolds.fa') else {}
                for g, comps in g2c.get(c, {}).items():
                    v = CD.get(g)
                    if not v: continue
                    atg, stop, tri = v
                    gi = GF.get(g.split('-T')[0]); edge = None
                    if gi and gi[0] in SL:
                        L = SL[gi[0]]; edge = (gi[1] <= EDGE_BP or (L - gi[2]) <= EDGE_BP)
                    for comp in comps:
                        C[comp]['n'] += 1
                        if not atg: C[comp]['no_start'] += 1
                        if not stop: C[comp]['no_stop'] += 1
                        if not tri: C[comp]['not_x3'] += 1
                        if edge: C[comp]['edge'] += 1
                        if (not atg) or (not stop) or (not tri) or (edge is True): C[comp]['incomplete'] += 1
            for comp in COMPARTMENTS:
                n = C[comp]['n']
                if n:
                    rows.append({'Species': sp, 'Compartment': comp, 'n_genes': n,
                                 'pct_no_start': 100*C[comp]['no_start']/n, 'pct_no_stop': 100*C[comp]['no_stop']/n,
                                 'pct_not_x3': 100*C[comp]['not_x3']/n, 'pct_contig_edge': 100*C[comp]['edge']/n,
                                 'pct_incomplete': 100*C[comp]['incomplete']/n})
            print(f'  {sp}: done', flush=True)
        orf_sp = pd.DataFrame(rows); orf_sp.to_csv(ORF_CACHE, index=False)
        print(f'Saved per-species ORF completeness -> {ORF_CACHE.name}')

    len_df = pd.read_csv(LEN_CACHE)
    _NA = {'', '-', 'nan', 'NaN', 'None'}
    def _clean(v): v = str(v).strip(); return '' if v in _NA else v

    # ---------- per (species, compartment) metrics ----------
    recs = []
    for sp in SP_KEYS:
        ogc = pd.read_csv(NB1_RESULTS / sp / f'{sp}_og_consensus.tsv', sep='\t')
        dc = pd.read_csv(NB4_RESULTS / sp / 'diamond_rare_classification.tsv', sep='\t')
        nohit = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'].astype(str))
        ld = len_df[len_df.Species == sp]
        for comp in COMPARTMENTS:
            if comp == 'Truly-rare':
                sub_ogc = ogc[ogc['Orthogroup'].astype(str).isin(nohit)]; len_ogs = ld[ld.TrulyRare]
            else:
                sub_ogc = ogc[ogc['Pangenome_Class'] == comp]; len_ogs = ld[ld.Class == comp]
            length = len_ogs['median_len'].mean()
            pfam = sub_ogc['PFAMs'].apply(lambda v: _clean(v) != '').mean() * 100 if len(sub_ogc) else np.nan
            m = orf_sp[(orf_sp.Species == sp) & (orf_sp.Compartment == comp)]
            inc = m['pct_incomplete'].iloc[0] if len(m) else np.nan
            recs.append({'Species': sp, 'Compartment': comp, 'length': length, 'pfam': pfam, 'incomplete': inc})
    per_sp = pd.DataFrame(recs)
    avg = per_sp.groupby('Compartment').agg(length=('length', 'mean'), length_sd=('length', 'std'),
                                            pfam=('pfam', 'mean'), pfam_sd=('pfam', 'std'),
                                            incomplete=('incomplete', 'mean'), incomplete_sd=('incomplete', 'std')).reindex(COMPARTMENTS)
    avg.to_csv(NB4_RESULTS / 'compartment_length_pfam_orf.csv')

    # ---------- figure: mean bar per compartment, with every species shown as a point ----------
    # Bars are means across the four species, whiskers the SD, points the individual species.
    GRID, AXIS, INK = '#e3e7eb', '#9aa4ad', '#1d2429'
    fig, axes = plt.subplots(1, 3, figsize=(FIG_W, FIG_H))
    x = np.arange(len(COMPARTMENTS))
    panels = [('length', 'Protein length', 'amino acids (mean per orthogroup)', '%.0f'),
              ('pfam', 'Pfam domain coverage', '% of orthogroups with a domain', '%.0f'),
              ('incomplete', 'Incomplete gene models', '% of genes', '%.2f')]
    for ax, (col, title, ylab, fmt) in zip(axes, panels):
        vals = avg[col].values
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID, linewidth=1.0)
        sd = np.nan_to_num(avg[col + '_sd'].values)
        ax.bar(x, vals, 0.60, color=[COMP_COLOR[c] for c in COMPARTMENTS],
               edgecolor='#33383d', linewidth=0.7, zorder=2,
               yerr=sd, error_kw=dict(ecolor='#33383d', elinewidth=1.1, capsize=3, capthick=1.1, zorder=3))
        tops = []
        for xi, comp in zip(x, COMPARTMENTS):
            pts = per_sp[per_sp.Compartment == comp]
            for _, r in pts.iterrows():
                ax.scatter(xi + SP_OFF[r.Species], r[col], marker=SP_MARK[r.Species], s=28,
                           facecolor='white', edgecolor=INK, linewidth=1.0, zorder=5)
            tops.append(max(vals[xi] + sd[xi], pts[col].max()))
        span = max(tops)
        for xi, v, t in zip(x, vals, tops):
            ax.text(xi, t + span * 0.035, fmt % v, ha='center', va='bottom',
                    fontsize=FS - 3, fontweight='semibold', color=INK)
        ax.set_xticks(x); ax.set_xticklabels(COMPARTMENTS, fontsize=FS - 2)
        ax.set_ylim(0, span * 1.13)
        ax.set_title(title, fontsize=FS, fontweight='bold', color=INK, pad=8)
        ax.set_ylabel(ylab, fontsize=FS - 3, color='#41505c')
        ax.tick_params(axis='both', labelsize=FS - 3, colors='#41505c', length=0)
        for spn in ('top', 'right', 'left'):
            ax.spines[spn].set_visible(False)
        ax.spines['bottom'].set_color(AXIS)
    axes[1].set_ylim(0, 100)

    handles = [plt.Line2D([0], [0], marker=SP_MARK[s], color='none', markerfacecolor='white',
                          markeredgecolor=INK, markeredgewidth=1.0, markersize=7,
                          label=r'$\mathit{A.\ %s}$' % s) for s in SP_KEYS]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.030), ncol=4,
               frameon=False, fontsize=FS - 3, handletextpad=0.4, columnspacing=1.8)
    fig.suptitle('Pangenome protein length, domain coverage and gene-model quality',
                 fontsize=FS + 2, fontweight='bold', color=INK, y=1.045)
    plt.tight_layout(rect=[0, 0.02, 1, 1.0], w_pad=1.4)
    fig.subplots_adjust(top=0.888)   # tight_layout leaves slack at the top; close it explicitly
    out = NB4_RESULTS / 'compartment_length_pfam_orf.png'
    print(f'\nSaved -> {out}  and compartment_length_pfam_orf.csv\n')
    print(avg[['length', 'pfam', 'incomplete']].round(2).to_string())

    return fig, avg


# =============================================================================
# NB4 INPUT LOADING, PREFLIGHT AND SUMMARY
# =============================================================================

def load_nb4_species_data(species_list, genus, base_path, nb1_results,
                          ani_excluded=None):
    """Per-species PAV matrix, orthogroup consensus table and class thresholds.

    Orthogroups carried only by ANI or contamination excluded genomes are
    dropped. Returns species_data keyed by species.
    """
    import glob
    import funpan_pangenome as fpp
    import funpan_utils as fpu

    GENUS = genus
    NB1_RESULTS = str(nb1_results)
    SPECIES_LIST = species_list
    if ani_excluded is None:
        ani_excluded = fpu.ANI_EXCLUDED

    _, pav_data, _ = fpp.load_all_orthofinder(SPECIES_LIST, str(base_path), NB1_RESULTS)

    species_data = {}
    for sp in SPECIES_LIST:
        og_path = os.path.join(NB1_RESULTS, sp, f'{sp}_og_consensus.tsv')
        if not os.path.exists(og_path):
            print(f'  [WARN] No OG consensus for {sp} at {og_path}')
            continue
        ogc = pd.read_csv(og_path, sep='\t')
        _excl = set().union(*ani_excluded.values())
        _priv = set()
        for _f in glob.glob(f'{base_path}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups.tsv') + \
                  glob.glob(f'{base_path}/{sp}/orthofinder_output/*/Orthogroups/Orthogroups_UnassignedGenes.tsv'):
            with open(_f) as fh:
                _hdr = fh.readline().rstrip('\n').split('\t')
                _ex = {i for i, h in enumerate(_hdr) if any(a in h for a in _excl)}
                if not _ex:
                    continue
                for _l in fh:
                    _r = _l.rstrip('\n').split('\t')
                    _p = [i for i in range(1, len(_r)) if _r[i].strip()]
                    if _p and all(i in _ex for i in _p):
                        _priv.add(_r[0])
        if _priv:
            ogc = ogc[~ogc['Orthogroup'].isin(_priv)].copy()
            print(f'  {sp}: dropped {len(_priv)} OGs private to excluded genomes')

        pav = pav_data.get(sp)
        if pav is not None:
            _, core_n, _, rare_n = fpp.determine_core_and_rare_thresholds(pav)
        else:
            core_n = rare_n = None

        species_data[sp] = {
            'pav':          pav,
            'og_consensus': ogc,
            'n_strains':    pav.shape[1] if pav is not None else 0,
            'core_n':       core_n,
            'rare_n':       rare_n,
        }
        n_rare = (ogc['Pangenome_Class'] == 'Rare').sum() if 'Pangenome_Class' in ogc.columns else 0
        print(f'  {sp}: PAV {pav.shape if pav is not None else "--"}, '
              f'OG consensus {ogc.shape}, {n_rare} rare OGs, '
              f'core>={core_n}, rare<={rare_n}')
    return species_data


def apply_truly_rare_filter(species_list, species_data, species_config, results_dir):
    """Reclassify Rare orthogroups outside the DIAMOND No-hit class.

    Edits species_data in place so that Pangenome_Class == 'Rare' means the
    truly-rare subset from here on. Returns the No-hit orthogroup set per
    species.
    """
    RESULTS_DIR = str(results_dir)
    SPECIES_LIST = species_list
    SPECIES_CONFIG = species_config

    truly_rare_by_sp = {}
    for sp in SPECIES_LIST:
        p = os.path.join(RESULTS_DIR, sp, 'diamond_rare_classification.tsv')
        if not os.path.exists(p):
            print(f'{SPECIES_CONFIG[sp]["label"]}: no DIAMOND table at {p}')
            continue
        dc = pd.read_csv(p, sep='\t')
        truly_rare_by_sp[sp] = set(dc.loc[dc['blast_class'] == 'No hit', 'Orthogroup'])

    for sp in SPECIES_LIST:
        if sp not in truly_rare_by_sp:
            continue
        label = SPECIES_CONFIG[sp]['label']
        ogc = species_data[sp]['og_consensus']
        keep = truly_rare_by_sp[sp]
        rare_mask = ogc['Pangenome_Class'] == 'Rare'
        n_rare = int(rare_mask.sum())
        excl_mask = rare_mask & (~ogc['Orthogroup'].isin(keep))
        n_excl = int(excl_mask.sum())
        n_kept = n_rare - n_excl
        ogc.loc[excl_mask, 'Pangenome_Class'] = 'Rare_excluded'
        print(f'  {label:14s}: {n_rare:>5d} Rare -> {n_kept:>5d} truly-rare kept | '
              f'{n_excl:>5d} reclassified as Rare_excluded '
              f'({100*n_kept/max(n_rare,1):.1f} % kept)')
    return truly_rare_by_sp


def nb4_preflight(species_list, species_root, nb0_results, nb1_results,
                  nb4_results, qc_table=None, verbose=True):
    """Report which NB4 inputs are present and create the output directories.

    Returns a DataFrame with one row per checked input.
    """
    import shutil

    nb0, nb1, nb4 = Path(nb0_results), Path(nb1_results), Path(nb4_results)
    root = Path(species_root)
    rows = [
        ('phenotype table', nb0 / 'phenotype_classified_for_gwas.csv'),
        ('assembly metadata', Path(qc_table) if qc_table else nb0 / 'qc_passed_metadata_enriched.csv'),
    ]
    for sp in species_list:
        rows += [
            (f'{sp} PAV matrix', nb1 / sp / f'{sp}_pav.tsv'),
            (f'{sp} OG consensus', nb1 / sp / f'{sp}_og_consensus.tsv'),
            (f'{sp} kinship (GRM)', nb1 / sp / f'{sp}_kinship.tsv'),
            (f'{sp} SNP PCs', nb1 / sp / f'{sp}_snp_pcs.tsv'),
            (f'{sp} genome directory', root / sp),
        ]
    check = pd.DataFrame([{'input': n, 'path': str(p), 'exists': p.exists()} for n, p in rows])

    nb4.mkdir(parents=True, exist_ok=True)
    for sp in species_list:
        (nb4 / sp).mkdir(exist_ok=True)

    if verbose:
        missing = check[~check['exists']]
        print(f'Inputs present: {int(check["exists"].sum())}/{len(check)}')
        if len(missing):
            print('Missing:')
            for _, r in missing.iterrows():
                print(f'  {r["input"]:26s} {r["path"]}')
        for tool in ('diamond', 'mash', 'blastp'):
            print(f'  {tool:8s} {"on PATH" if shutil.which(tool) else "not on PATH"}')
        print(f'Output directory: {nb4}')
    return check


def print_nb4_summary(species_list, nb4_results):
    """Rare and truly-rare counts per species, and every NB4 output file."""
    nb4 = Path(nb4_results)
    rows = []
    for sp in species_list:
        p = nb4 / sp / 'diamond_rare_classification.tsv'
        if not p.exists():
            continue
        dc = pd.read_csv(p, sep='\t')
        rows.append({
            'Species': f'A. {sp}',
            'Rare OGs': len(dc),
            'Truly-rare (No hit)': int((dc['blast_class'] == 'No hit').sum()),
            'Paralog': int((dc['blast_class'] == 'Full-length homolog').sum()),
            'Fragment': int((dc['blast_class'] == 'Fragment of longer gene').sum()),
            'No significant similarity': int((dc['blast_class'] == 'No significant similarity').sum()),
        })
    counts = pd.DataFrame(rows)
    if len(counts):
        counts['% truly-rare'] = (100 * counts['Truly-rare (No hit)'] / counts['Rare OGs']).round(1)

    expected = [
        'rare_genome_summary.csv',
        'diamond_classification_pies.png',
        'rare_gene_burden_all_genomes.csv',
        'rare_gene_characterisation_summary.csv',
        'og_length_by_class.csv',
        'orf_completeness_by_class.csv',
        'compartment_length_pfam_orf.csv',
        'compartment_length_pfam_orf.png',
        'cog_enrichment_rare_vs_nonrare.csv',
        'rare_xenolog_ncbi_blast_results.tsv',
        'rare_xenolog_kingdom_summary.csv',
        'rare_xenolog_confident_hits.csv',
        'rare_xenolog_taxonomy.png',
        'xenolog_neighbourhood.png',
        'burden_contrast_table.csv',
        'truly_rare_burden_by_phenotype.png',
        'pooled_kinship_nb2.tsv',
        'combined_tree_vcv.tsv',
        'rare_burden_kinship_control.csv',
        'burden_bidirectional_kinship_corrected.png',
        'burden_oryzae_vs_flavus.csv',
        'burden_oryzae_vs_flavus.png',
        'burden_manuscript_four_panel.csv',
        'burden_manuscript_four_panel.png',
    ]
    per_species = [
        '{sp}/diamond_rare_classification.tsv',
        '{sp}/og_gc_content.tsv',
        '{sp}/phylogenetic_signal.png',
        '{sp}/ancestral_niche_reconstruction.png',
        '{sp}/gain_loss_events.tsv',
        '{sp}/gain_loss_summary.tsv',
        'mash_clustering/{sp}/{sp}_mash_distance_matrix.tsv',
        'mash_clustering/{sp}/{sp}_primary_clusters.csv',
    ]
    names = list(expected)
    for sp in species_list:
        names += [t.format(sp=sp) for t in per_species]
    names += ['mash_clustering/all_species/all_species_mash_distance_matrix.csv',
              'mash_clustering/all_species/all_species_mash_heatmap.png']

    files = pd.DataFrame([{'path': str(nb4 / n), 'exists': (nb4 / n).exists()}
                          for n in names])
    print(f'Rare-genome outputs under {nb4}')
    print(f'  present: {int(files["exists"].sum())}/{len(files)}')
    for _, r in files.iterrows():
        print(f'  {"OK " if r["exists"] else "-- "}{r["path"]}')
    return counts, files
