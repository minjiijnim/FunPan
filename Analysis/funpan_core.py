#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
funpan_core.py
==============
Core-compartment analysis for the FunPan Aspergillus pangenome project:
cross-species conservation of literature-curated virulence and industrial
gene panels.

Used by NB5_CoreGenome.ipynb.
"""

import os
import re
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu
from Bio import SeqIO

import funpan_utils as fpu
import funpan_pangenome as fpp
import funpan_convergence as fpc

_fpu = fpu

import warnings
warnings.filterwarnings('ignore')

# --- shared constants ---------------------------------------------------------
NB5_SPECIES = ['fumigatus', 'flavus', 'niger', 'oryzae']
NB5_SP_DISPLAY = {'fumigatus': 'A. fumigatus', 'flavus': 'A. flavus',
                  'niger': 'A. niger', 'oryzae': 'A. oryzae'}
NB5_CLASSES = ['Core', 'Accessory', 'Rare', 'Absent']
NB5_CLASS_COLOR = {'Core': '#2ca02c', 'Accessory': '#e67e22',
                   'Rare': '#c0392b', 'Absent': '#ffffff'}
NB5_CLASS_COLORS = {'Core': '#2ca02c', 'Soft-core': '#7fbf7b',
                    'Accessory': '#e67e22', 'Rare': '#c0392b', 'Absent': '#d9d9d9'}
NB5_COG_DESC = {
    'J': 'Translation', 'A': 'RNA processing', 'K': 'Transcription',
    'L': 'Replication/repair', 'B': 'Chromatin', 'D': 'Cell cycle',
    'Y': 'Nuclear structure', 'V': 'Defense', 'T': 'Signal transduction',
    'M': 'Cell wall/membrane', 'N': 'Cell motility', 'Z': 'Cytoskeleton',
    'W': 'Extracellular structures', 'U': 'Intracellular trafficking',
    'O': 'PTM/chaperones', 'C': 'Energy production', 'G': 'Carbohydrate metab.',
    'E': 'Amino acid metab.', 'F': 'Nucleotide metab.', 'H': 'Coenzyme metab.',
    'I': 'Lipid metab.', 'P': 'Inorganic ion transport', 'Q': 'Secondary metabolites',
    'R': 'General function', 'S': 'Function unknown',
}


def _show(obj):
    """Print a DataFrame or object; stands in for the notebook's display()."""
    try:
        print(obj.to_string())
    except AttributeError:
        print(obj)


def load_genus_pangenome(species_list, combined_genecount):
    """Per-species class and prevalence for every genus orthogroup.

    Reads the combined OrthoFinder gene-count table, drops the ANI and
    contamination excluded genome columns, and derives Core, Accessory,
    Rare and Absent per species from the S-curve thresholds.
    """
    SPECIES_LIST = list(species_list)
    COMBINED_GENECOUNT = str(combined_genecount)

    gc           = fpc.load_combined_orthogroups_genecount(str(COMBINED_GENECOUNT))
    combined_ogs = gc['Orthogroup'].values
    count_mat    = gc.set_index('Orthogroup')

    # --- Drop ANI/contamination-excluded genomes (combined run predates that QC) ---
    _excl = set().union(*fpu.ANI_EXCLUDED.values())
    _drop = [c for c in count_mat.columns if any(a in str(c) for a in _excl)]
    count_mat = count_mat.drop(columns=_drop)
    print(f'  Dropped {len(_drop)} excluded genome column(s): {_drop}')


    def _col_species(col):
        sp = col.split('__')[0]
        for suf in ('_ani_filtered', '_filtered', '_old'):
            if sp.endswith(suf):
                return sp[:-len(suf)]
        return sp
    col_species = {c: _col_species(c) for c in count_mat.columns}
    sp_cols = {sp: [c for c in count_mat.columns if col_species[c] == sp]
               for sp in SPECIES_LIST}

    genus_class      = pd.DataFrame(index=combined_ogs)         # OG x species -> class
    genus_present    = pd.DataFrame(index=combined_ogs)         # OG x species -> prevalence
    genus_thresholds = {}
    for sp in SPECIES_LIST:
        cols = sp_cols[sp]
        sub  = (count_mat[cols] > 0).astype(int)
        n    = len(cols)
        present = sub.sum(axis=1).values
        _, core_n, _, rare_n = fpp.determine_core_and_rare_thresholds(sub)
        genus_thresholds[sp] = {'n': n, 'core_n': int(core_n), 'rare_n': int(rare_n)}
        cls = np.full(len(present), 'Accessory', dtype=object)
        cls[present == 0]                          = 'Absent'
        cls[(present >= 1) & (present <= rare_n)]  = 'Rare'
        cls[present >= core_n]                     = 'Core'
        genus_class[sp]   = cls
        genus_present[sp] = present / n
        print(f'  A. {sp:10s}: n={n:3d}  core>={int(core_n):3d}  rare<={int(rare_n):2d}  |  '
              f'Core={int((cls=="Core").sum()):5d}  Accessory={int((cls=="Accessory").sum()):5d}  '
              f'Rare={int((cls=="Rare").sum()):5d}  Absent={int((cls=="Absent").sum()):5d}')

    return {'genus_class': genus_class, 'genus_present': genus_present,
            'genus_thresholds': genus_thresholds, 'count_mat': count_mat,
            'sp_cols': sp_cols, 'combined_ogs': combined_ogs}


def core_compartment_summary(species_list, og_consensus, n_strains, results_dir):
    """Core, accessory and rare counts per species, with the strict-core split.
    """
    SPECIES_LIST = list(species_list)
    RESULTS_DIR = str(results_dir)

    rows = []
    for sp in SPECIES_LIST:
        df   = og_consensus[sp]
        N    = n_strains[sp]
        core = df[df['Pangenome_Class'] == 'Core']
        strict = int((core['n_genomes'] >= N).sum())
        rows.append({
            'Species': f'A. {sp}', 'n_strains': N, 'Total OGs': len(df),
            'Core': len(core), 'Core %': round(100*len(core)/len(df), 1),
            'Strict-core': strict, 'Soft-core': len(core)-strict,
            'Accessory': int((df['Pangenome_Class']=='Accessory').sum()),
            'Rare': int((df['Pangenome_Class']=='Rare').sum()),
            'Core PFAM-annot %': round(100*core['PFAMs'].notna().mean(), 1),
        })
    core_summary = pd.DataFrame(rows)
    core_summary.to_csv(os.path.join(RESULTS_DIR, 'core_compartment_summary.csv'), index=False)
    _show(core_summary)

    return core_summary


def load_curated_panels(panel_dir):
    """Curated literature panels and their reference protein sequences.

    Fetches the reference proteins on first use if the FASTA is missing.
    """
    PANEL_DIR = str(panel_dir)

    PANELS = {'fumigatus': 'fumigatus_virulence_panel.csv',
              'flavus':    'flavus_virulence_panel.csv',
              'niger':     'niger_industrial_panel.csv',
              'oryzae':    'oryzae_industrial_panel.csv'}
    PANEL_TYPE = {'fumigatus': 'virulence', 'flavus': 'virulence',
                  'niger': 'industrial', 'oryzae': 'industrial'}

    panels = {}
    for sp, fn in PANELS.items():
        p = pd.read_csv(os.path.join(PANEL_DIR, fn))
        p['species'] = sp
        panels[sp] = p
        print(f'  A. {sp:10s}: {len(p):2d} {PANEL_TYPE[sp]} genes  ({p["category"].nunique()} categories)')
    panel_all = pd.concat(panels.values(), ignore_index=True)

    SEQ_CACHE = os.path.join(PANEL_DIR, 'reference_proteins.faa')
    if not os.path.exists(SEQ_CACHE):
        print('\nFetching reference proteins ...')
        subprocess.run([sys.executable, os.path.join(PANEL_DIR, 'fetch_reference_proteins.py')],
                       check=True)
    ref_records = {}
    for rec in SeqIO.parse(SEQ_CACHE, 'fasta'):
        parts = rec.id.split('|')
        if len(parts) >= 2:
            ref_records[(parts[0], parts[1])] = rec
    n_seq = sum(1 for sp, g in zip(panel_all['species'], panel_all['gene'])
                if (sp, g) in ref_records)
    print(f'\nReference sequences available: {n_seq}/{len(panel_all)} '
          f'({100*n_seq/len(panel_all):.0f}%)')

    return {'panels': panels, 'panel_all': panel_all,
            'ref_records': ref_records, 'seq_cache': SEQ_CACHE,
            'panel_type': PANEL_TYPE}


def anchor_panels_to_orthogroups(species_list, panel_all, genus_class_df, panel_type, seq_cache, cache_dir, results_dir, combined_og_tsv, ref_prot_dir, reference_genomes, force_blast=False):
    """Anchor each panel reference protein to a genus orthogroup by blastp.

    Results are cached in ``{cache_dir}/panel_blastp.tsv``. Requires
    blastp and makeblastdb on PATH or in the interpreter's directory.
    """
    FORCE_BLAST = force_blast
    SPECIES_LIST = list(species_list)
    genus_class = genus_class_df
    PANEL_TYPE = panel_type
    SEQ_CACHE = str(seq_cache)
    CACHE_DIR = str(cache_dir)
    RESULTS_DIR = str(results_dir)
    COMBINED_OG_TSV = str(combined_og_tsv)
    REF_PROT_DIR = Path(ref_prot_dir)
    REFERENCE_GENOMES = reference_genomes

    import shutil
    BLAST_CACHE = os.path.join(CACHE_DIR, 'panel_blastp.tsv')
    NTHREADS    = max(1, (os.cpu_count() or 4))

    def _resolve(binary):
        # resolve a BLAST+ executable to an absolute path (PATH-independent)
        p = shutil.which(binary)
        if p:
            return p
        for cand in (os.path.dirname(sys.executable), os.path.expanduser('~/anaconda3/envs/funpan/bin')):
            fp = os.path.join(cand, binary)
            if os.path.exists(fp):
                return fp
        return binary
    MAKEBLASTDB, BLASTP = _resolve('makeblastdb'), _resolve('blastp')

    seqs_by_sp = defaultdict(list)
    for rec in SeqIO.parse(SEQ_CACHE, 'fasta'):
        seqs_by_sp[rec.id.split('|')[0]].append(rec)

    BCOLS = ['qseqid','sseqid','pident','length','qlen','slen','evalue','bitscore']
    if (not FORCE_BLAST) and os.path.exists(BLAST_CACHE):
        blast = pd.read_csv(BLAST_CACHE, sep='\t')
        print(f'Loaded cached blastp results ({len(blast)} hits).')
    else:
        frames = []
        tmp = tempfile.mkdtemp(prefix='nbx_blast_')
        try:
            for sp in SPECIES_LIST:
                if not seqs_by_sp[sp]:
                    continue
                ref_fa = str(REF_PROT_DIR / (REFERENCE_GENOMES[sp] + '.fa'))
                query  = os.path.join(tmp, f'{sp}_q.fa')
                db     = os.path.join(tmp, f'{sp}_db')
                out    = os.path.join(tmp, f'{sp}_o.tsv')
                SeqIO.write(seqs_by_sp[sp], query, 'fasta')
                subprocess.run([MAKEBLASTDB,'-in',ref_fa,'-dbtype','prot','-out',db],
                               check=True, capture_output=True)
                subprocess.run([BLASTP,'-query',query,'-db',db,'-out',out,
                                '-outfmt','6 '+' '.join(BCOLS),'-max_target_seqs','5',
                                '-evalue','1e-5','-num_threads',str(NTHREADS)],
                               check=True, capture_output=True)
                if os.path.getsize(out) > 0:
                    frames.append(pd.read_csv(out, sep='\t', names=BCOLS))
                print(f'  A. {sp:10s}: blastp {len(seqs_by_sp[sp])} queries vs reference genome')
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        blast = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=BCOLS)
        blast.to_csv(BLAST_CACHE, sep='\t', index=False)

    blast['qcov'] = 100 * blast['length'] / blast['qlen']
    best = (blast.sort_values('bitscore', ascending=False)
                 .drop_duplicates('qseqid', keep='first').copy())

    ref_cols = list(REFERENCE_GENOMES.values())
    og_ref = pd.read_csv(COMBINED_OG_TSV, sep='\t', usecols=['Orthogroup'] + ref_cols)
    prot2og = {sp: {} for sp in SPECIES_LIST}
    for _, r in og_ref.iterrows():
        og = r['Orthogroup']
        for sp, col in REFERENCE_GENOMES.items():
            v = r[col]
            if isinstance(v, str):
                for prot in v.split(', '):
                    prot2og[sp][prot.strip()] = og

    recs = []
    for _, row in best.iterrows():
        parts = str(row['qseqid']).split('|')
        sp, gene = parts[0], parts[1]
        recs.append({'species': sp, 'gene': gene,
                     'pident': round(row['pident'],1), 'qcov': round(row['qcov'],1),
                     'combined_og': prot2og[sp].get(row['sseqid'])})
    anchor = pd.DataFrame(recs).merge(
        panel_all[['species','gene','category','description']], on=['species','gene'], how='left')
    for s in SPECIES_LIST:
        anchor[f'class_{s}'] = anchor['combined_og'].map(genus_class[s])
    anchor['native_class'] = anchor.apply(lambda r: r.get(f'class_{r["species"]}'), axis=1)
    anchor['kind'] = anchor['species'].map(PANEL_TYPE)
    anchor.to_csv(os.path.join(RESULTS_DIR, 'panel_orthogroup_anchoring.csv'), index=False)
    n_anch = int(anchor['combined_og'].notna().sum())
    print(f'\nAnchored {n_anch}/{len(panel_all)} panel genes to a genus orthogroup.')

    return anchor


def derive_validated_panel(results_dir, filtered_dir):
    """Validated, positive-causation panel from the PMID audit.

    Keeps genes whose cited work demonstrates the trait at the gene level
    in that species, then drops silenced relics, suppressors and
    repressors. Writes panel_validated_poscausation.csv.
    """
    RES = Path(results_dir)
    FILT = Path(filtered_dir)

    import pandas as pd
    def _epi(s):
        s=str(s).lower().replace('aspergillus','').replace('a.','').strip()
        return s.split()[0] if s.split() else s
    _aud=pd.read_csv(RES/'panel_pmid_audit.csv'); _sup=pd.read_csv(FILT/'panel_cross_species_conservation_SUP.csv')
    _aud['g']=_aud['gene'].str.lower(); _aud['s']=_aud['species'].apply(_epi)
    _sup['g']=_sup['gene'].str.lower(); _sup['s']=_sup['species'].apply(_epi)
    _keep=set(zip(_aud.loc[_aud['final']=='supported','g'],_aud.loc[_aud['final']=='supported','s']))
    VAL=_sup[_sup.set_index(['g','s']).index.isin(_keep)].copy(); VAL['s']=VAL['species'].apply(_epi); VAL['gl']=VAL['gene'].str.lower()
    # Positive-causation filter: drop silenced relics, virulence suppressors, and CCR repressors; toxins count as virulence
    _excl={('oryzae','afld'),('oryzae','aflp'),('oryzae','aflr'),('oryzae','aflt'),('fumigatus','pes3'),('oryzae','crea'),('niger','crea'),('oryzae','creb')}
    _recl={('niger','fum21'),('niger','otaa')}
    VAL['kind_final']=VAL.apply(lambda r:'virulence' if (r['s'],r['gl']) in _recl else r['kind'],axis=1)
    VAL=VAL[~VAL.apply(lambda r:(r['s'],r['gl']) in _excl,axis=1)].reset_index(drop=True)
    VAL['category']=VAL.apply(lambda r:{'otaa':'Ochratoxin','fum21':'Fumonisin','laea':'Secondary metabolite regulator','hsba':'Polyester degradation','dmta':'Manganese transporter','agse':'Cell wall/glucan','rho1':'Cell wall/glucan','chsa':'Cell wall/chitin','mnt1':'Cell wall/mannan','glfa':'Cell wall/galactofuranose'}.get(r['gl'],r['category']),axis=1)  # name the mycotoxins instead of 'Other-SM'
    VAL['category']=VAL['category'].replace({'GAG':'Galactosaminogalactan','Secretion/UPR':'Unfolded protein response','UPR/secretion':'Unfolded protein response','pH/stress':'pH stress','Melanin/conidia':'Melanin biosynthesis','Velvet/development':'Velvet complex','Conidia/development':'Conidiation','Protease':'Aspartic protease','Amino acid biosynth':'Biosynthesis (auxotrophy)'})  # spell out abbreviations, unify UPR
    VAL.to_csv(FILT/'panel_validated_poscausation.csv',index=False)
    print(f"Validated positive-causation panel: {len(VAL)} genes, {VAL['combined_og'].nunique()} OGs")
    print("  kind_final:", VAL['kind_final'].value_counts().to_dict(), "| by species:", VAL['s'].value_counts().to_dict())
    print(f"  core in all 4: {(VAL['n_core']==4).sum()}/{len(VAL)} = {100*(VAL['n_core']==4).mean():.0f}%")

    return VAL


def plot_validated_conservation(VAL, genus_class_df, results_dir):
    """Number of species in which a validated gene is core, against background.

    Returns the figure and the one-sided Mann-Whitney comparison of the
    panel against the orthogroups that are not in it.
    """
    genus_class = genus_class_df
    RES = Path(results_dir)

    import matplotlib.pyplot as plt, numpy as np, pandas as pd
    core_in=pd.DataFrame({s:(genus_class[s]=='Core').astype(int) for s in ['fumigatus','flavus','niger','oryzae']})
    bg_ncore=core_in.sum(axis=1); bg_any=bg_ncore.loc[~bg_ncore.index.isin(set(VAL['combined_og'].dropna().unique()))]
    panel_any=VAL['n_core'].astype(int)
    BINS=[4,3,2,1,0]; LABELS=['All 4','3','2','1','None']; INK,MUTE,GRID='#1a1d21','#6b7280','#e6e8eb'; C_BG,C_PANEL='#c2c7cf','#2f9e44'
    xs=np.arange(len(BINS)); w=0.28
    fig,ax=plt.subplots(figsize=(4.6,6.7)); fig.patch.set_facecolor('white'); ax.set_facecolor('white')
    for k,(lab,data,color) in enumerate([('Genome background',bg_any,C_BG),('Trait panel',panel_any,C_PANEL)]):
        vals=pd.Series(list(data)).value_counts().reindex(BINS,fill_value=0); pct=100*vals/vals.sum(); off=(k-0.5)*w
        ax.bar(xs+off,pct.values,w,color=color,edgecolor='white',linewidth=1.3,zorder=3,label=f'{lab}  (n={len(data):,})')
        for x,v in zip(xs,pct.values):
            if v>=1: ax.text(x+off,v+1.4,f'{v:.0f}',ha='center',va='bottom',fontsize=8.5,color=MUTE)
    ax.set_axisbelow(True); ax.yaxis.grid(True,color=GRID,linewidth=1,zorder=0); ax.xaxis.grid(False)
    for _s in ('top','right'): ax.spines[_s].set_visible(False)
    for _s in ('left','bottom'): ax.spines[_s].set_color(MUTE); ax.spines[_s].set_linewidth(0.8)
    ax.tick_params(colors=MUTE,length=0,labelsize=10); ax.set_xticks(xs); ax.set_xticklabels(LABELS,fontsize=10.5,color=INK)
    ax.set_xlabel('Number of species in which the gene is core',fontsize=10.5,color=INK,labelpad=9)
    ax.set_ylabel('% of orthogroups',fontsize=10.5,color=INK,labelpad=9); ax.set_ylim(0,100); ax.set_yticks(range(0,101,20))
    for _t in ax.get_yticklabels(): _t.set_color(MUTE)
    ax.set_title('Cross-genus conservation of trait genes',fontsize=12.5,fontweight='bold',color=INK,loc='left',pad=16)
    ax.legend(loc='upper right',frameon=False,fontsize=9.5,handlelength=1.1,bbox_to_anchor=(1.0,1.0),labelcolor=INK)
    plt.tight_layout()
    outp = RES / 'fig5a_conservation.png'
    fig.savefig(outp, dpi=400, bbox_inches='tight')
    print(f"Fig 5.a saved (panel n={len(panel_any)}, background n={len(bg_any):,})")

    # One-sided Mann-Whitney U: is the panel more conserved than the background?
    U, p_cons = mannwhitneyu(panel_any, bg_any, alternative='greater')
    stats = {'n_panel': int(len(panel_any)), 'n_background': int(len(bg_any)),
             'panel_core_in_all4': int((panel_any == 4).sum()),
             'pct_panel_core_in_all4': round(100 * (panel_any == 4).mean(), 1),
             'pct_background_core_in_all4': round(100 * (bg_any == 4).mean(), 1),
             'mannwhitney_U': float(U), 'p_value': float(p_cons)}
    print(f"  panel core in all 4      : {stats['panel_core_in_all4']}/{stats['n_panel']} "
          f"({stats['pct_panel_core_in_all4']:.0f}%)")
    print(f"  background core in all 4 : {stats['pct_background_core_in_all4']:.0f}% "
          f"of {stats['n_background']:,}")
    print(f"  Mann-Whitney (panel > background): p = {stats['p_value']:.1e}")

    return fig, stats


def plot_validated_combined(VAL, results_dir):
    """Validated genes variable across the genus, over those shared between
    A. flavus and A. oryzae.
    """
    RES = Path(results_dir)

    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.patches import Patch as _Patch
    import matplotlib.gridspec as _gs
    cls_to_num={'Core':3,'Accessory':2,'Rare':1,'Absent':0}
    INK='#1a1d21'; PAL=['#f0f0f0','#c0392b','#e67e22','#2ca02c']; GAP='#b9bec5'; DISC='#1971c2'
    BOX_W=1.15; COL_DX=1.30; PAD_Y=0.14; _M=0.05
    def render(ax,d,spcols,splab,title,rowmax,ylfs):
        n=len(d); xs=[j*COL_DX for j in range(len(spcols))]
        M=np.array([[cls_to_num.get(r[f'class_{sp}'],0) for sp in spcols] for _,r in d.iterrows()],dtype=int) if n else np.zeros((0,len(spcols)),int)
        for i in range(n):
            for j in range(len(spcols)):
                ax.add_patch(plt.Rectangle((xs[j]-BOX_W/2,i-0.5+PAD_Y),BOX_W,1-2*PAD_Y,facecolor=PAL[M[i,j]],edgecolor=GAP,linewidth=0.6))
        for i,r in enumerate(d.itertuples()):          # blue box on the species the gene was curated from
            for sp in r.origins:
                if sp in spcols:
                    j=spcols.index(sp)
                    _x0,_y0=xs[j]-BOX_W/2,i-0.5+PAD_Y   # stroke on exact cell, clipped to it -> grows inward, outer edge aligned
                    _blue=plt.Rectangle((_x0,_y0),BOX_W,1-2*PAD_Y,fill=False,edgecolor=DISC,linewidth=7,zorder=6); ax.add_patch(_blue)
                    _blue.set_clip_path(plt.Rectangle((_x0,_y0),BOX_W,1-2*PAD_Y,transform=ax.transData))
        ax.set_xlim(xs[0]-BOX_W/2-0.2,xs[-1]+BOX_W/2+0.2); ax.set_ylim(rowmax-0.5,-0.5)
        ax.set_xticks(xs); ax.set_xticklabels(splab,fontsize=11,style='italic',color=INK,rotation=35,ha='right')
        ax.set_yticks(range(n)); ax.set_yticklabels([f'{r.gene} [{r.category}]' for r in d.itertuples()],fontsize=ylfs)
        ax.set_ylabel('Trait gene',fontsize=12,color=INK)
        ax.yaxis.set_label_coords(-0.65, 0.6)
        for _sp in ax.spines.values(): _sp.set_visible(False)
        ax.tick_params(length=0); ax.set_title(title,fontsize=13,fontweight='bold',color=INK,pad=6)

    subB=VAL[VAL['s'].isin(['flavus','oryzae'])].copy()
    Bg={k: subB[subB['kind_final']==k].sort_values(['category','gene']).reset_index(drop=True) for k in ['virulence','industrial']}
    for d in Bg.values(): d['origins']=d['s'].apply(lambda x:[x])
    _d=VAL.copy(); _d['gl']=_d['gene'].str.lower()
    SPC=['fumigatus','flavus','niger','oryzae']; SPL4=['A. fumigatus','A. flavus','A. niger','A. oryzae']
    for sp in SPC: _d[f'n_{sp}']=_d[f'class_{sp}'].map(cls_to_num)
    _agg={f'n_{sp}':'max' for sp in SPC}; _agg.update({'gene':'first','category':'first'})
    Dd=_d.groupby(['gl','kind_final']).agg(_agg).reset_index()
    _orig=_d.groupby(['gl','kind_final'])['s'].apply(lambda x:sorted(set(x))).reset_index(name='origins')
    Dd=Dd.merge(_orig,on=['gl','kind_final'])
    Dd['ncore4']=sum((Dd[f'n_{sp}']==3).astype(int) for sp in SPC)
    Dd=Dd[Dd['ncore4']!=4].copy()
    for sp in SPC: Dd[f'class_{sp}']=Dd[f'n_{sp}'].map({0:'Absent',1:'Rare',2:'Accessory',3:'Core'})
    Cg={k: Dd[Dd['kind_final']==k].sort_values(['ncore4','category','gene'],ascending=[False,True,True]).reset_index(drop=True) for k in ['virulence','industrial']}
    Cmax=max(len(d) for d in Cg.values()); nVir=len(Bg['virulence'])

    fig=plt.figure(figsize=(13,0.46*(2*Cmax)+2.6)); fig.patch.set_facecolor('white')
    gs=_gs.GridSpec(2,2,height_ratios=[1,1],width_ratios=[1,1],hspace=0.40,wspace=0.85)
    axAv=fig.add_subplot(gs[0,0]); render(axAv,Cg['virulence'],SPC,SPL4,f'Virulence  (n={len(Cg["virulence"])})',Cmax,11.5)
    axAi=fig.add_subplot(gs[0,1]); render(axAi,Cg['industrial'],SPC,SPL4,f'Industrial  (n={len(Cg["industrial"])})',Cmax,11.5)
    axBv=fig.add_subplot(gs[1,0]); render(axBv,Bg['virulence'],['flavus','oryzae'],['A. flavus','A. oryzae'],f'Virulence  (n={nVir})',nVir,11.5)
    axBi=fig.add_subplot(gs[1,1]); render(axBi,Bg['industrial'],['flavus','oryzae'],['A. flavus','A. oryzae'],f'Industrial  (n={len(Bg["industrial"])})',len(Bg['industrial']),11.5)
    fig.canvas.draw()
    _ya=max(axAv.get_position().y1,axAi.get_position().y1); _yb=max(axBv.get_position().y1,axBi.get_position().y1)
    fig.text(0.5,min(_ya+0.035,0.99),'Validated trait genes variable across all four species',ha='center',fontsize=14,fontweight='bold',color=INK)
    fig.text(0.5,_yb+0.035,'Validated trait genes shared between A. flavus and A. oryzae',ha='center',fontsize=14,fontweight='bold',color=INK)
    hs=[_Patch(facecolor=PAL[v],label=l,edgecolor=GAP,linewidth=0.6) for l,v in [('Core',3),('Accessory',2),('Rare',1),('Absent',0)]]
    hs.append(_Patch(facecolor='white',edgecolor=DISC,linewidth=1.7,label='Species of characterization'))
    fig.legend(handles=hs,loc='lower center',ncol=5,frameon=False,fontsize=11.5,bbox_to_anchor=(0.5,-0.005),labelcolor=INK)
    outp = RES / 'fig5bc_combined.png'
    fig.savefig(outp, dpi=350, bbox_inches='tight')
    print("Fig 5.bc combined saved (+ curated-species blue boxes)")

    return fig


# =============================================================================
# PREFLIGHT AND SUMMARY
# =============================================================================

def nb5_preflight(species_list, species_root, nb1_results, nb5_results,
                  combined_of, panel_dir, verbose=True):
    """Report which NB5 inputs are present and create the output directories.

    Returns a DataFrame with one row per checked input.
    """
    nb1, nb5 = Path(nb1_results), Path(nb5_results)
    comb, panels = Path(combined_of), Path(panel_dir)
    rows = [
        ('combined gene counts', comb / 'Orthogroups' / 'Orthogroups.GeneCount.tsv'),
        ('combined orthogroups', comb / 'Orthogroups' / 'Orthogroups.tsv'),
        ('reference proteomes', Path(species_root) / 'all_combined' /
         'orthofinder_output' / 'input_proteins'),
        ('reference protein FASTA', panels / 'reference_proteins.faa'),
    ]
    for sp in species_list:
        rows.append((f'{sp} OG consensus', nb1 / sp / f'{sp}_og_consensus.tsv'))
        rows.append((f'{sp} curated panel', panels /
                     (f'{sp}_virulence_panel.csv' if sp in ('fumigatus', 'flavus')
                      else f'{sp}_industrial_panel.csv')))
    check = pd.DataFrame([{'input': n, 'path': str(p), 'exists': p.exists()}
                          for n, p in rows])

    for d in (nb5, panels, nb5 / 'cache', nb5 / 'filtered'):
        d.mkdir(parents=True, exist_ok=True)

    if verbose:
        missing = check[~check['exists']]
        print(f'Inputs present: {int(check["exists"].sum())}/{len(check)}')
        if len(missing):
            print('Missing:')
            for _, r in missing.iterrows():
                print(f'  {r["input"]:26s} {r["path"]}')
        for tool in ('blastp', 'makeblastdb'):
            print(f'  {tool:12s} {"on PATH" if shutil.which(tool) else "not on PATH"}')
        print(f'Output directory: {nb5}')
    return check


def print_nb5_summary(nb5_results):
    """Panel and conservation counts, and every NB5 output file."""
    nb5 = Path(nb5_results)
    rows = []

    def _add(label, path, fn):
        p = nb5 / path
        if not p.exists():
            return
        try:
            rows.append({'table': label, **fn(pd.read_csv(p))})
        except Exception:
            pass

    _add('anchored panel', 'panel_cross_species_conservation.csv',
         lambda d: {'genes': len(d), 'orthogroups': d['combined_og'].nunique(),
                    'core in all 4': int((d['n_core'] == 4).sum()),
                    '% core in all 4': round(100 * (d['n_core'] == 4).mean(), 1)})
    _add('SUP-filtered panel', 'filtered/panel_cross_species_conservation_SUP.csv',
         lambda d: {'genes': len(d), 'orthogroups': d['combined_og'].nunique(),
                    'core in all 4': int((d['n_core'] == 4).sum()),
                    '% core in all 4': round(100 * (d['n_core'] == 4).mean(), 1)})
    _add('validated panel', 'filtered/panel_validated_poscausation.csv',
         lambda d: {'genes': len(d), 'orthogroups': d['combined_og'].nunique(),
                    'core in all 4': int((d['n_core'] == 4).sum()),
                    '% core in all 4': round(100 * (d['n_core'] == 4).mean(), 1)})
    counts = pd.DataFrame(rows)

    expected = [
        'core_compartment_summary.csv',
        'core_characterization.png',
        'panel_orthogroup_anchoring.csv',
        'panel_cross_species_conservation.csv',
        'panel_cross_species_heatmap.png',
        'conservation_spectrum.png',
        'panel_class_by_category_bw.png',
        'cross_species_conservation_matrix.csv',
        'cross_species_conservation_matrix.png',
        'fig_optionA_stacked_bars.png',
        'fig_optionB_donut_table.png',
        'fig_optionC_spectrum_exceptions.png',
        'fig_optionD_upset.png',
        'shared_toolkit_fumigatus_flavus.csv',
        'shared_toolkit_niger_oryzae.csv',
        'shared_toolkits.png',
        'flavus_oryzae_core_sharing.csv',
        'aflatoxin_cluster_oryzae.csv',
        'aflatoxin_cluster_oryzae.png',
        'flavus_oryzae_panel_pav.png',
        'panel_SUP_class_distribution.csv',
        'panel_SUP_frac_core_matrix.png',
        'panel_SUP_class_by_group.png',
        'panel_SUP_typeB_spotlight.png',
        'panel_SUP_integrated_frac_core.png',
        'panel_SUP_conservation_spectrum.png',
        'panel_SUP_flavus_oryzae_pav.png',
        'panel_SUP_class_by_category_bw.png',
        'filtered/panel_validated_poscausation.csv',
        'fig5a_conservation.png',
        'fig5bc_combined.png',
        'panels/reference_proteins.faa',
        'cache/panel_blastp.tsv',
    ]
    files = pd.DataFrame([{'path': str(nb5 / n), 'exists': (nb5 / n).exists()}
                          for n in expected])
    print(f'Core-genome outputs under {nb5}')
    print(f'  present: {int(files["exists"].sum())}/{len(files)}')
    for _, r in files.iterrows():
        print(f'  {"OK " if r["exists"] else "-- "}{r["path"]}')
    return counts, files
