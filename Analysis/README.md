# FunPan Analysis Notebooks

Six notebooks. Each reads the pipeline outputs under `Species/` and writes to its
own `NB*_Results/` folder.

See also: [project overview](../README.md), [shell pipeline](../Code/README.md),
[Python modules](MODULES.md).

## Running

```bash
conda activate funpan
cd /datadrive/Code
bash run_analysis_notebooks.sh            # all six
bash run_analysis_notebooks.sh --from NB2 # start at NB2
bash run_analysis_notebooks.sh --only NB3 # one notebook
```

Outputs are saved into the notebooks in place. To work interactively:

```bash
cd /datadrive/Analysis && jupyter notebook NB0_DataPrep.ipynb
```

## Order

```
NB0 -> NB1 -> NB2 -> NB3
               \
                -> NB4
                -> NB5
```

NB3 also needs the combined four-species OrthoFinder run
(`orthofinder_on_all_species.sh`).

## Section structure

Section 1 configures; the last section checks every expected output file.
Sections 2 onward open with a bootstrap cell that sets paths, imports the modules
and reloads that section's inputs, so a section runs from a fresh kernel without
the sections above it.

Every bootstrap cell starts with the path block:

```python
# --- Paths: edit these to point at your own data ---
ANALYSIS_DIR = Path('/datadrive/Analysis')
SPECIES_DIR  = Path('/datadrive/Species')
GENUS        = 'Aspergillus'
```

## The notebooks

### `NB0_DataPrep.ipynb`

| Section | |
|---|---|
| 2 | Metadata, quality control, ANI species verification |
| 3 | Two-tier phenotype classification from 8 metadata fields |
| 4 | Geographic and phenotype distribution figures |
| 5 | Per-species phenotype files for association testing |

### `NB1_Pangenome.ipynb`

| Section | |
|---|---|
| 2 | PAV and CNV matrices |
| 3 | Core/accessory/rare by S-curve inflection |
| 4 | Heap's law with permutation confidence intervals |
| 5 | Orthogroup annotation consensus tables |
| 6 | Functional enrichment per compartment |
| 7 | Kinship matrices from recombination-filtered SNPs |
| 8 | BGC and GCF matrices |

### `NB2_PanGWAS.ipynb`

| Section | |
|---|---|
| 2 | One-vs-rest phenotype contrasts |
| 3 | EMMA-style LMM per species and feature layer |
| 4 | Power analysis, minimum detectable OR at 80% power |
| 5 | Functional enrichment of significant orthogroups |
| 6 | BGC co-localisation |
| 7 | Cross-species comparison |

### `NB3_Convergence.ipynb`

| Section | |
|---|---|
| 2 | Species-to-genus orthogroup mapping |
| 3 | Convergence of the significant sets |
| 4 | Directional convergence across jointly tested orthogroups |
| 5 | Permutation null model |
| 6 | CAZy, protease and secretome comparisons |
| 7 | Functional convergence of enriched terms |

### `NB4_RareGenome.ipynb`

| Section | |
|---|---|
| 2 | DIAMOND reclassification of rare orthogroups |
| 3 | Restriction to the No-hit subset |
| 4 | Protein length, Pfam coverage, ORF completeness |
| 5 | COG enrichment |
| 6 | Xenolog screen on GC outliers |
| 7 | Phylogenetic signal and ancestral reconstruction |
| 8 | Mash clustering |
| 9 | Kinship-corrected burden models |

### `NB5_CoreGenome.ipynb`

| Section | |
|---|---|
| 2 | Curated panels anchored to genus orthogroups by blastp |
| 3 | PMID-validated subset, conservation against background |

## Result folders

```
NB0_Results/    Metadata, phenotype assignments, geographic map
NB1_Results/    PAV/CNV matrices, OG consensus, kinship, BGC/GCF matrices
NB2_Results/    Association results, enrichment, power analysis
NB3_Results/    Convergence results, null model, family comparisons
NB4_Results/    Rare-genome characterisation, phylogenetics, Mash clustering
NB5_Results/    Core-genome conservation figures and tables
```

`NB*_Results/` is generated and untracked, except the curated inputs under
`NB5_Results/`, which no code regenerates:

```
NB5_Results/panels/*_panel.csv                                # trait panels, one per species
NB5_Results/panel_pmid_audit.csv                              # per-gene PMID verification
NB5_Results/filtered/panel_cross_species_conservation_SUP.csv
```

## Caches

Expensive steps cache and reuse. Section 1 of each notebook reports which inputs
are present and creates the output directories.

| Step | Cache | Override |
|---|---|---|
| NB4 DIAMOND | `NB4_Results/{species}/diamond_rare_classification.tsv` | `run_diamond` |
| NB4 GC content | `NB4_Results/{species}/og_gc_content.tsv` | `FORCE_GC` |
| NB4 NCBI BLAST | `NB4_Results/rare_xenolog_ncbi_blast_results.tsv` | `FORCE_BLAST` |
| NB5 blastp | `NB5_Results/cache/panel_blastp.tsv` | `force_blast` |

## Binaries needed on PATH

`diamond` (NB4 §2), `mash` (NB4 §8), `blastp` and `makeblastdb` (NB5 §2). NB4 §6
needs network access for NCBI BLAST.
