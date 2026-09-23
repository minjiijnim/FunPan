# FunPan Analysis Notebooks

Six Jupyter notebooks covering the downstream analysis, from quality-filtered
genomes to the comparative results. Each notebook reads the pipeline outputs
under `Species/` and writes to its own `NB*_Results/` folder.

For project overview and environment setup see the [main README](../README.md).
For the shell pipeline that produces the inputs see
[Code/README.md](../Code/README.md). For the Python modules the notebooks import
see [MODULES.md](MODULES.md).

---

## Running the notebooks

All six in sequence, from the command line:

```bash
conda activate funpan
cd /datadrive/Code
bash run_analysis_notebooks.sh
```

Options:

- `bash run_analysis_notebooks.sh --from NB2` — skip NB0/NB1, start at NB2
- `bash run_analysis_notebooks.sh --only NB3` — run NB3 alone
- `bash run_analysis_notebooks.sh --help` — usage

The script shows a progress bar and saves executed outputs in place.

Or open one interactively:

```bash
conda activate funpan
cd /datadrive/Analysis
jupyter notebook NB0_DataPrep.ipynb
```

## Order and dependencies

Run NB0 first: everything else depends on its phenotype assignments. Then NB1,
which builds the matrices NB2 to NB5 consume. NB2 to NB5 can then run in any
order, except that NB3 needs NB2's association results.

```
NB0  ->  NB1  ->  NB2  ->  NB3
                   \
                    ->  NB4
                    ->  NB5
```

NB3 also needs the combined four-species OrthoFinder run from
`orthofinder_on_all_species.sh`.

## Section structure

Every notebook is organised the same way. Section 1 is configuration; the last
section is a summary that checks every expected output file. Sections 2 onward
each open with a **standalone bootstrap cell** that sets the paths, imports the
modules and reloads that section's inputs, so any section runs on its own from a
fresh kernel without executing the sections above it.

Paths live in an editable block at the top of every bootstrap cell:

```python
# --- Paths: edit these to point at your own data ---
ANALYSIS_DIR = Path('/datadrive/Analysis')
SPECIES_DIR  = Path('/datadrive/Species')
GENUS        = 'Aspergillus'
```

Change those lines to point the notebook at your own directories.

## The notebooks

### `NB0_DataPrep.ipynb`

Metadata assembly, quality control and phenotype assignment.

| Section | Contents |
|---|---|
| 2 | Metadata and quality control, ANI species verification |
| 3 | Two-tier phenotype classification from eight metadata fields |
| 4 | Geographic distribution and phenotype distribution figures |
| 5 | Per-species phenotype files for association testing |

### `NB1_Pangenome.ipynb`

Pangenome construction and per-species architecture.

| Section | Contents |
|---|---|
| 2 | PAV and CNV matrices from the per-species OrthoFinder runs |
| 3 | Core, accessory and rare classification by S-curve inflection |
| 4 | Heap's law openness with permutation confidence intervals |
| 5 | Orthogroup annotation consensus tables |
| 6 | Functional enrichment per compartment |
| 7 | Kinship matrices from the recombination-filtered SNP alignments |
| 8 | BGC and GCF presence/absence and copy-number matrices |

### `NB2_PanGWAS.ipynb`

Phenotype-labelled association testing, per species.

| Section | Contents |
|---|---|
| 2 | One-vs-rest phenotype contrasts |
| 3 | EMMA-style linear mixed-model association testing per feature layer |
| 4 | Power analysis, minimum detectable odds ratio at 80% power |
| 5 | Functional enrichment of the significant orthogroups |
| 6 | BGC co-localisation |
| 7 | Cross-species comparison of the results |

### `NB3_Convergence.ipynb`

Whether independent species use the same genes for the same lifestyle.

| Section | Contents |
|---|---|
| 2 | Mapping per-species orthogroups into the genus-level pangenome |
| 3 | Convergence of the significant sets, per layer and species pair |
| 4 | Directional convergence across orthogroups tested in both species |
| 5 | Permutation null model for the overlap |
| 6 | CAZy, protease and secretome family comparisons |
| 7 | Functional convergence of the enriched annotation terms |

### `NB4_RareGenome.ipynb`

The rare compartment, purified of paralogs and fragments.

| Section | Contents |
|---|---|
| 2 | DIAMOND reclassification of every rare orthogroup |
| 3 | Restriction to the No-hit, truly-rare subset |
| 4 | Compartment characterisation: protein length, Pfam coverage, ORF completeness |
| 5 | COG enrichment in the truly-rare compartment |
| 6 | Xenolog screen on GC-outlier orthogroups |
| 7 | Phylogenetic signal and ancestral niche reconstruction |
| 8 | Mash whole-genome distance clustering |
| 9 | Kinship-corrected per-genome burden models |

### `NB5_CoreGenome.ipynb`

The core-compartment counterpart of NB4: are canonical trait genes conserved?

| Section | Contents |
|---|---|
| 2 | Curated literature panels anchored to genus orthogroups by blastp |
| 3 | The experimentally validated subset, and its conservation against background |

## Result folders

Each notebook writes to its own folder:

```
Analysis/
├── NB0_Results/    Metadata, phenotype assignments, geographic map
├── NB1_Results/    PAV/CNV matrices, OG consensus, kinship, BGC/GCF matrices
├── NB2_Results/    Association results, enrichment, power analysis
├── NB3_Results/    Convergence results, null model, family comparisons
├── NB4_Results/    Rare-genome characterisation, phylogenetics, Mash clustering
└── NB5_Results/    Core-genome conservation figures and tables
```

`NB*_Results/` contents are generated and are not tracked, with one exception.
The hand-curated literature inputs under `NB5_Results/` are tracked, because no
code regenerates them:

```
NB5_Results/panels/*_panel.csv                                # curated trait panels, one per species
NB5_Results/panel_pmid_audit.csv                              # per-gene PMID verification
NB5_Results/filtered/panel_cross_species_conservation_SUP.csv
```

## Caches and first run

A new user starts with no cached results. Expensive steps cache their output and
reuse it on later runs, each with its own switch:

| Step | Cache | Switch |
|---|---|---|
| NB2 association testing | `NB2_Results/{species}/pangwas_results/` | rerun the section |
| NB4 DIAMOND classification | `NB4_Results/{species}/diamond_rare_classification.tsv` | `run_diamond` |
| NB4 GC content | `NB4_Results/{species}/og_gc_content.tsv` | `FORCE_GC` |
| NB4 NCBI BLAST | `NB4_Results/rare_xenolog_ncbi_blast_results.tsv` | `FORCE_BLAST` |
| NB5 blastp anchoring | `NB5_Results/cache/panel_blastp.tsv` | `force_blast` |

Section 1 of every notebook reports which inputs are present and creates the
output directories before anything runs.

## External tools

Beyond the conda environment, some sections need a binary on `PATH`:

- NB4 Section 2: `diamond`
- NB4 Section 6: network access for NCBI BLAST
- NB4 Section 8: `mash`
- NB5 Section 2: `blastp` and `makeblastdb`
