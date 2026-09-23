# FunPan Python Modules

Seven modules in `Analysis/` hold the analysis functions. The notebooks import
them at the top of every bootstrap cell and contain only configuration, function
calls and visualization.

For the notebooks themselves see [README.md](README.md). For the shell pipeline
see [Code/README.md](../Code/README.md).

---

## Overview

| Module | Purpose | Used by |
|---|---|---|
| `funpan_utils.py` | Shared constants, metadata loaders, annotation parsers, QC visualization | all |
| `funpan_classify.py` | Phenotype classification, NCBI BioProject enrichment, geocoding | NB0 |
| `funpan_pangenome.py` | Pangenome construction, classification, enrichment, Heap's law, kinship, BGC/GCF matrices | NB1 |
| `funpan_gwas.py` | Association testing, power analysis, functional enrichment | NB2 |
| `funpan_convergence.py` | Cross-species convergence in the genus-level pangenome | NB3 |
| `funpan_phylo.py` | Rare-genome characterisation, phylogenetic signal, Mash clustering | NB4 |
| `funpan_core.py` | Curated trait panels and core-compartment conservation | NB5 |

## Conventions

**Paths are parameters, not constants.** Functions take the directories they
read and write as arguments, passed in from the notebook's editable path block.
A reader changing where their data lives edits the notebook, not the module.

**Caches are explicit.** Any function that caches an expensive result takes a
`force` argument and says in its docstring where the cache lives.

**Figures return their handle.** Plotting functions save a PNG and return the
figure, so the notebook calls `plt.show()` rather than the module.

Every public function carries a docstring giving its purpose, its parameters and
what it returns. Those docstrings are the detailed reference; this file is the
map.

---

## `funpan_utils.py`

The base layer every other module builds on. Genus and species constants, the
shared colour palettes, the ANI and contamination exclusion sets, and the
loaders that read the pipeline's per-genome outputs: assembly metadata, BUSCO
scores, eggNOG annotations, InterProScan signatures, SignalP predictions and
dbCAN CAZyme calls. Also holds the QC and metadata figures used by NB0.

Imported by all six notebooks and by every other module.

## `funpan_classify.py`

Turns NCBI metadata into the phenotype labels the association tests use. A
weighted, rules-based classifier reads eight metadata fields and assigns each
isolate a provenance label (where the strain was obtained) and a phenotype label
(its biological association). Includes the NCBI BioProject and WGS-title
enrichment that fills gaps in the assembly records, strain-name rules for
isolates with no isolation source, and the geocoding that places strains on the
world map.

NCBI access takes `email` and `api_key` arguments; both default to `None`.

Used by NB0.

## `funpan_pangenome.py`

Builds the pangenome and describes its architecture. Reads the per-species
OrthoFinder output into presence/absence and copy-number matrices, derives the
core, accessory and rare boundaries per species by the S-curve inflection
method, and fits Heap's law to assess openness. Also assembles the orthogroup
annotation consensus tables, runs functional enrichment per compartment,
computes the SNP principal components and kinship matrices from the
recombination-filtered alignments, and builds the BGC and GCF matrices.

Used by NB1.

## `funpan_gwas.py`

Phenotype-labelled association testing. Fits an EMMA-style linear mixed model
per species, per contrast and per feature layer, with the kinship matrix as the
random-effect covariance, and applies Benjamini-Hochberg correction. Includes
the power analysis that reports the minimum detectable odds ratio at 80% power
under each contrast's case-control imbalance, the functional enrichment of the
significant orthogroups, and the BGC co-localisation test.

Used by NB2.

## `funpan_convergence.py`

Tests whether independent species reach the same lifestyle through the same
genes. Maps each species' results into the genus-level pangenome by protein
lookup, intersects the significant sets per layer and species pair, and
classifies the jointly tested orthogroups by whether the sign of effect agrees.
Includes the permutation null models, the CAZy, protease and secretome family
comparisons, and the Jaccard-based functional convergence testing.

Used by NB3.

## `funpan_phylo.py`

The rare compartment and the phylogenetic context. Reclassifies every rare
orthogroup by aligning it against its own species' core and accessory proteins,
separating genuine lineage-private genes from paralogs and gene fragments.
Characterises the resulting compartment by protein length, domain coverage and
gene-model completeness, screens GC-outlier orthogroups against NCBI for
horizontal transfer candidates, and fits the kinship-corrected per-genome burden
models. Also holds the Fritz and Purvis phylogenetic signal statistic, the
ancestral niche reconstruction, and the Mash whole-genome clustering.

Used by NB4.

## `funpan_core.py`

The core-compartment counterpart. Loads the curated literature trait panels,
anchors each reference protein to a genus orthogroup by blastp against that
species' RefSeq proteome, and derives the experimentally validated subset from
the per-gene PMID audit. Then measures how many species each validated gene is
core in, against a background of the orthogroups that are not in the panel.

Used by NB5.
