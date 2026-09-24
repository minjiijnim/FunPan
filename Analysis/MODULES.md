# FunPan Python Modules

Seven modules in `Analysis/` hold the analysis functions. Notebooks import them
and contain only configuration, calls and figures.

Docstrings are the reference. This file says which module does what.

See also: [notebooks](README.md), [shell pipeline](../Code/README.md).

## Modules

| Module | Contents | Used by |
|---|---|---|
| `funpan_utils.py` | Species constants, ANI exclusions, colour palettes. Loaders for assembly metadata, BUSCO, eggNOG, InterProScan, SignalP, dbCAN | all |
| `funpan_classify.py` | Weighted phenotype classifier over 8 metadata fields, NCBI BioProject and WGS-title enrichment, geocoding | NB0 |
| `funpan_pangenome.py` | PAV/CNV matrices, S-curve core/accessory/rare thresholds, Heap's law, orthogroup annotation tables, enrichment, SNP PCs and kinship, BGC/GCF matrices | NB1 |
| `funpan_gwas.py` | EMMA-style LMM association testing, BH correction, power analysis, enrichment of significant hits, BGC co-localisation | NB2 |
| `funpan_convergence.py` | Species-to-genus orthogroup mapping, set overlap and effect-direction agreement, permutation nulls, CAZy/protease/secretome comparisons, Jaccard functional convergence | NB3 |
| `funpan_phylo.py` | DIAMOND rare-orthogroup reclassification, compartment characterisation, GC-outlier xenolog screen, kinship-corrected burden models, Fritz and Purvis D, Fitch reconstruction, Mash clustering | NB4 |
| `funpan_core.py` | Curated trait panels, blastp anchoring to genus orthogroups, PMID-validated subset, conservation against background | NB5 |

## Conventions

- **Paths are parameters.** Functions take their input and output directories as
  arguments, supplied by the notebook's path block. No module hard-codes a data
  directory.
- **Defaults are relative to the repository.** `funpan_utils` derives
  `FUNPAN_ROOT` from its own file location, giving `SPECIES_ROOT`, `ANALYSIS_ROOT`
  and `DATA_ROOT`. Functions that take no explicit directory fall back to these,
  so the modules work wherever the repository is cloned. Set the `FUNPAN_ROOT`
  environment variable to override, or call `species_path(genus, species, root=...)`.
- **Caches take a `force` argument.** The docstring names the cache file.
- **Plotting functions save a PNG and return the figure.** The notebook calls
  `plt.show()`.
- `funpan_classify` takes `email` and `api_key` for NCBI access. Both default to
  `None`.
