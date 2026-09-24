# FunPan - Filamentous Fungi Pangenome Analysis Pipeline

A pipeline for downloading, quality-filtering, annotating and analysing
filamentous fungi pangenomes, with a focus on secondary metabolite biosynthesis
and comparative genomics.

## Contents

1. [Overview](#overview)
2. [Documentation](#documentation)
3. [Environment setup](#environment-setup)
4. [Quick start](#quick-start)
5. [Pipeline architecture](#pipeline-architecture)
6. [Repository layout](#repository-layout)
7. [Citation](#citation)

## Overview

Given a genus and species, FunPan downloads every NCBI assembly for that
species, filters them to a high-quality, non-redundant set, and annotates each
genome structurally and functionally. It then builds a pangenome, predicts
biosynthetic gene clusters, and infers a phylogeny from a recombination-filtered
core alignment. Six Jupyter notebooks analyse the results.

The pipeline is 13 shell scripts, one per step. Each finished step writes a
checkpoint, so an interrupted run resumes where it stopped. Scripts set their
thread count from the machine's core count. Where SRA has RNA-seq data for a
strain, step 1 downloads it as annotation evidence. Assembly metadata (44
columns, in `metadata.csv`) stays with each genome through to the notebooks.

## Documentation

| Document | Covers |
|---|---|
| [`Code/README.md`](Code/README.md) | The 13 pipeline scripts: parameters, inputs, outputs, logs |
| [`Analysis/README.md`](Analysis/README.md) | The six notebooks: running them, order, sections, caches |
| [`Analysis/MODULES.md`](Analysis/MODULES.md) | The seven Python modules |

## Environment setup

One conda environment, `funpan`, about 8 GB, holds everything except the three
exceptions below. Versions are pinned in the command and in `environment.yml`.

```bash
conda create -n funpan -c conda-forge -c bioconda \
  python=3.11 pip pandas biopython=1.79 \
  busco=6.0.0 bbmap fastani=1.34 mash=2.3 \
  funannotate=1.8.17 repeatmodeler=2.0.7 repeatmasker rmblast trf evidencemodeler \
  eggnog-mapper=2.1.13 orthofinder=3.1.0 diamond mafft \
  parsnp=2.1.5 hmmer=3.4 dbcan=5.2.1 iqtree=2.4.0 \
  openjdk=17 ncbi-datasets-cli=18.5.1 gawk bc
```

antiSMASH 8.0.2 and Gubbins conflict with the above, so their scripts build
`Code/.envs/antismash8/` and `Code/.envs/gubbins/` on first run. InterProScan is
downloaded by its own script and needs the OpenJDK above. SignalP 6 is
proprietary: install it by hand with `pip` from the DTU-licensed tarball.

Databases download on first use into `Data/`:

```
Data/
├── BUSCO/              # Lineage datasets
├── BiG-SCAPE/          # Pfam-A.hmm
├── dbCAN/              # CAZyme database
├── EggNOG/             # Orthology database
├── Funannotate/        # Gene prediction databases
└── InterProScan/       # Member-database signatures
```

## Quick start

```bash
conda activate funpan
cd Code
bash run_complete_pipeline.sh      # prompts: fresh, resume or one step; then genus, species
bash run_analysis_notebooks.sh     # once the pipeline finishes
```

Steps also run one at a time. See [`Code/README.md`](Code/README.md) for each
script's arguments, parameters and outputs, and
[`Analysis/README.md`](Analysis/README.md) for the notebooks.

## Pipeline architecture

```
1  download_genome_and_BUSCO       NCBI assemblies, RNA-seq, BUSCO
2  ani_and_filter_genome_QC        fastANI, Mash, quality filter
3  run_funannotate                 structural annotation
   |
   +-- 4  run_eggnog       7  run_interproscan   8  run_signalp   9  run_dbcan
   +-- 5  run_orthofinder
   +-- 6  run_antismash  ->  11  run_bigscape
   +-- 10 run_parsnp      ->  12  run_gubbins  ->  13  run_iqtree_gubbins
```

Steps 4 to 9 are independent of each other and all read the proteins from step 3.
Each completed step writes a marker to `Code/.pipeline_checkpoints/`, so a failed
run resumes where it stopped and the saved genus and species are offered back.

## Repository layout

```
FunPan/
├── environment.yml                          # Conda environment spec
├── Code/                                    # Pipeline scripts, one per step
│   ├── run_complete_pipeline.sh             # Controller
│   ├── orthofinder_on_all_species.sh        # Cross-species OrthoFinder, needed by NB3
│   ├── run_analysis_notebooks.sh            # Executes NB0 to NB5
│   ├── .pipeline_checkpoints/               # Resumption markers
│   ├── .pipeline_logfiles/                  # Execution logs
│   └── .envs/                               # antismash8/, gubbins/
│
├── Data/                                    # Reference databases
│
├── Analysis/
│   ├── funpan_utils.py                      # Shared constants and loaders
│   ├── funpan_classify.py                   # Phenotype classification (NB0)
│   ├── funpan_pangenome.py                  # Pangenome construction (NB1)
│   ├── funpan_gwas.py                       # Pan-GWAS and LMM testing (NB2)
│   ├── funpan_convergence.py                # Convergence testing (NB3)
│   ├── funpan_phylo.py                      # Rare genome and phylogenetics (NB4)
│   ├── funpan_core.py                       # Core genome trait panels (NB5)
│   ├── NB0_DataPrep.ipynb                   # Metadata, QC, ANI, phenotypes
│   ├── NB1_Pangenome.ipynb                  # Pangenome architecture
│   ├── NB2_PanGWAS.ipynb                    # Pan-GWAS and enrichment
│   ├── NB3_Convergence.ipynb                # Cross-species convergence
│   ├── NB4_RareGenome.ipynb                 # Rare genome analysis
│   ├── NB5_CoreGenome.ipynb                 # Core genome conservation
│   └── NB0_Results/ ... NB5_Results/        # Per-notebook outputs
│
└── Species/{genus}/{species}/               # Per-species pipeline outputs
    ├── genome/  protein/  genbank/  gff3/   # Downloaded assemblies
    ├── rna/{strain}/                        # RNA-seq reads
    ├── busco_output/                        # BUSCO per genome
    ├── filtered_genome/  filtered_protein/  # Quality-filtered set
    ├── funannotate_output/                  # Structural annotation
    ├── eggnog_output/                       # Functional annotation
    ├── interproscan_output/                 # Protein domains
    ├── signalp_output/                      # Signal peptides
    ├── dbcan_output/                        # CAZymes
    ├── orthofinder_output/                  # Orthogroups and species tree
    ├── antismash_output/                    # Biosynthetic gene clusters
    ├── bigscape_output/                     # BGC networks
    ├── parsnp_output/                       # Core genome alignment
    ├── gubbins_output/                      # Recombination-filtered alignment
    ├── iqtree_output_snv/                   # SNV phylogeny
    └── metadata.csv                         # Assembly metadata and BUSCO
```

## Citation

Cite the tools the pipeline calls:

- BUSCO: Manni et al., 2021. Molecular Biology and Evolution.
- Funannotate: Palmer & Stajich, 2020. Funannotate v1.8.1. Zenodo.
- EggNOG-mapper: Cantalapiedra et al., 2021. eLife.
- OrthoFinder: Emms & Kelly, 2019. Genome Biology.
- antiSMASH: Blin et al., 2023. Nucleic Acids Research.
- BiG-SCAPE: Navarro-Munoz et al., 2020. Nature Chemical Biology.
- InterProScan: Jones et al., 2014. Bioinformatics.
- SignalP: Teufel et al., 2022. Nature Biotechnology.
- dbCAN: Zheng et al., 2023. Nucleic Acids Research.
- Parsnp: Treangen et al., 2014. Genome Biology.
- Gubbins: Croucher et al., 2015. Nucleic Acids Research.
- IQ-TREE: Minh et al., 2020. Molecular Biology and Evolution.

Tested with *Aspergillus oryzae*, *A. niger*, *A. flavus* and *A. fumigatus*
pangenomes.
