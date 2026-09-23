# FunPan - Filamentous Fungi Pangenome Analysis Pipeline

A comprehensive pipeline for downloading, quality-filtering, annotating, and analyzing filamentous fungi pangenomes with a focus on secondary metabolite biosynthesis and comparative genomics.

---

## Table of Contents

1. [Overview](#overview)
2. [Documentation](#documentation)
3. [Environment Setup](#environment-setup)
4. [Quick Start](#quick-start)
5. [Pipeline Architecture](#pipeline-architecture)
6. [Repository Layout](#repository-layout)
7. [Citation](#citation)

---

## Overview

**FunPan** is an automated 13-step pipeline designed for:
- Downloading genome assemblies from NCBI
- Quality assessment using BUSCO
- ANI species verification, Mash deduplication, and quality-based genome filtering
- Structural annotation with Funannotate
- Functional annotation with EggNOG-mapper, InterProScan, SignalP, and dbCAN
- Comparative genomics with OrthoFinder
- Secondary metabolite cluster detection with antiSMASH
- Biosynthetic gene cluster (BGC) network analysis with BiG-SCAPE
- Core genome alignment and phylogenetics with Parsnp, Gubbins, and IQ-TREE

**Key Features:**
- Checkpoint-based resumable execution
- Parallel processing with automatic CPU detection
- Comprehensive metadata extraction and tracking
- ANI species verification and Mash deduplication
- Quality-based genome filtering
- RNA-seq evidence integration for annotation

---

## Documentation

This README covers the project as a whole. Component reference lives beside the
code it documents:

| Document | Covers |
|---|---|
| [`Code/README.md`](Code/README.md) | The 13-step shell pipeline:each  script's purpose, parameters, inputs and outputs, plus tool-specific settings and log locations |
| [`Analysis/README.md`](Analysis/README.md) | The six analysis notebooks: how to run them, their dependency order, section-by-section contents, result folders and caches |
| [`Analysis/MODULES.md`](Analysis/MODULES.md) | The seven Python modules the notebooks import: what each covers and which notebook uses it |

---

## Environment Setup

### Conda Environment

**Environment Name:** `funpan`

A minimal conda environment (~8 GB) containing all tools needed to run the full pipeline. Created with:

```bash
conda create -n funpan -c conda-forge -c bioconda \
  python=3.11 pip pandas biopython=1.79 \
  busco=6.0.0 bbmap fastani=1.34 mash=2.3 \
  funannotate=1.8.17 repeatmodeler=2.0.7 repeatmasker rmblast trf evidencemodeler \
  eggnog-mapper=2.1.13 orthofinder=3.1.0 diamond mafft \
  parsnp=2.1.5 hmmer=3.4 dbcan=5.2.1 iqtree=2.4.0 \
  openjdk=17 ncbi-datasets-cli=18.5.1 gawk bc
```

**Tools included in `funpan`:**
- BUSCO 6.0 + BBTools (genome quality assessment)
- fastANI + Mash (species verification and deduplication)
- Funannotate 1.8.17 + RepeatModeler/RepeatMasker (genome annotation)
- EggNOG-mapper 2.1.13 (functional annotation)
- OrthoFinder 3.1.0 + DIAMOND + MAFFT (orthology inference)
- HMMER 3.4 (HMM searches for dbCAN and BiG-SCAPE)
- dbCAN 5.2.1 (CAZyme annotation)
- Parsnp 2.1.5 (core genome alignment)
- IQ-TREE 2.4.0 (phylogenetic inference)
- NCBI Datasets CLI (genome download)
- OpenJDK 17 (for InterProScan)

**Managed separately (auto-created by their scripts):**
- antiSMASH 8.0.2 → `.envs/antismash8/`
- Gubbins → `.envs/gubbins/`

**Not included (proprietary):**
- SignalP 6 — requires manual `pip install` from DTU-licensed tarball

### Database Requirements

The pipeline automatically downloads and maintains databases in `/datadrive/Data/`:

```
Data/
├── BUSCO/              # BUSCO lineage datasets
├── BiG-SCAPE/          # Pfam-A.hmm for domain detection
├── dbCAN/              # CAZyme database
├── EggNOG/             # EggNOG orthology database
├── Funannotate/        # Gene prediction databases (UniProt, etc.)
├── InterProScan/       # InterProScan signatures (auto-installed)
└── GOATOOLS/           # GO term enrichment (tbc)
```

---

## Quick Start

### Run Complete Pipeline

```bash
conda activate funpan
cd /datadrive/Code
bash run_complete_pipeline.sh
```

**Interactive Prompts:**
1. Choose to start fresh, resume, or start from a specific step
2. Enter genus name (e.g., `Aspergillus`)
3. Enter species name (e.g., `oryzae`)

### Run Individual Steps

Each step can be run independently:

```bash
# Step 1: Download and BUSCO
bash download_genome_and_BUSCO.sh

# Step 2: ANI check, deduplication, and QC filtering
bash ani_and_filter_genome_QC.sh

# Step 3: Annotation
bash run_funannotate.sh

# Step 4: Functional annotation (EggNOG)
bash run_eggnog.sh

# Step 5: Comparative genomics (OrthoFinder)
bash run_orthofinder.sh

# Step 6: Secondary metabolites (antiSMASH)
bash run_antismash.sh

# Step 7: Protein domains (InterProScan)
bash run_interproscan.sh

# Step 8: Signal peptides (SignalP)
bash run_signalp.sh

# Step 9: CAZymes (dbCAN)
bash run_dbcan.sh

# Step 10: Core genome alignment (Parsnp)
bash run_parsnp.sh

# Step 11: BGC networks (BiG-SCAPE)
bash run_bigscape.sh

# Step 12: Recombination detection (Gubbins)
bash run_gubbins.sh

# Step 13: SNV-based phylogeny (IQ-TREE)
bash run_iqtree_gubbins.sh
```

Each script is documented in [`Code/README.md`](Code/README.md).

### Run the Analysis Notebooks

Once the pipeline has finished, the six downstream notebooks run with:

```bash
conda activate funpan
cd /datadrive/Code
bash run_analysis_notebooks.sh
```

See [`Analysis/README.md`](Analysis/README.md) for their order, contents and
caching behaviour.

---

## Pipeline Architecture

### Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│ Step 1: Download Genome & BUSCO                                 │
│  - Download from NCBI (genome, protein, GFF3, GenBank)         │
│  - Auto-discover and download RNA-seq data from SRA/ENA        │
│  - Run BUSCO quality assessment                                 │
│  - Remove GCA/GCF duplicates (keep RefSeq)                     │
│  - Generate comprehensive metadata CSV                          │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 2: ANI Check, Deduplication & QC Filtering                 │
│  - Verify species identity via fastANI against GCF reference   │
│  - Remove near-identical genomes via Mash distance clustering  │
│  - Filter by BUSCO completeness (>95%)                         │
│  - Filter by contig count (≤ median) and N50 (≥ 75th pctl)    │
│  - Ensure at least 1 GCF (RefSeq) genome included              │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 3: Funannotate Annotation                                  │
│  - Rename contigs to scaffold_N format                         │
│  - Clean assemblies (remove small/bad contigs)                 │
│  - Soft-mask repeats (RepeatModeler + RepeatMasker)            │
│  - Gene prediction with protein & RNA-seq evidence             │
│  - Structural annotation (CDS, mRNA, exons)                    │
└─────────────────────────────────────────────────────────────────┘
                            │
            ┌───────────────┼───────────────────────────┐
            │               │               │           │
            ▼               ▼               ▼           ▼
┌───────────────┐ ┌─────────────┐ ┌─────────────┐ ┌──────────┐
│ Step 4: EggNOG│ │ Step 5:     │ │ Step 6:     │ │Steps 7-9:│
│ Functional    │ │ OrthoFinder │ │ antiSMASH   │ │InterPro, │
│ Annotation    │ │ Comparative │ │ Secondary   │ │SignalP,  │
│               │ │ Genomics    │ │ Metabolites │ │dbCAN     │
└───────────────┘ └─────────────┘ └──────┬──────┘ └──────────┘
                                         │
                                         ▼
                              ┌─────────────────┐
                              │ Step 11:        │
                              │ BiG-SCAPE       │
                              │ BGC Networks    │
                              └─────────────────┘

            ┌─────────────────────────┐
            │ Step 10: Parsnp         │
            │ Core Genome Alignment   │
            └────────────┬────────────┘
                         │
                         ▼
            ┌─────────────────────────┐
            │ Step 12: Gubbins        │
            │ Recombination Detection │
            └────────────┬────────────┘
                         │
                         ▼
            ┌─────────────────────────┐
            │ Step 13: IQ-TREE        │
            │ SNV-based Phylogeny     │
            └─────────────────────────┘
```

### Checkpointing System

The pipeline uses a checkpoint system stored in `/datadrive/Code/.pipeline_checkpoints/`:
- Each completed step creates a checkpoint marker
- Pipeline can resume from last incomplete step
- User inputs (genus, species) are saved for resumption
- Prevents re-running completed steps accidentally

---

## Repository Layout

```
datadrive/
├── README.md                               # This file
├── environment.yml                         # Conda environment spec (funpan)
├── Code/                                    # Pipeline scripts
│   ├── run_complete_pipeline.sh             # Master pipeline controller
│   ├── download_genome_and_BUSCO.sh         # Step 1
│   ├── ani_and_filter_genome_QC.sh           # Step 2
│   ├── run_funannotate.sh                   # Step 3
│   ├── run_eggnog.sh                        # Step 4
│   ├── run_orthofinder.sh                   # Step 5
│   ├── run_antismash.sh                     # Step 6
│   ├── run_interproscan.sh                  # Step 7
│   ├── run_signalp.sh                       # Step 8
│   ├── run_dbcan.sh                         # Step 9
│   ├── run_parsnp.sh                        # Step 10
│   ├── run_bigscape.sh                      # Step 11
│   ├── run_gubbins.sh                       # Step 12
│   ├── run_iqtree_gubbins.sh                # Step 13
│   ├── orthofinder_on_all_species.sh        # Optional: cross-species OrthoFinder
│   ├── run_analysis_notebooks.sh            # Runs the NB0-NB5 notebooks
│   ├── .pipeline_checkpoints/               # Resumption markers
│   ├── .pipeline_logfiles/                  # Execution logs
│   └── .envs/                               # Isolated conda envs
│       ├── antismash8/
│       └── gubbins/
│
├── Data/                                    # Reference databases
│   ├── BUSCO/busco_downloads/
│   ├── BiG-SCAPE/Pfam-A.hmm
│   ├── dbCAN/
│   ├── EggNOG/
│   ├── Funannotate/
│   └── InterProScan/
│
├── Analysis/                                # Analysis notebooks, modules & results
│   ├── README.md                            # Notebook reference
│   ├── MODULES.md                           # Python module reference
│   ├── funpan_utils.py                      # Shared constants, loaders, QC viz
│   ├── funpan_classify.py                   # Phenotype classification (NB0)
│   ├── funpan_pangenome.py                  # Pangenome construction (NB1)
│   ├── funpan_gwas.py                       # Pan-GWAS / LMM testing (NB2)
│   ├── funpan_convergence.py                # Convergence testing (NB3)
│   ├── funpan_phylo.py                      # Rare-genome & phylogenetics (NB4)
│   ├── funpan_core.py                       # Core-genome trait panels (NB5)
│   ├── NB0_DataPrep.ipynb                   # Data prep, QC, ANI, phenotypes
│   ├── NB1_Pangenome.ipynb                  # Pangenome architecture
│   ├── NB2_PanGWAS.ipynb                    # Pan-GWAS & enrichment
│   ├── NB3_Convergence.ipynb                # Cross-species convergence
│   ├── NB4_RareGenome.ipynb                 # Rare genome analysis
│   ├── NB5_CoreGenome.ipynb                 # Core genome conservation
│   └── NB0_Results/ ... NB5_Results/        # Per-notebook output folders
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
    ├── orthofinder_output/                  # Orthogroups & species tree
    ├── antismash_output/                    # Biosynthetic gene clusters
    ├── bigscape_output/                     # BGC networks
    ├── parsnp_output/                       # Core genome alignment
    ├── gubbins_output/                      # Recombination-filtered alignment
    ├── iqtree_output_snv/                   # SNV phylogeny
    └── metadata.csv                         # Assembly metadata + BUSCO
```

Per-file contents of each pipeline output directory are documented in
[`Code/README.md`](Code/README.md); the analysis result folders in
[`Analysis/README.md`](Analysis/README.md).

---

## Citation

If you use FunPan in your research, please cite the individual tools:

- **BUSCO:** Manni et al., 2021. BUSCO update. Molecular Biology and Evolution.
- **Funannotate:** Palmer & Stajich, 2020. Funannotate v1.8.1. Zenodo.
- **EggNOG-mapper:** Cantalapiedra et al., 2021. eLife.
- **OrthoFinder:** Emms & Kelly, 2019. Genome Biology.
- **antiSMASH:** Blin et al., 2023. Nucleic Acids Research.
- **BiG-SCAPE:** Navarro-Munoz et al., 2020. Nature Chemical Biology.
- **InterProScan:** Jones et al., 2014. Bioinformatics.
- **SignalP:** Teufel et al., 2022. Nature Biotechnology.
- **dbCAN:** Zheng et al., 2023. Nucleic Acids Research.
- **Parsnp:** Treangen et al., 2014. Genome Biology.
- **Gubbins:** Croucher et al., 2015. Nucleic Acids Research.
- **IQ-TREE:** Minh et al., 2020. Molecular Biology and Evolution.

---

**Last Updated:** 2026-09-24
**Tested With:** Aspergillus oryzae, A. niger, A. flavus, A. fumigatus pangenomes
