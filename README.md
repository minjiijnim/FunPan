# FunPan - Filamentous Fungi Pangenome Analysis Pipeline

A comprehensive pipeline for downloading, quality-filtering, annotating, and analyzing filamentous fungi pangenomes with a focus on secondary metabolite biosynthesis and comparative genomics.

---

## Table of Contents

1. [Overview](#overview)
2. [Environment Setup](#environment-setup)
3. [Quick Start](#quick-start)
4. [Pipeline Architecture](#pipeline-architecture)
5. [Tools & Scripts Documentation](#tools--scripts-documentation)
6. [Tool-Specific Parameters](#tool-specific-parameters)
7. [Output Structure](#output-structure)
8. [Additional Standalone Scripts](#additional-standalone-scripts)
9. [Analysis Notebooks & Python Modules](#analysis-notebooks--python-modules)
10. [Log File Locations](#log-file-locations)

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

## Tools & Scripts Documentation

### 1. `run_complete_pipeline.sh` - Master Pipeline Controller

**Purpose:** Orchestrates all 13 pipeline steps with progress tracking and error handling.

**Key Features:**
- Interactive step selection (resume, fresh start, specific step)
- Real-time progress bar with colored output
- Automatic CPU detection (uses 75% of available cores)
- Checkpoint-based resumption
- Comprehensive logging to `.pipeline_logfiles/`

**Environment Variables:**
- `PIPELINE_THREADS`: Number of threads to use (default: auto-detected)

**Usage:**
```bash
bash run_complete_pipeline.sh
```

**Checkpoint Files:**
```
.pipeline_checkpoints/
├── step_1.completed
├── step_2.completed
├── ...
└── pipeline_inputs.txt  # Saved genus, species, lineage
```

---

### 2. `download_genome_and_BUSCO.sh` - Genome Acquisition & QC

**Purpose:** Download genomes from NCBI and assess quality with BUSCO.

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PIPELINE_THREADS` | auto (total-4) | CPU threads for BUSCO |
| `BUSCO_THREADS` | `PIPELINE_THREADS` | BUSCO-specific threads |
| `BUSCO_MODE` | `auto` | Lineage detection: `auto`, `auto_euk`, `auto_prok` |
| `BUSCO_OFFLINE` | `0` | Set to `1` for offline mode |
| `BUSCO_LINEAGE` | auto-detected | Force specific lineage (e.g., `eurotiomycetes_odb10`) |

**Process Flow:**
1. **Download:** Uses NCBI Datasets CLI to fetch:
   - Genome assemblies (`.fna`)
   - Protein sequences (`.faa`)
   - GenBank files (`.gbk`)
   - GFF3 annotations (`.gff3`)

2. **RNA-seq Discovery:**
   - Extracts BioSample IDs from assembly metadata
   - Queries SRA for linked RNA-seq runs
   - Downloads FASTQ files from ENA
   - Organizes by strain: `rna/{strain_name}/`

3. **Duplicate Handling:**
   - Identifies GCA/GCF pairs (same genome, GenBank vs RefSeq)
   - Keeps GCF (RefSeq), removes GCA
   - Cleans up orphaned BUSCO results

4. **BUSCO Analysis:**
   - Auto-detects lineage from first genome
   - Reuses lineage for all subsequent genomes
   - Saves lineage to `pipeline_inputs.txt`
   - Cleans up redundant lineage results

5. **Metadata Generation:**
   - Merges BUSCO scores into assembly metadata
   - Extracts isolation source, location, collection date
   - Outputs `metadata.csv` and `busco_updated.jsonl`

**Output:**
```
Species/{genus}/{species}/
├── genome/                      # Genome FASTA files
├── protein/                     # Protein FASTA files
├── genbank/                     # GenBank format files
├── gff3/                        # GFF3 annotation files
├── rna/{strain_name}/           # RNA-seq FASTQ files
│   ├── runs.txt                # List of SRA run IDs
│   └── {SRR}.fastq.gz          # Downloaded reads
├── busco_output/{genome_name}/  # BUSCO results per genome
├── metadata.csv                 # Comprehensive metadata table
├── busco_updated.jsonl          # JSONL with BUSCO scores
└── biosamples.tsv               # Assembly-BioSample mapping
```

---

### 3. `ani_and_filter_genome_QC.sh` - ANI Species Check, Deduplication & QC Filtering

**Purpose:** Verify species identity, remove near-identical duplicates, and filter genomes by assembly quality metrics.

**Usage:**
```bash
bash ani_and_filter_genome_QC.sh Aspergillus/flavus [ANI_THRESHOLD] [MASH_THRESHOLD]
```

**Process Flow:**

1. **ANI Species Verification (fastANI):**
   - Compares all genomes against the GCF reference genome
   - Removes genomes below the ANI threshold (default: 95%)
   - Genomes too divergent to compute ANI are also removed

2. **Mash Deduplication:**
   - Sketches ANI-passed genomes with Mash (sketch size 10,000)
   - Computes all-vs-all Mash distances
   - Clusters near-identical genomes (default threshold: 0.0001 = 99.99% identity)
   - Keeps the highest-quality representative per cluster

3. **QC Filtering (same criteria as before):**

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| **BUSCO Complete** | > 95% | Ensures genome completeness |
| **Contig Count** | ≤ median | Prefers more contiguous assemblies |
| **Contig N50** | ≥ 75th percentile | Prefers longer contigs |
| **GCF Requirement** | ≥ 1 | Always includes at least one RefSeq genome |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ANI_THRESHOLD` | 95 | Minimum ANI% to GCF reference |
| `MASH_THRESHOLD` | 0.0001 | Maximum Mash distance for duplicate detection (~99.99% identity) |

**Dependencies:** `fastANI`, `mash` — auto-installed via conda/mamba if not found on PATH.

**Input:**
- `metadata.csv` from Step 1
- Genome files in `genome/`
- Protein files in `protein/`

**Output:**
```
filtered_genome/     # High-quality genome assemblies
filtered_protein/    # Corresponding protein sequences
qc_work/             # Intermediate files (ANI results, Mash distances, dedup lists)
```

---

### 4. `run_funannotate.sh` - Structural & Functional Annotation

**Purpose:** Predict genes and annotate genomes using ab initio and evidence-based methods.

**Process Flow:**

1. **FASTA Header Renaming:** Converts contig names to `scaffold_1`, `scaffold_2`, etc.
2. **GFF3 Validation:** Checks for contig mismatches, mRNA features, duplicate IDs
3. **Assembly Cleaning:** Removes short contigs (`funannotate clean`)
4. **Repeat Masking:** De novo RepeatModeler library + RepeatMasker soft-masking
5. **Protein Evidence Selection:** Best GCF proteome as reference + strain-specific
6. **RNA-seq Evidence:** Auto-detects RNA-seq files per strain
7. **Gene Prediction:** EvidenceModeler integration with ab initio predictors

**Output:**
```
funannotate_output/{genome}/
├── {genome}.clean.fa            # Cleaned assembly
├── {genome}.softmasked.fa       # Repeat-masked assembly
├── predict_results/
│   ├── {genome}.gbk            # GenBank annotation
│   ├── {genome}.gff3           # GFF3 annotation
│   ├── {genome}.proteins.fa    # Predicted proteins
│   └── {genome}.transcripts.fa # Predicted transcripts
└── logfiles/
```

---

### 5. `run_eggnog.sh` - Functional Annotation

**Purpose:** Assign functional annotations (GO terms, KEGG pathways, domains) via orthology.

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
eggnog_output/
├── {genome}.emapper.annotations  # Tab-delimited functional annotations
├── {genome}.emapper.hits         # DIAMOND hits
└── {genome}.emapper.seed_orthologs
```

**Annotation Fields:** GO terms, KEGG KO numbers and pathways, COG functional categories

---

### 6. `run_orthofinder.sh` - Comparative Genomics

**Purpose:** Identify orthogroups and infer evolutionary relationships.

**Input:** `*.proteins.fa` from Funannotate predict results

**Options:** `-S diamond -M msa -A mafft`

**Output:**
```
orthofinder_output/Results_{date}/
├── Orthogroups/
│   ├── Orthogroups.tsv              # Gene-to-orthogroup mapping
│   ├── Orthogroups.GeneCount.tsv    # Counts per species
│   └── Orthogroups_SingleCopyOrthologues.txt
├── MultipleSequenceAlignments/      # MSA per orthogroup
├── Species_Tree/
│   └── SpeciesTree_rooted.txt       # Newick format
└── Gene_Trees/                      # Individual gene trees
```

---

### 7. `run_antismash.sh` - Secondary Metabolite Detection

**Purpose:** Identify and characterize biosynthetic gene clusters (BGCs).

**Isolated Environment:** Creates separate conda env at `.envs/antismash8` to avoid conflicts.

**antiSMASH Options:**
- `--taxon fungi`: Fungal-specific gene cluster detection
- `--genefinding-tool none`: Uses existing Funannotate annotations
- `--asf`, `--cb-knownclusters`, `--cb-general`, `--cb-subclusters`, `--pfam2go`

**Input:** GenBank files (`.gbk`) from Funannotate predict results

**Output:**
```
antismash_output/{genome}/
├── index.html                   # Interactive report
├── {genome}.gbk                 # Annotated GenBank with BGCs
├── *.region001.gbk              # Individual BGC GenBank files
└── knownclusterblast/           # MIBiG comparisons
```

---

### 8. `run_interproscan.sh` - Protein Domain Annotation

**Purpose:** Scan protein sequences against InterPro signature databases (Pfam, SMART, CDD, PANTHER, etc.).

**Features:**
- Auto-installs InterProScan if not found (downloads ~2.5 GB)
- Optimized tool set by default (excludes Hamap, ProSitePatterns, PRINTS for ~30-40% speedup)
- `TEST_MODE=1` to process only first file for testing

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
interproscan_output/{genome}/
├── {genome}.tsv                 # Tab-separated annotations
├── {genome}.gff3                # GFF3 format annotations
└── {genome}.xml                 # XML format annotations
```

---

### 9. `run_signalp.sh` - Signal Peptide Prediction

**Purpose:** Predict signal peptides in protein sequences using SignalP 6.

**Features:**
- Auto-tunes parallel jobs and torch threads based on CPU cores
- Supports both SignalP 6.x and legacy versions
- `TEST_MODE=1` to process only first file for testing

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
signalp_output/{genome}/
├── prediction_results.txt       # Signal peptide predictions
└── signalp.log                  # Run log
```

---

### 10. `run_dbcan.sh` - CAZyme Annotation

**Purpose:** Annotate carbohydrate-active enzymes (CAZymes) using the dbCAN database.

**Features:**
- Auto-downloads dbCAN database if not present
- Supports multiple detection methods (HMM, DIAMOND, dbCANsub)
- Skips completed samples automatically

**Input:** `*.proteins.fa` from Funannotate predict results

**Output:**
```
dbcan_output/{genome}/
└── overview.tsv                 # CAZyme annotations summary
```

---

### 11. `run_parsnp.sh` - Core Genome Alignment

**Purpose:** Perform rapid core genome alignment using Parsnp.

**Features:**
- Auto-selects best quality reference genome (prefers GCF/RefSeq)
- Requires at least 3 genomes
- Uses funannotate cleaned genomes (`*.clean.fa`)

**Input:** `*.clean.fa` from Funannotate output directories

**Output:**
```
parsnp_output/
├── parsnp.tree                  # Newick phylogenetic tree
├── parsnp.xmfa                  # Core genome alignment (XMFA)
└── parsnp.ggr                   # Gingr visualization file
```

---

### 12. `run_bigscape.sh` - BGC Network Analysis

**Purpose:** Compare BGCs across genomes and group into Gene Cluster Families (GCFs).

**Dependencies:** Auto-clones BiG-SCAPE from GitHub, downloads Pfam-A.hmm

**Process Flow:**
1. Collect `*.region*.gbk` from antiSMASH output
2. Rename with genome accession prefix into `antismash_gbk_unique/`
3. Run BiG-SCAPE clustering (GCF cutoff: 0.3)

**Output:**
```
bigscape_output/
├── index.html                   # Interactive network viewer
├── data_sqlite.db               # SQLite database of results
└── network_files/               # Cytoscape-compatible networks
```

---

### 13. `run_gubbins.sh` - Recombination Detection

**Purpose:** Detect and remove recombination regions from core genome alignment.

**Isolated Environment:** Creates separate conda env at `.envs/gubbins` with Gubbins, IQ-TREE, RAxML-NG, FastTree, and harvesttools.

**Features:**
- Converts Parsnp XMFA to FASTA via harvesttools
- Tries IQ-TREE first, falls back to RAxML-NG if it fails
- Requires at least 3 sequences

**Input:** `parsnp_output/parsnp.xmfa` from Step 10

**Output:**
```
gubbins_output/
├── gubbins.final_tree.tre                       # Recombination-corrected tree
├── gubbins.recombination_predictions.gff        # Detected recombination regions
├── gubbins.filtered_polymorphic_sites.fasta     # Filtered alignment
└── gubbins.node_labelled.final_tree.tre         # Node-labeled tree
```

---

### 14. `run_iqtree_gubbins.sh` - SNV-based Phylogeny

**Purpose:** Build phylogenetic tree from recombination-filtered SNVs using IQ-TREE 2.

**Features:**
- ModelFinder for optimal substitution model selection
- 1000 ultrafast bootstrap replicates with BNNI correction
- Requires IQ-TREE 2.x specifically

**Input:** `gubbins_output/gubbins.filtered_polymorphic_sites.fasta` from Step 12

**Output:**
```
iqtree_output_snv/
├── gubbins_tree.treefile        # Maximum likelihood tree
├── gubbins_tree.iqtree          # Full IQ-TREE report
└── gubbins_tree.log             # Run log
```

---

## Tool-Specific Parameters

**Funannotate:**
- Augustus training: automatic per species
- Protein evidence weight: best quality proteome + strain-specific
- GFF3 evidence weight: 10 (high confidence)
- Minimum protein length: 50 aa (Funannotate default)

**OrthoFinder:**
- Inflation parameter: 1.5 (MCL default)
- E-value threshold: 1e-3
- MSA tool: MAFFT
- Tree inference: FastTree (default) or IQ-TREE

**antiSMASH:**
- Detection strictness: relaxed (default)
- Minimum cluster size: follows antiSMASH defaults
- All detection modules enabled

**BiG-SCAPE:**
- Cutoff: 0.3 (balanced)
  - 0.1-0.2: very strict (many small GCFs)
  - 0.3-0.5: moderate (recommended)
  - 0.6-0.9: loose (few large GCFs)

**InterProScan:**
- Default applications: AntiFam, CDD, Coils, FunFam, Gene3D, MobiDBLite, NCBIfam, PANTHER, Pfam, PIRSF, PIRSR, ProSiteProfiles, SFLD, SMART, SUPERFAMILY
- Override with `INTERPRO_APPS="Pfam,SMART,CDD"`

**SignalP:**
- Organism: `euk` (eukaryote, default)
- Torch threads per job: 3 (default)
- Override parallel jobs with `SIGNALP_JOBS=30`

**dbCAN:**
- Methods: HMM, DIAMOND, dbCANsub (all by default)
- Override with `DBCAN_METHODS="hmm,diamond"`

---

## Output Structure

### Complete Directory Tree

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
├── Analysis/                                # Analysis modules, notebooks & results
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
└── Species/{genus}/{species}/               # Per-species outputs
    ├── genome/                              # Downloaded genomes
    ├── protein/                             # Downloaded proteins
    ├── genbank/                             # GenBank format
    ├── gff3/                                # GFF3 annotations
    │   ├── mismatched_gffs/                # Invalid GFFs
    │   └── non_evm_format/                 # Non-EVM compatible
    ├── rna/{strain}/                       # RNA-seq data
    ├── busco_output/{genome}/              # BUSCO per genome
    ├── filtered_genome/                    # Quality-filtered genomes
    ├── filtered_protein/                   # Filtered proteins
    ├── funannotate_output/{genome}/
    │   ├── *.clean.fa
    │   ├── *.softmasked.fa
    │   └── predict_results/
    │       ├── *.gbk
    │       ├── *.gff3
    │       ├── *.proteins.fa
    │       └── *.transcripts.fa
    ├── eggnog_output/
    │   └── *.emapper.annotations
    ├── orthofinder_output/Results_{date}/
    │   ├── Orthogroups/
    │   ├── MultipleSequenceAlignments/
    │   ├── Species_Tree/
    │   └── Gene_Trees/
    ├── antismash_output/{genome}/
    │   ├── index.html
    │   ├── *.region*.gbk
    │   └── antismash_gbk_unique/           # Renamed BGCs for BiG-SCAPE
    ├── interproscan_output/{genome}/
    │   ├── *.tsv
    │   ├── *.gff3
    │   └── *.xml
    ├── signalp_output/{genome}/
    │   └── prediction_results.txt
    ├── dbcan_output/{genome}/
    │   └── overview.tsv
    ├── parsnp_output/
    │   ├── parsnp.tree
    │   └── parsnp.xmfa
    ├── bigscape_output/
    │   ├── index.html
    │   └── data_sqlite.db
    ├── gubbins_output/
    │   ├── gubbins.final_tree.tre
    │   ├── gubbins.recombination_predictions.gff
    │   └── gubbins.filtered_polymorphic_sites.fasta
    ├── iqtree_output_snv/
    │   └── gubbins_tree.treefile
    ├── metadata.csv                        # Assembly metadata + BUSCO
    └── busco_updated.jsonl                 # JSONL format
```

## Additional Standalone Scripts

### `orthofinder_on_all_species.sh` - Cross-Species OrthoFinder

Runs OrthoFinder on protein files from multiple species simultaneously. Required before NB3 (convergence analysis).

```bash
bash orthofinder_on_all_species.sh
# Prompts for: genus, space-separated species list
# Output: Species/{genus}/all_combined/orthofinder_output/
```

---

## Analysis Notebooks & Python Modules

### Running All Notebooks

Run all six analysis notebooks in sequence from the command line:

```bash
conda activate funpan
cd /datadrive/Code
bash run_analysis_notebooks.sh
```

Options:
- `bash run_analysis_notebooks.sh --from NB2` — skip NB0/NB1, start from NB2
- `bash run_analysis_notebooks.sh --only NB3` — run only NB3
- `bash run_analysis_notebooks.sh --help` — show usage

The script shows a progress bar in the terminal and saves executed outputs in-place to each notebook. Results are also written to `NB*_Results/` folders.

### Notebooks

Six Jupyter notebooks in `/datadrive/Analysis/`, each covering one stage of the downstream analysis. Run sequentially: NB0 first (all others depend on its phenotype outputs), then NB1 (NB2-5 depend on its matrices), then NB2-5 in any order. NB5 is the core-compartment counterpart of NB4. All six run via `run_analysis_notebooks.sh`.

```bash
conda activate funpan
cd /datadrive/Analysis
jupyter notebook NB0_DataPrep.ipynb
```

| Notebook | Description |
|---|---|
| `NB0_DataPrep.ipynb` | QC filtering, ANI species verification, two-tier phenotype classification, geographic mapping |
| `NB1_Pangenome.ipynb` | PAV/CNV matrices, Core/Accessory/Rare classification by S-curve inflection, Heap's law openness, orthogroup annotation tables, functional enrichment, kinship matrices, BGC/GCF matrices |
| `NB2_PanGWAS.ipynb` | One-vs-rest phenotype contrasts, EMMA-style LMM association testing per species and feature layer, power analysis, functional enrichment of significant hits, BGC co-localisation |
| `NB3_Convergence.ipynb` | Species-to-genus orthogroup mapping, convergence of the significant sets, directional convergence across jointly tested orthogroups, permutation null model, CAZy/protease/secretome comparisons, functional convergence |
| `NB4_RareGenome.ipynb` | DIAMOND reclassification of the rare compartment, truly-rare filter, compartment characterisation (protein length, Pfam coverage, ORF completeness), COG enrichment, xenolog screen, phylogenetic signal, Mash clustering, kinship-corrected burden models |
| `NB5_CoreGenome.ipynb` | Curated literature trait panels anchored to genus orthogroups by blastp, then cross-species conservation of the experimentally validated subset |

### Python Modules

All analysis functions live in seven `.py` modules in `/datadrive/Analysis/`. Notebooks import these at the top and contain only configuration, function calls and visualization. Each notebook section opens with a standalone bootstrap cell that sets the paths, imports the modules and reloads that section's inputs, so any section runs on its own from a fresh kernel.

| Module | Functions | Used by |
|---|---|---|
| `funpan_utils.py` | Shared constants, metadata loaders, annotation parsers, QC visualization | All notebooks |
| `funpan_classify.py` | Weighted rules-based phenotype classification, NCBI BioProject enrichment, geocoding | NB0 |
| `funpan_pangenome.py` | Pangenome construction, classification, enrichment, Heap's law, SNP PCA/GRM, BGC/GCF matrices | NB1 |
| `funpan_gwas.py` | LMM association testing (EMMA-style), diagnostic plotting, functional enrichment, power analysis | NB2 |
| `funpan_convergence.py` | Combined pangenome loading, convergent hit identification, annotation transfer | NB3 |
| `funpan_phylo.py` | Rare genome characterisation, phylogenetic signal, ancestral reconstruction, gain/loss, Mash clustering | NB4 |
| `funpan_core.py` | Curated trait panels, orthogroup anchoring by blastp, cross-species core conservation | NB5 |

### Result Folders

Each notebook writes outputs to its own result folder:

```
Analysis/
├── NB0_Results/         # Metadata, phenotype classifications, world map
├── NB1_Results/         # PAV/CNV matrices, OG consensus, kinship, BGC/GCF matrices
├── NB2_Results/         # Pan-GWAS association results, enrichment, power analysis
├── NB3_Results/         # Cross-species convergence, null model, functional comparisons
├── NB4_Results/         # Rare genome characterisation, phylogenetics, Mash clustering
└── NB5_Results/         # Cross-species core-genome conservation (figures & tables)
```

`NB*_Results/` contents are generated and are not tracked, with one exception: the
hand-curated literature inputs under `NB5_Results/` are tracked, because no code
regenerates them.

```
NB5_Results/panels/*_panel.csv                           # curated trait panels, one per species
NB5_Results/panel_pmid_audit.csv                          # per-gene PMID verification
NB5_Results/filtered/panel_cross_species_conservation_SUP.csv
```

---

## Log File Locations

**Pipeline Logs:**
```bash
/datadrive/Code/.pipeline_logfiles/pipeline_YYYYMMDD_HHMMSS.log
```

**Tool-Specific Logs:**
```bash
# Funannotate
Species/{genus}/{species}/funannotate_output/{genome}/logfiles/

# BUSCO
Species/{genus}/{species}/busco_output/{genome}/logs/

# Gubbins
Species/{genus}/{species}/gubbins_output/gubbins.log
```

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
